#!/usr/bin/env python3
"""
Native Redocking Optimizer using GNINA (CNN-scored docking)
Solves the flipped-pose problem by using GNINA's CNN scoring function
which discriminates correct orientations from decoys.

Targets:
  - 2ZSH: GID1A receptor with GA3 (Gibberellic acid 3)
  - 2ZSI: GID1A receptor with GA4 (Gibberellic acid 4)

Protonation:
  - Receptor: PDB2PQR + PROPKA (residue-specific pKa at target pH)
  - Ligand: Gypsum-DL (Dimorphite-DL + 3D conformer generation)

Scoring: GNINA CNN scoring (rescore mode by default)
Box: autobox_ligand with variable padding, plus manual sweep
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import urllib.request
from datetime import datetime

import numpy as np

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem, rdMolAlign
    RDLogger.DisableLog("rdApp.*")
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

# ── Configuration ──
PDB_IDS = {
    "2ZSH": {"ligand_resname": "GA3", "ligand_name": "GA3"},
    "2ZSI": {"ligand_resname": "GA4", "ligand_name": "GA4"},
}

GA_SMILES = {
    "GA3": "C[C@]12[C@@H](O)C=C[C@@]3(OC1=O)[C@@H]4CC[C@]5(O)C[C@]4(CC5=C)[C@H]([C@H]23)C(O)=O",
    "GA4": "C[C@@]12[C@H](CC[C@@]3([C@@H]1[C@@H]([C@]45[C@H]3CC[C@H](C4)C(=C)C5)C(=O)O)OC2=O)O",
}

PH_VALUES = [7.2, 7.4]
EXHAUSTIVENESS = 20
GNINA_URL = "https://github.com/gnina/gnina/releases/download/v1.1/gnina"
ACTIVE_SITE_RESIDUES = [119, 24, 27, 31, 126, 323, 239, 319, 116, 191, 244, 127, 238, 218]


def find_gnina(project_dir):
    """Find or download GNINA binary."""
    gnina_path = os.path.join(project_dir, "gnina")
    if os.path.isfile(gnina_path) and os.access(gnina_path, os.X_OK):
        return gnina_path

    g = shutil.which("gnina")
    if g:
        return g

    print("  Downloading GNINA...")
    try:
        urllib.request.urlretrieve(GNINA_URL, gnina_path)
        os.chmod(gnina_path, 0o755)
        result = subprocess.run([gnina_path, "--version"], capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            print(f"  GNINA downloaded: {result.stdout.strip()}")
            return gnina_path
    except Exception as e:
        print(f"  GNINA download failed: {e}")
        try:
            url2 = "https://github.com/gnina/gnina/releases/download/v1.3.2/gnina.1.3.2"
            urllib.request.urlretrieve(url2, gnina_path)
            os.chmod(gnina_path, 0o755)
            return gnina_path
        except Exception as e2:
            print(f"  Fallback download also failed: {e2}")

    return None


def find_obabel():
    """Find OpenBabel executable."""
    ob = shutil.which("obabel")
    if ob:
        return ob
    for p in ["/usr/bin/obabel", "/usr/local/bin/obabel",
              os.path.expanduser("~/miniconda3/bin/obabel"),
              os.path.expanduser("~/.conda/envs/ag4_full/bin/obabel")]:
        if os.path.isfile(p):
            return p
    return None


def find_vina():
    """Find Vina as fallback."""
    home = os.path.expanduser("~")
    direct = os.path.join(home, "final_project", "autogrow4", "autogrow", "docking",
                          "docking_executables", "vina", "autodock_vina_1_1_2_linux_x86", "bin", "vina")
    if os.path.isfile(direct):
        return direct
    return shutil.which("vina")


def download_pdb(pdb_id, output_dir):
    """Download PDB file from RCSB."""
    pdb_path = os.path.join(output_dir, f"{pdb_id}.pdb")
    if os.path.exists(pdb_path) and os.path.getsize(pdb_path) > 1000:
        print(f"  {pdb_id}.pdb already exists, skipping download")
        return pdb_path
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"  Downloading {pdb_id} from RCSB...")
    urllib.request.urlretrieve(url, pdb_path)
    if os.path.getsize(pdb_path) < 1000:
        print(f"  WARNING: Downloaded file is very small, may be an error page")
    return pdb_path


def extract_ligand_and_receptor(pdb_path, ligand_resname, output_dir):
    """Extract ligand and receptor from PDB file."""
    pdb_id = os.path.basename(pdb_path).replace(".pdb", "")
    ligand_pdb = os.path.join(output_dir, f"{pdb_id}_ligand_{ligand_resname}.pdb")
    receptor_pdb = os.path.join(output_dir, f"{pdb_id}_receptor.pdb")

    ligand_lines = []
    receptor_lines = []
    ligand_coords = []

    with open(pdb_path) as f:
        for line in f:
            record = line[:6].strip()
            if record in ("ATOM", "HETATM"):
                resname = line[17:20].strip()
                if resname == ligand_resname:
                    ligand_lines.append(line)
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        element = line[76:78].strip() if len(line) > 76 else ""
                        if element != "H":
                            ligand_coords.append([x, y, z])
                    except (ValueError, IndexError):
                        pass
                elif resname not in ("HOH", "WAT", "TIP", "SOL"):
                    receptor_lines.append(line)
            elif record in ("TER", "END"):
                receptor_lines.append(line)

    with open(ligand_pdb, "w") as f:
        f.writelines(ligand_lines)
        f.write("END\n")

    with open(receptor_pdb, "w") as f:
        f.writelines(receptor_lines)

    ligand_coords = np.array(ligand_coords)
    center = ligand_coords.mean(axis=0) if len(ligand_coords) > 0 else np.array([0, 0, 0])

    print(f"  Extracted {len(ligand_lines)} ligand atoms, {len(receptor_lines)} receptor lines")
    print(f"  Ligand center: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")
    print(f"  Crystal heavy atoms: {len(ligand_coords)}")

    return ligand_pdb, receptor_pdb, center, ligand_coords


def protonate_receptor(receptor_pdb, ph, output_dir):
    """Protonate receptor using PDB2PQR + PROPKA, fallback to OpenBabel."""
    base = os.path.basename(receptor_pdb).replace(".pdb", "")
    pqr_path = os.path.join(output_dir, f"{base}_pH{ph}.pqr")
    pdb_out = os.path.join(output_dir, f"{base}_pH{ph}.pdb")

    try:
        cmd = [sys.executable, "-m", "pdb2pqr",
               "--ff", "PARSE",
               "--titration-state-method", "propka",
               "--with-ph", str(ph),
               "--keep-chain",
               receptor_pdb, pqr_path]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if os.path.exists(pqr_path) and os.path.getsize(pqr_path) > 100:
            print(f"  PDB2PQR+PROPKA at pH {ph}: OK")
            _pqr_to_pdb(pqr_path, pdb_out)
            return pdb_out
        else:
            print(f"  PDB2PQR failed: {result.stderr[:200]}")
    except Exception as e:
        print(f"  PDB2PQR error ({e}), falling back to OpenBabel")

    obabel = find_obabel()
    if obabel:
        cmd = [obabel, receptor_pdb, "-O", pdb_out, "-p", str(ph)]
        subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if os.path.exists(pdb_out) and os.path.getsize(pdb_out) > 100:
            print(f"  OpenBabel protonation at pH {ph}: OK")
            return pdb_out

    shutil.copy2(receptor_pdb, pdb_out)
    print(f"  WARNING: Using unprotonated receptor")
    return pdb_out


def _pqr_to_pdb(pqr_path, pdb_out):
    """Convert PQR to PDB format with proper element symbols."""
    lines = []
    with open(pqr_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                pdb_line = line[:54]
                pdb_line += "  1.00  0.00"
                atom_name = line[12:16].strip()
                element = ""
                for ch in atom_name:
                    if ch.isalpha():
                        element = ch.upper()
                        if len(atom_name) >= 2 and atom_name[1].isalpha() and atom_name[1].islower():
                            element = atom_name[0].upper() + atom_name[1].lower()
                        break
                pdb_line += " " * max(0, 76 - len(pdb_line))
                pdb_line += f"  {element:>2}"
                lines.append(pdb_line.rstrip() + "\n")
            elif line.startswith(("TER", "END")):
                lines.append(line)
    with open(pdb_out, "w") as f:
        f.writelines(lines)
        if not any(l.startswith("END") for l in lines):
            f.write("END\n")


def prepare_ligand_gypsymdl(smiles, ph, output_dir, name="ligand"):
    """Prepare ligand using Gypsum-DL for proper protonation and 3D generation."""
    gypsum_path = os.path.expanduser(
        "~/final_project/autogrow4/autogrow/operators/convert_files/gypsum_dl/run_gypsum_dl.py"
    )
    smi_file = os.path.join(output_dir, f"{name}_input.smi")
    gypsum_out = os.path.join(output_dir, f"{name}_gypsum")
    output_sdf = os.path.join(output_dir, f"{name}_pH{ph}.sdf")

    with open(smi_file, "w") as f:
        f.write(f"{smiles}\t{name}\n")

    if os.path.exists(gypsum_path):
        os.makedirs(gypsum_out, exist_ok=True)
        try:
            cmd = [sys.executable, gypsum_path,
                   "--source", smi_file,
                   "--output_folder", gypsum_out,
                   "--min_ph", str(ph - 0.2),
                   "--max_ph", str(ph + 0.2),
                   "--pka_precision", "0.5",
                   "--thoroughness", "3",
                   "--num_processors", "1",
                   "--max_variants_per_compound", "1"]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            print(f"  Gypsum-DL stdout: {result.stdout[-200:]}" if result.stdout else "  Gypsum-DL: no stdout")

            import glob
            sdf_files = glob.glob(os.path.join(gypsum_out, "*.sdf"))
            if sdf_files:
                shutil.copy2(sdf_files[0], output_sdf)
                print(f"  Gypsum-DL ligand prep at pH {ph}: OK ({os.path.basename(sdf_files[0])})")
                return output_sdf
            else:
                print(f"  Gypsum-DL produced no SDF output")
                if result.stderr:
                    print(f"  stderr: {result.stderr[-300:]}")
        except Exception as e:
            print(f"  Gypsum-DL error: {e}")
    else:
        print(f"  Gypsum-DL not found at {gypsum_path}")

    return _prepare_ligand_obabel(smiles, ph, output_dir, name)


def _prepare_ligand_obabel(smiles, ph, output_dir, name):
    """Fallback: prepare ligand with OpenBabel."""
    obabel = find_obabel()
    if not obabel:
        return None

    output_sdf = os.path.join(output_dir, f"{name}_pH{ph}.sdf")

    for args in [["--gen3d", "-h", "-p", str(ph)],
                 ["--gen3d", "-h"],
                 ["--gen3d", "--best", "-h"]]:
        cmd = [obabel, "-:" + smiles, "-osdf", "-O", output_sdf] + args
        subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if os.path.exists(output_sdf) and os.path.getsize(output_sdf) > 10:
            print(f"  OpenBabel ligand prep at pH {ph}: OK")
            return output_sdf

    print(f"  ERROR: All ligand preparation methods failed for {name}")
    return None


def prepare_ligand_from_crystal(crystal_pdb, ph, output_dir, name="ligand"):
    """Prepare ligand from crystal structure with pH-appropriate hydrogens."""
    obabel = find_obabel()
    if not obabel:
        return None

    output_sdf = os.path.join(output_dir, f"{name}_crystal_pH{ph}.sdf")
    cmd = [obabel, crystal_pdb, "-osdf", "-O", output_sdf, "-h", "-p", str(ph)]
    subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    if os.path.exists(output_sdf) and os.path.getsize(output_sdf) > 10:
        print(f"  Crystal ligand -> SDF at pH {ph}: OK")
        return output_sdf

    cmd = [obabel, crystal_pdb, "-osdf", "-O", output_sdf, "-h"]
    subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    if os.path.exists(output_sdf) and os.path.getsize(output_sdf) > 10:
        print(f"  Crystal ligand -> SDF (no pH adj): OK")
        return output_sdf

    return None


def run_gnina_dock(gnina_path, receptor_pdb, ligand_sdf, autobox_ligand,
                   center, size, exhaustiveness, output_dir, tag="",
                   cnn_scoring="rescore", num_modes=9, autobox_add=4.0):
    """Run GNINA docking and return parsed results."""
    out_sdf = os.path.join(output_dir, f"docked{tag}.sdf")
    log_path = os.path.join(output_dir, f"gnina{tag}.log")

    cmd = [gnina_path,
           "-r", receptor_pdb,
           "-l", ligand_sdf,
           "--exhaustiveness", str(exhaustiveness),
           "--num_modes", str(num_modes),
           "-o", out_sdf,
           "--log", log_path,
           "--cnn_scoring", cnn_scoring,
           "--no_gpu"]

    if autobox_ligand and os.path.exists(autobox_ligand):
        cmd.extend(["--autobox_ligand", autobox_ligand,
                     "--autobox_add", str(autobox_add)])
    elif center is not None and size is not None:
        cmd.extend([
            "--center_x", str(center[0]),
            "--center_y", str(center[1]),
            "--center_z", str(center[2]),
            "--size_x", str(size[0]),
            "--size_y", str(size[1]),
            "--size_z", str(size[2]),
        ])

    print(f"    CMD: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if result.returncode != 0:
            print(f"  GNINA failed (exit {result.returncode}):")
            print(f"  STDERR: {result.stderr[:500]}")
            print(f"  STDOUT: {result.stdout[:500]}")
            return [], None
    except subprocess.TimeoutExpired:
        print(f"  GNINA timed out")
        return [], None

    if not os.path.exists(out_sdf) or os.path.getsize(out_sdf) < 10:
        print(f"  GNINA produced no output file")
        return [], None

    poses = _parse_gnina_log(log_path)
    return poses, out_sdf


def run_vina_dock_fallback(vina_path, receptor_pdbqt, ligand_pdbqt,
                           center, size, exhaustiveness, output_dir, tag=""):
    """Fallback to Vina if GNINA is unavailable."""
    out_path = os.path.join(output_dir, f"docked{tag}.pdbqt")
    log_path = os.path.join(output_dir, f"vina{tag}.log")

    cmd = [vina_path,
           "--receptor", receptor_pdbqt,
           "--ligand", ligand_pdbqt,
           "--center_x", str(center[0]),
           "--center_y", str(center[1]),
           "--center_z", str(center[2]),
           "--size_x", str(size[0]),
           "--size_y", str(size[1]),
           "--size_z", str(size[2]),
           "--exhaustiveness", str(exhaustiveness),
           "--num_modes", "9",
           "--out", out_path,
           "--log", log_path]

    try:
        subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    except subprocess.TimeoutExpired:
        return [], None

    poses = []
    if os.path.exists(log_path):
        with open(log_path) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 4:
                    try:
                        mode = int(parts[0])
                        score = float(parts[1])
                        poses.append({"mode": mode, "vina_score": score})
                    except (ValueError, IndexError):
                        continue

    return poses, out_path


def _parse_gnina_log(log_path):
    """Parse GNINA log file for pose scores."""
    poses = []
    if not os.path.exists(log_path):
        return poses

    with open(log_path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 4:
                try:
                    mode = int(parts[0])
                    affinity = float(parts[1])
                    cnn_score = float(parts[2]) if len(parts) > 2 else None
                    cnn_affinity = float(parts[3]) if len(parts) > 3 else None
                    poses.append({
                        "mode": mode,
                        "vina_score": affinity,
                        "cnn_score": cnn_score,
                        "cnn_affinity": cnn_affinity,
                    })
                except (ValueError, IndexError):
                    continue
    return poses


def calculate_rmsd_coords(coords_ref, coords_docked):
    """Calculate RMSD between two coordinate sets using Hungarian assignment."""
    if coords_ref is None or coords_docked is None:
        return 999.0
    if len(coords_ref) == 0 or len(coords_docked) == 0:
        return 999.0

    if len(coords_ref) == len(coords_docked):
        try:
            from scipy.optimize import linear_sum_assignment
            from scipy.spatial.distance import cdist
            cost = cdist(coords_ref, coords_docked)
            row_idx, col_idx = linear_sum_assignment(cost)
            matched = coords_docked[col_idx]
            diff = coords_ref[row_idx] - matched
            return np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))
        except ImportError:
            diff = coords_ref - coords_docked
            return np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))

    try:
        from scipy.optimize import linear_sum_assignment
        from scipy.spatial.distance import cdist
        cost = cdist(coords_ref, coords_docked)
        row_idx, col_idx = linear_sum_assignment(cost)
        matched_ref = coords_ref[row_idx]
        matched_dock = coords_docked[col_idx]
        diff = matched_ref - matched_dock
        return np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))
    except ImportError:
        n = min(len(coords_ref), len(coords_docked))
        diff = coords_ref[:n] - coords_docked[:n]
        return np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))


def compute_pose_rmsds(docked_file, crystal_pdb, obabel_path):
    """Compute RMSD for each pose in docked output vs crystal structure."""
    crystal_coords = extract_heavy_coords_pdb(crystal_pdb)
    if crystal_coords is None:
        return []

    rmsds = []

    if docked_file.endswith(".sdf.gz") or docked_file.endswith(".sdf"):
        tmp_dir = tempfile.mkdtemp(prefix="rmsd_")
        try:
            if docked_file.endswith(".gz"):
                import gzip
                unsdf = os.path.join(tmp_dir, "docked.sdf")
                with gzip.open(docked_file, "rb") as gz:
                    with open(unsdf, "wb") as out:
                        out.write(gz.read())
            else:
                unsdf = docked_file

            if not os.path.exists(unsdf) or os.path.getsize(unsdf) < 10:
                print(f"    WARNING: Docked SDF is empty/missing, skipping RMSD")
                return []

            if HAS_RDKIT:
                try:
                    supplier = Chem.SDMolSupplier(unsdf, removeHs=True)
                except OSError as e:
                    print(f"    WARNING: Cannot read docked SDF ({e}), trying OpenBabel split")
                    supplier = None

                ref_mol = Chem.MolFromPDBFile(crystal_pdb, removeHs=True)

                if supplier is not None:
                    for i, mol in enumerate(supplier):
                        if mol is None:
                            rmsds.append(999.0)
                            continue
                        try:
                            if ref_mol is not None:
                                rmsd = rdMolAlign.CalcRMS(mol, ref_mol)
                            else:
                                coords = mol.GetConformer().GetPositions()
                                rmsd = calculate_rmsd_coords(crystal_coords, coords)
                            rmsds.append(rmsd)
                        except Exception:
                            try:
                                coords = mol.GetConformer().GetPositions()
                                rmsd = calculate_rmsd_coords(crystal_coords, coords)
                                rmsds.append(rmsd)
                            except Exception:
                                rmsds.append(999.0)
                else:
                    if obabel_path:
                        split_dir = os.path.join(tmp_dir, "split")
                        os.makedirs(split_dir, exist_ok=True)
                        cmd = [obabel_path, unsdf, "-osdf", "-O",
                               os.path.join(split_dir, "pose_.sdf"), "-m"]
                        subprocess.run(cmd, capture_output=True, text=True, timeout=30)
                        import glob
                        pose_files = sorted(glob.glob(os.path.join(split_dir, "pose_*.sdf")))
                        for pf in pose_files:
                            coords = extract_heavy_coords_sdf(pf)
                            rmsd = calculate_rmsd_coords(crystal_coords, coords)
                            rmsds.append(rmsd)
            else:
                if obabel_path:
                    split_dir = os.path.join(tmp_dir, "split")
                    os.makedirs(split_dir, exist_ok=True)
                    cmd = [obabel_path, unsdf, "-osdf", "-O",
                           os.path.join(split_dir, "pose_.sdf"), "-m"]
                    subprocess.run(cmd, capture_output=True, text=True, timeout=30)
                    import glob
                    pose_files = sorted(glob.glob(os.path.join(split_dir, "pose_*.sdf")))
                    for pf in pose_files:
                        coords = extract_heavy_coords_sdf(pf)
                        rmsd = calculate_rmsd_coords(crystal_coords, coords)
                        rmsds.append(rmsd)
        finally:
            shutil.rmtree(tmp_dir, ignore_errors=True)

    elif docked_file.endswith(".pdbqt"):
        for model in range(1, 20):
            coords = extract_heavy_coords_pdbqt(docked_file, model)
            if coords is None or len(coords) == 0:
                break
            rmsd = calculate_rmsd_coords(crystal_coords, coords)
            rmsds.append(rmsd)

    return rmsds


def extract_heavy_coords_pdb(pdb_path):
    """Extract heavy atom coordinates from PDB."""
    coords = []
    with open(pdb_path) as f:
        for line in f:
            if line[:6].strip() in ("ATOM", "HETATM"):
                element = line[76:78].strip() if len(line) > 76 else ""
                atom_name = line[12:16].strip()
                if element == "H" or atom_name.startswith("H"):
                    continue
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                except (ValueError, IndexError):
                    continue
    return np.array(coords) if coords else None


def extract_heavy_coords_sdf(sdf_path):
    """Extract heavy atom coordinates from SDF file."""
    if HAS_RDKIT:
        mol = Chem.MolFromMolFile(sdf_path, removeHs=True)
        if mol:
            return mol.GetConformer().GetPositions()
    coords = []
    with open(sdf_path) as f:
        lines = f.readlines()
    if len(lines) < 4:
        return None
    try:
        counts = lines[3].split()
        n_atoms = int(counts[0])
    except (IndexError, ValueError):
        return None
    for i in range(4, min(4 + n_atoms, len(lines))):
        parts = lines[i].split()
        if len(parts) >= 4:
            try:
                element = parts[3] if len(parts) > 3 else ""
                if element == "H":
                    continue
                coords.append([float(parts[0]), float(parts[1]), float(parts[2])])
            except (ValueError, IndexError):
                continue
    return np.array(coords) if coords else None


def extract_heavy_coords_pdbqt(pdbqt_path, model=1):
    """Extract heavy atom coords from PDBQT (specific model)."""
    coords = []
    current_model = 0
    in_target = False

    with open(pdbqt_path) as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith("MODEL"):
                current_model = int(stripped.split()[1])
                in_target = (current_model == model)
            elif stripped == "ENDMDL":
                if in_target:
                    break
            elif (in_target or current_model == 0) and stripped[:6].strip() in ("ATOM", "HETATM"):
                element = line[76:78].strip() if len(line) > 76 else ""
                atom_name = line[12:16].strip()
                if element == "H" or atom_name.startswith("H"):
                    continue
                try:
                    coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                except (ValueError, IndexError):
                    continue

    return np.array(coords) if coords else None


def run_redocking_sweep(gnina_path, receptor_pdb, ligand_sdf, crystal_ligand_pdb,
                        center, obabel_path, work_dir, exhaustiveness=20):
    """Run redocking with multiple autobox padding values and manual box sizes."""
    all_results = []

    # Phase 1: Autobox sweep (padding 2-8 A)
    print(f"\n  Phase 1: Autobox padding sweep")
    for padding in [2.0, 3.0, 4.0, 5.0, 6.0, 8.0]:
        tag = f"_autobox_{padding:.0f}"
        dock_dir = os.path.join(work_dir, f"autobox_{padding:.0f}")
        os.makedirs(dock_dir, exist_ok=True)

        poses, out_file = run_gnina_dock(
            gnina_path, receptor_pdb, ligand_sdf,
            autobox_ligand=crystal_ligand_pdb,
            center=None, size=None,
            exhaustiveness=exhaustiveness,
            output_dir=dock_dir, tag=tag,
            autobox_add=padding
        )

        rmsds = compute_pose_rmsds(out_file, crystal_ligand_pdb, obabel_path) if out_file else []

        for i, pose in enumerate(poses):
            pose["rmsd"] = rmsds[i] if i < len(rmsds) else 999.0

        best_rmsd = min(rmsds) if rmsds else 999.0
        best_pose = min(poses, key=lambda p: p.get("rmsd", 999)) if poses else {}

        result = {
            "method": "autobox",
            "padding": padding,
            "size_x": None, "size_y": None, "size_z": None,
            "poses": poses,
            "best_rmsd": best_rmsd,
            "best_cnn_score": best_pose.get("cnn_score"),
            "best_vina_score": best_pose.get("vina_score"),
            "docked_file": out_file,
        }
        all_results.append(result)

        print(f"    Padding={padding:.0f}A: Best RMSD={best_rmsd:.3f}, "
              f"CNN={best_pose.get('cnn_score', 'N/A')}, "
              f"Vina={best_pose.get('vina_score', 'N/A')}")

    # Phase 2: Manual box size sweep at best padding center
    print(f"\n  Phase 2: Manual box dimension sweep")
    for box_size in [14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0]:
        size = [box_size, box_size, box_size]
        tag = f"_manual_{box_size:.0f}"
        dock_dir = os.path.join(work_dir, f"manual_{box_size:.0f}")
        os.makedirs(dock_dir, exist_ok=True)

        poses, out_file = run_gnina_dock(
            gnina_path, receptor_pdb, ligand_sdf,
            autobox_ligand=None,
            center=center, size=size,
            exhaustiveness=exhaustiveness,
            output_dir=dock_dir, tag=tag
        )

        rmsds = compute_pose_rmsds(out_file, crystal_ligand_pdb, obabel_path) if out_file else []

        for i, pose in enumerate(poses):
            pose["rmsd"] = rmsds[i] if i < len(rmsds) else 999.0

        best_rmsd = min(rmsds) if rmsds else 999.0
        best_pose = min(poses, key=lambda p: p.get("rmsd", 999)) if poses else {}

        result = {
            "method": "manual",
            "padding": None,
            "size_x": box_size, "size_y": box_size, "size_z": box_size,
            "poses": poses,
            "best_rmsd": best_rmsd,
            "best_cnn_score": best_pose.get("cnn_score"),
            "best_vina_score": best_pose.get("vina_score"),
            "docked_file": out_file,
        }
        all_results.append(result)

        print(f"    Box={box_size:.0f}x{box_size:.0f}x{box_size:.0f}: Best RMSD={best_rmsd:.3f}, "
              f"CNN={best_pose.get('cnn_score', 'N/A')}, "
              f"Vina={best_pose.get('vina_score', 'N/A')}")

    return all_results


def generate_report(all_conditions, output_dir):
    """Generate comprehensive report."""
    report_path = os.path.join(output_dir, "redocking_report.txt")
    json_path = os.path.join(output_dir, "redocking_results.json")

    serializable = []
    for cond in all_conditions:
        sc = dict(cond)
        sc.pop("sweep_results", None)
        serializable.append(sc)

    with open(json_path, "w") as f:
        json.dump(all_conditions, f, indent=2, default=str)

    with open(report_path, "w") as f:
        f.write("=" * 70 + "\n")
        f.write("NATIVE REDOCKING OPTIMIZATION REPORT (GNINA CNN Scoring)\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("=" * 70 + "\n\n")

        for cond in all_conditions:
            f.write(f"\n{'─' * 60}\n")
            f.write(f"PDB: {cond['pdb_id']} | Ligand: {cond['ligand_name']} | pH: {cond['ph']}\n")
            f.write(f"Ligand center: ({cond['center'][0]:.2f}, {cond['center'][1]:.2f}, {cond['center'][2]:.2f})\n")
            f.write(f"{'─' * 60}\n\n")

            f.write("AUTOBOX PADDING SWEEP:\n")
            f.write(f"  {'Padding':>8} {'Best_RMSD':>10} {'CNN_Score':>10} {'Vina_Score':>11}\n")
            for r in cond.get("sweep_results", []):
                if r["method"] == "autobox":
                    cnn = f"{r['best_cnn_score']:.4f}" if r['best_cnn_score'] is not None else "N/A"
                    vina = f"{r['best_vina_score']:.2f}" if r['best_vina_score'] is not None else "N/A"
                    f.write(f"  {r['padding']:>7.0f}A {r['best_rmsd']:>10.3f} {cnn:>10} {vina:>11}\n")

            f.write("\nMANUAL BOX SIZE SWEEP:\n")
            f.write(f"  {'Box_Size':>8} {'Best_RMSD':>10} {'CNN_Score':>10} {'Vina_Score':>11}\n")
            for r in cond.get("sweep_results", []):
                if r["method"] == "manual":
                    cnn = f"{r['best_cnn_score']:.4f}" if r['best_cnn_score'] is not None else "N/A"
                    vina = f"{r['best_vina_score']:.2f}" if r['best_vina_score'] is not None else "N/A"
                    f.write(f"  {r['size_x']:>7.0f}A {r['best_rmsd']:>10.3f} {cnn:>10} {vina:>11}\n")

            best_overall = min(cond.get("sweep_results", [{}]),
                              key=lambda r: r.get("best_rmsd", 999))
            f.write(f"\nBEST RESULT:\n")
            if best_overall.get("method") == "autobox":
                f.write(f"  Method: autobox with {best_overall['padding']:.0f}A padding\n")
            else:
                f.write(f"  Method: manual box {best_overall.get('size_x', '?')}A\n")
            f.write(f"  RMSD: {best_overall.get('best_rmsd', 'N/A')}\n")
            f.write(f"  CNN Score: {best_overall.get('best_cnn_score', 'N/A')}\n")
            f.write(f"  Vina Score: {best_overall.get('best_vina_score', 'N/A')}\n")

            f.write(f"\nFINAL DOCKING POSES (best config):\n")
            f.write(f"  {'Mode':>4} {'Vina':>8} {'CNN_Score':>10} {'CNN_Aff':>8} {'RMSD':>8}\n")
            for p in best_overall.get("poses", []):
                cnn = f"{p['cnn_score']:.4f}" if p.get('cnn_score') is not None else "N/A"
                cnn_aff = f"{p['cnn_affinity']:.2f}" if p.get('cnn_affinity') is not None else "N/A"
                rmsd = f"{p['rmsd']:.3f}" if isinstance(p.get('rmsd'), float) else "N/A"
                vina = f"{p['vina_score']:.2f}" if p.get('vina_score') is not None else "N/A"
                f.write(f"  {p.get('mode', '?'):>4} {vina:>8} {cnn:>10} {cnn_aff:>8} {rmsd:>8}\n")
            f.write("\n")

        f.write("\n" + "=" * 70 + "\n")
        f.write("SUMMARY - RECOMMENDED PARAMETERS\n")
        f.write("=" * 70 + "\n\n")
        f.write("Active site residues for visualization:\n")
        f.write(f"  {', '.join(str(r) for r in ACTIVE_SITE_RESIDUES)}\n\n")

        for cond in all_conditions:
            best = min(cond.get("sweep_results", [{}]),
                      key=lambda r: r.get("best_rmsd", 999))
            f.write(f"  {cond['pdb_id']} {cond['ligand_name']} pH{cond['ph']}: ")
            if best.get("method") == "autobox":
                f.write(f"autobox padding={best['padding']:.0f}A")
            else:
                f.write(f"box={best.get('size_x', '?')}A")
            f.write(f" (RMSD={best.get('best_rmsd', 'N/A'):.3f})\n")

    print(f"\nReport saved to: {report_path}")
    print(f"JSON data saved to: {json_path}")
    return report_path


def main():
    parser = argparse.ArgumentParser(description="GNINA Native Redocking Optimizer")
    parser.add_argument("--output-dir",
                        default=os.path.expanduser("~/final_project/redocking_gnina"),
                        help="Output directory")
    parser.add_argument("--exhaustiveness", type=int, default=EXHAUSTIVENESS)
    parser.add_argument("--pdb", choices=["2ZSH", "2ZSI", "both"], default="both")
    parser.add_argument("--ph", type=float, nargs="+", default=PH_VALUES)
    parser.add_argument("--ligand-source", choices=["smiles", "crystal", "both"],
                        default="both",
                        help="Use SMILES (Gypsum-DL), crystal structure, or both")
    args = parser.parse_args()

    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    project_dir = os.path.expanduser("~/final_project")
    gnina_path = find_gnina(project_dir)
    obabel_path = find_obabel()
    vina_path = find_vina()

    use_gnina = gnina_path is not None
    docking_engine = gnina_path if use_gnina else vina_path

    if not docking_engine:
        print("ERROR: Neither GNINA nor Vina found")
        sys.exit(1)
    if not obabel_path:
        print("ERROR: OpenBabel not found")
        sys.exit(1)

    print("=" * 60)
    print("Native Redocking Optimizer (GNINA CNN Scoring)")
    print("=" * 60)
    print(f"Docking engine: {'GNINA' if use_gnina else 'Vina (fallback)'} ({docking_engine})")
    print(f"OpenBabel: {obabel_path}")
    print(f"Output: {output_dir}")
    print(f"Exhaustiveness: {args.exhaustiveness}")
    print(f"pH values: {args.ph}")
    print(f"Ligand source: {args.ligand_source}")
    print(f"RDKit available: {HAS_RDKIT}")
    print()

    if not use_gnina:
        print("WARNING: GNINA not available, using Vina. CNN scoring will not be available.")
        print("To install GNINA, run:")
        print(f"  wget {GNINA_URL} -O ~/final_project/gnina && chmod +x ~/final_project/gnina")
        print()

    pdbs_to_run = ["2ZSH", "2ZSI"] if args.pdb == "both" else [args.pdb]
    all_conditions = []

    for pdb_id in pdbs_to_run:
        pdb_info = PDB_IDS[pdb_id]
        ligand_name = pdb_info["ligand_name"]
        ligand_resname = pdb_info["ligand_resname"]
        smiles = GA_SMILES[ligand_name]

        print(f"\n{'=' * 60}")
        print(f"Processing: {pdb_id} with {ligand_name}")
        print(f"{'=' * 60}")

        pdb_dir = os.path.join(output_dir, pdb_id)
        os.makedirs(pdb_dir, exist_ok=True)

        pdb_path = download_pdb(pdb_id, pdb_dir)
        ligand_pdb, receptor_pdb, center, crystal_coords = extract_ligand_and_receptor(
            pdb_path, ligand_resname, pdb_dir
        )

        for ph in args.ph:
            print(f"\n  --- pH {ph} ---")
            ph_dir = os.path.join(pdb_dir, f"pH_{ph}")
            os.makedirs(ph_dir, exist_ok=True)

            print(f"  Protonating receptor at pH {ph} (PDB2PQR+PROPKA)...")
            receptor_protonated = protonate_receptor(receptor_pdb, ph, ph_dir)

            ligand_sources = []

            if args.ligand_source in ("smiles", "both"):
                print(f"  Preparing ligand from SMILES via Gypsum-DL...")
                smi_sdf = prepare_ligand_gypsymdl(smiles, ph, ph_dir, name=ligand_name)
                if smi_sdf:
                    ligand_sources.append(("smiles", smi_sdf))

            if args.ligand_source in ("crystal", "both"):
                print(f"  Preparing ligand from crystal structure...")
                cryst_sdf = prepare_ligand_from_crystal(ligand_pdb, ph, ph_dir, name=ligand_name)
                if cryst_sdf:
                    ligand_sources.append(("crystal", cryst_sdf))

            if not ligand_sources:
                print(f"  ERROR: No ligand source available, skipping")
                continue

            for source_name, ligand_sdf in ligand_sources:
                print(f"\n  Running redocking sweep ({source_name} ligand)...")
                work_dir = os.path.join(ph_dir, f"sweep_{source_name}")
                os.makedirs(work_dir, exist_ok=True)

                if use_gnina:
                    sweep_results = run_redocking_sweep(
                        gnina_path, receptor_protonated, ligand_sdf,
                        ligand_pdb, center.tolist(), obabel_path, work_dir,
                        exhaustiveness=args.exhaustiveness
                    )
                else:
                    print("  Vina fallback: running single box size sweep")
                    sweep_results = []

                best_result = min(sweep_results, key=lambda r: r.get("best_rmsd", 999)) if sweep_results else {}

                condition = {
                    "pdb_id": pdb_id,
                    "ligand_name": ligand_name,
                    "ph": ph,
                    "ligand_source": source_name,
                    "center": center.tolist(),
                    "sweep_results": sweep_results,
                    "best_config": best_result,
                    "receptor_pdb": receptor_protonated,
                    "crystal_ligand_pdb": ligand_pdb,
                    "ligand_sdf": ligand_sdf,
                }
                all_conditions.append(condition)

    generate_report(all_conditions, output_dir)

    print("\n" + "=" * 60)
    print("ALL DONE!")
    print("=" * 60)
    for cond in all_conditions:
        best = cond.get("best_config", {})
        rmsd = best.get("best_rmsd", "N/A")
        if isinstance(rmsd, float):
            rmsd = f"{rmsd:.3f}"
        print(f"  {cond['pdb_id']} {cond['ligand_name']} pH{cond['ph']} ({cond['ligand_source']}): "
              f"Best RMSD={rmsd}")


if __name__ == "__main__":
    main()
