#!/usr/bin/env python3
"""
Native Redocking Optimizer for AutoDock Vina
Finds optimal box dimensions by minimizing RMSD between docked poses
and crystallographic ligand coordinates.

Targets:
  - 2ZSH: GID1A receptor with GA3 (Gibberellic acid 3)
  - 2ZSI: GID1A receptor with GA4 (Gibberellic acid 4)

Runs at pH 7.2 and pH 7.4 protonation states.
Uses coordinate descent with golden-section search per axis.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.request
from datetime import datetime

import numpy as np

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem, rdMolAlign, Descriptors
    RDLogger.DisableLog("rdApp.*")
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    print("WARNING: RDKit not available. RMSD calculations will use coordinate-based method.")

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
BOX_MIN = 12.0
BOX_MAX = 30.0
BOX_STEP = 0.5
GOLDEN_RATIO = (np.sqrt(5) - 1) / 2
CONVERGENCE_TOL = 0.5
MAX_CYCLES = 3


def download_pdb(pdb_id, output_dir):
    """Download PDB file from RCSB."""
    pdb_path = os.path.join(output_dir, f"{pdb_id}.pdb")
    if os.path.exists(pdb_path):
        print(f"  {pdb_id}.pdb already exists, skipping download")
        return pdb_path
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"  Downloading {pdb_id} from RCSB...")
    urllib.request.urlretrieve(url, pdb_path)
    return pdb_path


def extract_ligand_and_receptor(pdb_path, ligand_resname, output_dir):
    """Extract ligand coordinates and receptor from PDB file."""
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
    print(f"  Ligand heavy atom coords: {len(ligand_coords)}")

    return ligand_pdb, receptor_pdb, center, ligand_coords


def protonate_receptor_pdb2pqr(receptor_pdb, ph, output_dir):
    """Protonate receptor using PDB2PQR + PROPKA."""
    base = os.path.basename(receptor_pdb).replace(".pdb", "")
    pqr_path = os.path.join(output_dir, f"{base}_pH{ph}.pqr")
    pdb_out = os.path.join(output_dir, f"{base}_pH{ph}.pdb")

    try:
        cmd = [
            sys.executable, "-m", "pdb2pqr",
            "--ff", "PARSE",
            "--titration-state-method", "propka",
            "--with-ph", str(ph),
            "--keep-chain",
            receptor_pdb, pqr_path
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if os.path.exists(pqr_path):
            print(f"  PDB2PQR protonation at pH {ph}: OK")
            _pqr_to_pdb(pqr_path, pdb_out)
            return pdb_out
        else:
            print(f"  PDB2PQR failed: {result.stderr[:200]}")
    except Exception as e:
        print(f"  PDB2PQR not available ({e}), falling back to OpenBabel")

    return _protonate_obabel(receptor_pdb, ph, pdb_out)


def _pqr_to_pdb(pqr_path, pdb_out):
    """Convert PQR back to PDB format."""
    lines = []
    with open(pqr_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                parts = line.split()
                pdb_line = line[:54]
                pdb_line += "  1.00  0.00"
                pdb_line += " " * max(0, 76 - len(pdb_line))
                element = parts[-1] if len(parts) > 0 else ""
                if len(element) <= 2:
                    pdb_line += f"  {element:>2}"
                lines.append(pdb_line.rstrip() + "\n")
            elif line.startswith(("TER", "END")):
                lines.append(line)
    with open(pdb_out, "w") as f:
        f.writelines(lines)
        if not any(l.startswith("END") for l in lines):
            f.write("END\n")


def _protonate_obabel(input_pdb, ph, output_pdb):
    """Fallback: protonate with OpenBabel."""
    obabel = _find_obabel()
    if not obabel:
        print("  WARNING: OpenBabel not found, using unprotonated structure")
        shutil.copy2(input_pdb, output_pdb)
        return output_pdb

    cmd = [obabel, input_pdb, "-O", output_pdb, "-p", str(ph)]
    subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    if os.path.exists(output_pdb) and os.path.getsize(output_pdb) > 100:
        print(f"  OpenBabel protonation at pH {ph}: OK")
        return output_pdb

    shutil.copy2(input_pdb, output_pdb)
    return output_pdb


def prepare_ligand_pdbqt(crystal_ligand_pdb, smiles, ph, obabel_path, output_dir, name="ligand"):
    """Prepare ligand PDBQT for docking.

    Primary: use the crystallographic ligand PDB (already has correct 3D coords),
    add hydrogens at target pH, convert to PDBQT.
    Fallback: generate 3D from SMILES if crystal approach fails.
    """
    pdbqt_path = os.path.join(output_dir, f"{name}_pH{ph}.pdbqt")
    protonated_pdb = os.path.join(output_dir, f"{name}_pH{ph}_H.pdb")

    # Primary approach: crystal ligand + OpenBabel protonation
    if crystal_ligand_pdb and os.path.exists(crystal_ligand_pdb):
        cmd = [obabel_path, crystal_ligand_pdb, "-O", protonated_pdb, "-h", "-p", str(ph)]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)

        if os.path.exists(protonated_pdb) and os.path.getsize(protonated_pdb) > 10:
            cmd = [obabel_path, protonated_pdb, "-opdbqt", "-O", pdbqt_path]
            subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            if os.path.exists(pdbqt_path) and os.path.getsize(pdbqt_path) > 10:
                print(f"  Crystal ligand -> PDBQT at pH {ph}: OK")
                return pdbqt_path

        cmd = [obabel_path, crystal_ligand_pdb, "-opdbqt", "-O", pdbqt_path, "-h"]
        subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        if os.path.exists(pdbqt_path) and os.path.getsize(pdbqt_path) > 10:
            print(f"  Crystal ligand -> PDBQT (no pH adjust): OK")
            return pdbqt_path

    # Fallback: generate 3D from SMILES
    print(f"  Crystal approach failed, trying SMILES -> 3D...")
    sdf_path = os.path.join(output_dir, f"{name}_pH{ph}.sdf")

    for gen3d_args in [["--gen3d", "-h", "-p", str(ph)],
                       ["--gen3d", "-h"],
                       ["--gen3d", "--best", "-h"]]:
        cmd = [obabel_path, "-:" + smiles, "-osdf", "-O", sdf_path] + gen3d_args
        subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if os.path.exists(sdf_path) and os.path.getsize(sdf_path) > 10:
            cmd = [obabel_path, sdf_path, "-opdbqt", "-O", pdbqt_path, "-h"]
            subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            if os.path.exists(pdbqt_path) and os.path.getsize(pdbqt_path) > 10:
                print(f"  SMILES -> 3D -> PDBQT at pH {ph}: OK")
                return pdbqt_path

    print(f"  ERROR: All ligand preparation methods failed for {name}")
    return None


def prepare_receptor_pdbqt(receptor_pdb, output_dir):
    """Convert receptor PDB to PDBQT using MGLTools or OpenBabel."""
    pdbqt_path = receptor_pdb.replace(".pdb", ".pdbqt")

    mgltools_prep = os.path.expanduser(
        "~/final_project/autogrow4/autogrow/docking/docking_executables/"
        "mgl_tools/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"
    )
    mgltools_python = os.path.expanduser(
        "~/final_project/autogrow4/autogrow/docking/docking_executables/"
        "mgl_tools/bin/pythonsh"
    )

    if os.path.exists(mgltools_prep) and os.path.exists(mgltools_python):
        cmd = [mgltools_python, mgltools_prep,
               "-r", receptor_pdb, "-o", pdbqt_path,
               "-A", "hydrogens"]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if os.path.exists(pdbqt_path) and os.path.getsize(pdbqt_path) > 100:
            print(f"  MGLTools receptor prep: OK")
            return pdbqt_path
        print(f"  MGLTools failed: {result.stderr[:200]}")

    obabel = _find_obabel()
    if obabel:
        cmd = [obabel, receptor_pdb, "-opdbqt", "-O", pdbqt_path, "-xr"]
        subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if os.path.exists(pdbqt_path) and os.path.getsize(pdbqt_path) > 100:
            print(f"  OpenBabel receptor prep: OK")
            return pdbqt_path

    return None


def _find_obabel():
    """Find OpenBabel executable."""
    ob = shutil.which("obabel")
    if ob:
        return ob
    for p in ["/usr/bin/obabel", "/usr/local/bin/obabel",
              os.path.expanduser("~/miniconda3/bin/obabel"),
              os.path.expanduser("~/.conda/envs/autogrow/bin/obabel"),
              os.path.expanduser("~/.conda/envs/ag4_full/bin/obabel")]:
        if os.path.isfile(p):
            return p
    return None


def _find_vina():
    """Find Vina executable."""
    home = os.path.expanduser("~")
    direct = os.path.join(home, "final_project", "autogrow4", "autogrow", "docking",
                          "docking_executables", "vina", "autodock_vina_1_1_2_linux_x86", "bin", "vina")
    if os.path.isfile(direct):
        return direct
    v = shutil.which("vina")
    if v:
        return v
    return None


def run_vina_dock(vina_path, receptor_pdbqt, ligand_pdbqt, center, size, exhaustiveness, tmp_dir):
    """Run Vina docking and return (score, docked_pdbqt_path)."""
    out_path = os.path.join(tmp_dir, "docked.pdbqt")
    log_path = os.path.join(tmp_dir, "vina.log")

    cmd = [
        vina_path,
        "--receptor", receptor_pdbqt,
        "--ligand", ligand_pdbqt,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(size[0]),
        "--size_y", str(size[1]),
        "--size_z", str(size[2]),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", "5",
        "--out", out_path,
        "--log", log_path,
    ]

    try:
        subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    except subprocess.TimeoutExpired:
        return None, None

    score = None
    if os.path.exists(log_path):
        with open(log_path) as f:
            for line in f:
                line = line.strip()
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        mode = int(parts[0])
                        if mode == 1:
                            score = float(parts[1])
                            break
                    except (ValueError, IndexError):
                        continue

    docked = out_path if os.path.exists(out_path) else None
    return score, docked


def extract_pose_coords(pdbqt_path, model=1):
    """Extract heavy atom coordinates from a PDBQT file (specific model)."""
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
                in_target = False
            elif (in_target or current_model == 0) and stripped[:6].strip() in ("ATOM", "HETATM"):
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

    if not coords and current_model == 0:
        with open(pdbqt_path) as f:
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


def calculate_rmsd(coords_ref, coords_docked):
    """Calculate RMSD between two coordinate sets.

    If atom counts differ, uses the minimum overlapping set with
    Hungarian assignment based on distance matrix.
    """
    if coords_ref is None or coords_docked is None:
        return 999.0
    if len(coords_ref) == 0 or len(coords_docked) == 0:
        return 999.0

    if len(coords_ref) == len(coords_docked):
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


def calculate_rmsd_rdkit(crystal_pdb, docked_pdbqt, obabel_path, tmp_dir):
    """Calculate RMSD using RDKit with symmetry handling."""
    if not HAS_RDKIT:
        crystal_coords = _extract_pdb_heavy_coords(crystal_pdb)
        docked_coords = extract_pose_coords(docked_pdbqt, model=1)
        return calculate_rmsd(crystal_coords, docked_coords)

    crystal_sdf = os.path.join(tmp_dir, "crystal.sdf")
    docked_sdf = os.path.join(tmp_dir, "docked.sdf")

    subprocess.run([obabel_path, crystal_pdb, "-osdf", "-O", crystal_sdf],
                   capture_output=True, text=True, timeout=30)
    docked_pdb_tmp = os.path.join(tmp_dir, "docked_pose1.pdb")
    _extract_model1_as_pdb(docked_pdbqt, docked_pdb_tmp)
    subprocess.run([obabel_path, docked_pdb_tmp, "-osdf", "-O", docked_sdf],
                   capture_output=True, text=True, timeout=30)

    try:
        ref_mol = Chem.MolFromMolFile(crystal_sdf, removeHs=True)
        probe_mol = Chem.MolFromMolFile(docked_sdf, removeHs=True)
        if ref_mol and probe_mol:
            rmsd = rdMolAlign.CalcRMS(probe_mol, ref_mol)
            return rmsd
    except Exception:
        pass

    crystal_coords = _extract_pdb_heavy_coords(crystal_pdb)
    docked_coords = extract_pose_coords(docked_pdbqt, model=1)
    return calculate_rmsd(crystal_coords, docked_coords)


def _extract_pdb_heavy_coords(pdb_path):
    """Extract heavy atom coordinates from PDB file."""
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


def _extract_model1_as_pdb(pdbqt_path, output_pdb):
    """Extract MODEL 1 from PDBQT and write as PDB."""
    lines = []
    in_model1 = False
    found_model = False

    with open(pdbqt_path) as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith("MODEL"):
                model_num = int(stripped.split()[1])
                if model_num == 1:
                    in_model1 = True
                    found_model = True
                elif found_model:
                    break
            elif stripped == "ENDMDL":
                if in_model1:
                    break
            elif in_model1 or not found_model:
                record = line[:6].strip()
                if record in ("ATOM", "HETATM"):
                    lines.append(line[:66].rstrip() + "\n")
                elif record not in ("ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF"):
                    lines.append(line)

    with open(output_pdb, "w") as f:
        f.writelines(lines)
        f.write("END\n")


def golden_section_search(func, a, b, tol=0.5):
    """Golden section search to find minimum of func in [a, b]."""
    evaluations = []
    c = b - GOLDEN_RATIO * (b - a)
    d = a + GOLDEN_RATIO * (b - a)

    c = round(c / 0.5) * 0.5
    d = round(d / 0.5) * 0.5

    fc = func(c)
    evaluations.append((c, fc))
    fd = func(d)
    evaluations.append((d, fd))

    while abs(b - a) > tol:
        if fc < fd:
            b = d
            d = c
            fd = fc
            c = b - GOLDEN_RATIO * (b - a)
            c = round(c / 0.5) * 0.5
            c = max(BOX_MIN, min(BOX_MAX, c))
            fc = func(c)
            evaluations.append((c, fc))
        else:
            a = c
            c = d
            fc = fd
            d = a + GOLDEN_RATIO * (b - a)
            d = round(d / 0.5) * 0.5
            d = max(BOX_MIN, min(BOX_MAX, d))
            fd = func(d)
            evaluations.append((d, fd))

    best_val, best_rmsd = min(evaluations, key=lambda x: x[1])
    return best_val, best_rmsd, evaluations


def optimize_box_dimensions(vina_path, receptor_pdbqt, ligand_pdbqt,
                            crystal_pdb, center, obabel_path, work_dir,
                            exhaustiveness=20):
    """Optimize box dimensions using coordinate descent + golden section."""
    print(f"\n  Starting box optimization...")
    print(f"  Center: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")
    print(f"  Search range: {BOX_MIN} - {BOX_MAX} A per axis")
    print(f"  Exhaustiveness: {exhaustiveness}")

    current_size = [21.0, 21.0, 21.0]
    all_results = []
    total_dockings = 0

    for cycle in range(MAX_CYCLES):
        print(f"\n  Cycle {cycle + 1}/{MAX_CYCLES}")
        prev_size = current_size.copy()

        for axis in range(3):
            axis_name = ["X", "Y", "Z"][axis]
            print(f"    Optimizing {axis_name} axis (current: {current_size})...")

            def eval_axis(val):
                nonlocal total_dockings
                size = current_size.copy()
                size[axis] = val
                tmp = tempfile.mkdtemp(prefix=f"redock_c{cycle}_{axis_name}_", dir=work_dir)
                try:
                    score, docked = run_vina_dock(
                        vina_path, receptor_pdbqt, ligand_pdbqt,
                        center, size, exhaustiveness, tmp
                    )
                    total_dockings += 1

                    if docked:
                        rmsd = calculate_rmsd_rdkit(crystal_pdb, docked, obabel_path, tmp)
                    else:
                        rmsd = 999.0

                    result = {
                        "cycle": cycle + 1, "axis": axis_name,
                        "size_x": size[0], "size_y": size[1], "size_z": size[2],
                        "vina_score": score, "rmsd": rmsd,
                        "docking_num": total_dockings
                    }
                    all_results.append(result)
                    print(f"      {axis_name}={val:.1f} -> RMSD={rmsd:.3f}, Score={score}")
                    return rmsd
                finally:
                    shutil.rmtree(tmp, ignore_errors=True)

            best_val, best_rmsd, evals = golden_section_search(
                eval_axis, BOX_MIN, BOX_MAX, tol=CONVERGENCE_TOL
            )
            current_size[axis] = best_val
            print(f"    Best {axis_name}: {best_val:.1f} A (RMSD={best_rmsd:.3f})")

        size_change = sum(abs(current_size[i] - prev_size[i]) for i in range(3))
        print(f"  Cycle {cycle + 1} result: size={current_size}, total change={size_change:.1f}")
        if size_change < CONVERGENCE_TOL:
            print(f"  Converged after {cycle + 1} cycles")
            break

    print(f"\n  Optimization complete: {total_dockings} total dockings")
    print(f"  Best box: {current_size[0]:.1f} x {current_size[1]:.1f} x {current_size[2]:.1f}")

    best_result = min(all_results, key=lambda r: r["rmsd"])
    return current_size, all_results, best_result


def run_final_docking(vina_path, receptor_pdbqt, ligand_pdbqt, crystal_pdb,
                      center, size, obabel_path, work_dir, exhaustiveness=20):
    """Run final docking with optimal parameters and save all poses."""
    print(f"\n  Final docking with optimal box: {size}")
    out_dir = os.path.join(work_dir, "final_poses")
    os.makedirs(out_dir, exist_ok=True)

    out_path = os.path.join(out_dir, "final_docked.pdbqt")
    log_path = os.path.join(out_dir, "final_vina.log")

    cmd = [
        vina_path,
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
        "--log", log_path,
    ]

    subprocess.run(cmd, capture_output=True, text=True, timeout=600)

    poses = []
    if os.path.exists(log_path):
        with open(log_path) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 4:
                    try:
                        mode = int(parts[0])
                        score = float(parts[1])
                        poses.append({"mode": mode, "score": score})
                    except (ValueError, IndexError):
                        continue

    rmsds = []
    if os.path.exists(out_path):
        for i, pose in enumerate(poses):
            tmp = tempfile.mkdtemp(prefix="rmsd_", dir=work_dir)
            try:
                docked_coords = extract_pose_coords(out_path, model=i + 1)
                crystal_coords = _extract_pdb_heavy_coords(crystal_pdb)
                rmsd = calculate_rmsd(crystal_coords, docked_coords)
                rmsds.append(rmsd)
                pose["rmsd"] = rmsd
            finally:
                shutil.rmtree(tmp, ignore_errors=True)

    print(f"  Final docking poses:")
    for p in poses:
        rmsd_str = f"{p.get('rmsd', 999):.3f}" if "rmsd" in p else "N/A"
        print(f"    Mode {p['mode']}: Score={p['score']:.2f}, RMSD={rmsd_str}")

    return poses, out_path


def generate_report(all_conditions, output_dir):
    """Generate comprehensive report."""
    report_path = os.path.join(output_dir, "redocking_report.txt")
    json_path = os.path.join(output_dir, "redocking_results.json")

    with open(json_path, "w") as f:
        json.dump(all_conditions, f, indent=2, default=str)

    with open(report_path, "w") as f:
        f.write("=" * 70 + "\n")
        f.write("NATIVE REDOCKING OPTIMIZATION REPORT\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("=" * 70 + "\n\n")

        for cond in all_conditions:
            f.write(f"\n{'─' * 60}\n")
            f.write(f"PDB: {cond['pdb_id']} | Ligand: {cond['ligand_name']} | pH: {cond['ph']}\n")
            f.write(f"{'─' * 60}\n")
            f.write(f"Ligand Center: ({cond['center'][0]:.2f}, {cond['center'][1]:.2f}, {cond['center'][2]:.2f})\n")
            f.write(f"Total Dockings: {cond['total_dockings']}\n\n")

            f.write(f"OPTIMAL BOX DIMENSIONS:\n")
            best = cond["best_box"]
            f.write(f"  X: {best[0]:.1f} A\n")
            f.write(f"  Y: {best[1]:.1f} A\n")
            f.write(f"  Z: {best[2]:.1f} A\n\n")

            f.write(f"BEST RESULT:\n")
            br = cond["best_result"]
            f.write(f"  RMSD: {br['rmsd']:.3f} A\n")
            f.write(f"  Vina Score: {br['vina_score']}\n\n")

            f.write(f"FINAL DOCKING POSES:\n")
            for p in cond.get("final_poses", []):
                rmsd_str = f"{p.get('rmsd', 'N/A')}"
                if isinstance(p.get('rmsd'), float):
                    rmsd_str = f"{p['rmsd']:.3f}"
                f.write(f"  Mode {p['mode']}: Score={p['score']:.2f} kcal/mol, RMSD={rmsd_str} A\n")

            f.write(f"\nOPTIMIZATION TRAJECTORY:\n")
            f.write(f"  {'#':>4} {'Cycle':>5} {'Axis':>4} {'Size_X':>7} {'Size_Y':>7} {'Size_Z':>7} {'RMSD':>8} {'Score':>8}\n")
            for i, r in enumerate(cond.get("trajectory", [])):
                f.write(f"  {i+1:>4} {r['cycle']:>5} {r['axis']:>4} "
                        f"{r['size_x']:>7.1f} {r['size_y']:>7.1f} {r['size_z']:>7.1f} "
                        f"{r['rmsd']:>8.3f} {str(r['vina_score']):>8}\n")
            f.write("\n")

        f.write("\n" + "=" * 70 + "\n")
        f.write("SUMMARY - RECOMMENDED BOX DIMENSIONS\n")
        f.write("=" * 70 + "\n\n")

        all_best = [c["best_box"] for c in all_conditions]
        avg_x = np.mean([b[0] for b in all_best])
        avg_y = np.mean([b[1] for b in all_best])
        avg_z = np.mean([b[2] for b in all_best])

        for cond in all_conditions:
            b = cond["best_box"]
            br = cond["best_result"]
            f.write(f"  {cond['pdb_id']} {cond['ligand_name']} pH{cond['ph']}: "
                    f"{b[0]:.1f} x {b[1]:.1f} x {b[2]:.1f} (RMSD={br['rmsd']:.3f})\n")

        f.write(f"\n  Average across conditions: {avg_x:.1f} x {avg_y:.1f} x {avg_z:.1f}\n")
        f.write(f"\n  Recommended for AutoGrow4 runs:\n")
        rec_x = round(avg_x / 0.5) * 0.5
        rec_y = round(avg_y / 0.5) * 0.5
        rec_z = round(avg_z / 0.5) * 0.5
        f.write(f"    --size_x {rec_x:.1f} --size_y {rec_y:.1f} --size_z {rec_z:.1f}\n")

    print(f"\nReport saved to: {report_path}")
    print(f"JSON data saved to: {json_path}")
    return report_path


def main():
    parser = argparse.ArgumentParser(description="Native Redocking Box Optimizer")
    parser.add_argument("--output-dir", default=os.path.expanduser("~/final_project/redocking_optimization"),
                        help="Output directory")
    parser.add_argument("--exhaustiveness", type=int, default=EXHAUSTIVENESS)
    parser.add_argument("--pdb", choices=["2ZSH", "2ZSI", "both"], default="both",
                        help="Which PDB to process")
    parser.add_argument("--ph", type=float, nargs="+", default=PH_VALUES,
                        help="pH values to test")
    args = parser.parse_args()

    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    vina_path = _find_vina()
    obabel_path = _find_obabel()

    if not vina_path:
        print("ERROR: Could not find Vina executable")
        sys.exit(1)
    if not obabel_path:
        print("ERROR: Could not find OpenBabel (obabel)")
        sys.exit(1)

    print("=" * 60)
    print("Native Redocking Box Optimizer")
    print("=" * 60)
    print(f"Vina: {vina_path}")
    print(f"OpenBabel: {obabel_path}")
    print(f"Output: {output_dir}")
    print(f"Exhaustiveness: {args.exhaustiveness}")
    print(f"pH values: {args.ph}")
    print(f"RDKit available: {HAS_RDKIT}")
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

            print(f"  Protonating receptor at pH {ph}...")
            receptor_protonated = protonate_receptor_pdb2pqr(receptor_pdb, ph, ph_dir)
            receptor_pdbqt = prepare_receptor_pdbqt(receptor_protonated, ph_dir)
            if not receptor_pdbqt:
                print(f"  ERROR: Failed to prepare receptor PDBQT")
                continue

            print(f"  Preparing ligand {ligand_name} at pH {ph}...")
            ligand_pdbqt = prepare_ligand_pdbqt(
                ligand_pdb, smiles, ph, obabel_path, ph_dir, name=ligand_name
            )
            if not ligand_pdbqt:
                print(f"  ERROR: Failed to prepare ligand PDBQT")
                continue

            work_dir = os.path.join(ph_dir, "optimization")
            os.makedirs(work_dir, exist_ok=True)

            best_box, trajectory, best_result = optimize_box_dimensions(
                vina_path, receptor_pdbqt, ligand_pdbqt,
                ligand_pdb, center.tolist(), obabel_path, work_dir,
                exhaustiveness=args.exhaustiveness
            )

            final_poses, final_pdbqt = run_final_docking(
                vina_path, receptor_pdbqt, ligand_pdbqt, ligand_pdb,
                center.tolist(), best_box, obabel_path, work_dir,
                exhaustiveness=args.exhaustiveness
            )

            shutil.copy2(final_pdbqt, os.path.join(ph_dir, f"final_docked_{ligand_name}_pH{ph}.pdbqt"))
            shutil.copy2(receptor_pdbqt, os.path.join(ph_dir, f"receptor_{pdb_id}_pH{ph}.pdbqt"))
            shutil.copy2(ligand_pdb, os.path.join(ph_dir, f"crystal_{ligand_name}.pdb"))

            condition = {
                "pdb_id": pdb_id,
                "ligand_name": ligand_name,
                "ph": ph,
                "center": center.tolist(),
                "best_box": best_box,
                "best_result": best_result,
                "trajectory": trajectory,
                "final_poses": final_poses,
                "total_dockings": len(trajectory),
                "receptor_pdbqt": receptor_pdbqt,
                "crystal_ligand_pdb": ligand_pdb,
                "final_docked_pdbqt": os.path.join(ph_dir, f"final_docked_{ligand_name}_pH{ph}.pdbqt"),
            }
            all_conditions.append(condition)

    report_path = generate_report(all_conditions, output_dir)

    print("\n" + "=" * 60)
    print("ALL DONE!")
    print("=" * 60)
    for cond in all_conditions:
        b = cond["best_box"]
        br = cond["best_result"]
        print(f"  {cond['pdb_id']} {cond['ligand_name']} pH{cond['ph']}: "
              f"{b[0]:.1f} x {b[1]:.1f} x {b[2]:.1f} (RMSD={br['rmsd']:.3f})")


if __name__ == "__main__":
    main()
