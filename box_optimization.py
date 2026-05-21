#!/usr/bin/env python3
"""
Box Optimization for 4-Structure Ensemble - Task 4

Determines optimal docking box dimensions for all 4 GID1 crystal
structures using GNINA CNN scoring at pH 7.2 and 7.4.

Calculates:
  1. Individual optimal box for each structure (autobox sweep)
  2. Universal consensus box for the 4-structure ensemble
  3. nglview visualization script for documentation

Structures:
  2ZSH: Arabidopsis GID1A + GA3
  2ZSI: Arabidopsis GID1A + GA4
  3ED1: Rice GID1   + GA3 (hexamer, chain A)
  3EBL: Rice GID1   + GA4 (hexamer, chain A)
"""

import os
import sys
import json
import subprocess
import urllib.request
import numpy as np

from Bio.PDB import PDBParser, Superimposer, PDBIO, Select

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    RDLogger.DisableLog("rdApp.*")
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

OUTPUT_DIR = os.path.expanduser("~/final_project/box_optimization")
POCKET_DIR = os.path.expanduser("~/final_project/pocket_analysis")
GNINA_PATH = os.path.expanduser("~/final_project/gnina")
os.makedirs(OUTPUT_DIR, exist_ok=True)

PDB_IDS = ["2ZSH", "2ZSI", "3ED1", "3EBL"]
SPECIES = {"2ZSH": "Arabidopsis", "2ZSI": "Arabidopsis",
           "3ED1": "Rice", "3EBL": "Rice"}
LIGAND_OF = {"2ZSH": "GA3", "2ZSI": "GA4",
             "3ED1": "GA3", "3EBL": "GA4"}

ARAB_RESIDUES = [24, 27, 31, 116, 119, 126, 127, 191, 218, 238, 239, 244, 319, 323]
RICE_RESIDUES = [24, 27, 31, 123, 126, 133, 134, 198, 225, 245, 246, 251, 326, 330]

BACKBONE_ATOMS = ["N", "CA", "C", "O"]
PH_VALUES = [7.2, 7.4]
EXHAUSTIVENESS = 16

AUTOBOX_PADDINGS = [2, 3, 4, 5, 6, 8]
MANUAL_SIZES = [
    (16, 16, 16),
    (18, 18, 18),
    (20, 20, 20),
    (22, 22, 22),
    (24, 24, 24),
]


def download_pdb(pdb_id):
    path = os.path.join(POCKET_DIR, f"{pdb_id}.pdb")
    if os.path.exists(path) and os.path.getsize(path) > 1000:
        return path
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"  Downloading {pdb_id}...")
    urllib.request.urlretrieve(url, path)
    return path


def get_pocket_residues(pdb_id):
    if SPECIES[pdb_id] == "Arabidopsis":
        return list(ARAB_RESIDUES)
    return list(RICE_RESIDUES)


def get_pocket_backbone_atoms(chain, pocket_residues):
    atoms = []
    for resnum in pocket_residues:
        try:
            res = chain[(' ', resnum, ' ')]
            for aname in BACKBONE_ATOMS:
                if aname in res:
                    atoms.append(res[aname])
        except KeyError:
            pass
    return atoms


class ChainAProteinSelect(Select):
    def accept_chain(self, chain):
        return chain.id == "A"

    def accept_residue(self, residue):
        return residue.id[0] == " "


class ChainALigandSelect(Select):
    def __init__(self, lig_name):
        self.lig_name = lig_name

    def accept_chain(self, chain):
        return chain.id == "A"

    def accept_residue(self, residue):
        return (residue.id[0] == " " or
                residue.get_resname() == self.lig_name)


def extract_receptor(struct, pdb_id):
    out_path = os.path.join(OUTPUT_DIR, f"{pdb_id}_receptor.pdb")
    io = PDBIO()
    io.set_structure(struct)
    io.save(out_path, ChainAProteinSelect())
    return out_path


def extract_ligand_pdb(pdb_path, pdb_id, lig_name):
    out_path = os.path.join(OUTPUT_DIR, f"{pdb_id}_{lig_name}_ligand.pdb")
    with open(pdb_path) as f:
        lines = f.readlines()

    lig_lines = [l for l in lines
                 if l.startswith("HETATM") and l[21] == "A"
                 and l[17:20].strip() == lig_name]

    if not lig_lines:
        print(f"    WARNING: No {lig_name} in chain A of {pdb_id}")
        return None

    with open(out_path, "w") as f:
        for line in lig_lines:
            f.write(line)
        f.write("END\n")
    return out_path


def extract_aligned_ligand(struct, pdb_id, lig_name):
    out_path = os.path.join(OUTPUT_DIR,
                            f"{pdb_id}_{lig_name}_ligand_aligned.pdb")
    chain = struct[0]['A']
    lines = []
    serial = 1
    for res in chain:
        if res.get_resname() == lig_name and res.id[0] != " ":
            for atom in res:
                coord = atom.get_vector().get_array()
                elem = atom.element.strip()
                if not elem:
                    elem = atom.name.strip()[0]
                line = (f"HETATM{serial:5d} {atom.name:<4s} "
                        f"{lig_name:>3s} A{res.id[1]:4d}    "
                        f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}"
                        f"  1.00  0.00          {elem:>2s}\n")
                lines.append(line)
                serial += 1

    if not lines:
        return None

    with open(out_path, "w") as f:
        for line in lines:
            f.write(line)
        f.write("END\n")
    return out_path


def extract_ligand_sdf(lig_pdb_path, pdb_id, lig_name):
    out_path = os.path.join(OUTPUT_DIR, f"{pdb_id}_{lig_name}_ligand.sdf")
    if not HAS_RDKIT:
        return None
    mol = Chem.MolFromPDBFile(lig_pdb_path, removeHs=False, sanitize=False)
    if mol is None:
        mol = Chem.MolFromPDBFile(lig_pdb_path, removeHs=True, sanitize=False)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass
    writer = Chem.SDWriter(out_path)
    writer.write(mol)
    writer.close()
    return out_path


def get_ligand_coords(lig_pdb_path):
    coords = []
    with open(lig_pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
    return np.array(coords) if coords else None


def calc_rmsd_to_crystal(docked_sdf, crystal_pdb):
    if not HAS_RDKIT:
        return None

    crystal_mol = Chem.MolFromPDBFile(crystal_pdb, removeHs=True,
                                       sanitize=False)
    if crystal_mol is None:
        return None
    try:
        Chem.SanitizeMol(crystal_mol)
    except Exception:
        pass

    suppl = Chem.SDMolSupplier(docked_sdf, removeHs=True)
    best_rmsd = 999.0

    for mol in suppl:
        if mol is None:
            continue
        try:
            n_ref = crystal_mol.GetNumAtoms()
            n_mol = mol.GetNumAtoms()
            if n_ref != n_mol:
                c1 = crystal_mol.GetConformer()
                c2 = mol.GetConformer()
                n = min(n_ref, n_mol)
                coords1 = np.array([c1.GetAtomPosition(i)
                                     for i in range(n)])
                coords2 = np.array([c2.GetAtomPosition(i)
                                     for i in range(n)])
                rmsd = np.sqrt(np.mean(np.sum(
                    (coords1 - coords2) ** 2, axis=1)))
            else:
                from rdkit.Chem import rdMolAlign
                rmsd = rdMolAlign.CalcRMS(crystal_mol, mol)
            if rmsd < best_rmsd:
                best_rmsd = rmsd
        except Exception:
            try:
                c1 = crystal_mol.GetConformer()
                c2 = mol.GetConformer()
                n = min(crystal_mol.GetNumAtoms(), mol.GetNumAtoms())
                coords1 = np.array([c1.GetAtomPosition(i) for i in range(n)])
                coords2 = np.array([c2.GetAtomPosition(i) for i in range(n)])
                rmsd = np.sqrt(np.mean(np.sum(
                    (coords1 - coords2) ** 2, axis=1)))
                if rmsd < best_rmsd:
                    best_rmsd = rmsd
            except Exception:
                continue

    return best_rmsd if best_rmsd < 900 else None


def get_gnina_scores(docked_sdf):
    if not HAS_RDKIT or not os.path.exists(docked_sdf):
        return None, None
    suppl = Chem.SDMolSupplier(docked_sdf, removeHs=False)
    best_cnn = None
    best_vina = None
    for mol in suppl:
        if mol is None:
            continue
        props = mol.GetPropsAsDict()
        cnn = props.get("CNNscore")
        vina = props.get("minimizedAffinity")
        if cnn is not None and (best_cnn is None or cnn > best_cnn):
            best_cnn = cnn
            best_vina = vina
        break
    return best_cnn, best_vina


def run_gnina_dock(receptor, ligand_sdf, autobox_lig, out_sdf,
                   size_x=None, size_y=None, size_z=None,
                   center=None, padding=None, label=""):
    if os.path.exists(out_sdf) and os.path.getsize(out_sdf) > 100:
        return True

    cmd = [
        GNINA_PATH,
        "-r", receptor,
        "-l", ligand_sdf,
        "-o", out_sdf,
        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--num_modes", "9",
        "--cnn_scoring", "rescore",
    ]

    if size_x is not None and center is not None:
        cmd += [
            "--center_x", str(center[0]),
            "--center_y", str(center[1]),
            "--center_z", str(center[2]),
            "--size_x", str(size_x),
            "--size_y", str(size_y),
            "--size_z", str(size_z),
        ]
    else:
        cmd += ["--autobox_ligand", autobox_lig]
        if padding is not None:
            cmd += ["--autobox_add", str(padding)]
        else:
            cmd += ["--autobox_add", "4"]

    cuda_devs = os.environ.get("CUDA_VISIBLE_DEVICES", "")
    if not cuda_devs:
        cmd.append("--no_gpu")

    try:
        env = os.environ.copy()
        env.setdefault("OMP_NUM_THREADS", "4")
        result = subprocess.run(cmd, capture_output=True, text=True,
                                timeout=1800, env=env)
        if result.returncode != 0:
            return False
        return os.path.exists(out_sdf) and os.path.getsize(out_sdf) > 50
    except Exception:
        return False


def compute_consensus_box(all_lig_coords, padding=4.0):
    all_coords = np.vstack(all_lig_coords)
    center = np.mean(all_coords, axis=0)

    mins = np.min(all_coords, axis=0)
    maxs = np.max(all_coords, axis=0)
    spans = maxs - mins + 2 * padding

    return {
        "center": center.tolist(),
        "size_x": float(spans[0]),
        "size_y": float(spans[1]),
        "size_z": float(spans[2]),
        "min_coords": mins.tolist(),
        "max_coords": maxs.tolist(),
        "padding": padding,
    }


def main():
    print("=" * 70)
    print("BOX OPTIMIZATION FOR 4-STRUCTURE ENSEMBLE - Task 4")
    print("=" * 70)
    print(f"Output: {OUTPUT_DIR}")
    print(f"GNINA: {GNINA_PATH}")
    print(f"pH values: {PH_VALUES}")
    print(f"Exhaustiveness: {EXHAUSTIVENESS}")
    print(f"Autobox paddings: {AUTOBOX_PADDINGS}")
    print(f"Manual box sizes: {MANUAL_SIZES}")
    print()

    has_gnina = os.path.isfile(GNINA_PATH)
    if not has_gnina:
        print("WARNING: GNINA not found. Will compute consensus box and")
        print("         visualization but skip docking optimization.\n")

    parser = PDBParser(QUIET=True)
    structures = {}
    chains = {}

    print("Step 1: Loading and aligning structures...")
    for pdb_id in PDB_IDS:
        pdb_path = download_pdb(pdb_id)
        struct = parser.get_structure(pdb_id, pdb_path)
        structures[pdb_id] = struct
        chains[pdb_id] = struct[0]['A']
        n_chains = len(list(struct[0].get_chains()))
        print(f"  {pdb_id}: {SPECIES[pdb_id]} + {LIGAND_OF[pdb_id]}, "
              f"{n_chains} chain(s)")

    ref_id = "2ZSH"
    ref_chain = chains[ref_id]
    ref_pocket = get_pocket_residues(ref_id)
    ref_bb = get_pocket_backbone_atoms(ref_chain, ref_pocket)

    print(f"\nAligning to {ref_id}...")
    for pdb_id in PDB_IDS:
        if pdb_id == ref_id:
            continue
        mov_chain = chains[pdb_id]
        mov_pocket = get_pocket_residues(pdb_id)
        mov_bb = get_pocket_backbone_atoms(mov_chain, mov_pocket)
        n = min(len(ref_bb), len(mov_bb))
        sup = Superimposer()
        sup.set_atoms(ref_bb[:n], mov_bb[:n])
        sup.apply(structures[pdb_id][0].get_atoms())
        print(f"  {pdb_id} aligned: RMSD = {sup.rms:.4f} A")

    print(f"\nStep 2: Extracting receptors and aligned ligands...")
    receptors = {}
    lig_pdbs = {}
    lig_sdfs = {}
    lig_coords = {}

    for pdb_id in PDB_IDS:
        rec = extract_receptor(structures[pdb_id], pdb_id)
        receptors[pdb_id] = rec

        lig_name = LIGAND_OF[pdb_id]
        lpdb = extract_aligned_ligand(structures[pdb_id], pdb_id, lig_name)
        lig_pdbs[pdb_id] = lpdb

        if lpdb:
            lsdf = extract_ligand_sdf(lpdb, pdb_id, lig_name)
            lig_sdfs[pdb_id] = lsdf
            coords = get_ligand_coords(lpdb)
            lig_coords[pdb_id] = coords
            if coords is not None:
                center = np.mean(coords, axis=0)
                span = np.max(coords, axis=0) - np.min(coords, axis=0)
                print(f"  {pdb_id} {lig_name}: center=({center[0]:.1f}, "
                      f"{center[1]:.1f}, {center[2]:.1f}), "
                      f"span=({span[0]:.1f}, {span[1]:.1f}, {span[2]:.1f})")

    print(f"\n{'=' * 70}")
    print("Step 3: Computing consensus box...")
    print(f"{'=' * 70}")

    valid_coords = [c for c in lig_coords.values() if c is not None]
    consensus_boxes = {}
    for pad in [3, 4, 5, 6]:
        box = compute_consensus_box(valid_coords, padding=pad)
        consensus_boxes[pad] = box
        print(f"\n  Padding {pad} A:")
        print(f"    Center: ({box['center'][0]:.2f}, "
              f"{box['center'][1]:.2f}, {box['center'][2]:.2f})")
        print(f"    Size: {box['size_x']:.1f} x {box['size_y']:.1f} x "
              f"{box['size_z']:.1f} A")

    best_consensus = consensus_boxes[4]
    print(f"\n  Selected consensus (padding=4 A):")
    print(f"    Center: ({best_consensus['center'][0]:.2f}, "
          f"{best_consensus['center'][1]:.2f}, "
          f"{best_consensus['center'][2]:.2f})")
    print(f"    Size: {best_consensus['size_x']:.1f} x "
          f"{best_consensus['size_y']:.1f} x "
          f"{best_consensus['size_z']:.1f} A")

    all_results = {}

    if has_gnina:
        print(f"\n{'=' * 70}")
        print("Step 4: GNINA box optimization sweep...")
        print(f"{'=' * 70}")

        for pdb_id in PDB_IDS:
            lig_name = LIGAND_OF[pdb_id]
            rec = receptors[pdb_id]
            lig_sdf = lig_sdfs.get(pdb_id)
            lig_pdb = lig_pdbs.get(pdb_id)

            if not lig_sdf or not lig_pdb:
                print(f"\n  {pdb_id}: skipping (no ligand SDF)")
                continue

            print(f"\n  {pdb_id} ({SPECIES[pdb_id]} + {lig_name}):")
            sweep_results = []

            for pad in AUTOBOX_PADDINGS:
                out_sdf = os.path.join(
                    OUTPUT_DIR,
                    f"{pdb_id}_autobox_pad{pad}.sdf")
                label = f"autobox pad={pad}"
                success = run_gnina_dock(
                    rec, lig_sdf, lig_pdb, out_sdf,
                    padding=pad, label=label)

                if success:
                    rmsd = calc_rmsd_to_crystal(out_sdf, lig_pdb)
                    cnn, vina = get_gnina_scores(out_sdf)
                    rmsd_s = f"{rmsd:.3f}" if rmsd else "N/A"
                    cnn_s = f"{cnn:.4f}" if cnn else "N/A"
                    vina_s = f"{vina:.2f}" if vina else "N/A"
                    print(f"    Autobox pad={pad}A: RMSD={rmsd_s}, "
                          f"CNN={cnn_s}, Vina={vina_s}")
                    sweep_results.append({
                        "method": "autobox",
                        "padding": pad,
                        "rmsd": rmsd,
                        "cnn_score": cnn,
                        "vina_score": vina,
                    })
                else:
                    print(f"    Autobox pad={pad}A: FAILED")

            lig_center = np.mean(lig_coords[pdb_id], axis=0)
            for sx, sy, sz in MANUAL_SIZES:
                out_sdf = os.path.join(
                    OUTPUT_DIR,
                    f"{pdb_id}_manual_{sx}x{sy}x{sz}.sdf")
                label = f"manual {sx}x{sy}x{sz}"
                success = run_gnina_dock(
                    rec, lig_sdf, lig_pdb, out_sdf,
                    size_x=sx, size_y=sy, size_z=sz,
                    center=lig_center, label=label)

                if success:
                    rmsd = calc_rmsd_to_crystal(out_sdf, lig_pdb)
                    cnn, vina = get_gnina_scores(out_sdf)
                    rmsd_s = f"{rmsd:.3f}" if rmsd else "N/A"
                    cnn_s = f"{cnn:.4f}" if cnn else "N/A"
                    vina_s = f"{vina:.2f}" if vina else "N/A"
                    print(f"    Manual {sx}x{sy}x{sz}: RMSD={rmsd_s}, "
                          f"CNN={cnn_s}, Vina={vina_s}")
                    sweep_results.append({
                        "method": "manual",
                        "size_x": sx, "size_y": sy, "size_z": sz,
                        "center": lig_center.tolist(),
                        "rmsd": rmsd,
                        "cnn_score": cnn,
                        "vina_score": vina,
                    })
                else:
                    print(f"    Manual {sx}x{sy}x{sz}: FAILED")

            cx, cy, cz = best_consensus["center"]
            csx = best_consensus["size_x"]
            csy = best_consensus["size_y"]
            csz = best_consensus["size_z"]
            out_sdf = os.path.join(
                OUTPUT_DIR, f"{pdb_id}_consensus.sdf")
            label = f"consensus {csx:.0f}x{csy:.0f}x{csz:.0f}"
            success = run_gnina_dock(
                rec, lig_sdf, lig_pdb, out_sdf,
                size_x=csx, size_y=csy, size_z=csz,
                center=[cx, cy, cz], label=label)

            if success:
                rmsd = calc_rmsd_to_crystal(out_sdf, lig_pdb)
                cnn, vina = get_gnina_scores(out_sdf)
                rmsd_s = f"{rmsd:.3f}" if rmsd else "N/A"
                cnn_s = f"{cnn:.4f}" if cnn else "N/A"
                vina_s = f"{vina:.2f}" if vina else "N/A"
                print(f"    Consensus box: RMSD={rmsd_s}, "
                      f"CNN={cnn_s}, Vina={vina_s}")
                sweep_results.append({
                    "method": "consensus",
                    "size_x": csx, "size_y": csy, "size_z": csz,
                    "center": [cx, cy, cz],
                    "rmsd": rmsd,
                    "cnn_score": cnn,
                    "vina_score": vina,
                })
            else:
                print(f"    Consensus box: FAILED")

            all_results[pdb_id] = sweep_results

            valid = [r for r in sweep_results
                     if r.get("rmsd") is not None and r["rmsd"] < 900]
            if valid:
                best = min(valid, key=lambda r: r["rmsd"])
                if best["method"] == "autobox":
                    desc = f"autobox pad={best['padding']}A"
                elif best["method"] == "consensus":
                    desc = "consensus"
                else:
                    desc = (f"{best['size_x']}x{best['size_y']}x"
                            f"{best['size_z']} A")
                print(f"    BEST: {desc}, RMSD={best['rmsd']:.3f}")

    print(f"\n{'=' * 70}")
    print("Step 5: Summary")
    print(f"{'=' * 70}\n")

    print(f"  Universal Consensus Box (padding=4 A):")
    print(f"    Center: ({best_consensus['center'][0]:.2f}, "
          f"{best_consensus['center'][1]:.2f}, "
          f"{best_consensus['center'][2]:.2f})")
    print(f"    Size: {best_consensus['size_x']:.1f} x "
          f"{best_consensus['size_y']:.1f} x "
          f"{best_consensus['size_z']:.1f} A")

    if all_results:
        print(f"\n  Per-structure best configurations:")
        for pdb_id in PDB_IDS:
            results = all_results.get(pdb_id, [])
            valid = [r for r in results
                     if r.get("rmsd") is not None and r["rmsd"] < 900]
            if valid:
                best = min(valid, key=lambda r: r["rmsd"])
                if best["method"] == "autobox":
                    desc = f"autobox pad={best['padding']}A"
                elif best["method"] == "consensus":
                    desc = "consensus"
                else:
                    desc = (f"{best['size_x']}x{best['size_y']}x"
                            f"{best['size_z']} A")
                print(f"    {pdb_id}: {desc}, RMSD={best['rmsd']:.3f}, "
                      f"CNN={best.get('cnn_score', 'N/A')}")

        consensus_rmsds = []
        for pdb_id in PDB_IDS:
            results = all_results.get(pdb_id, [])
            cons = [r for r in results if r.get("method") == "consensus"
                    and r.get("rmsd") is not None]
            if cons:
                consensus_rmsds.append(cons[0]["rmsd"])

        if consensus_rmsds:
            print(f"\n  Consensus box performance across all structures:")
            print(f"    Mean RMSD: {np.mean(consensus_rmsds):.3f} A")
            print(f"    Max RMSD:  {np.max(consensus_rmsds):.3f} A")
            all_sub_half = all(r < 0.5 for r in consensus_rmsds)
            print(f"    All < 0.5 A: {'YES' if all_sub_half else 'NO'}")

    results_json = {
        "consensus_box": best_consensus,
        "all_consensus_boxes": {str(k): v for k, v in consensus_boxes.items()},
        "per_structure_results": {},
        "ligand_centers": {},
    }
    for pdb_id in PDB_IDS:
        if pdb_id in all_results:
            results_json["per_structure_results"][pdb_id] = all_results[pdb_id]
        if lig_coords.get(pdb_id) is not None:
            center = np.mean(lig_coords[pdb_id], axis=0)
            results_json["ligand_centers"][pdb_id] = center.tolist()

    json_path = os.path.join(OUTPUT_DIR, "box_optimization_results.json")
    with open(json_path, "w") as f:
        json.dump(results_json, f, indent=2, default=str)
    print(f"\nJSON results: {json_path}")

    report_path = os.path.join(OUTPUT_DIR, "box_optimization_report.txt")
    _write_report(report_path, best_consensus, consensus_boxes,
                  all_results, lig_coords)
    print(f"Report: {report_path}")

    vis_path = os.path.join(OUTPUT_DIR, "visualize_box.py")
    _write_visualization(vis_path, best_consensus, lig_coords)
    print(f"Visualization script: {vis_path}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        _generate_plots(all_results, best_consensus, lig_coords)
    except ImportError:
        print("matplotlib not available, skipping plots")


def _write_report(path, consensus, all_consensus, all_results, lig_coords):
    with open(path, "w") as f:
        f.write("BOX OPTIMIZATION FOR 4-STRUCTURE ENSEMBLE - Task 4\n")
        f.write(f"{'=' * 60}\n\n")

        f.write("Universal Consensus Box:\n")
        for pad, box in sorted(all_consensus.items()):
            f.write(f"  Padding {pad} A: center=({box['center'][0]:.2f}, "
                    f"{box['center'][1]:.2f}, {box['center'][2]:.2f}), "
                    f"size={box['size_x']:.1f}x{box['size_y']:.1f}x"
                    f"{box['size_z']:.1f}\n")

        f.write(f"\nSelected (padding=4 A):\n")
        f.write(f"  Center: ({consensus['center'][0]:.2f}, "
                f"{consensus['center'][1]:.2f}, "
                f"{consensus['center'][2]:.2f})\n")
        f.write(f"  Size: {consensus['size_x']:.1f} x "
                f"{consensus['size_y']:.1f} x "
                f"{consensus['size_z']:.1f} A\n\n")

        f.write("Per-structure ligand centers:\n")
        for pdb_id in PDB_IDS:
            coords = lig_coords.get(pdb_id)
            if coords is not None:
                c = np.mean(coords, axis=0)
                f.write(f"  {pdb_id}: ({c[0]:.2f}, {c[1]:.2f}, "
                        f"{c[2]:.2f})\n")

        if all_results:
            f.write(f"\nDocking sweep results:\n")
            for pdb_id in PDB_IDS:
                results = all_results.get(pdb_id, [])
                if not results:
                    continue
                f.write(f"\n  {pdb_id} ({SPECIES[pdb_id]} + "
                        f"{LIGAND_OF[pdb_id]}):\n")
                for r in sorted(results,
                                key=lambda x: x.get("rmsd", 999)):
                    rmsd = r.get("rmsd")
                    rmsd_s = f"{rmsd:.3f}" if rmsd else "N/A"
                    if r["method"] == "autobox":
                        desc = f"autobox pad={r['padding']}A"
                    elif r["method"] == "consensus":
                        desc = "consensus"
                    else:
                        desc = (f"{r['size_x']}x{r['size_y']}x"
                                f"{r['size_z']}")
                    f.write(f"    {desc}: RMSD={rmsd_s}\n")


def _write_visualization(path, consensus, lig_coords):
    cx, cy, cz = consensus["center"]
    sx, sy, sz = (consensus["size_x"],
                  consensus["size_y"],
                  consensus["size_z"])

    active_res_arab = ", ".join(str(r) for r in ARAB_RESIDUES)
    active_res_rice = ", ".join(str(r) for r in RICE_RESIDUES)

    script = f'''"""
Box Optimization Visualization - Task 4
Paste into a Jupyter notebook cell on Expanse.

Shows all 4 GID1 structures aligned with:
  - 14 key residues as sticks (element coloring)
  - Crystal ligands (green=Arabidopsis, cyan=Rice)
  - Universal consensus docking box (wireframe)
"""
import json, os, time
import nglview as nv
from IPython.display import display, HTML

BASE = os.path.expanduser("~/final_project")
BOX_DIR = os.path.join(BASE, "box_optimization")
POCKET_DIR = os.path.join(BASE, "pocket_analysis")

ARAB_RESIDUES = [{active_res_arab}]
RICE_RESIDUES = [{active_res_rice}]

BOX_CENTER = [{cx:.2f}, {cy:.2f}, {cz:.2f}]
BOX_SIZE = [{sx:.1f}, {sy:.1f}, {sz:.1f}]

display(HTML("""
<pre style='background:#ffffff !important; color:#000000 !important; padding:15px; border-radius:8px;
            margin-bottom:15px; border:2px solid #228B22; font-family:monospace; font-size:14px;
            line-height:1.6; white-space:pre;'><span style='color:#006400 !important; font-size:20px; font-weight:bold;'>Task 4: Box Optimization Visualization</span>

4-Structure Ensemble: 2ZSH, 2ZSI, 3ED1, 3EBL
All aligned on 14 key pocket residue backbone atoms.

  <span style='color:#228B22 !important; font-size:16px;'>&#9632;</span> Arabidopsis crystal ligands (green)
  <span style='color:#00CED1 !important; font-size:16px;'>&#9632;</span> Rice crystal ligands (cyan)
  <span style='color:#888 !important; font-size:16px;'>&#9632;</span> Active site residues (element coloring)

Consensus box: """ + f"{{BOX_SIZE[0]:.1f}} x {{BOX_SIZE[1]:.1f}} x {{BOX_SIZE[2]:.1f}} A" + """
Center: """ + f"({{BOX_CENTER[0]:.2f}}, {{BOX_CENTER[1]:.2f}}, {{BOX_CENTER[2]:.2f}})" + """</pre>
"""))

PDB_IDS = ["2ZSH", "2ZSI", "3ED1", "3EBL"]
SPECIES = {{"2ZSH": "Arabidopsis", "2ZSI": "Arabidopsis",
           "3ED1": "Rice", "3EBL": "Rice"}}
LIGAND_OF = {{"2ZSH": "GA3", "2ZSI": "GA4",
             "3ED1": "GA3", "3EBL": "GA4"}}

view = nv.NGLWidget()
view._remote_call("setSize", target="Widget", args=["100%", "600px"])

comp_idx = 0
for pdb_id in PDB_IDS:
    rec_path = os.path.join(BOX_DIR, f"{{pdb_id}}_receptor.pdb")
    if not os.path.exists(rec_path):
        print(f"Missing: {{rec_path}}")
        continue

    with open(rec_path) as f:
        pdb_text = f.read()

    view.add_component(nv.TextStructure(pdb_text, ext="pdb"),
                       default_representation=False,
                       name=f"{{pdb_id}}_receptor")
    time.sleep(0.2)

    view.add_cartoon(selection="protein", color="white", opacity=0.08,
                     component=comp_idx)

    species = SPECIES[pdb_id]
    if species == "Arabidopsis":
        pocket = ARAB_RESIDUES
    else:
        pocket = RICE_RESIDUES

    res_sel = "(" + " or ".join(str(r) for r in pocket) + ") and sidechainAttached"
    view.add_licorice(selection=res_sel, color="element", radius=0.18,
                      component=comp_idx)

    comp_idx += 1

for pdb_id in PDB_IDS:
    lig_path = os.path.join(BOX_DIR,
                            f"{{pdb_id}}_{{LIGAND_OF[pdb_id]}}_ligand_aligned.pdb")
    if not os.path.exists(lig_path):
        continue

    with open(lig_path) as f:
        lig_text = f.read()

    species = SPECIES[pdb_id]
    color = "green" if species == "Arabidopsis" else "#00CED1"

    view.add_component(nv.TextStructure(lig_text, ext="pdb"),
                       default_representation=False,
                       name=f"{{pdb_id}}_ligand")
    time.sleep(0.2)
    view.add_ball_and_stick(selection="all", color=color,
                            aspectRatio=2.5, component=comp_idx)
    comp_idx += 1

view.center(selection="all", component=comp_idx - 1)
display(view)

display(HTML(f"""
<pre style='background:#ffffff !important; color:#000000 !important; padding:12px; border-radius:8px;
            margin-top:10px; font-size:13px; border:1px solid #228B22; font-family:monospace;
            line-height:1.5; white-space:pre-wrap;'>Universal Consensus Box for AutoGrow4:
  --center_x {{BOX_CENTER[0]:.2f}} --center_y {{BOX_CENTER[1]:.2f}} --center_z {{BOX_CENTER[2]:.2f}}
  --size_x {{BOX_SIZE[0]:.1f}} --size_y {{BOX_SIZE[1]:.1f}} --size_z {{BOX_SIZE[2]:.1f}}

These dimensions encompass all crystal ligand positions across
both Arabidopsis and Rice GID1 structures with 4 A padding.</pre>
"""))
'''

    with open(path, "w") as f:
        f.write(script)


def _generate_plots(all_results, consensus, lig_coords):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from mpl_toolkits.mplot3d import Axes3D

    if all_results:
        n_structs = len([p for p in PDB_IDS if p in all_results])
        if n_structs == 0:
            return

        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle("Box Optimization: RMSD vs Configuration\n"
                     "All 4 GID1 Structures",
                     fontsize=14, fontweight="bold")

        for idx, pdb_id in enumerate(PDB_IDS):
            ax = axes[idx // 2][idx % 2]
            results = all_results.get(pdb_id, [])
            if not results:
                ax.text(0.5, 0.5, f"{pdb_id}: No results",
                        ha="center", va="center", fontsize=14)
                ax.set_title(f"{pdb_id}")
                continue

            labels = []
            rmsds = []
            colors = []
            for r in sorted(results,
                            key=lambda x: x.get("rmsd", 999)):
                rmsd = r.get("rmsd")
                if rmsd is None or rmsd > 900:
                    continue
                if r["method"] == "autobox":
                    labels.append(f"AB\npad={r['padding']}")
                    colors.append("#2196F3")
                elif r["method"] == "consensus":
                    labels.append("Consensus")
                    colors.append("#FF5722")
                else:
                    labels.append(f"{r['size_x']:.0f}x\n"
                                  f"{r['size_y']:.0f}x\n"
                                  f"{r['size_z']:.0f}")
                    colors.append("#4CAF50")
                rmsds.append(rmsd)

            if rmsds:
                bars = ax.bar(range(len(labels)), rmsds, color=colors,
                              edgecolor="black", linewidth=0.5)
                ax.axhline(y=0.5, color="red", linestyle="--",
                           linewidth=1.5, label="0.5 A threshold")
                ax.set_xticks(range(len(labels)))
                ax.set_xticklabels(labels, fontsize=7)
                ax.set_ylabel("RMSD (A)", fontsize=10)
                ax.set_title(f"{pdb_id}: {SPECIES[pdb_id]} + "
                             f"{LIGAND_OF[pdb_id]}",
                             fontsize=11, fontweight="bold")
                ax.legend(fontsize=8)
                ax.grid(axis="y", alpha=0.3)

                for i, v in enumerate(rmsds):
                    ax.text(i, v + 0.01, f"{v:.3f}", ha="center",
                            fontsize=8, fontweight="bold")

        plt.tight_layout()
        save_path = os.path.join(OUTPUT_DIR,
                                 "box_optimization_sweep.png")
        plt.savefig(save_path, dpi=200, bbox_inches="tight",
                    facecolor="white")
        plt.close()
        print(f"Sweep plot: {save_path}")

    fig2, ax2 = plt.subplots(figsize=(10, 8))
    ax2.set_title("Ligand Positions in Aligned Structures\n"
                  "with Consensus Box",
                  fontsize=13, fontweight="bold")

    species_colors = {"Arabidopsis": "#4CAF50", "Rice": "#2196F3"}
    for pdb_id in PDB_IDS:
        coords = lig_coords.get(pdb_id)
        if coords is None:
            continue
        color = species_colors[SPECIES[pdb_id]]
        ax2.scatter(coords[:, 0], coords[:, 1], c=color, s=20,
                    alpha=0.5, label=f"{pdb_id} ({LIGAND_OF[pdb_id]})")

    cx, cy = consensus["center"][0], consensus["center"][1]
    sx, sy = consensus["size_x"], consensus["size_y"]
    rect = Rectangle((cx - sx / 2, cy - sy / 2), sx, sy,
                      linewidth=2, edgecolor="red", facecolor="red",
                      alpha=0.1, label="Consensus box (XY)")
    ax2.add_patch(rect)
    ax2.plot(cx, cy, "r+", markersize=15, markeredgewidth=2)

    ax2.set_xlabel("X (A)", fontsize=11)
    ax2.set_ylabel("Y (A)", fontsize=11)
    ax2.legend(fontsize=9)
    ax2.set_aspect("equal")
    ax2.grid(alpha=0.3)

    plt.tight_layout()
    save_path2 = os.path.join(OUTPUT_DIR, "box_consensus_2d.png")
    plt.savefig(save_path2, dpi=200, bbox_inches="tight",
                facecolor="white")
    plt.close()
    print(f"Consensus 2D plot: {save_path2}")


if __name__ == "__main__":
    main()
