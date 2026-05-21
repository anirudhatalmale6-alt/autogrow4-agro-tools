#!/usr/bin/env python3
"""
Cross-Docking IFP Matrix - Task 3

Rigid cross-docking validation for universal ligand feasibility.
Docks GA3 and GA4 across all 4 GID1 structures using GNINA.
Computes Interaction Fingerprints (IFP) and Tanimoto similarity.

Docking matrix (8 runs):
  Self-dock:  GA3->2ZSH, GA4->2ZSI, GA3->3ED1, GA4->3EBL
  Cross-dock: GA4->2ZSH, GA3->2ZSI, GA4->3ED1, GA3->3EBL

Decision rule:
  Tc >= 0.85 -> rigid ensemble model acceptable
  Tc <  0.85 -> induced fit model needed
"""

import os
import sys
import subprocess
import json
import numpy as np
from collections import OrderedDict

from Bio.PDB import PDBParser, PDBIO, Select
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolfiles

OUTPUT_DIR = os.path.expanduser("~/final_project/cross_docking")
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

DOCKING_MATRIX = [
    ("2ZSH", "GA3", True,  "Self-dock: GA3 in Arabidopsis+GA3"),
    ("2ZSH", "GA4", False, "Cross-dock: GA4 in Arabidopsis+GA3 receptor"),
    ("2ZSI", "GA4", True,  "Self-dock: GA4 in Arabidopsis+GA4"),
    ("2ZSI", "GA3", False, "Cross-dock: GA3 in Arabidopsis+GA4 receptor"),
    ("3ED1", "GA3", True,  "Self-dock: GA3 in Rice+GA3"),
    ("3ED1", "GA4", False, "Cross-dock: GA4 in Rice+GA3 receptor"),
    ("3EBL", "GA4", True,  "Self-dock: GA4 in Rice+GA4"),
    ("3EBL", "GA3", False, "Cross-dock: GA3 in Rice+GA4 receptor"),
]

EXHAUSTIVENESS = 16
NUM_MODES = 9

IFP_INTERACTIONS = [
    "hbond_don",
    "hbond_acc",
    "hydrophobic",
    "aromatic",
    "ionic",
]

HBOND_CUTOFF = 3.5
HYDROPHOBIC_CUTOFF = 4.5
AROMATIC_CUTOFF = 5.5
IONIC_CUTOFF = 4.0

HBOND_ELEMENTS = {"N", "O", "S"}
HYDROPHOBIC_RESIDUES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"}
AROMATIC_RESIDUES = {"PHE", "TYR", "TRP", "HIS"}
POSITIVE_RESIDUES = {"ARG", "LYS", "HIS"}
NEGATIVE_RESIDUES = {"ASP", "GLU"}


def download_pdb(pdb_id):
    import urllib.request
    path = os.path.join(POCKET_DIR, f"{pdb_id}.pdb")
    if os.path.exists(path) and os.path.getsize(path) > 1000:
        return path
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"  Downloading {pdb_id}...")
    urllib.request.urlretrieve(url, path)
    return path


class ProteinOnlySelect(Select):
    def accept_chain(self, chain):
        return chain.id == "A"

    def accept_residue(self, residue):
        return residue.id[0] == " "


def extract_receptor(pdb_path, pdb_id):
    out_path = os.path.join(OUTPUT_DIR, f"{pdb_id}_receptor.pdb")
    if os.path.exists(out_path) and os.path.getsize(out_path) > 500:
        return out_path
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure(pdb_id, pdb_path)
    io = PDBIO()
    io.set_structure(struct)
    io.save(out_path, ProteinOnlySelect())
    print(f"    Receptor saved: {out_path}")
    return out_path


def extract_crystal_ligand(pdb_path, pdb_id, lig_name):
    out_pdb = os.path.join(OUTPUT_DIR, f"{pdb_id}_{lig_name}_crystal.pdb")
    out_sdf = os.path.join(OUTPUT_DIR, f"{pdb_id}_{lig_name}_crystal.sdf")

    with open(pdb_path) as f:
        lines = f.readlines()

    lig_lines = []
    for line in lines:
        if line.startswith("HETATM"):
            chain = line[21]
            resname = line[17:20].strip()
            if chain == "A" and resname == lig_name:
                lig_lines.append(line)

    if not lig_lines:
        print(f"    WARNING: No {lig_name} found in chain A of {pdb_id}")
        return None, None

    with open(out_pdb, "w") as f:
        for line in lig_lines:
            f.write(line)
        f.write("END\n")

    mol = Chem.MolFromPDBFile(out_pdb, removeHs=False, sanitize=False)
    if mol is None:
        mol = Chem.MolFromPDBFile(out_pdb, removeHs=True, sanitize=False)
    if mol is not None:
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            pass
        writer = Chem.SDWriter(out_sdf)
        writer.write(mol)
        writer.close()
        print(f"    Ligand {lig_name}: {mol.GetNumAtoms()} atoms -> {out_sdf}")
    else:
        print(f"    WARNING: Could not parse {lig_name} with RDKit, "
              f"using PDB directly")
        out_sdf = None

    return out_pdb, out_sdf


def get_ligand_center(lig_pdb_path):
    coords = []
    with open(lig_pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
    if not coords:
        return None
    return np.mean(coords, axis=0)


def run_gnina(receptor_pdb, ligand_sdf, autobox_lig, out_sdf, label):
    if os.path.exists(out_sdf) and os.path.getsize(out_sdf) > 100:
        print(f"    {label}: already docked, skipping")
        return True

    cmd = [
        GNINA_PATH,
        "-r", receptor_pdb,
        "-l", ligand_sdf,
        "--autobox_ligand", autobox_lig,
        "--autobox_add", "4",
        "-o", out_sdf,
        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--num_modes", str(NUM_MODES),
        "--cnn_scoring", "rescore",
    ]

    cuda_devs = os.environ.get("CUDA_VISIBLE_DEVICES", "")
    if not cuda_devs:
        cmd.append("--no_gpu")

    try:
        env = os.environ.copy()
        env.setdefault("OMP_NUM_THREADS", "4")

        print(f"    {label}: running GNINA...")
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=3600, env=env)

        if result.returncode != 0:
            print(f"    {label}: GNINA failed (exit {result.returncode})")
            stderr_lines = result.stderr.strip().split("\n")[-5:]
            for line in stderr_lines:
                print(f"      {line}")
            return False

        if not os.path.exists(out_sdf) or os.path.getsize(out_sdf) < 50:
            print(f"    {label}: no output SDF produced")
            return False

        print(f"    {label}: docking complete")
        return True

    except subprocess.TimeoutExpired:
        print(f"    {label}: GNINA timed out after 1 hour")
        return False
    except Exception as e:
        print(f"    {label}: error - {e}")
        return False


def parse_gnina_scores(sdf_path):
    scores = []
    if not os.path.exists(sdf_path):
        return scores
    suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
    for i, mol in enumerate(suppl):
        if mol is None:
            continue
        entry = {"pose": i}
        for prop in mol.GetPropsAsDict():
            entry[prop] = mol.GetPropsAsDict()[prop]
        if "minimizedAffinity" in entry:
            entry["vina_score"] = entry["minimizedAffinity"]
        if "CNNscore" in entry:
            entry["cnn_score"] = entry["CNNscore"]
        if "CNNaffinity" in entry:
            entry["cnn_affinity"] = entry["CNNaffinity"]
        scores.append(entry)
    return scores


def get_pocket_residues(pdb_id):
    if SPECIES[pdb_id] == "Arabidopsis":
        return list(ARAB_RESIDUES)
    return list(RICE_RESIDUES)


def get_residue_atoms(chain, pocket_residues):
    res_data = {}
    for resnum in pocket_residues:
        try:
            res = chain[(' ', resnum, ' ')]
            atoms = []
            for atom in res:
                if atom.element.strip() == "H":
                    continue
                atoms.append({
                    "name": atom.name,
                    "element": atom.element.strip().upper(),
                    "coord": atom.get_vector().get_array(),
                })
            res_data[resnum] = {
                "resname": res.get_resname(),
                "resnum": resnum,
                "atoms": atoms,
            }
        except KeyError:
            pass
    return res_data


def get_ligand_atoms_from_sdf(sdf_path, pose_idx=0):
    suppl = Chem.SDMolSupplier(sdf_path, removeHs=True)
    atoms = []
    for i, mol in enumerate(suppl):
        if i != pose_idx or mol is None:
            continue
        conf = mol.GetConformer()
        for j in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(j)
            pos = conf.GetAtomPosition(j)
            atoms.append({
                "element": atom.GetSymbol().upper(),
                "coord": np.array([pos.x, pos.y, pos.z]),
                "is_aromatic": atom.GetIsAromatic(),
            })
        break
    return atoms


def get_ligand_atoms_from_pdb(pdb_path):
    atoms = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                elem = line[76:78].strip().upper()
                if not elem:
                    name = line[12:16].strip()
                    elem = name[0]
                if elem == "H":
                    continue
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms.append({
                    "element": elem,
                    "coord": np.array([x, y, z]),
                    "is_aromatic": False,
                })
    return atoms


def compute_ifp(res_data, lig_atoms, pocket_residues):
    n_res = len(pocket_residues)
    n_types = len(IFP_INTERACTIONS)
    fp = np.zeros(n_res * n_types, dtype=int)

    lig_coords = np.array([a["coord"] for a in lig_atoms])
    lig_elements = [a["element"] for a in lig_atoms]
    lig_donor = [e in HBOND_ELEMENTS for e in lig_elements]
    lig_carbon = [e == "C" for e in lig_elements]

    for i, resnum in enumerate(pocket_residues):
        if resnum not in res_data:
            continue

        rd = res_data[resnum]
        resname = rd["resname"]

        for atom in rd["atoms"]:
            dists = np.linalg.norm(lig_coords - atom["coord"], axis=1)
            min_dist = np.min(dists)

            if atom["element"] in HBOND_ELEMENTS:
                close_lig = dists < HBOND_CUTOFF
                for j, is_close in enumerate(close_lig):
                    if is_close and lig_donor[j]:
                        fp[i * n_types + 0] = 1
                        break

            if atom["element"] in HBOND_ELEMENTS:
                close_lig = dists < HBOND_CUTOFF
                for j, is_close in enumerate(close_lig):
                    if is_close and lig_elements[j] in HBOND_ELEMENTS:
                        fp[i * n_types + 1] = 1
                        break

            if (resname in HYDROPHOBIC_RESIDUES and
                    atom["element"] == "C"):
                close_lig = dists < HYDROPHOBIC_CUTOFF
                for j, is_close in enumerate(close_lig):
                    if is_close and lig_carbon[j]:
                        fp[i * n_types + 2] = 1
                        break

            if resname in AROMATIC_RESIDUES:
                if min_dist < AROMATIC_CUTOFF:
                    fp[i * n_types + 3] = 1

            if resname in POSITIVE_RESIDUES or resname in NEGATIVE_RESIDUES:
                if min_dist < IONIC_CUTOFF:
                    if any(lig_elements[j] in HBOND_ELEMENTS
                           for j in range(len(lig_elements))
                           if dists[j] < IONIC_CUTOFF):
                        fp[i * n_types + 4] = 1

    return fp


def tanimoto(fp1, fp2):
    a = fp1.astype(bool)
    b = fp2.astype(bool)
    intersection = np.sum(a & b)
    union = np.sum(a | b)
    if union == 0:
        return 1.0
    return float(intersection) / float(union)


def ifp_to_string(fp, pocket_residues, res_data):
    n_types = len(IFP_INTERACTIONS)
    parts = []
    for i, resnum in enumerate(pocket_residues):
        bits = fp[i * n_types:(i + 1) * n_types]
        resname = res_data.get(resnum, {}).get("resname", "???")
        bit_str = "".join(str(b) for b in bits)
        parts.append(f"{resname}{resnum}:{bit_str}")
    return " | ".join(parts)


def main():
    print("=" * 70)
    print("CROSS-DOCKING IFP MATRIX - Task 3")
    print("Rigid Cross-Docking Validation")
    print("=" * 70)
    print(f"Output: {OUTPUT_DIR}")
    print(f"GNINA: {GNINA_PATH}")
    print(f"Exhaustiveness: {EXHAUSTIVENESS}")
    print(f"Decision: Tc >= 0.85 = rigid OK, < 0.85 = need induced fit")
    print()

    if not os.path.isfile(GNINA_PATH):
        print(f"ERROR: GNINA not found at {GNINA_PATH}")
        print("Download with: wget https://github.com/gnina/gnina/releases/"
              "download/v1.1/gnina -O ~/final_project/gnina && "
              "chmod +x ~/final_project/gnina")
        sys.exit(1)

    parser = PDBParser(QUIET=True)

    print("Step 1: Preparing structures...")
    receptors = {}
    crystal_ligs_pdb = {}
    crystal_ligs_sdf = {}

    for pdb_id in PDB_IDS:
        print(f"\n  {pdb_id} ({SPECIES[pdb_id]} + {LIGAND_OF[pdb_id]}):")
        pdb_path = download_pdb(pdb_id)
        receptors[pdb_id] = extract_receptor(pdb_path, pdb_id)
        lig_name = LIGAND_OF[pdb_id]
        pdb_lig, sdf_lig = extract_crystal_ligand(pdb_path, pdb_id, lig_name)
        crystal_ligs_pdb[pdb_id] = pdb_lig
        crystal_ligs_sdf[pdb_id] = sdf_lig

    ga3_sources = [pid for pid in PDB_IDS if LIGAND_OF[pid] == "GA3"]
    ga4_sources = [pid for pid in PDB_IDS if LIGAND_OF[pid] == "GA4"]

    print(f"\n  GA3 ligand sources: {ga3_sources}")
    print(f"  GA4 ligand sources: {ga4_sources}")

    print(f"\n{'=' * 70}")
    print("Step 2: Running GNINA cross-docking matrix...")
    print(f"{'=' * 70}")

    docking_results = {}
    for receptor_id, lig_name, is_self, desc in DOCKING_MATRIX:
        if lig_name == "GA3":
            source_id = ga3_sources[0]
        else:
            source_id = ga4_sources[0]

        lig_sdf = crystal_ligs_sdf[source_id]
        if lig_sdf is None:
            lig_sdf_path = crystal_ligs_pdb[source_id]
        else:
            lig_sdf_path = lig_sdf

        autobox_lig = crystal_ligs_pdb[receptor_id]

        run_label = f"{lig_name}->{receptor_id}"
        out_sdf = os.path.join(OUTPUT_DIR,
                               f"docked_{lig_name}_in_{receptor_id}.sdf")

        print(f"\n  [{run_label}] {desc}")
        success = run_gnina(
            receptors[receptor_id], lig_sdf_path, autobox_lig, out_sdf,
            run_label)

        scores = parse_gnina_scores(out_sdf) if success else []
        best_cnn = None
        best_vina = None
        if scores:
            best = max(scores, key=lambda s: s.get("cnn_score", 0))
            best_cnn = best.get("cnn_score")
            best_vina = best.get("vina_score")
            print(f"    Best CNN score: {best_cnn:.4f}, "
                  f"Vina: {best_vina:.2f} kcal/mol")

        docking_results[(receptor_id, lig_name)] = {
            "receptor": receptor_id,
            "ligand": lig_name,
            "is_self": is_self,
            "description": desc,
            "success": success,
            "output_sdf": out_sdf if success else None,
            "scores": scores,
            "best_cnn": best_cnn,
            "best_vina": best_vina,
        }

    print(f"\n{'=' * 70}")
    print("Step 3: Computing Interaction Fingerprints...")
    print(f"{'=' * 70}")

    structures = {}
    for pdb_id in PDB_IDS:
        pdb_path = download_pdb(pdb_id)
        structures[pdb_id] = parser.get_structure(pdb_id, pdb_path)

    crystal_ifps = {}
    for pdb_id in PDB_IDS:
        chain = structures[pdb_id][0]['A']
        pocket_res = get_pocket_residues(pdb_id)
        res_data = get_residue_atoms(chain, pocket_res)

        lig_atoms = get_ligand_atoms_from_pdb(crystal_ligs_pdb[pdb_id])
        ifp = compute_ifp(res_data, lig_atoms, pocket_res)

        crystal_ifps[pdb_id] = {
            "ifp": ifp,
            "res_data": res_data,
            "pocket_res": pocket_res,
            "n_bits_on": int(np.sum(ifp)),
        }

        print(f"\n  {pdb_id} crystal ({LIGAND_OF[pdb_id]}): "
              f"{np.sum(ifp)}/{len(ifp)} bits on")
        print(f"    {ifp_to_string(ifp, pocket_res, res_data)}")

    docked_ifps = {}
    for (receptor_id, lig_name), result in docking_results.items():
        if not result["success"] or result["output_sdf"] is None:
            continue

        chain = structures[receptor_id][0]['A']
        pocket_res = get_pocket_residues(receptor_id)
        res_data = get_residue_atoms(chain, pocket_res)

        scores = result["scores"]
        if not scores:
            continue

        best_idx = 0
        best_cnn = 0
        for s in scores:
            if s.get("cnn_score", 0) > best_cnn:
                best_cnn = s.get("cnn_score", 0)
                best_idx = s["pose"]

        lig_atoms = get_ligand_atoms_from_sdf(result["output_sdf"], best_idx)
        if not lig_atoms:
            lig_atoms = get_ligand_atoms_from_sdf(result["output_sdf"], 0)

        if not lig_atoms:
            print(f"\n  {lig_name}->{receptor_id}: could not read docked pose")
            continue

        ifp = compute_ifp(res_data, lig_atoms, pocket_res)

        docked_ifps[(receptor_id, lig_name)] = {
            "ifp": ifp,
            "res_data": res_data,
            "pocket_res": pocket_res,
            "is_self": result["is_self"],
            "n_bits_on": int(np.sum(ifp)),
        }

        dock_type = "self" if result["is_self"] else "cross"
        print(f"\n  {lig_name}->{receptor_id} [{dock_type}]: "
              f"{np.sum(ifp)}/{len(ifp)} bits on")
        print(f"    {ifp_to_string(ifp, pocket_res, res_data)}")

    print(f"\n{'=' * 70}")
    print("Step 4: Tanimoto Similarity Matrix")
    print(f"{'=' * 70}")

    tc_results = []
    for (receptor_id, lig_name), docked in docked_ifps.items():
        crystal = crystal_ifps[receptor_id]

        pocket_res_d = docked["pocket_res"]
        pocket_res_c = crystal["pocket_res"]

        if len(docked["ifp"]) == len(crystal["ifp"]):
            tc = tanimoto(docked["ifp"], crystal["ifp"])
        else:
            tc = tanimoto(docked["ifp"], crystal["ifp"][:len(docked["ifp"])])

        dock_type = "SELF" if docked["is_self"] else "CROSS"
        status = "PASS" if tc >= 0.85 else "FAIL"

        result = docking_results[(receptor_id, lig_name)]
        tc_results.append({
            "receptor": receptor_id,
            "ligand": lig_name,
            "is_self": docked["is_self"],
            "dock_type": dock_type,
            "tanimoto": tc,
            "status": status,
            "crystal_bits": crystal["n_bits_on"],
            "docked_bits": docked["n_bits_on"],
            "cnn_score": result.get("best_cnn"),
            "vina_score": result.get("best_vina"),
        })

        print(f"\n  {lig_name} -> {receptor_id} [{dock_type}]:")
        print(f"    Tanimoto: {tc:.4f}  [{status}]")
        print(f"    Crystal bits: {crystal['n_bits_on']}, "
              f"Docked bits: {docked['n_bits_on']}")

    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print(f"{'=' * 70}\n")

    print(f"  {'Run':<25} {'Type':<7} {'Tc':>6} {'CNN':>8} "
          f"{'Vina':>8} {'Status':>7}")
    print(f"  {'─' * 65}")
    for r in tc_results:
        cnn_s = f"{r['cnn_score']:.4f}" if r['cnn_score'] is not None else "N/A"
        vina_s = f"{r['vina_score']:.2f}" if r['vina_score'] is not None else "N/A"
        label = f"{r['ligand']}->{r['receptor']}"
        print(f"  {label:<25} {r['dock_type']:<7} {r['tanimoto']:>6.4f} "
              f"{cnn_s:>8} {vina_s:>8} {r['status']:>7}")

    self_results = [r for r in tc_results if r["is_self"]]
    cross_results = [r for r in tc_results if not r["is_self"]]

    if self_results:
        mean_self_tc = np.mean([r["tanimoto"] for r in self_results])
        print(f"\n  Self-docking mean Tc: {mean_self_tc:.4f}")

    if cross_results:
        mean_cross_tc = np.mean([r["tanimoto"] for r in cross_results])
        min_cross_tc = np.min([r["tanimoto"] for r in cross_results])
        all_pass = all(r["status"] == "PASS" for r in cross_results)

        print(f"  Cross-docking mean Tc: {mean_cross_tc:.4f}")
        print(f"  Cross-docking min Tc:  {min_cross_tc:.4f}")
        print(f"  All cross >= 0.85:     {'YES' if all_pass else 'NO'}")

    print(f"\n{'=' * 70}")
    if cross_results and all(r["status"] == "PASS" for r in cross_results):
        print("  DECISION: RIGID ENSEMBLE MODEL is ACCEPTABLE.")
        print("  All cross-docking Tanimoto coefficients >= 0.85.")
        print("  The 4 PDB files can be used as a rigid ensemble for")
        print("  universal ligand design without induced fit modeling.")
    elif cross_results:
        failing = [r for r in cross_results if r["status"] == "FAIL"]
        print("  DECISION: INDUCED FIT MODEL may be NECESSARY.")
        print(f"  {len(failing)} cross-docking run(s) below 0.85 threshold:")
        for f in failing:
            print(f"    {f['ligand']}->{f['receptor']}: Tc = {f['tanimoto']:.4f}")
    else:
        print("  DECISION: Insufficient data for conclusion.")
    print(f"{'=' * 70}")

    report_path = os.path.join(OUTPUT_DIR, "cross_docking_ifp_report.txt")
    _write_report(report_path, tc_results, crystal_ifps, docked_ifps)
    print(f"\nReport saved to: {report_path}")

    results_json = os.path.join(OUTPUT_DIR, "cross_docking_results.json")
    json_data = []
    for r in tc_results:
        json_data.append({
            "receptor": r["receptor"],
            "ligand": r["ligand"],
            "is_self": r["is_self"],
            "tanimoto": r["tanimoto"],
            "cnn_score": r["cnn_score"],
            "vina_score": r["vina_score"],
            "status": r["status"],
        })
    with open(results_json, "w") as f:
        json.dump(json_data, f, indent=2)
    print(f"JSON results saved to: {results_json}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        _generate_plots(tc_results)
    except ImportError:
        print("matplotlib not available, skipping plots")


def _write_report(path, tc_results, crystal_ifps, docked_ifps):
    with open(path, "w") as f:
        f.write("CROSS-DOCKING IFP MATRIX - Task 3\n")
        f.write("Rigid Cross-Docking Validation\n")
        f.write(f"{'=' * 60}\n\n")

        f.write("Parameters:\n")
        f.write(f"  GNINA exhaustiveness: {EXHAUSTIVENESS}\n")
        f.write(f"  CNN scoring: rescore\n")
        f.write(f"  IFP interactions: {IFP_INTERACTIONS}\n")
        f.write(f"  H-bond cutoff: {HBOND_CUTOFF} A\n")
        f.write(f"  Hydrophobic cutoff: {HYDROPHOBIC_CUTOFF} A\n")
        f.write(f"  Aromatic cutoff: {AROMATIC_CUTOFF} A\n")
        f.write(f"  Ionic cutoff: {IONIC_CUTOFF} A\n")
        f.write(f"  Decision threshold: Tc >= 0.85\n\n")

        f.write("Crystal Reference IFPs:\n")
        for pdb_id in PDB_IDS:
            crystal = crystal_ifps[pdb_id]
            f.write(f"  {pdb_id} ({LIGAND_OF[pdb_id]}): "
                    f"{crystal['n_bits_on']}/{len(crystal['ifp'])} bits\n")
            f.write(f"    {ifp_to_string(crystal['ifp'], crystal['pocket_res'], crystal['res_data'])}\n")
        f.write("\n")

        f.write("Tanimoto Results:\n")
        for r in tc_results:
            label = f"{r['ligand']}->{r['receptor']}"
            f.write(f"  {label}: Tc={r['tanimoto']:.4f} "
                    f"[{r['dock_type']}] [{r['status']}]\n")

        cross_results = [r for r in tc_results if not r["is_self"]]
        if cross_results:
            all_pass = all(r["status"] == "PASS" for r in cross_results)
            f.write(f"\nConclusion: Rigid ensemble is "
                    f"{'ACCEPTABLE' if all_pass else 'INSUFFICIENT'}\n")


def _generate_plots(tc_results):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = [f"{r['ligand']}->\n{r['receptor']}" for r in tc_results]
    tcs = [r["tanimoto"] for r in tc_results]
    colors = ["#4CAF50" if r["is_self"] else "#2196F3" for r in tc_results]
    edge_colors = ["green" if r["status"] == "PASS" else "red"
                   for r in tc_results]

    fig, ax = plt.subplots(figsize=(14, 7))
    bars = ax.bar(range(len(labels)), tcs, color=colors,
                  edgecolor=edge_colors, linewidth=2)

    ax.axhline(y=0.85, color="red", linestyle="--", linewidth=2,
               label="0.85 threshold")

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel("Tanimoto Coefficient (IFP)", fontsize=12)
    ax.set_title("Cross-Docking IFP Matrix: Rigid Ensemble Validation\n"
                 "Green = Self-dock, Blue = Cross-dock, "
                 "Red border = FAIL",
                 fontsize=13, fontweight="bold")
    ax.set_ylim(0, 1.1)
    ax.grid(axis="y", alpha=0.3)

    for i, v in enumerate(tcs):
        ax.text(i, v + 0.02, f"{v:.3f}", ha="center", fontsize=10,
                fontweight="bold")

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#4CAF50", label="Self-dock"),
        Patch(facecolor="#2196F3", label="Cross-dock"),
        plt.Line2D([0], [0], color="red", linestyle="--",
                    label="0.85 threshold"),
    ]
    ax.legend(handles=legend_elements, fontsize=11, loc="lower right")

    plt.tight_layout()
    save_path = os.path.join(OUTPUT_DIR, "cross_docking_ifp_matrix.png")
    plt.savefig(save_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"IFP matrix plot saved to: {save_path}")

    if len(tc_results) >= 4:
        fig2, axes2 = plt.subplots(1, 2, figsize=(14, 6))
        fig2.suptitle("Cross-Docking Scores: CNN vs Vina",
                      fontsize=13, fontweight="bold")

        cnn_scores = [r["cnn_score"] for r in tc_results
                      if r["cnn_score"] is not None]
        vina_scores = [r["vina_score"] for r in tc_results
                       if r["vina_score"] is not None]
        tc_vals = [r["tanimoto"] for r in tc_results
                   if r["cnn_score"] is not None]
        plot_colors = ["#4CAF50" if r["is_self"] else "#2196F3"
                       for r in tc_results if r["cnn_score"] is not None]
        plot_labels = [f"{r['ligand']}->{r['receptor']}" for r in tc_results
                       if r["cnn_score"] is not None]

        if cnn_scores and tc_vals:
            ax1 = axes2[0]
            ax1.scatter(cnn_scores, tc_vals[:len(cnn_scores)],
                        c=plot_colors[:len(cnn_scores)], s=100,
                        edgecolors="black", zorder=5)
            for i, label in enumerate(plot_labels[:len(cnn_scores)]):
                ax1.annotate(label, (cnn_scores[i], tc_vals[i]),
                             fontsize=8, ha="center",
                             textcoords="offset points", xytext=(0, 10))
            ax1.axhline(y=0.85, color="red", linestyle="--", alpha=0.7)
            ax1.set_xlabel("CNN Score", fontsize=11)
            ax1.set_ylabel("IFP Tanimoto", fontsize=11)
            ax1.set_title("CNN Score vs IFP Tanimoto", fontsize=11)
            ax1.grid(alpha=0.3)

        if vina_scores and tc_vals:
            ax2 = axes2[1]
            ax2.scatter(vina_scores[:len(tc_vals)], tc_vals[:len(vina_scores)],
                        c=plot_colors[:len(vina_scores)], s=100,
                        edgecolors="black", zorder=5)
            for i, label in enumerate(plot_labels[:len(vina_scores)]):
                ax2.annotate(label, (vina_scores[i], tc_vals[i]),
                             fontsize=8, ha="center",
                             textcoords="offset points", xytext=(0, 10))
            ax2.axhline(y=0.85, color="red", linestyle="--", alpha=0.7)
            ax2.set_xlabel("Vina Score (kcal/mol)", fontsize=11)
            ax2.set_ylabel("IFP Tanimoto", fontsize=11)
            ax2.set_title("Vina Score vs IFP Tanimoto", fontsize=11)
            ax2.grid(alpha=0.3)

        plt.tight_layout()
        save_path2 = os.path.join(OUTPUT_DIR,
                                  "cross_docking_scores_vs_ifp.png")
        plt.savefig(save_path2, dpi=200, bbox_inches="tight",
                    facecolor="white")
        plt.close()
        print(f"Scores vs IFP plot saved to: {save_path2}")


if __name__ == "__main__":
    main()
