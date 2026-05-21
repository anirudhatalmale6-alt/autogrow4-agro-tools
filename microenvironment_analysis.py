#!/usr/bin/env python3
"""
Microenvironment Analysis - Task 2: Sidechain Centroid Integration

For each of the 14 key GID1 binding pocket residues (+/- 3 flanking):
  - Void volume in 4.5 A sphere around sidechain centroid
  - Integrated electrostatic potential on 3.5-5.0 A radial shell
  - Solvent-accessible surface area (SASA)
  - Hydrophobicity (Kyte-Doolittle)

Compares within-species (induced fit) and cross-species (feasibility).
Flags electrostatic inversions that would inhibit universal ligand design.

Structures:
  2ZSH: Arabidopsis GID1A + GA3
  2ZSI: Arabidopsis GID1A + GA4
  3ED1: Rice GID1   + GA3 (hexamer, chain A)
  3EBL: Rice GID1   + GA4 (hexamer, chain A)
"""

import os
import sys
import copy
import urllib.request
import numpy as np
from collections import OrderedDict

from Bio.PDB import PDBParser, Superimposer, PDBIO, Select

OUTPUT_DIR = os.path.expanduser("~/final_project/pocket_analysis")
os.makedirs(OUTPUT_DIR, exist_ok=True)

PDB_IDS = ["2ZSH", "2ZSI", "3ED1", "3EBL"]
SPECIES = {"2ZSH": "Arabidopsis", "2ZSI": "Arabidopsis",
           "3ED1": "Rice", "3EBL": "Rice"}
LIGANDS = {"2ZSH": "GA3", "2ZSI": "GA4",
           "3ED1": "GA3", "3EBL": "GA4"}

ARAB_RESIDUES = [24, 27, 31, 116, 119, 126, 127, 191, 218, 238, 239, 244, 319, 323]
RICE_RESIDUES = [24, 27, 31, 123, 126, 133, 134, 198, 225, 245, 246, 251, 326, 330]

BACKBONE_ATOMS = ["N", "CA", "C", "O"]
BACKBONE_SET = {"N", "CA", "C", "O", "H", "HA", "HA2", "HA3", "OXT"}

VDW_RADII = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80,
    "P": 1.80, "SE": 1.90, "FE": 1.94, "ZN": 1.39, "MG": 1.73,
}

HYDROPHOBICITY = {
    "ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5,
    "MET": 1.9, "ALA": 1.8, "GLY": -0.4, "THR": -0.7, "SER": -0.8,
    "TRP": -0.9, "TYR": -1.3, "PRO": -1.6, "HIS": -3.2, "GLU": -3.5,
    "GLN": -3.5, "ASP": -3.5, "ASN": -3.5, "LYS": -3.9, "ARG": -4.5,
}

PARTIAL_CHARGES = {
    ("ARG", "NH1"): 0.34, ("ARG", "NH2"): 0.34, ("ARG", "NE"): 0.16,
    ("ARG", "CZ"): 0.16,
    ("LYS", "NZ"): 1.0,
    ("ASP", "OD1"): -0.5, ("ASP", "OD2"): -0.5,
    ("GLU", "OE1"): -0.5, ("GLU", "OE2"): -0.5,
    ("SER", "OG"): -0.25, ("THR", "OG1"): -0.25,
    ("TYR", "OH"): -0.20,
    ("ASN", "OD1"): -0.30, ("ASN", "ND2"): 0.30,
    ("GLN", "OE1"): -0.30, ("GLN", "NE2"): 0.30,
    ("HIS", "ND1"): -0.10, ("HIS", "NE2"): 0.10,
    ("CYS", "SG"): -0.10,
    ("TRP", "NE1"): 0.10,
}

BACKBONE_CHARGES = {"N": 0.15, "O": -0.40, "C": 0.25}

FLANKING = 3
VOID_RADIUS = 4.5
VOID_GRID_SPACING = 0.3
SHELL_INNER = 3.5
SHELL_OUTER = 5.0
SHELL_N_POINTS = 500
SHELL_N_LAYERS = 4
ELECTROSTATIC_CUTOFF = 15.0


def download_pdb(pdb_id, output_dir):
    path = os.path.join(output_dir, f"{pdb_id}.pdb")
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


def get_expanded_residues(pocket_residues):
    expanded = set()
    for r in pocket_residues:
        for offset in range(-FLANKING, FLANKING + 1):
            expanded.add(r + offset)
    return sorted(expanded)


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


def save_chain_a_pdb(struct, pdb_id, output_dir):
    out_path = os.path.join(output_dir, f"{pdb_id}_chainA_aligned.pdb")
    io = PDBIO()
    io.set_structure(struct)
    io.save(out_path, ChainAProteinSelect())
    return out_path


def get_atom_charge(resname, atomname):
    key = (resname, atomname)
    if key in PARTIAL_CHARGES:
        return PARTIAL_CHARGES[key]
    if atomname in BACKBONE_CHARGES:
        return BACKBONE_CHARGES[atomname]
    return 0.0


def get_atom_radius(atom):
    elem = atom.element.strip().upper()
    if not elem:
        name = atom.name.strip()
        elem = name[0] if name else "C"
    return VDW_RADII.get(elem, 1.70)


def sidechain_centroid(residue):
    coords = []
    for atom in residue:
        if atom.element.strip() != "H" and atom.name not in BACKBONE_SET:
            coords.append(atom.get_vector().get_array())
    if not coords:
        if "CA" in residue:
            return residue["CA"].get_vector().get_array()
        return None
    return np.mean(coords, axis=0)


def fibonacci_sphere(n_points, radius, center):
    golden = (1 + np.sqrt(5)) / 2
    indices = np.arange(n_points)
    theta = np.arccos(1 - 2 * (indices + 0.5) / n_points)
    phi = 2 * np.pi * indices / golden
    x = center[0] + radius * np.sin(theta) * np.cos(phi)
    y = center[1] + radius * np.sin(theta) * np.sin(phi)
    z = center[2] + radius * np.cos(theta)
    return np.column_stack([x, y, z])


def calc_void_volume(centroid, atom_coords, atom_radii):
    r = VOID_RADIUS
    n = int(2 * r / VOID_GRID_SPACING) + 1
    lin = np.linspace(-r, r, n)
    gx, gy, gz = np.meshgrid(lin, lin, lin, indexing='ij')
    grid = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])

    dist_from_center = np.linalg.norm(grid, axis=1)
    in_sphere = grid[dist_from_center <= r]
    if len(in_sphere) == 0:
        return 0.0

    in_sphere_abs = in_sphere + centroid
    sphere_vol = (4.0 / 3.0) * np.pi * r ** 3

    nearby_mask = np.linalg.norm(atom_coords - centroid, axis=1) < (r + 3.0)
    nearby_coords = atom_coords[nearby_mask]
    nearby_radii = atom_radii[nearby_mask]

    if len(nearby_coords) == 0:
        return sphere_vol

    occupied = np.zeros(len(in_sphere_abs), dtype=bool)
    for ac, ar in zip(nearby_coords, nearby_radii):
        dist = np.linalg.norm(in_sphere_abs - ac, axis=1)
        occupied |= (dist <= ar)

    void_frac = 1.0 - np.sum(occupied) / len(in_sphere_abs)
    return void_frac * sphere_vol


def calc_electrostatic_shell(centroid, atom_coords, atom_charges):
    nearby_mask = np.linalg.norm(atom_coords - centroid, axis=1) < ELECTROSTATIC_CUTOFF
    near_coords = atom_coords[nearby_mask]
    near_charges = atom_charges[nearby_mask]

    charged_mask = np.abs(near_charges) > 1e-6
    near_coords = near_coords[charged_mask]
    near_charges = near_charges[charged_mask]

    if len(near_charges) == 0:
        return 0.0

    total = 0.0
    for shell_r in np.linspace(SHELL_INNER, SHELL_OUTER, SHELL_N_LAYERS):
        pts = fibonacci_sphere(SHELL_N_POINTS, shell_r, centroid)
        for pt in pts:
            dists = np.linalg.norm(near_coords - pt, axis=1)
            dists = np.maximum(dists, 0.5)
            total += np.sum(near_charges / dists)

    n_total = SHELL_N_POINTS * SHELL_N_LAYERS
    return total / n_total


def calc_sasa_residue(structure, chain_id, resnum):
    try:
        from Bio.PDB.SASA import ShrakeRupley
        sr = ShrakeRupley()
        sr.compute(structure[0], level="R")
        chain = structure[0][chain_id]
        res = chain[(' ', resnum, ' ')]
        return res.sasa
    except Exception:
        return None


def collect_protein_atoms(chain):
    coords, charges, radii = [], [], []
    for res in chain:
        if res.id[0] != " ":
            continue
        resname = res.get_resname()
        for atom in res:
            if atom.element.strip() == "H":
                continue
            coords.append(atom.get_vector().get_array())
            charges.append(get_atom_charge(resname, atom.name))
            radii.append(get_atom_radius(atom))
    return np.array(coords), np.array(charges), np.array(radii)


def analyze_structure(pdb_id, chain, structure):
    pocket_res = get_pocket_residues(pdb_id)
    expanded = get_expanded_residues(pocket_res)

    atom_coords, atom_charges, atom_radii = collect_protein_atoms(chain)

    results = OrderedDict()
    for resnum in expanded:
        try:
            res = chain[(' ', resnum, ' ')]
        except KeyError:
            continue

        centroid = sidechain_centroid(res)
        if centroid is None:
            continue

        void_vol = calc_void_volume(centroid, atom_coords, atom_radii)
        elec_pot = calc_electrostatic_shell(centroid, atom_coords, atom_charges)
        sasa = calc_sasa_residue(structure, "A", resnum)
        hydro = HYDROPHOBICITY.get(res.get_resname(), 0.0)

        is_key = resnum in pocket_res
        if is_key:
            arab_num = ARAB_RESIDUES[pocket_res.index(resnum)]
        else:
            arab_num = None

        results[resnum] = {
            "resnum": resnum,
            "resname": res.get_resname(),
            "is_key": is_key,
            "arab_num": arab_num,
            "centroid": centroid,
            "void_volume": void_vol,
            "electrostatic": elec_pot,
            "sasa": sasa,
            "hydrophobicity": hydro,
        }

    return results


def compute_pairwise(res1_data, res2_data, pocket1, pocket2):
    diffs = []
    for i, (r1num, r2num, arab_num) in enumerate(
            zip(pocket1, pocket2, ARAB_RESIDUES)):
        if r1num not in res1_data or r2num not in res2_data:
            continue

        d1 = res1_data[r1num]
        d2 = res2_data[r2num]

        dv = d2["void_volume"] - d1["void_volume"]
        de = d2["electrostatic"] - d1["electrostatic"]
        ds = (d2["sasa"] - d1["sasa"]) if (
            d1["sasa"] is not None and d2["sasa"] is not None) else None

        e1 = d1["electrostatic"]
        e2 = d2["electrostatic"]
        inversion = False
        if abs(e1) > 0.01 and abs(e2) > 0.01:
            inversion = (e1 > 0) != (e2 > 0)

        diffs.append({
            "arab_resnum": arab_num,
            "resname": d1["resname"],
            "res1_num": r1num,
            "res2_num": r2num,
            "delta_void": dv,
            "delta_elec": de,
            "delta_sasa": ds,
            "elec1": e1,
            "elec2": e2,
            "void1": d1["void_volume"],
            "void2": d2["void_volume"],
            "inversion": inversion,
        })

    return diffs


def compute_flanking_pairwise(res1_data, res2_data, pocket1, pocket2):
    all_diffs = OrderedDict()

    for i, (key1, key2, arab_num) in enumerate(
            zip(pocket1, pocket2, ARAB_RESIDUES)):
        offset = key2 - key1
        for flank in range(-FLANKING, FLANKING + 1):
            r1 = key1 + flank
            r2 = key2 + flank
            if r1 not in res1_data or r2 not in res2_data:
                continue
            d1 = res1_data[r1]
            d2 = res2_data[r2]

            label = f"{d1['resname']}{arab_num}{'+' if flank > 0 else ''}{flank}" if flank != 0 else f"{d1['resname']}{arab_num}"

            e1 = d1["electrostatic"]
            e2 = d2["electrostatic"]
            inversion = False
            if abs(e1) > 0.01 and abs(e2) > 0.01:
                inversion = (e1 > 0) != (e2 > 0)

            all_diffs[(arab_num, flank)] = {
                "label": label,
                "arab_key": arab_num,
                "flank_offset": flank,
                "is_key": flank == 0,
                "resname1": d1["resname"],
                "resname2": d2["resname"],
                "delta_void": d2["void_volume"] - d1["void_volume"],
                "delta_elec": d2["electrostatic"] - d1["electrostatic"],
                "elec1": e1, "elec2": e2,
                "void1": d1["void_volume"], "void2": d2["void_volume"],
                "inversion": inversion,
            }

    return all_diffs


def main():
    print("=" * 70)
    print("MICROENVIRONMENT ANALYSIS - Task 2")
    print("Sidechain Centroid Integration (4.5 A Shell)")
    print("=" * 70)
    print(f"Output: {OUTPUT_DIR}")
    print(f"Void sphere: {VOID_RADIUS} A, grid {VOID_GRID_SPACING} A")
    print(f"Electrostatic shell: {SHELL_INNER}-{SHELL_OUTER} A, "
          f"{SHELL_N_POINTS} pts x {SHELL_N_LAYERS} layers")
    print()

    parser = PDBParser(QUIET=True)
    structures = {}
    chains = {}

    print("Loading and aligning structures...")
    for pdb_id in PDB_IDS:
        pdb_path = download_pdb(pdb_id, OUTPUT_DIR)
        struct = parser.get_structure(pdb_id, pdb_path)
        structures[pdb_id] = struct
        chains[pdb_id] = struct[0]['A']
        n_chains = len(list(struct[0].get_chains()))
        print(f"  {pdb_id}: {SPECIES[pdb_id]} + {LIGANDS[pdb_id]}, "
              f"{n_chains} chain(s)")

    ref_id = "2ZSH"
    ref_chain = chains[ref_id]
    ref_pocket = get_pocket_residues(ref_id)
    ref_bb = get_pocket_backbone_atoms(ref_chain, ref_pocket)

    print(f"\nAligning all structures to {ref_id} on pocket backbone...")
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

    print("\nComputing microenvironment for each structure...")
    all_data = {}
    for pdb_id in PDB_IDS:
        print(f"\n  {pdb_id} ({SPECIES[pdb_id]} + {LIGANDS[pdb_id]}):")
        pocket = get_pocket_residues(pdb_id)
        expanded = get_expanded_residues(pocket)
        print(f"    Key residues: {pocket}")
        print(f"    Expanded (+/-{FLANKING}): {len(expanded)} residues")

        data = analyze_structure(pdb_id, chains[pdb_id], structures[pdb_id])
        all_data[pdb_id] = data

        key_data = {k: v for k, v in data.items() if v["is_key"]}
        vols = [v["void_volume"] for v in key_data.values()]
        elecs = [v["electrostatic"] for v in key_data.values()]
        print(f"    Total pocket void volume: {sum(vols):.1f} A^3")
        print(f"    Mean electrostatic potential: {np.mean(elecs):.4f}")
        print(f"    Void volume range: {min(vols):.1f} - {max(vols):.1f} A^3")

    pairs = [
        ("2ZSH", "2ZSI", "Same species, different ligand (induced fit)"),
        ("3ED1", "3EBL", "Same species, different ligand (induced fit)"),
        ("2ZSH", "3ED1", "Cross-species, same ligand (GA3)"),
        ("2ZSH", "3EBL", "Cross-species, different ligand"),
        ("2ZSI", "3ED1", "Cross-species, different ligand"),
        ("2ZSI", "3EBL", "Cross-species, same ligand (GA4)"),
    ]

    print(f"\n{'=' * 70}")
    print("PAIRWISE MICROENVIRONMENT COMPARISON")
    print(f"{'=' * 70}")

    all_pair_results = []
    all_inversions = []

    for pdb1, pdb2, desc in pairs:
        pocket1 = get_pocket_residues(pdb1)
        pocket2 = get_pocket_residues(pdb2)

        diffs = compute_pairwise(
            all_data[pdb1], all_data[pdb2], pocket1, pocket2)
        flank_diffs = compute_flanking_pairwise(
            all_data[pdb1], all_data[pdb2], pocket1, pocket2)

        dvs = [d["delta_void"] for d in diffs]
        des = [d["delta_elec"] for d in diffs]
        inversions = [d for d in diffs if d["inversion"]]

        is_cross = SPECIES[pdb1] != SPECIES[pdb2]

        result = {
            "pdb1": pdb1, "pdb2": pdb2,
            "description": desc,
            "is_cross_species": is_cross,
            "key_diffs": diffs,
            "flank_diffs": flank_diffs,
            "mean_delta_void": np.mean(np.abs(dvs)),
            "max_delta_void": np.max(np.abs(dvs)),
            "mean_delta_elec": np.mean(np.abs(des)),
            "max_delta_elec": np.max(np.abs(des)),
            "inversions": inversions,
        }
        all_pair_results.append(result)

        print(f"\n  {pdb1} vs {pdb2}: {desc}")
        print(f"    Mean |DeltaV|: {result['mean_delta_void']:.2f} A^3")
        print(f"    Max  |DeltaV|: {result['max_delta_void']:.2f} A^3")
        print(f"    Mean |DeltaE|: {result['mean_delta_elec']:.4f}")
        print(f"    Max  |DeltaE|: {result['max_delta_elec']:.4f}")

        if inversions:
            print(f"    ELECTROSTATIC INVERSIONS: {len(inversions)}")
            for inv in inversions:
                print(f"      {inv['resname']}{inv['arab_resnum']}: "
                      f"{inv['elec1']:+.4f} -> {inv['elec2']:+.4f}")
                all_inversions.append((pdb1, pdb2, inv))
        else:
            print(f"    Electrostatic inversions: none")

        print(f"\n    {'Res':>6} {'Name':>4} {'V1':>8} {'V2':>8} "
              f"{'DV':>8} {'E1':>8} {'E2':>8} {'DE':>8} {'Flag':>6}")
        print(f"    {'─' * 70}")
        for d in diffs:
            flag = "INVERT" if d["inversion"] else ""
            print(f"    {d['arab_resnum']:>6} {d['resname']:>4} "
                  f"{d['void1']:>8.1f} {d['void2']:>8.1f} "
                  f"{d['delta_void']:>+8.1f} "
                  f"{d['elec1']:>8.4f} {d['elec2']:>8.4f} "
                  f"{d['delta_elec']:>+8.4f} {flag:>6}")

    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print(f"{'=' * 70}\n")

    cross = [r for r in all_pair_results if r["is_cross_species"]]
    within = [r for r in all_pair_results if not r["is_cross_species"]]

    if within:
        print("  Within-species (ligand-induced changes):")
        for r in within:
            print(f"    {r['pdb1']} vs {r['pdb2']}: "
                  f"mean|DV|={r['mean_delta_void']:.2f} A^3, "
                  f"mean|DE|={r['mean_delta_elec']:.4f}")

    if cross:
        print("\n  Cross-species:")
        for r in cross:
            inv_count = len(r["inversions"])
            inv_str = f", {inv_count} inversions" if inv_count > 0 else ""
            print(f"    {r['pdb1']} vs {r['pdb2']}: "
                  f"mean|DV|={r['mean_delta_void']:.2f} A^3, "
                  f"mean|DE|={r['mean_delta_elec']:.4f}{inv_str}")

    cross_inversions = [(r["pdb1"], r["pdb2"], inv)
                         for r in cross for inv in r["inversions"]]
    print(f"\n  Total cross-species electrostatic inversions: "
          f"{len(cross_inversions)}")

    if cross_inversions:
        print("  WARNING: Electrostatic inversions detected!")
        print("  These residues have opposite charge character between species")
        print("  and may require special attention in universal ligand design:")
        for pdb1, pdb2, inv in cross_inversions:
            print(f"    {inv['resname']}{inv['arab_resnum']} "
                  f"({pdb1}: {inv['elec1']:+.4f}, "
                  f"{pdb2}: {inv['elec2']:+.4f})")
    else:
        print("  No electrostatic inversions in cross-species comparisons.")
        print("  Electrostatic landscape is conserved across species.")

    if cross:
        mean_cross_dv = np.mean([r["mean_delta_void"] for r in cross])
        mean_within_dv = np.mean([r["mean_delta_void"] for r in within]) if within else 0
        print(f"\n  Cross-species mean |DeltaV|: {mean_cross_dv:.2f} A^3")
        print(f"  Within-species mean |DeltaV|: {mean_within_dv:.2f} A^3")
        ratio = mean_cross_dv / mean_within_dv if mean_within_dv > 0 else float('inf')
        print(f"  Cross/Within ratio: {ratio:.2f}x")
        if ratio < 2.0:
            print("  Cross-species steric changes are comparable to "
                  "ligand-induced changes.")
            print("  This supports universal ligand feasibility.")
        else:
            print("  Cross-species steric changes are significantly larger "
                  "than ligand-induced changes.")

    print(f"\n{'=' * 70}")
    if cross and len(cross_inversions) == 0:
        print("  CONCLUSION: Microenvironment is CONSERVED across species.")
        print("  No electrostatic inversions. Steric differences are modest.")
        print("  Universal ligand design is supported by Task 2 analysis.")
    elif cross and len(cross_inversions) <= 2:
        print("  CONCLUSION: Microenvironment is MOSTLY CONSERVED.")
        print(f"  {len(cross_inversions)} electrostatic inversion(s) detected -")
        print("  these residues need attention but do not necessarily")
        print("  preclude universal ligand design.")
    else:
        print("  CONCLUSION: Significant microenvironment differences detected.")
        print("  Universal ligand design faces electrostatic challenges.")
    print(f"{'=' * 70}")

    report_path = os.path.join(OUTPUT_DIR, "microenvironment_report.txt")
    _write_report(report_path, all_data, all_pair_results, cross_inversions)
    print(f"\nReport saved to: {report_path}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        _generate_plots(all_pair_results, all_data)
    except ImportError:
        print("matplotlib not available, skipping plots")


def _write_report(path, all_data, pair_results, cross_inversions):
    with open(path, "w") as f:
        f.write("MICROENVIRONMENT ANALYSIS - Task 2\n")
        f.write("Sidechain Centroid Integration (4.5 A Shell)\n")
        f.write(f"{'=' * 60}\n\n")

        f.write("Parameters:\n")
        f.write(f"  Void volume sphere: {VOID_RADIUS} A radius, "
                f"{VOID_GRID_SPACING} A grid\n")
        f.write(f"  Electrostatic shell: {SHELL_INNER}-{SHELL_OUTER} A, "
                f"{SHELL_N_POINTS} sample points x {SHELL_N_LAYERS} layers\n")
        f.write(f"  Flanking residues: +/- {FLANKING}\n\n")

        for pdb_id in PDB_IDS:
            data = all_data[pdb_id]
            key_data = {k: v for k, v in data.items() if v["is_key"]}
            vols = [v["void_volume"] for v in key_data.values()]
            f.write(f"{pdb_id} ({SPECIES[pdb_id]} + {LIGANDS[pdb_id]}):\n")
            f.write(f"  Total pocket void volume: {sum(vols):.1f} A^3\n")
            f.write(f"  Per-residue:\n")
            for resnum, d in key_data.items():
                f.write(f"    {d['resname']}{resnum}: "
                        f"V={d['void_volume']:.1f} A^3, "
                        f"E={d['electrostatic']:.4f}, "
                        f"Hydro={d['hydrophobicity']:.1f}")
                if d['sasa'] is not None:
                    f.write(f", SASA={d['sasa']:.1f} A^2")
                f.write("\n")
            f.write("\n")

        for r in pair_results:
            f.write(f"\n{r['pdb1']} vs {r['pdb2']}: {r['description']}\n")
            f.write(f"  Mean |DeltaV|: {r['mean_delta_void']:.2f} A^3\n")
            f.write(f"  Max  |DeltaV|: {r['max_delta_void']:.2f} A^3\n")
            f.write(f"  Mean |DeltaE|: {r['mean_delta_elec']:.4f}\n")
            f.write(f"  Max  |DeltaE|: {r['max_delta_elec']:.4f}\n")
            f.write(f"  Inversions: {len(r['inversions'])}\n\n")

            for d in r['key_diffs']:
                flag = " INVERSION" if d['inversion'] else ""
                f.write(f"    {d['resname']}{d['arab_resnum']}: "
                        f"DV={d['delta_void']:+.1f} A^3, "
                        f"DE={d['delta_elec']:+.4f}{flag}\n")

        f.write(f"\n{'=' * 60}\n")
        f.write(f"Cross-species electrostatic inversions: "
                f"{len(cross_inversions)}\n")
        for pdb1, pdb2, inv in cross_inversions:
            f.write(f"  {inv['resname']}{inv['arab_resnum']} "
                    f"({pdb1}: {inv['elec1']:+.4f}, "
                    f"{pdb2}: {inv['elec2']:+.4f})\n")


def _generate_plots(pair_results, all_data):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=(22, 14))
    fig.suptitle(
        "Microenvironment Analysis: "
        "Delta Void Volume & Delta Electrostatic Potential\n"
        "14 Key GID1 Binding Pocket Residues",
        fontsize=15, fontweight="bold", y=0.99)

    for idx, result in enumerate(pair_results):
        ax = axes[idx // 3][idx % 3]
        diffs = result["key_diffs"]

        labels = [f"{d['resname']}\n{d['arab_resnum']}" for d in diffs]
        dvs = [d["delta_void"] for d in diffs]
        des_scaled = [d["delta_elec"] * 50 for d in diffs]

        x = np.arange(len(labels))
        w = 0.35

        bars1 = ax.bar(x - w / 2, dvs, w, label="Delta Void Vol (A^3)",
                        color="#2196F3", edgecolor="black", linewidth=0.5,
                        alpha=0.8)
        bars2 = ax.bar(x + w / 2, des_scaled, w,
                        label="Delta Elec. Pot. (x50)",
                        color="#FF9800", edgecolor="black", linewidth=0.5,
                        alpha=0.8)

        for i, d in enumerate(diffs):
            if d["inversion"]:
                ax.plot(x[i], 0, marker='v', color='red', markersize=12,
                        zorder=5)
                ax.annotate("INV", (x[i], 0), textcoords="offset points",
                            xytext=(0, -18), ha='center', fontsize=7,
                            color='red', fontweight='bold')

        ax.axhline(y=0, color="black", linewidth=1)
        ax.set_xlabel("Residue", fontsize=9)
        ax.set_ylabel("Delta Value", fontsize=9)

        cross_tag = "CROSS-SPECIES" if result["is_cross_species"] else "WITHIN-SPECIES"
        inv_tag = f" | {len(result['inversions'])} inversions" if result['inversions'] else ""
        ax.set_title(f"{result['pdb1']} vs {result['pdb2']}\n"
                     f"{result['description']}\n"
                     f"[{cross_tag}{inv_tag}]",
                     fontsize=9, fontweight="bold")

        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=7)
        ax.legend(fontsize=7, loc="best")
        ax.grid(axis="y", alpha=0.3)

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    save_path = os.path.join(OUTPUT_DIR, "microenvironment_deltav_deltae.png")
    plt.savefig(save_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"DeltaV/DeltaE plot saved to: {save_path}")

    fig2, axes2 = plt.subplots(2, 3, figsize=(24, 14))
    fig2.suptitle(
        "Microenvironment Landscape: Key Residues +/- 3 Flanking\n"
        "Delta Void Volume (blue) and Delta Electrostatic Potential (orange)",
        fontsize=15, fontweight="bold", y=0.99)

    for idx, result in enumerate(pair_results):
        ax = axes2[idx // 3][idx % 3]
        flank_diffs = result["flank_diffs"]

        if not flank_diffs:
            continue

        sorted_keys = sorted(flank_diffs.keys())
        labels = []
        dvs = []
        des_scaled = []
        is_key_flags = []

        for key in sorted_keys:
            d = flank_diffs[key]
            labels.append(d["label"])
            dvs.append(d["delta_void"])
            des_scaled.append(d["delta_elec"] * 50)
            is_key_flags.append(d["is_key"])

        x = np.arange(len(labels))
        colors_v = ["#1565C0" if ik else "#90CAF9" for ik in is_key_flags]
        colors_e = ["#E65100" if ik else "#FFB74D" for ik in is_key_flags]

        ax.bar(x - 0.2, dvs, 0.35, color=colors_v,
               edgecolor="black", linewidth=0.3, alpha=0.8)
        ax.bar(x + 0.2, des_scaled, 0.35, color=colors_e,
               edgecolor="black", linewidth=0.3, alpha=0.8)

        ax.axhline(y=0, color="black", linewidth=1)
        ax.set_title(f"{result['pdb1']} vs {result['pdb2']}\n"
                     f"{result['description']}", fontsize=9,
                     fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=90, fontsize=5)
        ax.grid(axis="y", alpha=0.3)

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    save_path2 = os.path.join(OUTPUT_DIR, "microenvironment_landscape.png")
    plt.savefig(save_path2, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Landscape plot saved to: {save_path2}")

    fig3, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
    fig3.suptitle("Microenvironment Summary: All Pairwise Comparisons",
                  fontsize=14, fontweight="bold")

    labels = [f"{r['pdb1']} vs\n{r['pdb2']}" for r in pair_results]
    mean_dvs = [r["mean_delta_void"] for r in pair_results]
    mean_des = [r["mean_delta_elec"] for r in pair_results]
    colors = ["#4CAF50" if not r["is_cross_species"] else "#2196F3"
              for r in pair_results]

    bars1 = ax1.bar(range(len(labels)), mean_dvs, color=colors,
                    edgecolor="black", linewidth=0.5)
    ax1.set_xticks(range(len(labels)))
    ax1.set_xticklabels(labels, fontsize=9)
    ax1.set_ylabel("Mean |Delta Void Volume| (A^3)", fontsize=11)
    ax1.set_title("Steric Changes", fontsize=12, fontweight="bold")
    ax1.grid(axis="y", alpha=0.3)
    for i, v in enumerate(mean_dvs):
        ax1.text(i, v + 0.1, f"{v:.1f}", ha="center", fontsize=9,
                 fontweight="bold")

    bars2 = ax2.bar(range(len(labels)), mean_des, color=colors,
                    edgecolor="black", linewidth=0.5)
    ax2.set_xticks(range(len(labels)))
    ax2.set_xticklabels(labels, fontsize=9)
    ax2.set_ylabel("Mean |Delta Electrostatic Potential|", fontsize=11)
    ax2.set_title("Electrostatic Changes", fontsize=12, fontweight="bold")
    ax2.grid(axis="y", alpha=0.3)
    for i, v in enumerate(mean_des):
        ax2.text(i, v + 0.001, f"{v:.4f}", ha="center", fontsize=9,
                 fontweight="bold")

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#4CAF50", label="Within species"),
        Patch(facecolor="#2196F3", label="Cross species"),
    ]
    ax1.legend(handles=legend_elements, fontsize=9)
    ax2.legend(handles=legend_elements, fontsize=9)

    plt.tight_layout()
    save_path3 = os.path.join(OUTPUT_DIR, "microenvironment_summary.png")
    plt.savefig(save_path3, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Summary plot saved to: {save_path3}")

    fig4, axes4 = plt.subplots(2, 2, figsize=(16, 12))
    fig4.suptitle("Absolute Microenvironment: Per-Structure Values\n"
                  "14 Key Binding Pocket Residues",
                  fontsize=14, fontweight="bold")

    for i, pdb_id in enumerate(PDB_IDS):
        ax = axes4[i // 2][i % 2]
        data = all_data[pdb_id]
        pocket = get_pocket_residues(pdb_id)

        labels_abs = []
        vols = []
        elecs = []
        for resnum in pocket:
            if resnum in data:
                d = data[resnum]
                labels_abs.append(f"{d['resname']}\n{d['arab_num']}")
                vols.append(d["void_volume"])
                elecs.append(d["electrostatic"])

        x = np.arange(len(labels_abs))
        ax_twin = ax.twinx()
        ax.bar(x - 0.2, vols, 0.35, color="#2196F3", alpha=0.7,
               label="Void Volume (A^3)")
        ax_twin.bar(x + 0.2, elecs, 0.35, color="#FF9800", alpha=0.7,
                    label="Electrostatic")

        ax.set_xlabel("Residue", fontsize=9)
        ax.set_ylabel("Void Volume (A^3)", fontsize=10, color="#2196F3")
        ax_twin.set_ylabel("Electrostatic Potential", fontsize=10,
                           color="#FF9800")
        ax.set_title(f"{pdb_id}: {SPECIES[pdb_id]} + {LIGANDS[pdb_id]}",
                     fontsize=11, fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(labels_abs, fontsize=7)
        ax.grid(axis="y", alpha=0.3)

        lines1, labels1 = ax.get_legend_handles_labels()
        lines2, labels2 = ax_twin.get_legend_handles_labels()
        ax.legend(lines1 + lines2, labels1 + labels2, fontsize=7,
                  loc="upper right")

    plt.tight_layout()
    save_path4 = os.path.join(OUTPUT_DIR, "microenvironment_absolute.png")
    plt.savefig(save_path4, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Absolute values plot saved to: {save_path4}")


if __name__ == "__main__":
    main()
