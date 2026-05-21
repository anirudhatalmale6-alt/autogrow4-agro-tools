#!/usr/bin/env python3
"""
Conserved Water Mapping - Task 5

Maps crystallographic water molecules in the GID1 binding pocket,
identifies conserved water positions across all 4 structures, and
analyzes 4 specific water-mediated hydrogen bonding networks:

  1. Ser116--HOH--Ser191  (structural water bridge)
  2. Tyr127--HOH--GA 3-OH (ligand recognition water)
  3. Phe238--HOH--GA3 13-OH (GA3-specific, weak)
  4. Asp289 hydration shell (6 waters, bottom of pocket)

Structures:
  2ZSH: Arabidopsis GID1A + GA3
  2ZSI: Arabidopsis GID1A + GA4
  3ED1: Rice GID1   + GA3 (hexamer, chain A)
  3EBL: Rice GID1   + GA4 (hexamer, chain A)
"""

import os
import sys
import json
import urllib.request
import numpy as np
from collections import defaultdict

from Bio.PDB import PDBParser, Superimposer, NeighborSearch

OUTPUT_DIR = os.path.expanduser("~/final_project/water_mapping")
POCKET_DIR = os.path.expanduser("~/final_project/pocket_analysis")
os.makedirs(OUTPUT_DIR, exist_ok=True)

PDB_IDS = ["2ZSH", "2ZSI", "3ED1", "3EBL"]
SPECIES = {"2ZSH": "Arabidopsis", "2ZSI": "Arabidopsis",
           "3ED1": "Rice", "3EBL": "Rice"}
LIGAND_OF = {"2ZSH": "GA3", "2ZSI": "GA4",
             "3ED1": "GA3", "3EBL": "GA4"}

ARAB_RESIDUES = [24, 27, 31, 116, 119, 126, 127, 191, 218, 238, 239, 244, 319, 323]
RICE_RESIDUES = [24, 27, 31, 123, 126, 133, 134, 198, 225, 245, 246, 251, 326, 330]

BACKBONE_ATOMS = ["N", "CA", "C", "O"]

WATER_BRIDGE_RESIDUES_ARAB = {
    "bridge1_ser116": 116,
    "bridge1_ser191": 191,
    "bridge2_tyr127": 127,
    "bridge3_phe238": 238,
    "bridge4_asp289": 289,
}

ARAB_TO_RICE_OFFSET = {
    24: 24, 27: 27, 31: 31,
    116: 123, 119: 126, 126: 133, 127: 134,
    191: 198, 218: 225, 238: 245, 239: 246,
    244: 251, 289: 296, 319: 326, 323: 330,
}

HBOND_DIST = 3.5
WATER_SEARCH_RADIUS = 4.0
CLUSTER_CUTOFF = 2.0
CONSERVED_MIN = 3


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


def map_residue(arab_resnum, pdb_id):
    if SPECIES[pdb_id] == "Arabidopsis":
        return arab_resnum
    return ARAB_TO_RICE_OFFSET.get(arab_resnum, arab_resnum)


def get_residue(chain, resnum):
    try:
        return chain[(' ', resnum, ' ')]
    except KeyError:
        return None


def get_waters_near_residue(chain, resnum, all_waters, cutoff=WATER_SEARCH_RADIUS):
    res = get_residue(chain, resnum)
    if res is None:
        return []

    res_atoms = list(res.get_atoms())
    nearby = []
    for wat in all_waters:
        wat_o = None
        for a in wat.get_atoms():
            if a.element == "O":
                wat_o = a
                break
        if wat_o is None:
            continue

        min_dist = min(
            np.linalg.norm(wat_o.get_vector().get_array() -
                           a.get_vector().get_array())
            for a in res_atoms
        )
        if min_dist <= cutoff:
            nearby.append({
                "water": wat,
                "oxygen": wat_o,
                "coord": wat_o.get_vector().get_array(),
                "dist_to_residue": min_dist,
                "resid": wat.id[1],
                "bfactor": wat_o.get_bfactor(),
            })

    return sorted(nearby, key=lambda w: w["dist_to_residue"])


def get_waters_near_ligand(chain, lig_name, all_waters, cutoff=WATER_SEARCH_RADIUS):
    lig_atoms = []
    for res in chain:
        if res.get_resname() == lig_name and res.id[0] != " ":
            lig_atoms.extend(list(res.get_atoms()))

    if not lig_atoms:
        return []

    nearby = []
    for wat in all_waters:
        wat_o = None
        for a in wat.get_atoms():
            if a.element == "O":
                wat_o = a
                break
        if wat_o is None:
            continue

        min_dist = min(
            np.linalg.norm(wat_o.get_vector().get_array() -
                           a.get_vector().get_array())
            for a in lig_atoms
        )
        if min_dist <= cutoff:
            nearby.append({
                "water": wat,
                "oxygen": wat_o,
                "coord": wat_o.get_vector().get_array(),
                "dist_to_ligand": min_dist,
                "resid": wat.id[1],
                "bfactor": wat_o.get_bfactor(),
            })

    return sorted(nearby, key=lambda w: w["dist_to_ligand"])


def get_all_chain_a_waters(struct):
    waters = []
    for chain in struct[0]:
        if chain.id == "A":
            for res in chain:
                if res.get_resname() == "HOH":
                    waters.append(res)
    return waters


def get_pocket_waters(chain, pocket_residues, all_waters, cutoff=5.0):
    pocket_atoms = []
    for resnum in pocket_residues:
        res = get_residue(chain, resnum)
        if res:
            pocket_atoms.extend(list(res.get_atoms()))

    nearby = []
    seen = set()
    for wat in all_waters:
        wat_o = None
        for a in wat.get_atoms():
            if a.element == "O":
                wat_o = a
                break
        if wat_o is None:
            continue

        min_dist = min(
            np.linalg.norm(wat_o.get_vector().get_array() -
                           a.get_vector().get_array())
            for a in pocket_atoms
        )
        if min_dist <= cutoff and wat.id[1] not in seen:
            seen.add(wat.id[1])
            nearby.append({
                "coord": wat_o.get_vector().get_array(),
                "dist": min_dist,
                "resid": wat.id[1],
                "bfactor": wat_o.get_bfactor(),
            })

    return nearby


def cluster_waters(all_water_positions, cutoff=CLUSTER_CUTOFF):
    clusters = []
    assigned = set()

    positions = []
    labels = []
    for pdb_id, waters in all_water_positions.items():
        for i, w in enumerate(waters):
            positions.append(w["coord"])
            labels.append((pdb_id, i))

    if not positions:
        return []

    positions = np.array(positions)

    for i in range(len(positions)):
        if i in assigned:
            continue

        cluster_members = [(labels[i][0], positions[i])]
        assigned.add(i)

        for j in range(i + 1, len(positions)):
            if j in assigned:
                continue
            if labels[j][0] == labels[i][0]:
                continue
            dist = np.linalg.norm(positions[i] - positions[j])
            if dist <= cutoff:
                cluster_members.append((labels[j][0], positions[j]))
                assigned.add(j)

        if len(cluster_members) >= 2:
            coords = np.array([m[1] for m in cluster_members])
            structs = set(m[0] for m in cluster_members)
            clusters.append({
                "center": np.mean(coords, axis=0),
                "n_structures": len(structs),
                "structures": sorted(structs),
                "spread": np.max(np.std(coords, axis=0)),
                "members": cluster_members,
            })

    clusters.sort(key=lambda c: -c["n_structures"])
    return clusters


def hierarchical_cluster_waters(all_water_positions, cutoff=CLUSTER_CUTOFF):
    all_coords = []
    all_labels = []
    for pdb_id, waters in all_water_positions.items():
        for w in waters:
            all_coords.append(w["coord"])
            all_labels.append(pdb_id)

    if not all_coords:
        return []

    all_coords = np.array(all_coords)
    n = len(all_coords)
    used = [False] * n
    clusters = []

    for i in range(n):
        if used[i]:
            continue
        members = [(all_labels[i], all_coords[i])]
        used[i] = True
        center = all_coords[i].copy()

        for j in range(i + 1, n):
            if used[j]:
                continue
            dist = np.linalg.norm(center - all_coords[j])
            if dist <= cutoff:
                members.append((all_labels[j], all_coords[j]))
                used[j] = True
                coords = np.array([m[1] for m in members])
                center = np.mean(coords, axis=0)

        structs = set(m[0] for m in members)
        coords = np.array([m[1] for m in members])
        clusters.append({
            "center": np.mean(coords, axis=0),
            "n_structures": len(structs),
            "structures": sorted(structs),
            "spread": float(np.max(np.std(coords, axis=0))) if len(members) > 1 else 0.0,
            "n_waters": len(members),
            "members": members,
        })

    clusters.sort(key=lambda c: (-c["n_structures"], -c["n_waters"]))
    return clusters


def analyze_bridge1(structures, chains):
    print("\n  Bridge 1: Ser116--HOH--Ser191")
    print("  " + "-" * 50)
    results = {}

    for pdb_id in PDB_IDS:
        chain = chains[pdb_id]
        all_waters = get_all_chain_a_waters(structures[pdb_id])
        ser116 = map_residue(116, pdb_id)
        ser191 = map_residue(191, pdb_id)

        waters_116 = get_waters_near_residue(chain, ser116, all_waters, HBOND_DIST)
        waters_191 = get_waters_near_residue(chain, ser191, all_waters, HBOND_DIST)

        ids_116 = set(w["resid"] for w in waters_116)
        ids_191 = set(w["resid"] for w in waters_191)
        bridging = ids_116 & ids_191

        bridging_details = []
        for w in waters_116:
            if w["resid"] in bridging:
                d191 = None
                for w2 in waters_191:
                    if w2["resid"] == w["resid"]:
                        d191 = w2["dist_to_residue"]
                        break
                bridging_details.append({
                    "resid": w["resid"],
                    "coord": w["coord"].tolist(),
                    "dist_ser116": w["dist_to_residue"],
                    "dist_ser191": d191,
                    "bfactor": w["bfactor"],
                })

        species_label = SPECIES[pdb_id]
        print(f"    {pdb_id} ({species_label}): "
              f"waters near Ser{ser116}={len(waters_116)}, "
              f"near Ser{ser191}={len(waters_191)}, "
              f"bridging={len(bridging)}")
        for bd in bridging_details:
            print(f"      HOH {bd['resid']}: d(Ser{ser116})={bd['dist_ser116']:.2f}, "
                  f"d(Ser{ser191})={bd['dist_ser191']:.2f}, "
                  f"B={bd['bfactor']:.1f}")

        results[pdb_id] = {
            "ser_a": ser116,
            "ser_b": ser191,
            "n_bridging": len(bridging),
            "details": bridging_details,
            "conserved": len(bridging) > 0,
        }

    conserved_count = sum(1 for r in results.values() if r["conserved"])
    print(f"    CONSERVED: {conserved_count}/4 structures")
    return results


def analyze_bridge2(structures, chains):
    print("\n  Bridge 2: Tyr127--HOH--GA 3-OH")
    print("  " + "-" * 50)
    results = {}

    for pdb_id in PDB_IDS:
        chain = chains[pdb_id]
        all_waters = get_all_chain_a_waters(structures[pdb_id])
        tyr127 = map_residue(127, pdb_id)
        lig_name = LIGAND_OF[pdb_id]

        waters_tyr = get_waters_near_residue(chain, tyr127, all_waters, HBOND_DIST)
        waters_lig = get_waters_near_ligand(chain, lig_name, all_waters, HBOND_DIST)

        ids_tyr = set(w["resid"] for w in waters_tyr)
        ids_lig = set(w["resid"] for w in waters_lig)
        bridging = ids_tyr & ids_lig

        bridging_details = []
        for w in waters_tyr:
            if w["resid"] in bridging:
                d_lig = None
                for w2 in waters_lig:
                    if w2["resid"] == w["resid"]:
                        d_lig = w2["dist_to_ligand"]
                        break
                bridging_details.append({
                    "resid": w["resid"],
                    "coord": w["coord"].tolist(),
                    "dist_tyr127": w["dist_to_residue"],
                    "dist_ligand": d_lig,
                    "bfactor": w["bfactor"],
                })

        species_label = SPECIES[pdb_id]
        print(f"    {pdb_id} ({species_label}, {lig_name}): "
              f"waters near Tyr{tyr127}={len(waters_tyr)}, "
              f"near {lig_name}={len(waters_lig)}, "
              f"bridging={len(bridging)}")
        for bd in bridging_details:
            print(f"      HOH {bd['resid']}: d(Tyr{tyr127})={bd['dist_tyr127']:.2f}, "
                  f"d({lig_name})={bd['dist_ligand']:.2f}, "
                  f"B={bd['bfactor']:.1f}")

        results[pdb_id] = {
            "tyr": tyr127,
            "ligand": lig_name,
            "n_bridging": len(bridging),
            "details": bridging_details,
            "conserved": len(bridging) > 0,
        }

    conserved_count = sum(1 for r in results.values() if r["conserved"])
    print(f"    CONSERVED: {conserved_count}/4 structures")
    return results


def analyze_bridge3(structures, chains):
    print("\n  Bridge 3: Phe238--HOH--GA3 13-OH (GA3-specific)")
    print("  " + "-" * 50)
    results = {}

    for pdb_id in PDB_IDS:
        chain = chains[pdb_id]
        all_waters = get_all_chain_a_waters(structures[pdb_id])
        phe238 = map_residue(238, pdb_id)
        lig_name = LIGAND_OF[pdb_id]

        waters_phe = get_waters_near_residue(chain, phe238, all_waters, WATER_SEARCH_RADIUS)
        waters_lig = get_waters_near_ligand(chain, lig_name, all_waters, WATER_SEARCH_RADIUS)

        ids_phe = set(w["resid"] for w in waters_phe)
        ids_lig = set(w["resid"] for w in waters_lig)
        bridging = ids_phe & ids_lig

        bridging_details = []
        for w in waters_phe:
            if w["resid"] in bridging:
                d_lig = None
                for w2 in waters_lig:
                    if w2["resid"] == w["resid"]:
                        d_lig = w2["dist_to_ligand"]
                        break
                bridging_details.append({
                    "resid": w["resid"],
                    "coord": w["coord"].tolist(),
                    "dist_phe238": w["dist_to_residue"],
                    "dist_ligand": d_lig,
                    "bfactor": w["bfactor"],
                })

        is_ga3 = lig_name == "GA3"
        species_label = SPECIES[pdb_id]
        marker = " (expected)" if is_ga3 else " (GA4, bridge unlikely)"
        print(f"    {pdb_id} ({species_label}, {lig_name}){marker}: "
              f"waters near Phe{phe238}={len(waters_phe)}, "
              f"bridging={len(bridging)}")
        for bd in bridging_details:
            print(f"      HOH {bd['resid']}: d(Phe{phe238})={bd['dist_phe238']:.2f}, "
                  f"d({lig_name})={bd['dist_ligand']:.2f}, "
                  f"B={bd['bfactor']:.1f}")

        results[pdb_id] = {
            "phe": phe238,
            "ligand": lig_name,
            "is_ga3": is_ga3,
            "n_bridging": len(bridging),
            "details": bridging_details,
            "conserved": len(bridging) > 0,
        }

    ga3_conserved = sum(1 for r in results.values()
                        if r["is_ga3"] and r["conserved"])
    ga4_present = sum(1 for r in results.values()
                      if not r["is_ga3"] and r["conserved"])
    print(f"    GA3 structures with bridge: {ga3_conserved}/2")
    print(f"    GA4 structures with bridge: {ga4_present}/2 "
          f"(expected 0 - no 13-OH)")
    return results


def analyze_bridge4(structures, chains):
    print("\n  Bridge 4: Asp289 hydration shell (6 waters)")
    print("  " + "-" * 50)
    results = {}

    for pdb_id in PDB_IDS:
        chain = chains[pdb_id]
        all_waters = get_all_chain_a_waters(structures[pdb_id])
        asp289 = map_residue(289, pdb_id)

        waters = get_waters_near_residue(chain, asp289, all_waters, WATER_SEARCH_RADIUS)

        hbond_waters = [w for w in waters if w["dist_to_residue"] <= HBOND_DIST]

        species_label = SPECIES[pdb_id]
        print(f"    {pdb_id} ({species_label}): "
              f"waters near Asp{asp289}={len(waters)}, "
              f"H-bond distance (<{HBOND_DIST}A)={len(hbond_waters)}")
        for w in hbond_waters:
            print(f"      HOH {w['resid']}: d={w['dist_to_residue']:.2f}, "
                  f"B={w['bfactor']:.1f}")

        results[pdb_id] = {
            "asp": asp289,
            "n_total_nearby": len(waters),
            "n_hbond": len(hbond_waters),
            "details": [{
                "resid": w["resid"],
                "coord": w["coord"].tolist(),
                "dist": w["dist_to_residue"],
                "bfactor": w["bfactor"],
            } for w in hbond_waters],
            "matches_literature": len(hbond_waters) >= 4,
        }

    counts = [r["n_hbond"] for r in results.values()]
    print(f"    Water count range: {min(counts)}-{max(counts)} "
          f"(literature: 6)")
    return results


def main():
    print("=" * 70)
    print("CONSERVED WATER MAPPING - Task 5")
    print("=" * 70)
    print(f"Output: {OUTPUT_DIR}")
    print(f"H-bond cutoff: {HBOND_DIST} A")
    print(f"Water search radius: {WATER_SEARCH_RADIUS} A")
    print(f"Cluster cutoff: {CLUSTER_CUTOFF} A")
    print(f"Conserved = present in >= {CONSERVED_MIN}/4 structures")
    print()

    parser = PDBParser(QUIET=True)
    structures = {}
    chains = {}

    print("Step 1: Loading and aligning structures...")
    for pdb_id in PDB_IDS:
        pdb_path = download_pdb(pdb_id)
        struct = parser.get_structure(pdb_id, pdb_path)
        structures[pdb_id] = struct
        chains[pdb_id] = struct[0]['A']

    ref_id = "2ZSH"
    ref_chain = chains[ref_id]
    ref_pocket = get_pocket_residues(ref_id)
    ref_bb = get_pocket_backbone_atoms(ref_chain, ref_pocket)

    print(f"  Aligning to {ref_id}...")
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
        print(f"    {pdb_id}: RMSD = {sup.rms:.4f} A")

    print(f"\nStep 2: Extracting pocket waters...")
    pocket_waters = {}
    for pdb_id in PDB_IDS:
        all_waters = get_all_chain_a_waters(structures[pdb_id])
        pocket_res = get_pocket_residues(pdb_id)
        pw = get_pocket_waters(chains[pdb_id], pocket_res, all_waters, cutoff=5.0)
        pocket_waters[pdb_id] = pw
        print(f"  {pdb_id}: {len(all_waters)} total waters, "
              f"{len(pw)} within 5 A of pocket")

    print(f"\n{'=' * 70}")
    print("Step 3: Analyzing specific water bridges...")
    print(f"{'=' * 70}")

    bridge1 = analyze_bridge1(structures, chains)
    bridge2 = analyze_bridge2(structures, chains)
    bridge3 = analyze_bridge3(structures, chains)
    bridge4 = analyze_bridge4(structures, chains)

    print(f"\n{'=' * 70}")
    print("Step 4: Global conserved water clustering...")
    print(f"{'=' * 70}")

    clusters = hierarchical_cluster_waters(pocket_waters, cutoff=CLUSTER_CUTOFF)

    conserved = [c for c in clusters if c["n_structures"] >= CONSERVED_MIN]
    fully_conserved = [c for c in clusters if c["n_structures"] == 4]

    print(f"\n  Total water clusters: {len(clusters)}")
    print(f"  Conserved (>={CONSERVED_MIN} structures): {len(conserved)}")
    print(f"  Fully conserved (4/4): {len(fully_conserved)}")

    for i, c in enumerate(conserved):
        center = c["center"]
        print(f"\n  Cluster {i+1}: {c['n_structures']}/4 structures "
              f"({', '.join(c['structures'])})")
        print(f"    Center: ({center[0]:.2f}, {center[1]:.2f}, "
              f"{center[2]:.2f})")
        print(f"    Spread: {c['spread']:.2f} A, "
              f"Total waters: {c['n_waters']}")

    print(f"\n{'=' * 70}")
    print("Step 5: Summary and implications")
    print(f"{'=' * 70}\n")

    b1_cons = sum(1 for r in bridge1.values() if r["conserved"])
    b2_cons = sum(1 for r in bridge2.values() if r["conserved"])
    b3_ga3 = sum(1 for r in bridge3.values() if r["is_ga3"] and r["conserved"])
    b4_match = sum(1 for r in bridge4.values() if r["matches_literature"])

    print("  Water bridge conservation:")
    print(f"    Bridge 1 (Ser116--HOH--Ser191): {b1_cons}/4 structures")
    print(f"    Bridge 2 (Tyr127--HOH--GA 3-OH): {b2_cons}/4 structures")
    print(f"    Bridge 3 (Phe238--HOH--GA3 13-OH): {b3_ga3}/2 GA3 structures")
    print(f"    Bridge 4 (Asp289 hydration): {b4_match}/4 structures "
          f"(>=4 waters)")

    print(f"\n  Global pocket water conservation:")
    print(f"    Fully conserved waters (4/4): {len(fully_conserved)}")
    print(f"    Mostly conserved (3/4): "
          f"{len([c for c in conserved if c['n_structures'] == 3])}")

    if b1_cons >= 3:
        print("\n  Bridge 1 IMPLICATION: The Ser116-Ser191 water network is")
        print("    conserved across species. Universal ligand should NOT")
        print("    displace this structural water - design around it.")
    else:
        print("\n  Bridge 1 IMPLICATION: The Ser116-Ser191 water bridge is")
        print("    NOT fully conserved. May be dispensable for binding.")

    if b2_cons >= 3:
        print("\n  Bridge 2 IMPLICATION: The Tyr127-ligand water bridge is")
        print("    conserved. The 3-OH pharmacophore is critical - universal")
        print("    ligand MUST retain a 3-OH equivalent for this interaction.")
    else:
        print("\n  Bridge 2 IMPLICATION: The Tyr127-ligand water is variable.")
        print("    The 3-OH pharmacophore may be flexible in design.")

    if b3_ga3 == 2:
        print("\n  Bridge 3 IMPLICATION: GA3-specific 13-OH water bridge is")
        print("    confirmed. Since GA4 lacks 13-OH and still binds well,")
        print("    the universal ligand does NOT need a 13-OH equivalent.")
    else:
        print("\n  Bridge 3 IMPLICATION: GA3-specific water bridge is")
        print("    inconsistent in the crystal structures.")

    if b4_match >= 3:
        print("\n  Bridge 4 IMPLICATION: Asp289 hydration shell is conserved.")
        print("    These waters form a structural network at the pocket base.")
        print("    Include these waters as part of the receptor in docking.")
    else:
        print("\n  Bridge 4 IMPLICATION: Asp289 hydration is variable.")
        print("    Consider flexible water treatment in this region.")

    print(f"\n  AutoGrow4 recommendation:")
    if len(fully_conserved) >= 3:
        print(f"    Include {len(fully_conserved)} fully conserved waters")
        print(f"    as explicit waters in the receptor PDB for docking.")
        print(f"    These waters are structural and should NOT be displaced.")
    else:
        print(f"    Few fully conserved waters found. Consider implicit")
        print(f"    solvation or flexible water treatment in docking.")

    results_json = {
        "parameters": {
            "hbond_cutoff": HBOND_DIST,
            "search_radius": WATER_SEARCH_RADIUS,
            "cluster_cutoff": CLUSTER_CUTOFF,
            "conserved_min": CONSERVED_MIN,
        },
        "pocket_water_counts": {
            pdb_id: len(pw) for pdb_id, pw in pocket_waters.items()
        },
        "bridge1_ser116_ser191": {
            k: {kk: vv for kk, vv in v.items() if kk != "water"}
            for k, v in bridge1.items()
        },
        "bridge2_tyr127_ligand": {
            k: {kk: vv for kk, vv in v.items() if kk != "water"}
            for k, v in bridge2.items()
        },
        "bridge3_phe238_ga3": {
            k: {kk: vv for kk, vv in v.items() if kk != "water"}
            for k, v in bridge3.items()
        },
        "bridge4_asp289": bridge4,
        "conserved_clusters": [{
            "center": c["center"].tolist(),
            "n_structures": c["n_structures"],
            "structures": c["structures"],
            "spread": c["spread"],
            "n_waters": c["n_waters"],
        } for c in conserved],
        "fully_conserved_clusters": [{
            "center": c["center"].tolist(),
            "n_structures": c["n_structures"],
            "structures": c["structures"],
            "spread": c["spread"],
        } for c in fully_conserved],
        "summary": {
            "bridge1_conserved": f"{b1_cons}/4",
            "bridge2_conserved": f"{b2_cons}/4",
            "bridge3_ga3_conserved": f"{b3_ga3}/2",
            "bridge4_matches_lit": f"{b4_match}/4",
            "fully_conserved_waters": len(fully_conserved),
            "mostly_conserved_waters": len(conserved),
        },
    }

    json_path = os.path.join(OUTPUT_DIR, "water_mapping_results.json")
    with open(json_path, "w") as f:
        json.dump(results_json, f, indent=2, default=str)
    print(f"\nJSON results: {json_path}")

    report_path = os.path.join(OUTPUT_DIR, "water_mapping_report.txt")
    _write_report(report_path, bridge1, bridge2, bridge3, bridge4,
                  conserved, fully_conserved, pocket_waters)
    print(f"Report: {report_path}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        _generate_plots(bridge1, bridge2, bridge3, bridge4,
                        conserved, fully_conserved, pocket_waters)
    except ImportError:
        print("matplotlib not available, skipping plots")

    vis_path = os.path.join(OUTPUT_DIR, "visualize_waters.py")
    _write_visualization(vis_path, conserved, fully_conserved,
                         bridge1, bridge4)
    print(f"Visualization: {vis_path}")


def _write_report(path, b1, b2, b3, b4, conserved, fully_conserved,
                  pocket_waters):
    with open(path, "w") as f:
        f.write("CONSERVED WATER MAPPING - Task 5\n")
        f.write(f"{'=' * 60}\n\n")

        f.write("Water Bridge Analysis:\n\n")

        f.write("Bridge 1: Ser116--HOH--Ser191\n")
        for pdb_id in PDB_IDS:
            r = b1[pdb_id]
            status = "PRESENT" if r["conserved"] else "ABSENT"
            f.write(f"  {pdb_id}: {status}, {r['n_bridging']} bridging "
                    f"water(s)\n")
            for d in r["details"]:
                f.write(f"    HOH {d['resid']}: d(Ser{r['ser_a']})="
                        f"{d['dist_ser116']:.2f}, d(Ser{r['ser_b']})="
                        f"{d['dist_ser191']:.2f}, B={d['bfactor']:.1f}\n")

        f.write(f"\nBridge 2: Tyr127--HOH--GA 3-OH\n")
        for pdb_id in PDB_IDS:
            r = b2[pdb_id]
            status = "PRESENT" if r["conserved"] else "ABSENT"
            f.write(f"  {pdb_id} ({r['ligand']}): {status}, "
                    f"{r['n_bridging']} bridging water(s)\n")
            for d in r["details"]:
                f.write(f"    HOH {d['resid']}: d(Tyr{r['tyr']})="
                        f"{d['dist_tyr127']:.2f}, d({r['ligand']})="
                        f"{d['dist_ligand']:.2f}, B={d['bfactor']:.1f}\n")

        f.write(f"\nBridge 3: Phe238--HOH--GA3 13-OH (GA3-specific)\n")
        for pdb_id in PDB_IDS:
            r = b3[pdb_id]
            ga_label = " (GA3)" if r["is_ga3"] else " (GA4, no 13-OH)"
            status = "PRESENT" if r["conserved"] else "ABSENT"
            f.write(f"  {pdb_id}{ga_label}: {status}, "
                    f"{r['n_bridging']} bridging water(s)\n")

        f.write(f"\nBridge 4: Asp289 hydration shell\n")
        for pdb_id in PDB_IDS:
            r = b4[pdb_id]
            f.write(f"  {pdb_id}: {r['n_hbond']} waters within H-bond "
                    f"distance of Asp{r['asp']}\n")
            for d in r["details"]:
                f.write(f"    HOH {d['resid']}: d={d['dist']:.2f}, "
                        f"B={d['bfactor']:.1f}\n")

        f.write(f"\n\nGlobal Conserved Water Positions:\n")
        f.write(f"{'=' * 60}\n\n")
        f.write(f"Total pocket waters per structure:\n")
        for pdb_id in PDB_IDS:
            f.write(f"  {pdb_id}: {len(pocket_waters[pdb_id])}\n")

        f.write(f"\nFully conserved (4/4 structures):\n")
        for i, c in enumerate(fully_conserved):
            ctr = c["center"]
            f.write(f"  Water {i+1}: ({ctr[0]:.2f}, {ctr[1]:.2f}, "
                    f"{ctr[2]:.2f}), spread={c['spread']:.2f} A\n")

        f.write(f"\nMostly conserved (3/4 structures):\n")
        mostly = [c for c in conserved if c["n_structures"] == 3]
        for i, c in enumerate(mostly):
            ctr = c["center"]
            f.write(f"  Water {i+1}: ({ctr[0]:.2f}, {ctr[1]:.2f}, "
                    f"{ctr[2]:.2f}), in {', '.join(c['structures'])}, "
                    f"spread={c['spread']:.2f} A\n")

        f.write(f"\n\nImplications for Universal Ligand Design:\n")
        f.write(f"{'=' * 60}\n\n")

        b1_cons = sum(1 for r in b1.values() if r["conserved"])
        b2_cons = sum(1 for r in b2.values() if r["conserved"])

        if b1_cons >= 3:
            f.write("- Ser116-Ser191 water bridge is CONSERVED. Design\n")
            f.write("  ligand to work WITH this water, not displace it.\n\n")

        if b2_cons >= 3:
            f.write("- Tyr127-ligand water bridge is CONSERVED. The 3-OH\n")
            f.write("  group (or H-bond equivalent) is essential.\n\n")

        f.write("- GA3's 13-OH water bridge (Phe238) is GA3-SPECIFIC.\n")
        f.write("  Not required for binding (GA4 binds without it).\n")
        f.write("  Universal ligand does NOT need 13-OH equivalent.\n\n")

        b4_match = sum(1 for r in b4.values() if r["matches_literature"])
        if b4_match >= 3:
            f.write("- Asp289 hydration shell is CONSERVED. Include these\n")
            f.write("  waters as explicit receptor waters in AutoGrow4.\n\n")

        f.write(f"- {len(fully_conserved)} fully conserved water position(s)\n")
        f.write(f"  should be included in the receptor for docking.\n")


def _generate_plots(b1, b2, b3, b4, conserved, fully_conserved,
                    pocket_waters):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle("Conserved Water Mapping - Task 5\n"
                 "Water Bridge Analysis Across 4 GID1 Structures",
                 fontsize=14, fontweight="bold")

    ax1 = axes[0][0]
    ax1.set_title("Bridge 1: Ser116--HOH--Ser191", fontweight="bold")
    pdb_labels = PDB_IDS
    counts = [b1[p]["n_bridging"] for p in PDB_IDS]
    colors = ["#4CAF50" if c > 0 else "#F44336" for c in counts]
    bars = ax1.bar(pdb_labels, counts, color=colors, edgecolor="black")
    ax1.set_ylabel("Bridging waters")
    ax1.set_ylim(0, max(counts + [1]) + 1)
    for i, v in enumerate(counts):
        ax1.text(i, v + 0.1, str(v), ha="center", fontweight="bold")
    ax1.axhline(y=0.5, color="gray", linestyle="--", alpha=0.3)

    ax2 = axes[0][1]
    ax2.set_title("Bridge 2: Tyr127--HOH--GA 3-OH", fontweight="bold")
    counts2 = [b2[p]["n_bridging"] for p in PDB_IDS]
    colors2 = ["#4CAF50" if c > 0 else "#F44336" for c in counts2]
    ax2.bar(pdb_labels, counts2, color=colors2, edgecolor="black")
    ax2.set_ylabel("Bridging waters")
    ax2.set_ylim(0, max(counts2 + [1]) + 1)
    for i, v in enumerate(counts2):
        lbl = f"{v}\n({LIGAND_OF[PDB_IDS[i]]})"
        ax2.text(i, v + 0.1, lbl, ha="center", fontweight="bold",
                 fontsize=9)

    ax3 = axes[1][0]
    ax3.set_title("Bridge 4: Asp289 Hydration Shell", fontweight="bold")
    counts4 = [b4[p]["n_hbond"] for p in PDB_IDS]
    colors4 = ["#4CAF50" if c >= 4 else "#FFC107" if c >= 2
               else "#F44336" for c in counts4]
    ax3.bar(pdb_labels, counts4, color=colors4, edgecolor="black")
    ax3.axhline(y=6, color="blue", linestyle="--", linewidth=1.5,
                label="Literature (6 waters)")
    ax3.set_ylabel("Waters within H-bond distance")
    ax3.set_ylim(0, max(counts4 + [7]) + 1)
    for i, v in enumerate(counts4):
        ax3.text(i, v + 0.15, str(v), ha="center", fontweight="bold")
    ax3.legend(fontsize=9)

    ax4 = axes[1][1]
    ax4.set_title("Global Conserved Water Clusters", fontweight="bold")
    n_by_cons = defaultdict(int)
    for c in conserved:
        n_by_cons[c["n_structures"]] += 1

    cons_labels = ["4/4\n(fully)", "3/4\n(mostly)"]
    cons_vals = [n_by_cons.get(4, 0), n_by_cons.get(3, 0)]
    cons_colors = ["#2196F3", "#90CAF9"]
    ax4.bar(cons_labels, cons_vals, color=cons_colors, edgecolor="black",
            width=0.5)
    ax4.set_ylabel("Number of water clusters")
    for i, v in enumerate(cons_vals):
        ax4.text(i, v + 0.15, str(v), ha="center", fontweight="bold",
                 fontsize=12)

    total_pocket = {p: len(pw) for p, pw in pocket_waters.items()}
    info = "Pocket waters: " + ", ".join(
        f"{p}={total_pocket[p]}" for p in PDB_IDS)
    ax4.text(0.5, 0.95, info, transform=ax4.transAxes, fontsize=8,
             ha="center", va="top",
             bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    for ax in axes.flat:
        ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    save_path = os.path.join(OUTPUT_DIR, "water_mapping_analysis.png")
    plt.savefig(save_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Analysis plot: {save_path}")

    fig2, ax = plt.subplots(figsize=(12, 8))
    ax.set_title("Water Bridge Conservation Summary\n"
                 "Across Arabidopsis and Rice GID1 Structures",
                 fontsize=13, fontweight="bold")

    bridge_labels = [
        "Bridge 1\nSer116-HOH-Ser191\n(structural)",
        "Bridge 2\nTyr127-HOH-GA 3-OH\n(ligand recognition)",
        "Bridge 3\nPhe238-HOH-GA3 13-OH\n(GA3-specific)",
        "Bridge 4\nAsp289 hydration\n(pocket base)",
    ]

    b3_ga3 = sum(1 for r in b3.values() if r["is_ga3"] and r["conserved"])
    b4_match = sum(1 for r in b4.values() if r["matches_literature"])

    conservation = [
        sum(1 for r in b1.values() if r["conserved"]),
        sum(1 for r in b2.values() if r["conserved"]),
        b3_ga3,
        b4_match,
    ]
    max_possible = [4, 4, 2, 4]

    x = np.arange(len(bridge_labels))
    width = 0.35

    bars1 = ax.bar(x - width / 2, conservation, width, label="Conserved",
                   color="#4CAF50", edgecolor="black")
    bars2 = ax.bar(x + width / 2, max_possible, width, label="Max possible",
                   color="#E0E0E0", edgecolor="black")

    ax.set_ylabel("Structures", fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(bridge_labels, fontsize=9)
    ax.legend(fontsize=11)
    ax.set_ylim(0, 5.5)

    for i, (c, m) in enumerate(zip(conservation, max_possible)):
        pct = c / m * 100
        ax.text(i - width / 2, c + 0.1, f"{c}/{m}\n({pct:.0f}%)",
                ha="center", fontweight="bold", fontsize=10)

    plt.tight_layout()
    save_path2 = os.path.join(OUTPUT_DIR, "water_bridge_summary.png")
    plt.savefig(save_path2, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Summary plot: {save_path2}")


def _write_visualization(path, conserved, fully_conserved, b1, b4):
    fc_coords = []
    for c in fully_conserved:
        ctr = c["center"]
        fc_coords.append(f"[{ctr[0]:.2f}, {ctr[1]:.2f}, {ctr[2]:.2f}]")

    fc_list = ", ".join(fc_coords) if fc_coords else "[]"

    arab_res = ", ".join(str(r) for r in ARAB_RESIDUES)
    rice_res = ", ".join(str(r) for r in RICE_RESIDUES)

    script = f'''"""
Conserved Water Visualization - Task 5
Paste into Jupyter notebook (ag4_full kernel).

Shows all 4 aligned GID1 structures with:
  - Active site residues (sticks)
  - Crystal ligands (green=Arabidopsis, cyan=Rice)
  - Conserved water positions (red spheres)
"""
import os, time
import nglview as nv
from IPython.display import display, HTML

BASE = os.path.expanduser("~/final_project")
WATER_DIR = os.path.join(BASE, "water_mapping")
BOX_DIR = os.path.join(BASE, "box_optimization")

ARAB_RESIDUES = [{arab_res}]
RICE_RESIDUES = [{rice_res}]

CONSERVED_WATERS = [{fc_list}]

display(HTML("""
<pre style='background:#ffffff !important; color:#000000 !important; padding:15px; border-radius:8px;
            margin-bottom:15px; border:2px solid #1565C0; font-family:monospace; font-size:14px;
            line-height:1.6; white-space:pre;'><span style='color:#1565C0 !important; font-size:20px; font-weight:bold;'>Task 5: Conserved Water Mapping</span>

4-Structure Ensemble: 2ZSH, 2ZSI, 3ED1, 3EBL

  <span style='color:#228B22 !important; font-size:16px;'>&#9632;</span> Arabidopsis ligands (green)
  <span style='color:#00CED1 !important; font-size:16px;'>&#9632;</span> Rice ligands (cyan)
  <span style='color:#FF0000 !important; font-size:16px;'>&#9679;</span> Conserved water positions (red spheres)
  <span style='color:#888 !important; font-size:16px;'>&#9632;</span> Active site residues (element coloring)

Conserved waters: """ + str(len(CONSERVED_WATERS)) + """ positions (present in 4/4 structures)</pre>
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
    pocket = ARAB_RESIDUES if species == "Arabidopsis" else RICE_RESIDUES
    res_sel = "(" + " or ".join(str(r) for r in pocket) + ") and sidechainAttached"
    view.add_licorice(selection=res_sel, color="element", radius=0.18,
                      component=comp_idx)
    comp_idx += 1

for pdb_id in PDB_IDS:
    lig_name = LIGAND_OF[pdb_id]
    lig_path = os.path.join(BOX_DIR,
                            f"{{pdb_id}}_{{lig_name}}_ligand_aligned.pdb")
    if not os.path.exists(lig_path):
        lig_path = os.path.join(BOX_DIR,
                                f"{{pdb_id}}_{{lig_name}}_ligand.pdb")
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

if CONSERVED_WATERS:
    water_pdb_lines = []
    for i, coord in enumerate(CONSERVED_WATERS):
        x, y, z = coord
        water_pdb_lines.append(
            f"HETATM{{i+1:5d}}  O   HOH A{{i+1:4d}}    "
            f"{{x:8.3f}}{{y:8.3f}}{{z:8.3f}}  1.00 10.00           O"
        )
    water_pdb = "\\n".join(water_pdb_lines) + "\\nEND\\n"

    view.add_component(nv.TextStructure(water_pdb, ext="pdb"),
                       default_representation=False,
                       name="conserved_waters")
    time.sleep(0.2)
    view.add_spacefill(selection="all", color="red", radius=0.5,
                       opacity=0.7, component=comp_idx)
    comp_idx += 1

view.center()
display(view)
'''

    with open(path, "w") as f:
        f.write(script)


if __name__ == "__main__":
    main()
