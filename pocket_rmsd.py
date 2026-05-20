#!/usr/bin/env python3
"""
Pocket RMSD Analysis - Universal Ligand Feasibility (Task 1)

Aligns 4 GID1 crystal structures on 14 key binding pocket residues
and calculates per-residue and overall pocket RMSD.

Structures:
  2ZSH: Arabidopsis GID1A + GA3
  2ZSI: Arabidopsis GID1A + GA4
  3ED1: Rice GID1   + GA3 (hexamer, chain A used)
  3EBL: Rice GID1   + GA4 (hexamer, chain A used)

Decision rule: Pocket RMSD < 0.5 A => universal ligand feasible
"""
import os
import sys
import urllib.request
import numpy as np

from Bio.PDB import PDBParser, Superimposer, PDBIO, Select

OUTPUT_DIR = os.path.expanduser("~/final_project/pocket_analysis")
os.makedirs(OUTPUT_DIR, exist_ok=True)

PDB_IDS = ["2ZSH", "2ZSI", "3ED1", "3EBL"]

SPECIES = {
    "2ZSH": "Arabidopsis", "2ZSI": "Arabidopsis",
    "3ED1": "Rice", "3EBL": "Rice",
}
LIGANDS = {
    "2ZSH": "GA3", "2ZSI": "GA4",
    "3ED1": "GA3", "3EBL": "GA4",
}

ARAB_RESIDUES = [24, 27, 31, 116, 119, 126, 127, 191, 218, 238, 239, 244, 319, 323]

RICE_RESIDUES = [24, 27, 31, 123, 126, 133, 134, 198, 225, 245, 246, 251, 326, 330]

RESIDUE_MAP = {}
for a, r in zip(ARAB_RESIDUES, RICE_RESIDUES):
    RESIDUE_MAP[("Arabidopsis", a)] = a
    RESIDUE_MAP[("Rice", a)] = r

BACKBONE_ATOMS = ["N", "CA", "C", "O"]


def download_pdb(pdb_id, output_dir):
    path = os.path.join(output_dir, f"{pdb_id}.pdb")
    if os.path.exists(path) and os.path.getsize(path) > 1000:
        return path
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"  Downloading {pdb_id}...")
    urllib.request.urlretrieve(url, path)
    return path


def get_pocket_residues(pdb_id):
    species = SPECIES[pdb_id]
    if species == "Arabidopsis":
        return ARAB_RESIDUES
    else:
        return RICE_RESIDUES


def get_pocket_atoms(chain, pocket_residues, atom_names=None):
    if atom_names is None:
        atom_names = BACKBONE_ATOMS
    atoms = []
    for resnum in pocket_residues:
        try:
            res = chain[(' ', resnum, ' ')]
            for aname in atom_names:
                if aname in res:
                    atoms.append(res[aname])
        except KeyError:
            pass
    return atoms


def get_sidechain_atoms(chain, pocket_residues):
    atoms = []
    for resnum in pocket_residues:
        try:
            res = chain[(' ', resnum, ' ')]
            for atom in res:
                if atom.name not in ["N", "CA", "C", "O", "H", "HA"]:
                    atoms.append(atom)
        except KeyError:
            pass
    return atoms


def calc_per_residue_rmsd(chain1, res_list1, chain2, res_list2):
    results = []
    for r1, r2, arab_num in zip(res_list1, res_list2, ARAB_RESIDUES):
        try:
            res1 = chain1[(' ', r1, ' ')]
            res2 = chain2[(' ', r2, ' ')]
        except KeyError:
            results.append({
                "arab_resnum": arab_num,
                "resname": "???",
                "res1_num": r1, "res2_num": r2,
                "backbone_rmsd": None,
                "all_heavy_rmsd": None,
            })
            continue

        bb_coords1, bb_coords2 = [], []
        for aname in BACKBONE_ATOMS:
            if aname in res1 and aname in res2:
                bb_coords1.append(res1[aname].get_vector().get_array())
                bb_coords2.append(res2[aname].get_vector().get_array())

        bb_rmsd = None
        if bb_coords1:
            c1 = np.array(bb_coords1)
            c2 = np.array(bb_coords2)
            bb_rmsd = np.sqrt(np.mean(np.sum((c1 - c2) ** 2, axis=1)))

        heavy1, heavy2 = [], []
        matched_names = set()
        for atom in res1:
            if atom.element != "H" and atom.name in [a.name for a in res2]:
                matched_names.add(atom.name)
        for aname in sorted(matched_names):
            heavy1.append(res1[aname].get_vector().get_array())
            heavy2.append(res2[aname].get_vector().get_array())

        all_rmsd = None
        if heavy1:
            c1 = np.array(heavy1)
            c2 = np.array(heavy2)
            all_rmsd = np.sqrt(np.mean(np.sum((c1 - c2) ** 2, axis=1)))

        results.append({
            "arab_resnum": arab_num,
            "resname": res1.get_resname(),
            "res1_num": r1, "res2_num": r2,
            "backbone_rmsd": bb_rmsd,
            "all_heavy_rmsd": all_rmsd,
        })

    return results


class ChainASelect(Select):
    def accept_chain(self, chain):
        return chain.id == "A"

    def accept_residue(self, residue):
        return True


def extract_chain_a(struct, pdb_id, output_dir):
    out_path = os.path.join(output_dir, f"{pdb_id}_chainA.pdb")
    io = PDBIO()
    io.set_structure(struct)
    io.save(out_path, ChainASelect())
    return out_path


def main():
    print("=" * 70)
    print("POCKET RMSD ANALYSIS - Universal Ligand Feasibility (Task 1)")
    print("=" * 70)
    print(f"Output: {OUTPUT_DIR}")
    print(f"Decision rule: Pocket RMSD < 0.5 A => feasible\n")

    parser = PDBParser(QUIET=True)
    structures = {}
    chains = {}

    print("Loading structures...")
    for pdb_id in PDB_IDS:
        pdb_path = download_pdb(pdb_id, OUTPUT_DIR)
        struct = parser.get_structure(pdb_id, pdb_path)
        structures[pdb_id] = struct
        chains[pdb_id] = struct[0]['A']

        n_chains = len(list(struct[0].get_chains()))
        pocket_res = get_pocket_residues(pdb_id)
        n_found = sum(1 for r in pocket_res
                      if (' ', r, ' ') in chains[pdb_id])

        print(f"  {pdb_id}: {SPECIES[pdb_id]} + {LIGANDS[pdb_id]}, "
              f"{n_chains} chain(s), {n_found}/14 pocket residues found")

        if n_chains > 2:
            extract_chain_a(struct, pdb_id, OUTPUT_DIR)
            print(f"    Extracted chain A monomer")

    print(f"\nResidue mapping verified:")
    print(f"  Arabidopsis: {ARAB_RESIDUES}")
    print(f"  Rice (+7):   {RICE_RESIDUES}")

    print(f"\n  Amino acid conservation check:")
    for i, arab_num in enumerate(ARAB_RESIDUES):
        row = f"    Arab {arab_num:>3}: "
        types = []
        for pdb_id in PDB_IDS:
            pocket_res = get_pocket_residues(pdb_id)
            resnum = pocket_res[i]
            try:
                res = chains[pdb_id][(' ', resnum, ' ')]
                types.append(f"{pdb_id}={res.get_resname()}{resnum}")
            except KeyError:
                types.append(f"{pdb_id}=???")
        all_same = len(set(t.split("=")[1][:3] for t in types)) == 1
        conserved = "CONSERVED" if all_same else "VARIABLE"
        row += " | ".join(types) + f"  [{conserved}]"
        print(row)

    ref_pdb = "2ZSH"
    ref_chain = chains[ref_pdb]
    ref_pocket = get_pocket_residues(ref_pdb)
    ref_atoms = get_pocket_atoms(ref_chain, ref_pocket)

    print(f"\n{'=' * 70}")
    print(f"PAIRWISE POCKET ALIGNMENT & RMSD")
    print(f"Reference: {ref_pdb} ({SPECIES[ref_pdb]} + {LIGANDS[ref_pdb]})")
    print(f"Alignment atoms: backbone (N, CA, C, O) of 14 pocket residues")
    print(f"{'=' * 70}")

    all_pair_results = []

    pairs = [
        ("2ZSH", "2ZSI", "Same species, different ligand"),
        ("3ED1", "3EBL", "Same species, different ligand"),
        ("2ZSH", "3ED1", "Cross-species, same ligand (GA3)"),
        ("2ZSH", "3EBL", "Cross-species, different ligand"),
        ("2ZSI", "3ED1", "Cross-species, different ligand"),
        ("2ZSI", "3EBL", "Cross-species, same ligand (GA4)"),
    ]

    for pdb1, pdb2, description in pairs:
        chain1 = chains[pdb1]
        chain2 = chains[pdb2]
        pocket1 = get_pocket_residues(pdb1)
        pocket2 = get_pocket_residues(pdb2)

        atoms1 = get_pocket_atoms(chain1, pocket1)
        atoms2 = get_pocket_atoms(chain2, pocket2)

        n_atoms = min(len(atoms1), len(atoms2))
        atoms1 = atoms1[:n_atoms]
        atoms2 = atoms2[:n_atoms]

        sup = Superimposer()
        sup.set_atoms(atoms1, atoms2)
        sup.apply(chain2.get_atoms())

        overall_rmsd = sup.rms

        per_res = calc_per_residue_rmsd(chain1, pocket1, chain2, pocket2)

        bb_rmsds = [r["backbone_rmsd"] for r in per_res if r["backbone_rmsd"] is not None]
        heavy_rmsds = [r["all_heavy_rmsd"] for r in per_res if r["all_heavy_rmsd"] is not None]

        result = {
            "pdb1": pdb1, "pdb2": pdb2,
            "description": description,
            "overall_rmsd": overall_rmsd,
            "mean_bb_rmsd": np.mean(bb_rmsds) if bb_rmsds else None,
            "max_bb_rmsd": np.max(bb_rmsds) if bb_rmsds else None,
            "mean_heavy_rmsd": np.mean(heavy_rmsds) if heavy_rmsds else None,
            "max_heavy_rmsd": np.max(heavy_rmsds) if heavy_rmsds else None,
            "per_residue": per_res,
            "feasible": overall_rmsd < 0.5,
        }
        all_pair_results.append(result)

        status = "PASS" if overall_rmsd < 0.5 else "FAIL"
        print(f"\n  {pdb1} vs {pdb2}: {description}")
        print(f"    Overall pocket RMSD (backbone): {overall_rmsd:.4f} A  [{status}]")
        print(f"    Mean per-residue backbone RMSD:  {result['mean_bb_rmsd']:.4f} A")
        print(f"    Max per-residue backbone RMSD:   {result['max_bb_rmsd']:.4f} A")
        print(f"    Mean per-residue all-heavy RMSD: {result['mean_heavy_rmsd']:.4f} A")
        print(f"    Max per-residue all-heavy RMSD:  {result['max_heavy_rmsd']:.4f} A")

        print(f"\n    {'Res':>6} {'Name':>4} "
              f"{'BB_RMSD':>8} {'Heavy_RMSD':>11} {'Status':>7}")
        print(f"    {'─' * 42}")
        for r in per_res:
            bb = f"{r['backbone_rmsd']:.4f}" if r['backbone_rmsd'] is not None else "N/A"
            hv = f"{r['all_heavy_rmsd']:.4f}" if r['all_heavy_rmsd'] is not None else "N/A"
            flag = ""
            if r['all_heavy_rmsd'] is not None and r['all_heavy_rmsd'] > 0.5:
                flag = " *HIGH*"
            print(f"    {r['arab_resnum']:>6} {r['resname']:>4} "
                  f"{bb:>8} {hv:>11}{flag}")

    print(f"\n{'=' * 70}")
    print(f"SUMMARY")
    print(f"{'=' * 70}\n")

    print(f"  {'Comparison':<35} {'RMSD (A)':>10} {'Status':>8}")
    print(f"  {'─' * 55}")
    for r in all_pair_results:
        status = "PASS" if r['feasible'] else "FAIL"
        label = f"{r['pdb1']} vs {r['pdb2']}"
        print(f"  {label:<35} {r['overall_rmsd']:>10.4f} {status:>8}")

    cross_species = [r for r in all_pair_results
                     if SPECIES[r['pdb1']] != SPECIES[r['pdb2']]]
    same_species = [r for r in all_pair_results
                    if SPECIES[r['pdb1']] == SPECIES[r['pdb2']]]

    if cross_species:
        mean_cross = np.mean([r['overall_rmsd'] for r in cross_species])
        max_cross = np.max([r['overall_rmsd'] for r in cross_species])
        all_pass = all(r['feasible'] for r in cross_species)

        print(f"\n  Cross-species pocket RMSD:")
        print(f"    Mean: {mean_cross:.4f} A")
        print(f"    Max:  {max_cross:.4f} A")
        print(f"    All < 0.5 A: {'YES' if all_pass else 'NO'}")

    if same_species:
        mean_same = np.mean([r['overall_rmsd'] for r in same_species])
        print(f"\n  Within-species pocket RMSD (ligand-induced changes):")
        print(f"    Mean: {mean_same:.4f} A")

    print(f"\n{'=' * 70}")
    if cross_species and all(r['feasible'] for r in cross_species):
        print("  CONCLUSION: Universal ligand is FEASIBLE by pocket geometry.")
        print("  All cross-species pocket RMSDs are below 0.5 A threshold.")
        print("  The 14 key binding residues are structurally conserved")
        print("  between Arabidopsis thaliana and Oryza sativa (rice) GID1.")
    elif cross_species:
        failing = [r for r in cross_species if not r['feasible']]
        print("  CONCLUSION: Universal ligand faces CHALLENGES.")
        print(f"  {len(failing)} cross-species comparison(s) exceed 0.5 A:")
        for f in failing:
            print(f"    {f['pdb1']} vs {f['pdb2']}: {f['overall_rmsd']:.4f} A")
        high_res = []
        for r in cross_species:
            for pr in r['per_residue']:
                if pr['all_heavy_rmsd'] is not None and pr['all_heavy_rmsd'] > 0.5:
                    high_res.append((r['pdb1'], r['pdb2'],
                                    pr['arab_resnum'], pr['resname'],
                                    pr['all_heavy_rmsd']))
        if high_res:
            print(f"\n  Residues with highest deviation:")
            for p1, p2, resnum, resname, rmsd in sorted(set(high_res)):
                print(f"    {resname}{resnum} in {p1} vs {p2}: {rmsd:.4f} A")
    print(f"{'=' * 70}")

    report_path = os.path.join(OUTPUT_DIR, "pocket_rmsd_report.txt")
    with open(report_path, "w") as f:
        f.write("POCKET RMSD ANALYSIS - Universal Ligand Feasibility\n")
        f.write(f"{'=' * 60}\n\n")
        f.write("Structures analyzed:\n")
        for pdb_id in PDB_IDS:
            f.write(f"  {pdb_id}: {SPECIES[pdb_id]} GID1 + {LIGANDS[pdb_id]}\n")
        f.write(f"\nArabidopsis pocket residues: {ARAB_RESIDUES}\n")
        f.write(f"Rice pocket residues:        {RICE_RESIDUES}\n")
        f.write(f"Decision threshold: 0.5 A\n\n")

        for r in all_pair_results:
            status = "PASS" if r['feasible'] else "FAIL"
            f.write(f"{r['pdb1']} vs {r['pdb2']}: {r['description']}\n")
            f.write(f"  Overall backbone RMSD: {r['overall_rmsd']:.4f} A [{status}]\n")
            f.write(f"  Mean per-residue BB:   {r['mean_bb_rmsd']:.4f} A\n")
            f.write(f"  Max per-residue BB:    {r['max_bb_rmsd']:.4f} A\n")
            f.write(f"  Mean all-heavy:        {r['mean_heavy_rmsd']:.4f} A\n")
            f.write(f"  Max all-heavy:         {r['max_heavy_rmsd']:.4f} A\n\n")

            f.write(f"  {'Res':>6} {'Name':>4} {'BB_RMSD':>8} {'Heavy_RMSD':>11}\n")
            for pr in r['per_residue']:
                bb = f"{pr['backbone_rmsd']:.4f}" if pr['backbone_rmsd'] else "N/A"
                hv = f"{pr['all_heavy_rmsd']:.4f}" if pr['all_heavy_rmsd'] else "N/A"
                f.write(f"  {pr['arab_resnum']:>6} {pr['resname']:>4} {bb:>8} {hv:>11}\n")
            f.write("\n")

        if cross_species:
            all_pass = all(r['feasible'] for r in cross_species)
            f.write(f"\nCONCLUSION: Universal ligand is "
                    f"{'FEASIBLE' if all_pass else 'CHALLENGED'} "
                    f"by pocket geometry.\n")

    print(f"\nReport saved to: {report_path}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        _generate_plots(all_pair_results)
    except ImportError:
        print("matplotlib not available, skipping plots")


def _generate_plots(all_pair_results):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=(20, 12))
    fig.suptitle("Pocket RMSD Analysis: 14 Key GID1 Binding Residues",
                 fontsize=16, fontweight="bold", y=0.98)

    for idx, result in enumerate(all_pair_results):
        ax = axes[idx // 3][idx % 3]
        per_res = result["per_residue"]

        labels = [f"{r['resname']}{r['arab_resnum']}" for r in per_res]
        bb_vals = [r['backbone_rmsd'] if r['backbone_rmsd'] else 0
                   for r in per_res]
        hv_vals = [r['all_heavy_rmsd'] if r['all_heavy_rmsd'] else 0
                   for r in per_res]

        x = np.arange(len(labels))
        w = 0.35
        bars1 = ax.bar(x - w / 2, bb_vals, w, label="Backbone",
                        color="#2196F3", edgecolor="black", linewidth=0.5)
        bars2 = ax.bar(x + w / 2, hv_vals, w, label="All Heavy",
                        color="#FF9800", edgecolor="black", linewidth=0.5)

        ax.axhline(y=0.5, color="red", linestyle="--", linewidth=1.5,
                   label="0.5 A threshold")

        ax.set_xlabel("Residue", fontsize=9)
        ax.set_ylabel("RMSD (A)", fontsize=9)
        ax.set_title(f"{result['pdb1']} vs {result['pdb2']}\n"
                     f"{result['description']}\n"
                     f"Overall: {result['overall_rmsd']:.4f} A",
                     fontsize=10, fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
        ax.legend(fontsize=7, loc="upper right")
        ax.set_ylim(0, max(max(hv_vals) * 1.3, 0.6))
        ax.grid(axis="y", alpha=0.3)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    save_path = os.path.join(OUTPUT_DIR, "pocket_rmsd_analysis.png")
    plt.savefig(save_path, dpi=200, bbox_inches="tight",
                facecolor="white")
    plt.close()
    print(f"Plot saved to: {save_path}")

    fig2, ax = plt.subplots(figsize=(10, 6))
    labels = [f"{r['pdb1']} vs {r['pdb2']}" for r in all_pair_results]
    rmsds = [r['overall_rmsd'] for r in all_pair_results]
    colors = []
    for r in all_pair_results:
        if SPECIES[r['pdb1']] == SPECIES[r['pdb2']]:
            colors.append("#4CAF50")
        else:
            colors.append("#2196F3")

    bars = ax.bar(range(len(labels)), rmsds, color=colors,
                  edgecolor="black", linewidth=0.5)
    ax.axhline(y=0.5, color="red", linestyle="--", linewidth=2,
               label="0.5 A threshold")

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=10)
    ax.set_ylabel("Pocket RMSD (A)", fontsize=12)
    ax.set_title("Overall Pocket RMSD: All Pairwise Comparisons",
                 fontsize=14, fontweight="bold")
    ax.legend(fontsize=10)
    ax.grid(axis="y", alpha=0.3)

    for i, v in enumerate(rmsds):
        ax.text(i, v + 0.01, f"{v:.4f}", ha="center", fontsize=9,
                fontweight="bold")

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#4CAF50", label="Within species"),
        Patch(facecolor="#2196F3", label="Cross species"),
        plt.Line2D([0], [0], color="red", linestyle="--", label="0.5 A threshold"),
    ]
    ax.legend(handles=legend_elements, fontsize=10)

    plt.tight_layout()
    save_path2 = os.path.join(OUTPUT_DIR, "pocket_rmsd_summary.png")
    plt.savefig(save_path2, dpi=200, bbox_inches="tight",
                facecolor="white")
    plt.close()
    print(f"Summary plot saved to: {save_path2}")


if __name__ == "__main__":
    main()
