#!/usr/bin/env python3
"""
Collect DiffGUI + GNINA results and score through the plant mobility pipeline.
Produces a ranked CSV combining binding affinity and plant systemicity scores.
"""

import csv
import glob
import os
import sys
from pathlib import Path

from rdkit import Chem

# Add project root to path for plant_scoring imports
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from plant_scoring.pipeline import process_compound


OUTPUT_DIR = "diffgui_run"


def read_sdf_mols(sdf_path):
    """Yield (mol, props_dict) for each molecule in an SDF file."""
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=True)
    for mol in suppl:
        if mol is None:
            continue
        props = {}
        for name in mol.GetPropsAsDict():
            props[name] = mol.GetPropsAsDict()[name]
        yield mol, props


def collect_gnina_scores(rescore_dir):
    """Map (scaffold, mol_id) -> {cnn_score, cnn_affinity, vina_score}."""
    scores = {}
    for sdf_path in glob.glob(os.path.join(rescore_dir, "*.sdf")):
        basename = Path(sdf_path).stem
        parts = basename.split("_")
        scaffold = parts[0] if parts else "unknown"

        for mol, props in read_sdf_mols(sdf_path):
            smi = Chem.MolToSmiles(mol)
            key = (scaffold, smi)
            scores[key] = {
                "cnn_score": props.get("CNNscore", ""),
                "cnn_affinity": props.get("CNNaffinity", ""),
                "vina_score": props.get("minimizedAffinity",
                              props.get("Vina", "")),
            }
    return scores


def main():
    results = []
    mol_counter = 0

    # Collect GNINA rescores if available
    rescore_dir = os.path.join(OUTPUT_DIR, "gnina_rescores")
    gnina_scores = {}
    if os.path.isdir(rescore_dir):
        gnina_scores = collect_gnina_scores(rescore_dir)
        print(f"Loaded {len(gnina_scores)} GNINA rescores")

    # Process each scaffold's output
    outputs_dir = os.path.join(OUTPUT_DIR, "outputs")
    scaffold_dirs = sorted(glob.glob(os.path.join(outputs_dir, "*")))

    for scaffold_dir in scaffold_dirs:
        if not os.path.isdir(scaffold_dir):
            continue
        scaffold = os.path.basename(scaffold_dir)
        sdf_files = glob.glob(os.path.join(scaffold_dir, "*.sdf"))

        scaffold_count = 0
        for sdf_path in sdf_files:
            for mol, props in read_sdf_mols(sdf_path):
                smi = Chem.MolToSmiles(mol)
                mol_counter += 1
                scaffold_count += 1
                mol_id = f"{scaffold}_{scaffold_count:04d}"

                # Plant scoring
                try:
                    r = process_compound(smi, compound_id=mol_id,
                                         use_pubchem=False, curator="diffgui")
                    scoring = r.get("scoring", {})
                    comp = scoring.get("composite_result", {})
                    gs = comp.get("gate_scores", {})
                    row_data = r.get("row", {})

                    result = {
                        "scaffold": scaffold,
                        "mol_id": mol_id,
                        "smiles": smi,
                        "cnn_affinity": "",
                        "vina_score": "",
                        "composite_phloem_score": comp.get("composite", ""),
                        "gate1": gs.get("gate1", ""),
                        "gate2": gs.get("gate2", ""),
                        "gate3": gs.get("gate3", ""),
                        "gate4": gs.get("gate4", ""),
                        "gate5": gs.get("gate5", ""),
                        "bottleneck": comp.get("bottleneck", ""),
                        "design_flags": "; ".join(comp.get("design_flags", [])),
                        "module7_class": row_data.get("module7_class", ""),
                    }
                except Exception as e:
                    result = {
                        "scaffold": scaffold, "mol_id": mol_id,
                        "smiles": smi, "cnn_affinity": "",
                        "vina_score": "",
                        "composite_phloem_score": "ERROR",
                        "gate1": "", "gate2": "", "gate3": "",
                        "gate4": "", "gate5": "",
                        "bottleneck": str(e)[:60],
                        "design_flags": "", "module7_class": "",
                    }

                # Merge GNINA scores
                key = (scaffold, smi)
                if key in gnina_scores:
                    result["cnn_affinity"] = gnina_scores[key].get("cnn_affinity", "")
                    result["vina_score"] = gnina_scores[key].get("vina_score", "")

                results.append(result)

        print(f"  {scaffold}: {scaffold_count} molecules")

    # Sort by composite phloem score (descending)
    def _sort_key(r):
        try:
            return -float(r["composite_phloem_score"])
        except (ValueError, TypeError):
            return 0
    results.sort(key=_sort_key)

    # Write CSV
    csv_path = os.path.join(OUTPUT_DIR, "diffgui_ranked_results.csv")
    fieldnames = [
        "scaffold", "mol_id", "smiles", "cnn_affinity", "vina_score",
        "composite_phloem_score", "gate1", "gate2", "gate3", "gate4", "gate5",
        "bottleneck", "design_flags", "module7_class",
    ]
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"\nTotal molecules: {mol_counter}")
    print(f"Ranked results: {csv_path}")

    # Summary stats per scaffold
    print(f"\n{'Scaffold':<20} {'Count':>6} {'Avg Phloem':>11} {'Best Phloem':>12}")
    print("-" * 55)
    from collections import defaultdict
    by_scaffold = defaultdict(list)
    for r in results:
        try:
            score = float(r["composite_phloem_score"])
            by_scaffold[r["scaffold"]].append(score)
        except (ValueError, TypeError):
            pass
    for sc in sorted(by_scaffold):
        vals = by_scaffold[sc]
        print(f"{sc:<20} {len(vals):>6} {sum(vals)/len(vals):>11.4f} {max(vals):>12.4f}")


if __name__ == "__main__":
    main()
