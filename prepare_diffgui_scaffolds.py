#!/usr/bin/env python3
"""
DiffGUI Scaffold Placement & Submission Pipeline
=================================================
Prepares 15 MCS scaffolds for fragment-conditioned de novo ligand
generation in the GID1 gibberellin receptor (PDB 2ZSH, 1.8 A).

Outputs:
  - 3D-placed scaffold SDFs at the GA3 binding centroid
  - Fragment index files (atoms to freeze during diffusion)
  - Per-scaffold DiffGUI YAML configs (frag_cond mode)
  - Slurm array job for SDSC Expanse (gpu-shared, account iit135)
  - GNINA CNN rescoring script
  - Results collector that merges binding + plant mobility scores
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms, rdmolops

# ── 15 MCS Scaffolds ─────────────────────────────────────────────────

SCAFFOLDS = {
    "quinoline":      "c1ccc2ncccc2c1",
    "isoquinoline":   "c1ccc2cnccc2c1",
    "naphthalene":    "c1ccc2ccccc2c1",
    "benzothiazole":  "c1ccc2ncsc2c1",
    "benzoxazole":    "c1ccc2ncoc2c1",
    "phthalimide":    "O=C1NC(=O)c2ccccc21",
    "saccharin":      "O=C1NS(=O)(=O)c2ccccc21",
    "quinazolinone":  "O=c1[nH]cnc2ccccc12",
    "coumarin":       "O=c1ccc2ccccc2o1",
    "phthalazinone":  "O=c1[nH]ncc2ccccc12",
    "oxadiazole_124": "c1nonn1",
    "oxadiazole_134": "c1nnoc1",
    "pyrimidine":     "c1ccncn1",
    "triazine":       "c1ncncn1",
    "carbazole":      "c1ccc2c(c1)[nH]c1ccccc12",
}

GA3_CENTROID = np.array([51.033, 59.452, 37.370])


def embed_and_optimise(mol):
    mol_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    status = AllChem.EmbedMolecule(mol_h, params)
    if status != 0:
        params.useRandomCoords = True
        status = AllChem.EmbedMolecule(mol_h, params)
    if status != 0:
        return None
    try:
        AllChem.MMFFOptimizeMolecule(mol_h, maxIters=500)
    except Exception:
        pass
    return Chem.RemoveHs(mol_h)


def get_centroid(mol):
    conf = mol.GetConformer()
    coords = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
    return coords.mean(axis=0)


def translate_to(mol, target):
    current = get_centroid(mol)
    shift = target - current
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (pos.x + shift[0],
                                  pos.y + shift[1],
                                  pos.z + shift[2]))
    return mol


def place_scaffold(name, smiles, centroid):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"  WARNING: invalid SMILES for {name}: {smiles}")
        return None
    mol = embed_and_optimise(mol)
    if mol is None:
        print(f"  WARNING: embedding failed for {name}")
        return None
    mol = translate_to(mol, centroid)
    mol.SetProp("_Name", name)
    return mol


def write_sdf(mol, path):
    writer = Chem.SDWriter(str(path))
    writer.write(mol)
    writer.close()


def write_frag_indices(mol, path):
    indices = list(range(mol.GetNumAtoms()))
    with open(path, "w") as f:
        f.write(" ".join(str(i) for i in indices) + "\n")


def write_diffgui_config(name, scaffold_dir, pocket_path, num_mols, config_path):
    rel_pocket = pocket_path
    rel_sdf = f"{scaffold_dir}/{name}_placed.sdf"
    rel_frag = f"{scaffold_dir}/{name}_frag_indices.txt"

    config = f"""model:
  checkpoint: ckpt/trained.pt
  target: {rel_pocket}
  ligand: {rel_sdf}
  frag: {rel_frag}
  gen_mode: frag_cond
  logp: 2.00
  tpsa: 80
  sa: 3.00
  qed: 0.60
  aff: 10.00

bond_predictor: ckpt/bond_trained.pt

sample:
  seed: 2024
  batch_size: 2
  num_mols: {num_mols}
  save_traj_prob: 0.0
  sample: True
  sample_method: priori
  mode: pocket
  add_edge: None
  gui_strength: 3.0
  guidance:
    - uncertainty
    - 1.e-4

data:
  name: protein_ligand
  path: data/PDBbind_v2020_pocket10
  split: data/PDBbind_pocket10_split.pt
  protein_root: data/PDBbind_v2020
  dataset: pdbbind
  transform:
    ligand_atom_mode: aromatic
    random_rot: False
    sample: False
"""
    with open(config_path, "w") as f:
        f.write(config)


def write_slurm_array(output_dir, scaffold_names, account="iit135"):
    n = len(scaffold_names)
    names_file = Path(output_dir) / "scaffold_names.txt"
    with open(names_file, "w") as f:
        for name in scaffold_names:
            f.write(name + "\n")

    script = f"""#!/bin/bash
#SBATCH --job-name=diffgui_frag
#SBATCH --account={account}
#SBATCH --partition=gpu-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=1-{n}
#SBATCH --output={output_dir}/logs/diffgui_%A_%a.out
#SBATCH --error={output_dir}/logs/diffgui_%A_%a.err

module purge
module load cpu/0.17.3b
module load gpu/0.17.3b
module load anaconda3/2021.05
conda activate diffgui

DIFFGUI_DIR="/expanse/lustre/scratch/$USER/temp_project/DiffGui"
cd "$DIFFGUI_DIR"

SCAFFOLD=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {names_file})
CONFIG="{output_dir}/configs/${{SCAFFOLD}}_sample.yml"
OUTDIR="{output_dir}/outputs/${{SCAFFOLD}}"

echo "=== DiffGUI frag_cond: $SCAFFOLD ==="
echo "Config: $CONFIG"
echo "Output: $OUTDIR"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null)"
echo ""

mkdir -p "$OUTDIR"

python scripts/sample.py \\
    --outdir "$OUTDIR" \\
    --config "$CONFIG" \\
    --device cuda:0

echo ""
echo "=== Done: $SCAFFOLD ==="
ls -la "$OUTDIR/" 2>/dev/null
"""
    script_path = Path(output_dir) / "submit_diffgui_all.sh"
    with open(script_path, "w") as f:
        f.write(script)
    os.chmod(script_path, 0o755)
    return script_path


def write_gnina_rescore(output_dir, pocket_pdb, protein_pdb="2ZSH.pdb"):
    script = f"""#!/bin/bash
#SBATCH --job-name=gnina_rescore
#SBATCH --account=iit135
#SBATCH --partition=gpu-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output={output_dir}/logs/gnina_rescore_%j.out
#SBATCH --error={output_dir}/logs/gnina_rescore_%j.err

# GNINA CNN + Vina rescoring of DiffGUI outputs against 2ZSH
# Requires gnina installed in the conda environment or available as module

module purge
module load cpu/0.17.3b
module load gpu/0.17.3b
module load anaconda3/2021.05
conda activate diffgui

OUTDIR="{output_dir}"
RECEPTOR="{protein_pdb}"

echo "=== GNINA CNN Rescoring ==="
echo "Receptor: $RECEPTOR"
echo ""

RESCORE_DIR="$OUTDIR/gnina_rescores"
mkdir -p "$RESCORE_DIR"

for SCAFFOLD_DIR in "$OUTDIR"/outputs/*/; do
    SCAFFOLD=$(basename "$SCAFFOLD_DIR")
    echo "Rescoring: $SCAFFOLD"

    # Find all generated SDF files in the scaffold output
    for SDF in "$SCAFFOLD_DIR"/*.sdf; do
        [ -f "$SDF" ] || continue
        BASENAME=$(basename "$SDF" .sdf)
        OUT_SDF="$RESCORE_DIR/${{SCAFFOLD}}_${{BASENAME}}_rescored.sdf"

        gnina \\
            -r "$RECEPTOR" \\
            -l "$SDF" \\
            --score_only \\
            --cnn_scoring \\
            -o "$OUT_SDF" \\
            --seed 42 \\
            2>/dev/null

        if [ -f "$OUT_SDF" ]; then
            echo "  $BASENAME: done"
        else
            echo "  $BASENAME: FAILED"
        fi
    done
done

echo ""
echo "=== Rescoring Complete ==="
echo "Results in: $RESCORE_DIR/"
ls "$RESCORE_DIR/"*.sdf 2>/dev/null | wc -l
echo "SDF files generated"
"""
    script_path = Path(output_dir) / "rescore_gnina.sh"
    with open(script_path, "w") as f:
        f.write(script)
    os.chmod(script_path, 0o755)
    return script_path


def write_results_collector(output_dir):
    script = '''#!/usr/bin/env python3
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


OUTPUT_DIR = "''' + output_dir + '''"


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

    print(f"\\nTotal molecules: {mol_counter}")
    print(f"Ranked results: {csv_path}")

    # Summary stats per scaffold
    print(f"\\n{'Scaffold':<20} {'Count':>6} {'Avg Phloem':>11} {'Best Phloem':>12}")
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
'''
    script_path = Path(output_dir) / "collect_results.py"
    with open(script_path, "w") as f:
        f.write(script)
    os.chmod(script_path, 0o755)
    return script_path


# ── Main ──────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Prepare DiffGUI scaffold placement and submission pipeline"
    )
    parser.add_argument("--pocket", default="2ZSH_pocket10.pdb",
                        help="Path to pocket PDB")
    parser.add_argument("--centroid", nargs=3, type=float,
                        default=[51.033, 59.452, 37.370],
                        help="Binding site centroid x y z")
    parser.add_argument("--num_mols", type=int, default=200,
                        help="Number of molecules per scaffold")
    parser.add_argument("--output_dir", default="diffgui_run",
                        help="Base output directory")
    parser.add_argument("--protein_pdb", default="2ZSH.pdb",
                        help="Full protein PDB for GNINA rescoring")
    args = parser.parse_args()

    centroid = np.array(args.centroid)
    output_dir = args.output_dir

    scaffold_dir = os.path.join(output_dir, "scaffolds")
    config_dir = os.path.join(output_dir, "configs")
    log_dir = os.path.join(output_dir, "logs")
    outputs_dir = os.path.join(output_dir, "outputs")

    for d in [scaffold_dir, config_dir, log_dir, outputs_dir]:
        os.makedirs(d, exist_ok=True)

    print("=" * 70)
    print("DiffGUI Scaffold Placement Pipeline")
    print(f"  Target: GID1 (PDB 2ZSH)")
    print(f"  Pocket: {args.pocket}")
    print(f"  Centroid: ({centroid[0]:.3f}, {centroid[1]:.3f}, {centroid[2]:.3f})")
    print(f"  Scaffolds: {len(SCAFFOLDS)}")
    print(f"  Mols/scaffold: {args.num_mols}")
    print(f"  Output: {output_dir}/")
    print("=" * 70)

    # 1. Generate and place scaffolds
    print("\n--- Placing scaffolds at GA3 centroid ---")
    placed = []
    for name, smiles in SCAFFOLDS.items():
        mol = place_scaffold(name, smiles, centroid)
        if mol is None:
            continue

        sdf_path = os.path.join(scaffold_dir, f"{name}_placed.sdf")
        frag_path = os.path.join(scaffold_dir, f"{name}_frag_indices.txt")

        write_sdf(mol, sdf_path)
        write_frag_indices(mol, frag_path)

        n_atoms = mol.GetNumAtoms()
        final_centroid = get_centroid(mol)
        print(f"  {name:<18} {n_atoms:>3} atoms  "
              f"center=({final_centroid[0]:.1f}, {final_centroid[1]:.1f}, "
              f"{final_centroid[2]:.1f})  -> {sdf_path}")
        placed.append(name)

    # 2. Write DiffGUI configs
    print(f"\n--- Writing DiffGUI configs ({len(placed)}) ---")
    for name in placed:
        cfg_path = os.path.join(config_dir, f"{name}_sample.yml")
        write_diffgui_config(name, scaffold_dir, args.pocket,
                             args.num_mols, cfg_path)
        print(f"  {cfg_path}")

    # 3. Write Slurm array script
    print("\n--- Writing Slurm submission script ---")
    slurm_path = write_slurm_array(output_dir, placed)
    print(f"  {slurm_path}")

    # 4. Write GNINA rescore script
    print("\n--- Writing GNINA rescoring script ---")
    gnina_path = write_gnina_rescore(output_dir, args.pocket, args.protein_pdb)
    print(f"  {gnina_path}")

    # 5. Write results collector
    print("\n--- Writing results collector ---")
    collector_path = write_results_collector(output_dir)
    print(f"  {collector_path}")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print(f"  Scaffolds placed: {len(placed)}/{len(SCAFFOLDS)}")
    print(f"  Expected ligands: {len(placed) * args.num_mols}")
    print()
    print("TO RUN ON EXPANSE:")
    print(f"  1. Copy {output_dir}/ to Expanse DiffGui/data/ directory")
    print(f"  2. Ensure checkpoints are at DiffGui/ckpt/")
    print(f"  3. sbatch {slurm_path}")
    print(f"  4. After generation: sbatch {gnina_path}")
    print(f"  5. After rescoring: python {collector_path}")
    print("=" * 70)


if __name__ == "__main__":
    main()
