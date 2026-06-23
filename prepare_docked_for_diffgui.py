#!/usr/bin/env python3
"""
Prepare GNINA-docked scaffold poses for DiffGUI frag_cond generation.

Takes the docked SDF files from GNINA and generates:
1. Placed SDF files (one per pose)
2. Fragment index files (atoms to freeze)
3. DiffGUI YAML configs with agrochemical-optimized constraints
4. Slurm submission script

Usage:
  python prepare_docked_for_diffgui.py
  python prepare_docked_for_diffgui.py --poses_dir scaffold_docking/docked_poses
  python prepare_docked_for_diffgui.py --num_poses 2 --num_mols 100
"""

import argparse
import os
import yaml
from pathlib import Path

from rdkit import Chem


def extract_poses(sdf_path, max_poses=3):
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    poses = []
    for i, mol in enumerate(suppl):
        if mol is None:
            continue
        if i >= max_poses:
            break
        cnn_score = mol.GetPropsAsDict().get("CNNscore", None)
        cnn_affinity = mol.GetPropsAsDict().get("CNNaffinity", None)
        vina_score = mol.GetPropsAsDict().get("minimizedAffinity", None)
        poses.append({
            "mol": mol,
            "pose_idx": i + 1,
            "cnn_score": cnn_score,
            "cnn_affinity": cnn_affinity,
            "vina_score": vina_score,
        })
    return poses


def write_pose_sdf(mol, path):
    w = Chem.SDWriter(str(path))
    mol_noH = Chem.RemoveHs(mol)
    w.write(mol_noH)
    w.close()
    return mol_noH.GetNumAtoms()


def write_frag_indices(num_atoms, path):
    indices = list(range(num_atoms))
    with open(path, "w") as f:
        f.write(" ".join(str(i) for i in indices) + "\n")


def make_config(scaffold_name, pose_idx, num_mols, pocket_path,
                sdf_path, frag_path):
    cfg = {
        "model": {
            "checkpoint": "ckpt/trained.pt",
            "target": str(pocket_path),
            "ligand": str(sdf_path),
            "frag": str(frag_path),
            "gen_mode": "frag_cond",
            "logp": 1.50,
            "tpsa": 70,
            "sa": 3.00,
            "qed": 0.60,
            "aff": 10.00,
        },
        "bond_predictor": "ckpt/bond_trained.pt",
        "sample": {
            "seed": 2024 + pose_idx,
            "batch_size": 2,
            "num_mols": num_mols,
            "save_traj_prob": 0.0,
            "sample": True,
            "sample_method": "priori",
            "mode": "pocket",
            "add_edge": None,
            "gui_strength": 3.0,
            "guidance": ["uncertainty", 1.0e-4],
        },
        "data": {
            "name": "protein_ligand",
            "path": "data/PDBbind_v2020_pocket10",
            "split": "data/PDBbind_pocket10_split.pt",
            "protein_root": "data/PDBbind_v2020",
            "dataset": "pdbbind",
            "transform": {
                "ligand_atom_mode": "aromatic",
                "random_rot": False,
                "sample": False,
            },
        },
    }
    return cfg


def make_slurm_script(run_names, run_dir):
    n = len(run_names)
    script = f"""#!/bin/bash
#SBATCH --job-name=diffgui_cnn
#SBATCH --account=iit135
#SBATCH --partition=gpu-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=1-{n}
#SBATCH --output={run_dir}/logs/diffgui_%A_%a.out
#SBATCH --error={run_dir}/logs/diffgui_%A_%a.err

module purge
module load cpu/0.17.3b
module load gpu/0.17.3b
module load anaconda3/2021.05
conda activate diffgui

WORK_DIR="/expanse/lustre/scratch/$USER/temp_project/DiffGui"
cd "$WORK_DIR"

NAMES=({' '.join(run_names)})
NAME=${{NAMES[$SLURM_ARRAY_TASK_ID-1]}}
CONFIG="{run_dir}/configs/${{NAME}}_sample.yml"
OUTDIR="{run_dir}/outputs/${{NAME}}"

mkdir -p "$OUTDIR"

echo "=== DiffGUI frag_cond: $NAME ==="
echo "Config: $CONFIG"
echo "Output: $OUTDIR"

python scripts/sample.py \\
    --outdir "$OUTDIR" \\
    --config "$CONFIG" \\
    --device cuda:0

echo "=== Done: $NAME ==="
"""
    return script


def main():
    parser = argparse.ArgumentParser(
        description="Prepare GNINA-docked scaffolds for DiffGUI")
    parser.add_argument("--poses_dir", default="scaffold_docking/docked_poses",
                        help="Directory with GNINA docked SDF files")
    parser.add_argument("--pocket", default="2ZSH_pocket10.pdb",
                        help="Pocket PDB file")
    parser.add_argument("--num_poses", type=int, default=3,
                        help="Max poses per scaffold to use")
    parser.add_argument("--num_mols", type=int, default=200,
                        help="Molecules to generate per pose")
    parser.add_argument("--output_dir", default="diffgui_cnn_run",
                        help="Output directory for DiffGUI files")
    args = parser.parse_args()

    poses_dir = Path(args.poses_dir)
    out_dir = Path(args.output_dir)

    for subdir in ["scaffolds", "configs", "logs", "outputs"]:
        (out_dir / subdir).mkdir(parents=True, exist_ok=True)

    if not poses_dir.exists():
        print(f"ERROR: {poses_dir} not found.")
        print("Run dock_scaffolds_gnina.sh first to generate docked poses.")
        return

    all_run_names = []
    all_results = []

    sdf_files = sorted(poses_dir.glob("*_docked.sdf"))
    if not sdf_files:
        print(f"No *_docked.sdf files found in {poses_dir}")
        return

    print(f"{'Scaffold':<18} {'Pose':>4} {'CNN':>6} {'Affinity':>8} "
          f"{'Vina':>6} {'Atoms':>5}")
    print("=" * 60)

    for sdf_path in sdf_files:
        scaffold_name = sdf_path.stem.replace("_docked", "")
        poses = extract_poses(sdf_path, args.num_poses)

        for pose in poses:
            pidx = pose["pose_idx"]
            run_name = f"{scaffold_name}_pose{pidx}"

            sdf_out = out_dir / "scaffolds" / f"{run_name}.sdf"
            frag_out = out_dir / "scaffolds" / f"{run_name}_frag.txt"
            cfg_out = out_dir / "configs" / f"{run_name}_sample.yml"

            num_atoms = write_pose_sdf(pose["mol"], sdf_out)
            write_frag_indices(num_atoms, frag_out)

            cfg = make_config(scaffold_name, pidx, args.num_mols,
                              args.pocket, sdf_out, frag_out)
            with open(cfg_out, "w") as f:
                yaml.dump(cfg, f, default_flow_style=False, sort_keys=False)

            all_run_names.append(run_name)

            cnn = f"{pose['cnn_score']:.3f}" if pose['cnn_score'] else "---"
            aff = f"{pose['cnn_affinity']:.2f}" if pose['cnn_affinity'] else "---"
            vin = f"{pose['vina_score']:.2f}" if pose['vina_score'] else "---"
            print(f"{scaffold_name:<18} {pidx:>4} {cnn:>6} {aff:>8} "
                  f"{vin:>6} {num_atoms:>5}")

            all_results.append({
                "scaffold": scaffold_name, "pose": pidx,
                "run_name": run_name, "atoms": num_atoms,
            })

    with open(out_dir / "run_names.txt", "w") as f:
        for rn in all_run_names:
            f.write(rn + "\n")

    slurm = make_slurm_script(all_run_names, str(out_dir))
    slurm_path = out_dir / "submit_all.sh"
    with open(slurm_path, "w") as f:
        f.write(slurm)
    os.chmod(slurm_path, 0o755)

    print(f"\n{'=' * 60}")
    print(f"Prepared {len(all_run_names)} DiffGUI runs")
    print(f"  Scaffolds: {len(sdf_files)}")
    print(f"  Poses per scaffold: up to {args.num_poses}")
    print(f"  Molecules per pose: {args.num_mols}")
    print(f"  Total candidates: {len(all_run_names) * args.num_mols}")
    print(f"\nOutput: {out_dir}/")
    print(f"  Configs: {out_dir}/configs/")
    print(f"  Scaffolds: {out_dir}/scaffolds/")
    print(f"  Submit: sbatch {slurm_path}")
    print(f"\nProperty constraints (agrochemical-optimized):")
    print(f"  logP target: 1.5 (window: 0.5-2.5)")
    print(f"  TPSA target: 70 (window: 40-100)")
    print(f"  SA target: 3.0")
    print(f"  QED target: 0.6")
    print(f"  Affinity target: 10.0")
    print(f"  MW filter: <=420 Da (applied post-generation)")


if __name__ == "__main__":
    main()
