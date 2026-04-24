#!/usr/bin/env python3
"""
Submit AutoGrow4 for the next generation.

Usage:
    python submit_next_gen.py                          # Continue with current config, one more gen
    python submit_next_gen.py --crossovers 2500 --mutations 2500 --elite 50
    python submit_next_gen.py --target-gen 5           # Run up to generation 5
    python submit_next_gen.py --walltime 48:00:00      # Set walltime
    python submit_next_gen.py --dry-run                # Preview without submitting

This script:
1. Detects the last completed generation in the output folder
2. Updates the config JSON with start_a_new_run=False and new parameters
3. Cleans up any incomplete generation directories from timed-out jobs
4. Generates and submits a SLURM script for the next generation(s)
"""

import argparse
import json
import os
import glob
import subprocess
import sys


PROJECT_DIR = os.path.expanduser("~/final_project")
CONFIG_PATH = os.path.join(PROJECT_DIR, "autogrow_config.json")
AUTOGROW_DIR = os.path.join(PROJECT_DIR, "autogrow4")


def find_last_completed_gen(output_dir):
    """Find the last generation with a ranked.smi file."""
    runs = sorted(glob.glob(os.path.join(output_dir, "Run_*")))
    if not runs:
        return None, None
    # Filter out Failed runs
    valid_runs = [r for r in runs if "Failed" not in os.path.basename(r)]
    if not valid_runs:
        return None, None
    last_run = valid_runs[-1]
    last_completed = None
    gens = sorted(glob.glob(os.path.join(last_run, "generation_*")))
    for g in gens:
        gen_name = os.path.basename(g)
        if "Failed" in gen_name:
            continue
        ranked = glob.glob(os.path.join(g, "*ranked*.smi"))
        if ranked and os.path.getsize(ranked[0]) > 0:
            gen_num = int(gen_name.replace("generation_", ""))
            last_completed = gen_num
    return last_run, last_completed


def clean_incomplete_generations(run_dir, last_completed):
    """Rename incomplete generation directories (no ranked.smi) to *_Failed_N."""
    if last_completed is None:
        return
    gens = sorted(glob.glob(os.path.join(run_dir, "generation_*")))
    for g in gens:
        gen_name = os.path.basename(g)
        if "Failed" in gen_name:
            continue
        try:
            gen_num = int(gen_name.replace("generation_", ""))
        except ValueError:
            continue
        if gen_num <= last_completed:
            continue
        ranked = glob.glob(os.path.join(g, "*ranked*.smi"))
        has_ranked = ranked and os.path.getsize(ranked[0]) > 0
        if not has_ranked:
            failed_name = g + "_Failed_0"
            i = 0
            while os.path.exists(failed_name):
                i += 1
                failed_name = g + f"_Failed_{i}"
            print(f"  Renaming incomplete {gen_name} -> {os.path.basename(failed_name)}")
            os.rename(g, failed_name)


def main():
    parser = argparse.ArgumentParser(
        description="Submit AutoGrow4 for next generation(s)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Continue one more generation with 2500/2500/50 population:
  python submit_next_gen.py --crossovers 2500 --mutations 2500 --elite 50

  # Run generations 3 through 5 in one job:
  python submit_next_gen.py --target-gen 5 --crossovers 2500 --mutations 2500 --elite 50

  # Preview what would happen (no submission):
  python submit_next_gen.py --dry-run --target-gen 5
        """
    )
    parser.add_argument("--crossovers", type=int, help="Number of crossovers")
    parser.add_argument("--mutations", type=int, help="Number of mutations")
    parser.add_argument("--elite", type=int, help="Number of elite to advance")
    parser.add_argument("--top-seed", type=int, help="Top mols to seed next gen")
    parser.add_argument("--target-gen", type=int, help="Run up to this generation number")
    parser.add_argument("--walltime", default="48:00:00", help="SLURM walltime (default: 48:00:00)")
    parser.add_argument("--nodes", type=int, default=1, help="Number of nodes (default: 1)")
    parser.add_argument("--ntasks", type=int, default=90, help="Tasks per node (default: 90)")
    parser.add_argument("--mem", default="180G", help="Memory (default: 180G)")
    parser.add_argument("--account", default="iit135", help="SLURM account")
    parser.add_argument("--partition", default="shared", help="SLURM partition")
    parser.add_argument("--dry-run", action="store_true", help="Preview without submitting")
    args = parser.parse_args()

    # Load current config
    if not os.path.exists(CONFIG_PATH):
        print(f"ERROR: Config not found: {CONFIG_PATH}")
        print("Run the notebook first to generate a config.")
        sys.exit(1)

    with open(CONFIG_PATH) as f:
        config = json.load(f)

    output_dir = config["root_output_folder"]
    print(f"Output directory: {output_dir}")

    # Find last completed generation
    run_dir, last_gen = find_last_completed_gen(output_dir)
    if last_gen is not None:
        print(f"Last completed generation: {last_gen} (in {os.path.basename(run_dir)})")
        next_gen = last_gen + 1
        # Clean up incomplete generation dirs from timed-out jobs
        clean_incomplete_generations(run_dir, last_gen)
    else:
        print("No completed generations found. Starting from generation 0.")
        next_gen = 0

    # Determine target generation
    if args.target_gen is not None:
        target_gen = args.target_gen
    else:
        target_gen = next_gen  # Just one more generation

    if last_gen is not None and target_gen <= last_gen:
        print(f"ERROR: Target generation {target_gen} already completed (last={last_gen}).")
        sys.exit(1)

    print(f"Will run generation(s) {next_gen} through {target_gen}")

    # Update config for continuation
    config["start_a_new_run"] = False
    config["number_of_generations"] = target_gen

    if args.crossovers is not None:
        config["number_of_crossovers"] = args.crossovers
    if args.mutations is not None:
        config["number_of_mutants"] = args.mutations
    if args.elite is not None:
        config["number_elitism_advance_from_previous_gen"] = args.elite
    if args.top_seed is not None:
        config["top_mols_to_seed_next_generation"] = args.top_seed

    print(f"  Crossovers: {config['number_of_crossovers']}")
    print(f"  Mutations: {config['number_of_mutants']}")
    print(f"  Elite: {config['number_elitism_advance_from_previous_gen']}")
    print(f"  Top seed: {config['top_mols_to_seed_next_generation']}")

    # Save updated config
    with open(CONFIG_PATH, "w") as f:
        json.dump(config, f, indent=4)
    print(f"Config updated: {CONFIG_PATH}")

    # Generate SLURM script
    slurm_lines = [
        "#!/bin/bash",
        f"#SBATCH --account={args.account}",
        f"#SBATCH --partition={args.partition}",
        f"#SBATCH --nodes={args.nodes}",
        f"#SBATCH --ntasks-per-node={args.ntasks}",
        f"#SBATCH --mem={args.mem}",
        f"#SBATCH --time={args.walltime}",
        f'#SBATCH --job-name="ag4_gen{next_gen}-{target_gen}"',
        '#SBATCH --output="autogrow4.%j.%N.out"',
        "#SBATCH --export=ALL",
        "#SBATCH --requeue",
        "#SBATCH --mail-user=mschnur@hawk.illinoistech",
        "#SBATCH --mail-type=ALL",
        "",
        "# Environment setup",
        "module purge",
        "module load cpu",
        "",
        "source ~/miniconda3/etc/profile.d/conda.sh",
        "conda activate ag4_full",
        "",
        f"cd {AUTOGROW_DIR}",
        "",
        "# Cache prerun (single processor, prevents race conditions)",
        f"python run_autogrow.py -c",
        "",
        "# Run AutoGrow4",
        f"python run_autogrow.py -j {CONFIG_PATH}",
        "PYTHON_EXIT=$?",
        "",
        "if [ $PYTHON_EXIT -ne 0 ]; then",
        '    echo "AutoGrow4 FAILED with exit code $PYTHON_EXIT"',
        "    exit $PYTHON_EXIT",
        "fi",
        f'echo "AutoGrow4 completed successfully (gen {next_gen}-{target_gen})."',
    ]
    slurm_script = "\n".join(slurm_lines) + "\n"

    slurm_path = os.path.join(PROJECT_DIR, "autogrow_slurm.sh")
    with open(slurm_path, "w") as f:
        f.write(slurm_script)
    print(f"SLURM script: {slurm_path}")

    if args.dry_run:
        print("\n=== DRY RUN - SLURM script preview ===")
        print(slurm_script)
        print("(Not submitted. Remove --dry-run to submit.)")
    else:
        result = subprocess.run(["sbatch", slurm_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        if result.returncode == 0:
            print(f"\n{result.stdout.strip()}")
            print("Job submitted! Monitor with: sacct -j <JOBID> --format=JobID,State,Elapsed")
        else:
            print(f"ERROR submitting job: {result.stderr.strip()}")
            sys.exit(1)


if __name__ == "__main__":
    main()
