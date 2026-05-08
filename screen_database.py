#!/usr/bin/env python3
"""
Virtual screening: dock a chunk of SMILES against a receptor using AutoDock Vina.

Usage (standalone):
  python3 screen_database.py --smi database.smi --receptor receptor.pdbqt \
    --center 55.0 59.5 35.5 --size 21.0 21.0 22.5 --exhaustiveness 20 \
    --start 0 --count 500 --output results_chunk.tsv

Usage (Slurm array - called by submit_screening.sh):
  Each array task processes a chunk defined by SLURM_ARRAY_TASK_ID.

Workflow per compound:
  1. SMILES -> 3D SDF (OpenBabel)
  2. SDF -> PDBQT (OpenBabel)
  3. PDBQT -> Vina docking -> score
  4. Write result to TSV
"""
import argparse
import os
import subprocess
import sys
import tempfile
import shutil


def find_vina():
    """Find the Vina executable."""
    home = os.path.expanduser("~")
    candidates = [
        os.path.join(home, "final_project", "autogrow4", "autogrow", "docking", "docking_executables", "vina"),
        shutil.which("vina"),
        shutil.which("vina_1.2.5_linux_x86_64"),
        "/usr/local/bin/vina",
    ]
    for c in candidates:
        if c and os.path.isfile(c):
            if not os.access(c, os.X_OK):
                try:
                    os.chmod(c, 0o755)
                    print(f"  Fixed permissions on {c}")
                except OSError:
                    pass
            if os.access(c, os.X_OK):
                return c
    for c in candidates:
        if c and os.path.isfile(c):
            return c
    return None


def find_obabel():
    """Find OpenBabel executable."""
    ob = shutil.which("obabel")
    if ob:
        return ob
    for p in ["/usr/bin/obabel", "/usr/local/bin/obabel",
              os.path.expanduser("~/miniconda3/bin/obabel"),
              os.path.expanduser("~/.conda/envs/autogrow/bin/obabel")]:
        if os.path.isfile(p):
            return p
    return None


def smiles_to_pdbqt(smiles, obabel_path, tmp_dir):
    """Convert SMILES to 3D PDBQT using OpenBabel."""
    sdf_path = os.path.join(tmp_dir, "mol.sdf")
    pdbqt_path = os.path.join(tmp_dir, "mol.pdbqt")

    result = subprocess.run(
        [obabel_path, "-:" + smiles, "-osdf", "-O", sdf_path,
         "--gen3d", "--best", "-h"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, timeout=60
    )
    if not os.path.exists(sdf_path) or os.path.getsize(sdf_path) < 10:
        return None

    result = subprocess.run(
        [obabel_path, sdf_path, "-opdbqt", "-O", pdbqt_path, "-h"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, timeout=30
    )
    if not os.path.exists(pdbqt_path) or os.path.getsize(pdbqt_path) < 10:
        return None

    return pdbqt_path


def run_vina(vina_path, receptor_path, ligand_pdbqt, center, size, exhaustiveness, tmp_dir):
    """Run Vina docking and return best score."""
    out_path = os.path.join(tmp_dir, "docked.pdbqt")
    log_path = os.path.join(tmp_dir, "vina.log")

    cmd = [
        vina_path,
        "--receptor", receptor_path,
        "--ligand", ligand_pdbqt,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(size[0]),
        "--size_y", str(size[1]),
        "--size_z", str(size[2]),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", "1",
        "--out", out_path,
        "--log", log_path,
    ]

    try:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, timeout=300)
    except subprocess.TimeoutExpired:
        return None, None

    score = None
    if os.path.exists(log_path):
        with open(log_path) as f:
            for line in f:
                line = line.strip()
                if line.startswith("1") or (len(line) > 3 and line[0:4].strip().isdigit()):
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            mode = int(parts[0])
                            if mode == 1:
                                score = float(parts[1])
                                break
                        except (ValueError, IndexError):
                            continue

    docked_path = out_path if os.path.exists(out_path) else None
    return score, docked_path


def load_smiles(smi_path, start, count):
    """Load a chunk of SMILES from a .smi file."""
    compounds = []
    with open(smi_path) as f:
        for i, line in enumerate(f):
            if i < start:
                continue
            if len(compounds) >= count:
                break
            parts = line.strip().split()
            if len(parts) >= 2:
                compounds.append((parts[0], parts[1]))
            elif len(parts) == 1:
                compounds.append((parts[0], f"mol_{i}"))
    return compounds


def main():
    parser = argparse.ArgumentParser(description="Virtual screening with AutoDock Vina")
    parser.add_argument("--smi", required=True, help="Input SMILES file")
    parser.add_argument("--receptor", required=True, help="Receptor PDBQT file")
    parser.add_argument("--center", nargs=3, type=float, required=True, help="Box center (x y z)")
    parser.add_argument("--size", nargs=3, type=float, required=True, help="Box size (x y z)")
    parser.add_argument("--exhaustiveness", type=int, default=20)
    parser.add_argument("--start", type=int, default=0, help="Start index in SMILES file")
    parser.add_argument("--count", type=int, default=500, help="Number of compounds to process")
    parser.add_argument("--output", required=True, help="Output TSV file")
    parser.add_argument("--save-poses", help="Directory to save docked poses (top scorers)")
    parser.add_argument("--save-threshold", type=float, default=-7.0, help="Save poses below this score")
    args = parser.parse_args()

    vina_path = find_vina()
    obabel_path = find_obabel()

    if not vina_path:
        print("ERROR: Could not find Vina executable")
        sys.exit(1)
    if not obabel_path:
        print("ERROR: Could not find OpenBabel (obabel)")
        sys.exit(1)

    print(f"Vina: {vina_path}")
    print(f"OpenBabel: {obabel_path}")
    print(f"Receptor: {args.receptor}")
    print(f"Box center: {args.center}")
    print(f"Box size: {args.size}")
    print(f"Exhaustiveness: {args.exhaustiveness}")

    compounds = load_smiles(args.smi, args.start, args.count)
    print(f"Loaded {len(compounds)} compounds (start={args.start})")

    if args.save_poses:
        os.makedirs(args.save_poses, exist_ok=True)

    with open(args.output, "w") as fout:
        fout.write("name\tsmiles\tscore\tstatus\n")

        for idx, (smiles, name) in enumerate(compounds):
            tmp_dir = tempfile.mkdtemp(prefix="vscreen_")
            try:
                pdbqt = smiles_to_pdbqt(smiles, obabel_path, tmp_dir)
                if pdbqt is None:
                    fout.write(f"{name}\t{smiles}\t\tconversion_failed\n")
                    if (idx + 1) % 50 == 0:
                        print(f"  [{idx+1}/{len(compounds)}] {name}: conversion failed")
                    continue

                score, docked = run_vina(
                    vina_path, args.receptor, pdbqt,
                    args.center, args.size, args.exhaustiveness, tmp_dir
                )

                if score is not None:
                    fout.write(f"{name}\t{smiles}\t{score:.2f}\tok\n")
                    if args.save_poses and docked and score < args.save_threshold:
                        pose_dest = os.path.join(args.save_poses, f"{name}_docked.pdbqt")
                        shutil.copy2(docked, pose_dest)
                else:
                    fout.write(f"{name}\t{smiles}\t\tdocking_failed\n")

                if (idx + 1) % 50 == 0:
                    print(f"  [{idx+1}/{len(compounds)}] {name}: {score if score else 'failed'}")
                fout.flush()

            except Exception as e:
                fout.write(f"{name}\t{smiles}\t\terror: {str(e)[:50]}\n")
            finally:
                shutil.rmtree(tmp_dir, ignore_errors=True)

    print(f"Done. Results written to {args.output}")


if __name__ == "__main__":
    main()
