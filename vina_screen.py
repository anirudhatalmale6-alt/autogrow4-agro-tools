#!/usr/bin/env python3
"""
Virtual Screen: Dock a SMILES database against a receptor with QuickVina2/Vina.
Parallelized with multiprocessing for HPC use.

Uses RDKit for 3D generation and built-in PDBQT writer (no obabel/MGLTools needed).

Usage:
  python vina_screen.py -i library.smi -r receptor.pdbqt -o results
  python vina_screen.py -i library.smi -r receptor.pdbqt -o results --resume
"""

import argparse
import multiprocessing
import os
import shutil
import subprocess
import sys
import time

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
except ImportError:
    print("ERROR: RDKit not found. Activate conda environment first.")
    sys.exit(1)

AUTOGROW_DIR = os.path.expanduser("~/final_project/autogrow4")

QVINA2_PATH = os.path.join(
    AUTOGROW_DIR, "autogrow", "docking", "docking_executables",
    "q_vina_2", "q_vina_2_1_linux", "qvina2.1")
VINA_PATH = os.path.join(
    AUTOGROW_DIR, "autogrow", "docking", "docking_executables",
    "vina", "autodock_vina_1_1_2_linux_x86", "bin", "vina")


AD4_ATOM_TYPES = {
    (6, False): "C",
    (6, True): "A",
    (7, False): "NA",
    (7, True): "NA",
    (8, False): "OA",
    (8, True): "OA",
    (9, False): "F",
    (15, False): "P",
    (16, False): "SA",
    (16, True): "SA",
    (17, False): "Cl",
    (35, False): "Br",
    (53, False): "I",
    (14, False): "Si",
    (5, False): "B",
}


def get_ad4_type(atom):
    """Get AutoDock4 atom type for an RDKit atom."""
    num = atom.GetAtomicNum()
    if num == 1:
        for nb in atom.GetNeighbors():
            if nb.GetAtomicNum() in (7, 8, 16):
                return "HD"
        return "H"
    arom = atom.GetIsAromatic()
    return AD4_ATOM_TYPES.get((num, arom), atom.GetSymbol())


def mol_to_pdbqt(mol, path):
    """Write RDKit mol to PDBQT format with Gasteiger charges."""
    AllChem.ComputeGasteigerCharges(mol)
    conf = mol.GetConformer()

    lines = ["ROOT"]
    atom_idx = 0
    for i, atom in enumerate(mol.GetAtoms()):
        pos = conf.GetAtomPosition(i)
        ad4_type = get_ad4_type(atom)
        try:
            charge = float(atom.GetDoubleProp("_GasteigerCharge"))
            if charge != charge:  # NaN check
                charge = 0.0
        except Exception:
            charge = 0.0

        atom_idx += 1
        name = atom.GetSymbol()
        if len(name) == 1:
            name = f" {name}  "
        elif len(name) == 2:
            name = f" {name} "
        else:
            name = f"{name:4s}"

        line = (f"ATOM  {atom_idx:5d} {name}"
                f" LIG     1    "
                f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}"
                f"  1.00  0.00    "
                f"{charge:+6.3f} {ad4_type:<2s}")
        lines.append(line)

    lines.append("ENDROOT")
    lines.append("TORSDOF 0")

    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")

    return True


def find_docking_exe():
    if os.path.isfile(QVINA2_PATH) and os.access(QVINA2_PATH, os.X_OK):
        return QVINA2_PATH
    if os.path.isfile(VINA_PATH) and os.access(VINA_PATH, os.X_OK):
        return VINA_PATH
    for name in ["qvina2", "qvina2.1", "qvina02", "vina"]:
        p = shutil.which(name)
        if p:
            return p
    return None


_WORKER_CONFIG = {}


def init_worker(config):
    global _WORKER_CONFIG
    _WORKER_CONFIG = config
    os.makedirs(config["tmp_dir"], exist_ok=True)


def dock_one(item):
    smiles, name = item
    config = _WORKER_CONFIG
    pid = os.getpid()
    work_dir = os.path.join(config["tmp_dir"], f"worker_{pid}")
    os.makedirs(work_dir, exist_ok=True)

    pdbqt_path = os.path.join(work_dir, "lig.pdbqt")
    out_path = os.path.join(work_dir, "lig_out.pdbqt")

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return name, smiles, None, "invalid_smiles"

        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = hash(smiles) & 0xFFFFFFFF
        result = AllChem.EmbedMolecule(mol, params)
        if result != 0:
            params.useRandomCoords = True
            result = AllChem.EmbedMolecule(mol, params)
            if result != 0:
                return name, smiles, None, "embed_failed"

        try:
            AllChem.UFFOptimizeMolecule(mol, maxIters=200)
        except Exception:
            pass

        if not mol_to_pdbqt(mol, pdbqt_path):
            return name, smiles, None, "pdbqt_failed"

        if not os.path.exists(pdbqt_path) or os.path.getsize(pdbqt_path) == 0:
            return name, smiles, None, "conversion_failed"

        cmd = [
            config["docking_exe"],
            "--receptor", config["receptor"],
            "--ligand", pdbqt_path,
            "--center_x", str(config["center_x"]),
            "--center_y", str(config["center_y"]),
            "--center_z", str(config["center_z"]),
            "--size_x", str(config["size_x"]),
            "--size_y", str(config["size_y"]),
            "--size_z", str(config["size_z"]),
            "--exhaustiveness", str(config["exhaustiveness"]),
            "--num_modes", "1",
            "--out", out_path,
            "--cpu", "1",
        ]

        proc = subprocess.run(cmd, capture_output=True, text=True,
                              timeout=config["timeout"])

        score = None
        if os.path.exists(out_path):
            with open(out_path) as f:
                for line in f:
                    if "VINA RESULT" in line:
                        parts = line.split()
                        for i, p in enumerate(parts):
                            if p == "RESULT:":
                                score = float(parts[i + 1])
                                break
                        if score is not None:
                            break

        if score is None and proc.returncode != 0:
            stderr = proc.stderr[:80] if proc.stderr else "unknown"
            return name, smiles, None, stderr

        return name, smiles, score, "ok" if score is not None else "no_score"

    except subprocess.TimeoutExpired:
        return name, smiles, None, "timeout"
    except Exception as e:
        return name, smiles, None, str(e)[:80]
    finally:
        for f in [pdbqt_path, out_path]:
            try:
                os.unlink(f)
            except OSError:
                pass


def main():
    parser = argparse.ArgumentParser(
        description="Virtual screen: dock SMILES database with Vina/QuickVina2")
    parser.add_argument("--input", "-i", required=True,
                        help="Input .smi file (SMILES<tab>Name)")
    parser.add_argument("--receptor", "-r", required=True,
                        help="Receptor .pdbqt file")
    parser.add_argument("--output", "-o", required=True,
                        help="Output prefix (writes _scores.tsv and _top10.smi)")
    parser.add_argument("--center_x", type=float, default=55.0)
    parser.add_argument("--center_y", type=float, default=59.5)
    parser.add_argument("--center_z", type=float, default=35.5)
    parser.add_argument("--size_x", type=float, default=19.0)
    parser.add_argument("--size_y", type=float, default=19.0)
    parser.add_argument("--size_z", type=float, default=19.0)
    parser.add_argument("--exhaustiveness", type=int, default=1,
                        help="Vina exhaustiveness (default: 1 for fast screening)")
    parser.add_argument("--timeout", type=int, default=120,
                        help="Per-compound docking timeout in seconds")
    parser.add_argument("--nprocs", type=int, default=-1,
                        help="Number of workers (-1 = all available CPUs)")
    parser.add_argument("--top_percent", type=float, default=10.0,
                        help="Keep top N%% by docking score (default: 10)")
    parser.add_argument("--resume", action="store_true",
                        help="Resume from checkpoint")
    parser.add_argument("--chunk_size", type=int, default=500,
                        help="Write checkpoint every N compounds")

    args = parser.parse_args()

    scores_file = args.output + "_scores.tsv"
    top_file = args.output + f"_top{int(args.top_percent)}.smi"
    tmp_dir = args.output + "_tmp"

    docking_exe = find_docking_exe()
    if not docking_exe:
        print("ERROR: No docking executable found (QuickVina2 or Vina)")
        sys.exit(1)

    print(f"PDBQT conversion: built-in RDKit (Gasteiger charges + AD4 atom types)")
    print(f"Docking engine: {docking_exe}")
    print(f"Receptor: {args.receptor}")
    print(f"Center: ({args.center_x}, {args.center_y}, {args.center_z})")
    print(f"Box size: ({args.size_x}, {args.size_y}, {args.size_z})")
    print(f"Exhaustiveness: {args.exhaustiveness}")
    print(f"Top percent to keep: {args.top_percent}%")

    compounds = []
    with open(args.input) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                compounds.append((parts[0], parts[1]))
            else:
                compounds.append((parts[0], f"mol_{len(compounds)}"))

    print(f"\nLoaded {len(compounds):,} compounds from {args.input}")

    done_set = set()
    if args.resume and os.path.exists(scores_file):
        with open(scores_file) as f:
            for line in f:
                if line.startswith("name\t"):
                    continue
                parts = line.strip().split("\t")
                if parts:
                    done_set.add(parts[0])
        print(f"Resuming: {len(done_set):,} already scored, "
              f"{len(compounds) - len(done_set):,} remaining")

    todo = [(smi, name) for smi, name in compounds if name not in done_set]

    if not todo:
        print("All compounds already scored. Generating top list...")
    else:
        nprocs = args.nprocs if args.nprocs > 0 else multiprocessing.cpu_count()
        nprocs = min(nprocs, len(todo))
        print(f"Using {nprocs} workers")

        est_time = len(todo) * 5 / nprocs
        print(f"Estimated time: {est_time/3600:.1f} hours "
              f"(~5 sec/compound, {nprocs} workers)")
        print()

        config = {
            "receptor": os.path.abspath(args.receptor),
            "center_x": args.center_x,
            "center_y": args.center_y,
            "center_z": args.center_z,
            "size_x": args.size_x,
            "size_y": args.size_y,
            "size_z": args.size_z,
            "exhaustiveness": args.exhaustiveness,
            "timeout": args.timeout,
            "docking_exe": docking_exe,
            "tmp_dir": os.path.abspath(tmp_dir),
        }

        os.makedirs(tmp_dir, exist_ok=True)

        write_header = not os.path.exists(scores_file)
        scores_f = open(scores_file, "a")
        if write_header:
            scores_f.write("name\tsmiles\tscore\tstatus\n")
            scores_f.flush()

        scored = 0
        failed = 0
        best_score = 0
        start_time = time.time()

        pool = multiprocessing.Pool(nprocs, initializer=init_worker,
                                    initargs=(config,))

        try:
            for result in pool.imap_unordered(dock_one, todo,
                                              chunksize=1):
                name, smiles, score, status = result
                score_str = f"{score:.2f}" if score is not None else "NA"
                scores_f.write(f"{name}\t{smiles}\t{score_str}\t{status}\n")

                if score is not None:
                    scored += 1
                    if score < best_score:
                        best_score = score
                else:
                    failed += 1

                total = scored + failed
                if total % args.chunk_size == 0:
                    scores_f.flush()
                    elapsed = time.time() - start_time
                    rate = total / elapsed if elapsed > 0 else 0
                    remaining = (len(todo) - total) / rate if rate > 0 else 0
                    print(f"  Progress: {total:>10,}/{len(todo):,} | "
                          f"Scored: {scored:,} | Failed: {failed:,} | "
                          f"Best: {best_score:.1f} | "
                          f"Rate: {rate:.1f}/sec | "
                          f"ETA: {remaining/3600:.1f}h")

        except KeyboardInterrupt:
            print(f"\nInterrupted. {scored + failed:,} compounds processed.")
            print("Use --resume to continue.")
        finally:
            pool.terminate()
            pool.join()
            scores_f.close()

        elapsed = time.time() - start_time
        print(f"\nScreening complete!")
        print(f"  Scored: {scored:,}")
        print(f"  Failed: {failed:,}")
        print(f"  Best score: {best_score:.2f}")
        print(f"  Time: {elapsed/3600:.1f} hours")

    print(f"\nGenerating top {args.top_percent}% list...")

    all_scores = []
    with open(scores_file) as f:
        for line in f:
            if line.startswith("name\t"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 3 and parts[2] != "NA":
                try:
                    all_scores.append((parts[0], parts[1], float(parts[2])))
                except ValueError:
                    pass

    all_scores.sort(key=lambda x: x[2])
    n_keep = max(1, int(len(all_scores) * args.top_percent / 100))

    top = all_scores[:n_keep]

    with open(top_file, "w") as f:
        for name, smiles, score in top:
            f.write(f"{smiles}\t{name}\n")

    print(f"  Total scored: {len(all_scores):,}")
    print(f"  Keeping top {args.top_percent}%: {n_keep:,} compounds")
    if top:
        print(f"  Score range: {top[0][2]:.2f} to {top[-1][2]:.2f}")
    print(f"  Output: {top_file}")

    try:
        shutil.rmtree(tmp_dir, ignore_errors=True)
    except Exception:
        pass

    print("\nDone!")


if __name__ == "__main__":
    main()
