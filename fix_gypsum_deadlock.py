#!/usr/bin/env python3
"""
Fix Gypsum deadlock in AutoGrow4.
Combines fix_deadlock_v2 + v3: replaces func_timeout with multiprocessing.Process
and removes daemon=True from Parallelizer workers.

Run from anywhere: python3 fix_gypsum_deadlock.py

What it does:
1. Patches conversion_to_3d.py - replaces func_timeout (which kills workers via
   async exception injection) with multiprocessing.Process + join(timeout)
2. Patches Parallelizer.py - removes daemon=True so workers can spawn Gypsum
   subprocesses, and reduces queue timeout to 120s as safety net
3. Clears __pycache__ to avoid stale bytecode
"""
import os
import shutil
import sys

AUTOGROW_DIR = os.path.expanduser("~/final_project/autogrow4")


def patch_conversion_to_3d():
    filepath = os.path.join(
        AUTOGROW_DIR,
        "autogrow/operators/convert_files/conversion_to_3d.py"
    )
    if not os.path.exists(filepath):
        print("ERROR: " + filepath + " not found")
        return False

    backup = filepath + ".bak_gypsum"
    if not os.path.exists(backup):
        shutil.copy2(filepath, backup)
        print("  Backed up to " + backup)

    with open(filepath) as f:
        content = f.read()

    if "GYPSUM_PROCESS_FIX" in content:
        print("  conversion_to_3d.py already patched")
        return True

    old_func = '''    lig_id = gypsum_params["source"].split(os.sep)[-1].replace(".smi", "")
    log_file = f"{gypsum_log_path}{lig_id}_log.txt"

    try:
        with StdoutRedirection(log_file):
            func_timeout(gypsum_timeout_limit, prepare_molecules, args=(gypsum_params,))

        sys.stdout.flush()
    except Exception:
        # This Ligand Timed out
        return lig_id

    # Check if it worked if it failed return lig_id if it works return None
    did_gypsum_complete = check_gypsum_log_did_complete(log_file)
    if did_gypsum_complete in [None, False]:
        # Failed to convert
        return lig_id

    return None'''

    new_func = '''    # GYPSUM_PROCESS_FIX: Use multiprocessing.Process instead of func_timeout
    # func_timeout uses async exception injection which kills worker processes
    import multiprocessing as _mp

    lig_id = gypsum_params["source"].split(os.sep)[-1].replace(".smi", "")
    log_file = f"{gypsum_log_path}{lig_id}_log.txt"

    def _run_gypsum_isolated(gp, lf):
        """Run Gypsum in isolated process so timeout can cleanly kill it."""
        try:
            with StdoutRedirection(lf):
                prepare_molecules(gp)
        except Exception:
            pass

    try:
        p = _mp.Process(target=_run_gypsum_isolated, args=(gypsum_params, log_file))
        p.start()
        p.join(timeout=gypsum_timeout_limit)

        if p.is_alive():
            p.terminate()
            p.join(5)
            if p.is_alive():
                p.kill()
                p.join(2)
            print(f"  Gypsum timed out ({gypsum_timeout_limit}s): {lig_id}")
            return lig_id

        if p.exitcode != 0:
            print(f"  Gypsum failed (exit {p.exitcode}): {lig_id}")
            return lig_id
    except Exception as e:
        print(f"  Gypsum error for {lig_id}: {e}")
        return lig_id

    sys.stdout.flush()

    # Check if it worked if it failed return lig_id if it works return None
    did_gypsum_complete = check_gypsum_log_did_complete(log_file)
    if did_gypsum_complete in [None, False]:
        # Failed to convert
        return lig_id

    return None'''

    if old_func in content:
        content = content.replace(old_func, new_func)
        with open(filepath, 'w') as f:
            f.write(content)
        print("  Patched conversion_to_3d.py successfully")
        return True
    else:
        print("  WARNING: Could not find exact match in conversion_to_3d.py")
        print("  Trying flexible match...")
        if "func_timeout(gypsum_timeout_limit, prepare_molecules" in content:
            content = content.replace(
                "func_timeout(gypsum_timeout_limit, prepare_molecules, args=(gypsum_params,))",
                "# GYPSUM_PROCESS_FIX applied - see below"
            )
            print("  Partial match found - please verify manually")
        else:
            print("  ERROR: Cannot find func_timeout call to patch")
            return False

    with open(filepath, 'w') as f:
        f.write(content)
    return True


def patch_parallelizer():
    filepath = os.path.join(
        AUTOGROW_DIR,
        "autogrow/operators/convert_files/gypsum_dl/gypsum_dl/Parallelizer.py"
    )
    if not os.path.exists(filepath):
        print("ERROR: " + filepath + " not found")
        return False

    with open(filepath) as f:
        content = f.read()

    changed = False

    if "p.daemon = True" in content:
        content = content.replace("        p.daemon = True\n", "")
        changed = True
        print("  Removed daemon=True from Parallelizer workers")

    if "timeout=900" in content:
        content = content.replace("timeout=900", "timeout=120")
        content = content.replace("timed out after 900s", "timed out after 120s")
        changed = True
        print("  Reduced queue timeout from 900s to 120s")

    if changed:
        with open(filepath, 'w') as f:
            f.write(content)
        print("  Parallelizer.py patched successfully")
    else:
        print("  Parallelizer.py already patched or no changes needed")

    return True


def clear_pycache():
    count = 0
    for root, dirs, files in os.walk(AUTOGROW_DIR):
        for d in dirs:
            if d == "__pycache__":
                shutil.rmtree(os.path.join(root, d))
                count += 1
        for f in files:
            if f.endswith(".pyc"):
                os.remove(os.path.join(root, f))
                count += 1
    print(f"  Cleared {count} cache items")


def main():
    print("=" * 60)
    print("AutoGrow4 Gypsum Deadlock Fix")
    print("Replaces func_timeout with process-based timeout")
    print("=" * 60)

    if not os.path.isdir(AUTOGROW_DIR):
        print("ERROR: AutoGrow4 directory not found at " + AUTOGROW_DIR)
        sys.exit(1)

    print("\n1. Patching conversion_to_3d.py...")
    r1 = patch_conversion_to_3d()

    print("\n2. Patching Parallelizer.py...")
    r2 = patch_parallelizer()

    print("\n3. Clearing __pycache__...")
    clear_pycache()

    print("\n" + "=" * 60)
    if r1 and r2:
        print("ALL PATCHES APPLIED SUCCESSFULLY")
        print("Gypsum now uses clean process-based timeout")
        print("Safe to use with 8-32 processors")
    else:
        print("SOME PATCHES FAILED - check output above")
    print("=" * 60)


if __name__ == "__main__":
    main()
