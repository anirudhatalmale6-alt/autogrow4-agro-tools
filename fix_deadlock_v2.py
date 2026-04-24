#!/usr/bin/env python3
"""
Fix AutoGrow4 deadlock v2 - Fix the actual root cause.

The func_timeout library uses async exception injection (PyThreadState_SetAsyncExc)
which kills worker processes entirely, bypassing try/except blocks.
This causes the multiprocessing queue to deadlock.

Fix: Replace func_timeout in Gypsum conversion with multiprocessing.Process + join(timeout)
which is a clean, reliable timeout mechanism.

Also reduces queue timeout from 900 to 120 seconds as additional safety.
"""

import os
import shutil

AUTOGROW_DIR = os.path.expanduser("~/final_project/autogrow4")


def patch_conversion_to_3d():
    """Replace func_timeout with multiprocessing.Process timeout in Gypsum conversion"""
    filepath = os.path.join(
        AUTOGROW_DIR,
        "autogrow/operators/convert_files/conversion_to_3d.py"
    )

    if not os.path.exists(filepath):
        print(f"ERROR: {filepath} not found")
        return False

    backup = filepath + ".bak2"
    if not os.path.exists(backup):
        shutil.copy2(filepath, backup)
        print(f"  Backed up: {backup}")

    with open(filepath) as f:
        content = f.read()

    if "GYPSUM_PROCESS_FIX" in content:
        print("  conversion_to_3d.py already patched (v2)")
        return True

    # Replace the run_gypsum_multiprocessing function
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
        print("  Patched run_gypsum_multiprocessing() with process-based timeout")
    else:
        print("  WARNING: Could not find exact match, trying flexible patch")
        # Try matching key part
        content = content.replace(
            "        with StdoutRedirection(log_file):\n            func_timeout(gypsum_timeout_limit, prepare_molecules, args=(gypsum_params,))",
            "        # GYPSUM_PROCESS_FIX: see fix_deadlock_v2.py\n        import multiprocessing as _mp\n        def _run(gp, lf):\n            try:\n                with StdoutRedirection(lf):\n                    prepare_molecules(gp)\n            except Exception:\n                pass\n        p = _mp.Process(target=_run, args=(gypsum_params, log_file))\n        p.start()\n        p.join(timeout=gypsum_timeout_limit)\n        if p.is_alive():\n            p.terminate()\n            p.join(5)\n            if p.is_alive(): p.kill()\n            return lig_id"
        )

    with open(filepath, 'w') as f:
        f.write(content)

    print("  conversion_to_3d.py patched successfully")
    return True


def patch_parallelizer_timeout():
    """Reduce queue timeout from 900 to 120 seconds"""
    filepath = os.path.join(
        AUTOGROW_DIR,
        "autogrow/operators/convert_files/gypsum_dl/gypsum_dl/Parallelizer.py"
    )

    if not os.path.exists(filepath):
        print(f"ERROR: {filepath} not found")
        return False

    with open(filepath) as f:
        content = f.read()

    if "timeout=120" in content:
        print("  Parallelizer.py queue timeout already set to 120s")
        return True

    content = content.replace("timeout=900", "timeout=120")
    content = content.replace("timed out after 900s", "timed out after 120s")

    with open(filepath, 'w') as f:
        f.write(content)

    print("  Reduced queue timeout from 900 to 120 seconds")
    return True


def patch_gypsum_timeout_config():
    """Reduce gypsum_timeout_limit from 60 to 20 seconds"""
    import json
    config_path = os.path.expanduser("~/final_project/autogrow_config.json")

    with open(config_path) as f:
        config = json.load(f)

    if config.get("gypsum_timeout_limit", 60) > 20:
        config["gypsum_timeout_limit"] = 20
        with open(config_path, 'w') as f:
            json.dump(config, f, indent=4)
        print("  Reduced gypsum_timeout_limit from 60 to 20 seconds")
    else:
        print("  gypsum_timeout_limit already reasonable")

    return True


def main():
    print("=" * 60)
    print("AutoGrow4 Deadlock Fix v2")
    print("Root cause: func_timeout kills worker processes")
    print("=" * 60)

    print("\n1. Patching conversion_to_3d.py (replace func_timeout)...")
    patch_conversion_to_3d()

    print("\n2. Reducing Parallelizer.py queue timeout...")
    patch_parallelizer_timeout()

    print("\n3. Updating gypsum_timeout_limit in config...")
    patch_gypsum_timeout_config()

    print("\n" + "=" * 60)
    print("All v2 patches applied!")
    print("Gypsum now uses clean process-based timeout instead of func_timeout")
    print("Queue timeout reduced to 120s as safety net")
    print("=" * 60)


if __name__ == "__main__":
    main()
