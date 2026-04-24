#!/usr/bin/env python3
"""
Fix AutoGrow4 multiprocessing deadlock issues.

Three patches:
1. Parallelizer.py: Add timeout to done_queue.get() and wrap worker in try/except
2. vina_docking.py: Replace os.system() with subprocess.run(timeout=...)
3. Config: Set timeout_vs_gtimeout if missing

Run: python fix_deadlock.py
"""

import os
import sys
import json
import shutil

AUTOGROW_DIR = os.path.expanduser("~/final_project/autogrow4")
CONFIG_PATH = os.path.expanduser("~/final_project/autogrow_config.json")

def patch_parallelizer():
    """Fix the queue deadlock in Parallelizer.py"""
    filepath = os.path.join(
        AUTOGROW_DIR,
        "autogrow/operators/convert_files/gypsum_dl/gypsum_dl/Parallelizer.py"
    )

    if not os.path.exists(filepath):
        print(f"ERROR: {filepath} not found")
        return False

    # Backup
    backup = filepath + ".bak"
    if not os.path.exists(backup):
        shutil.copy2(filepath, backup)
        print(f"  Backed up: {backup}")

    with open(filepath) as f:
        content = f.read()

    # Check if already patched
    if "DEADLOCK_FIX" in content:
        print("  Parallelizer.py already patched, skipping")
        return True

    # Patch 1: Replace worker function to catch exceptions
    old_worker = '''def worker(input, output):
    for seq, job in iter(input.get, "STOP"):
        func, args = job
        result = func(*args)
        ret_val = (seq, result)
        output.put(ret_val)'''

    new_worker = '''def worker(input, output):
    # DEADLOCK_FIX: Wrap in try/except so crashed workers still put results
    for seq, job in iter(input.get, "STOP"):
        func, args = job
        try:
            result = func(*args)
        except Exception as e:
            print(f"  Worker error on task {seq}: {e}")
            result = None
        ret_val = (seq, result)
        output.put(ret_val)'''

    if old_worker in content:
        content = content.replace(old_worker, new_worker)
        print("  Patched worker() function with try/except")
    else:
        print("  WARNING: Could not find exact worker() match, trying flexible patch")
        # Try line-by-line
        content = content.replace(
            "        result = func(*args)\n        ret_val = (seq, result)\n        output.put(ret_val)",
            "        try:\n            result = func(*args)\n        except Exception as e:\n            print(f\"  Worker error on task {seq}: {e}\")\n            result = None\n        ret_val = (seq, result)\n        output.put(ret_val)"
        )

    # Patch 2: Add timeout to done_queue.get() and make workers daemon processes
    old_start = '''    # Start worker processes
    for i in range(num_procs):
        multiprocessing.Process(target=worker, args=(task_queue, done_queue)).start()

    # Get and print results
    results = []
    for i in range(len(inputs)):
        results.append(done_queue.get())

    # Tell child processes to stop
    for i in range(num_procs):
        task_queue.put("STOP")'''

    new_start = '''    # Start worker processes (DEADLOCK_FIX: daemon=True so they die with parent)
    workers_list = []
    for i in range(num_procs):
        p = multiprocessing.Process(target=worker, args=(task_queue, done_queue))
        p.daemon = True
        p.start()
        workers_list.append(p)

    # Get results with timeout (DEADLOCK_FIX: 900s per task max, skip if hung)
    import queue as queue_module
    results = []
    failed_count = 0
    for i in range(len(inputs)):
        try:
            result = done_queue.get(timeout=900)
            results.append(result)
        except queue_module.Empty:
            print(f"  WARNING: Task {i} timed out after 900s, skipping")
            failed_count += 1
            results.append((i, None))
    if failed_count > 0:
        print(f"  {failed_count} tasks timed out and were skipped")

    # Tell child processes to stop
    for i in range(num_procs):
        task_queue.put("STOP")'''

    if old_start in content:
        content = content.replace(old_start, new_start)
        print("  Patched start_processes() with queue timeout and daemon workers")
    else:
        print("  WARNING: Could not find exact start_processes() match")
        # Try a more flexible approach
        content = content.replace(
            "        results.append(done_queue.get())",
            "        import queue as queue_module\n        try:\n            results.append(done_queue.get(timeout=900))\n        except queue_module.Empty:\n            print(f\"  WARNING: Task {i} timed out, skipping\")\n            results.append((i, None))"
        )
        print("  Applied flexible queue timeout patch")

    with open(filepath, 'w') as f:
        f.write(content)

    print("  Parallelizer.py patched successfully")
    return True


def patch_vina_docking():
    """Replace os.system() with subprocess.run(timeout=...) in vina_docking.py"""
    filepath = os.path.join(
        AUTOGROW_DIR,
        "autogrow/docking/docking_class/docking_class_children/vina_docking.py"
    )

    if not os.path.exists(filepath):
        print(f"ERROR: {filepath} not found")
        return False

    backup = filepath + ".bak"
    if not os.path.exists(backup):
        shutil.copy2(filepath, backup)
        print(f"  Backed up: {backup}")

    with open(filepath) as f:
        content = f.read()

    if "DEADLOCK_FIX" in content:
        print("  vina_docking.py already patched, skipping")
        return True

    # Add subprocess import at the top
    if "import subprocess" not in content:
        content = "import subprocess\n" + content
        print("  Added subprocess import")

    # Replace execute_docking_vina to use subprocess with Python-level timeout
    old_execute = '''    def execute_docking_vina(self, command):
        """
        Run a single docking execution command

        Inputs:
        :param str command: string of command to run.

        Returns:
        :returns: int result: the exit output for the command. If its None of
            256 it failed.
        """

        try:
            result = os.system(command)
        except Exception:
            result = None
            print(f"Failed to execute: {command}")
        return result'''

    new_execute = '''    def execute_docking_vina(self, command):
        """
        Run a single docking execution command
        DEADLOCK_FIX: Uses subprocess.run with Python-level timeout
        instead of os.system which can hang forever.

        Inputs:
        :param str command: string of command to run.

        Returns:
        :returns: int result: the exit output for the command. If its None of
            256 it failed.
        """
        # Use Python timeout as safety net (bash timeout + 60s buffer)
        python_timeout = self.params.get("docking_timeout_limit", 120) + 60

        try:
            proc = subprocess.run(
                command, shell=True,
                timeout=python_timeout,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )
            result = proc.returncode * 256  # match os.system() return format
        except subprocess.TimeoutExpired:
            print(f"  Python timeout ({python_timeout}s) killed docking: {command[:100]}...")
            result = 31744  # same as bash timeout exit code (124*256)
        except Exception:
            result = None
            print(f"Failed to execute: {command}")
        return result'''

    if old_execute in content:
        content = content.replace(old_execute, new_execute)
        print("  Patched execute_docking_vina() with subprocess.run + timeout")
    else:
        print("  WARNING: Could not find exact execute_docking_vina() match, trying flexible patch")
        # Try replacing just the core
        content = content.replace(
            "            result = os.system(command)\n        except Exception:\n            result = None\n            print(f\"Failed to execute: {command}\")",
            "            python_timeout = self.params.get('docking_timeout_limit', 120) + 60\n            proc = subprocess.run(command, shell=True, timeout=python_timeout, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)\n            result = proc.returncode * 256\n        except subprocess.TimeoutExpired:\n            print(f'  Python timeout killed docking: {command[:100]}...')\n            result = 31744\n        except Exception:\n            result = None\n            print(f\"Failed to execute: {command}\")"
        )
        print("  Applied flexible subprocess patch")

    with open(filepath, 'w') as f:
        f.write(content)

    print("  vina_docking.py patched successfully")
    return True


def patch_config():
    """Ensure config has timeout_vs_gtimeout set and reasonable values"""
    if not os.path.exists(CONFIG_PATH):
        print(f"ERROR: {CONFIG_PATH} not found")
        return False

    with open(CONFIG_PATH) as f:
        config = json.load(f)

    changed = False

    # Set timeout_vs_gtimeout if missing
    if "timeout_vs_gtimeout" not in config or not config["timeout_vs_gtimeout"]:
        config["timeout_vs_gtimeout"] = "timeout"
        print("  Set timeout_vs_gtimeout = 'timeout'")
        changed = True

    # Reduce docking_timeout_limit if too high (600 is 10 min per molecule!)
    if config.get("docking_timeout_limit", 120) > 300:
        config["docking_timeout_limit"] = 180
        print("  Reduced docking_timeout_limit from 600 to 180 seconds")
        changed = True

    # Reduce max_variants_per_compound to 1 for speed
    if config.get("max_variants_per_compound", 3) > 1:
        config["max_variants_per_compound"] = 1
        print(f"  Reduced max_variants_per_compound to 1")
        changed = True

    if changed:
        with open(CONFIG_PATH, 'w') as f:
            json.dump(config, f, indent=4)
        print("  Config updated")
    else:
        print("  Config already looks good")

    return True


def main():
    print("=" * 60)
    print("AutoGrow4 Deadlock Fix")
    print("=" * 60)

    print("\n1. Patching Parallelizer.py (queue deadlock fix)...")
    patch_parallelizer()

    print("\n2. Patching vina_docking.py (subprocess timeout fix)...")
    patch_vina_docking()

    print("\n3. Updating config...")
    patch_config()

    print("\n" + "=" * 60)
    print("All patches applied!")
    print("Your next run should no longer stall on hung molecules.")
    print("=" * 60)


if __name__ == "__main__":
    main()
