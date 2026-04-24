#!/usr/bin/env python3
"""
Fix v3: Remove daemon=True from Parallelizer workers.
daemon processes can't spawn children, which breaks the v2 Gypsum process-based timeout.
"""
import os

filepath = os.path.expanduser(
    "~/final_project/autogrow4/autogrow/operators/convert_files/gypsum_dl/gypsum_dl/Parallelizer.py"
)

with open(filepath) as f:
    content = f.read()

# Remove daemon=True line
old = """    # Start worker processes (DEADLOCK_FIX: daemon=True so they die with parent)
    workers_list = []
    for i in range(num_procs):
        p = multiprocessing.Process(target=worker, args=(task_queue, done_queue))
        p.daemon = True
        p.start()
        workers_list.append(p)"""

new = """    # Start worker processes (DEADLOCK_FIX: non-daemon so they can spawn Gypsum subprocesses)
    workers_list = []
    for i in range(num_procs):
        p = multiprocessing.Process(target=worker, args=(task_queue, done_queue))
        p.start()
        workers_list.append(p)"""

if old in content:
    content = content.replace(old, new)
    with open(filepath, 'w') as f:
        f.write(content)
    print("Removed daemon=True from Parallelizer workers")
else:
    print("Could not find exact match - checking if already fixed")
    if "p.daemon = True" in content:
        content = content.replace("        p.daemon = True\n", "")
        with open(filepath, 'w') as f:
            f.write(content)
        print("Removed p.daemon = True line")
    else:
        print("daemon=True already removed")
