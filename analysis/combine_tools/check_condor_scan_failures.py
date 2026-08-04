#!/usr/bin/env python3
"""
Scan condor_jobs/logs/*.out (matched up against the generated job scripts) for
the "negative event yield" EFT-validity failure signature seen in cbBRe_cbe -
a job can exit cleanly (return value 0) but still come back with far fewer
scanned points than requested because combine's grid loop stopped early after
hitting an unphysical (negative) predicted yield near the edge of the box.

This does NOT try to predict which operators will hit it ahead of time - it
just tells you, after a batch of condor jobs finishes, which operator pairs
actually did, and how many of their jobs were affected, so those can be
targeted with a narrower range instead of the whole thing being re-run blind.

Usage:
    python3 check_condor_scan_failures.py --condor-dir condor_jobs
"""

import argparse
import os
import re

FAIL_SIGNATURE = "Number of events is negative or error"


def get_pair_from_script(script_path):
    with open(script_path) as f:
        content = f.read()
    m = re.search(r"--redefineSignalPOIs=k_(\w+),k_(\w+)", content)
    if m:
        return m.group(1), m.group(2)
    return None, None


def get_expected_points(script_path):
    with open(script_path) as f:
        content = f.read()
    m = re.search(r"--firstPoint (\d+) --lastPoint (\d+)", content)
    if m:
        first, last = int(m.group(1)), int(m.group(2))
        return last - first + 1
    return None


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--condor-dir", default="condor_jobs")
    args = p.parse_args()

    joblist_path = os.path.join(args.condor_dir, "joblist.txt")
    with open(joblist_path) as f:
        scripts = [line.strip() for line in f if line.strip()]

    logs_dir = os.path.join(args.condor_dir, "logs")

    # Process IDs are assigned in the same order the queue file was read, so
    # line N of joblist.txt corresponds to logs/N.out
    results = {}
    missing_logs = 0

    for proc_id, script_path in enumerate(scripts):
        op1, op2 = get_pair_from_script(script_path)
        if op1 is None:
            continue
        pair = f"{op1}_{op2}"
        expected = get_expected_points(script_path) or 0

        out_path = os.path.join(logs_dir, f"{proc_id}.out")
        failed = False
        if os.path.isfile(out_path):
            with open(out_path, errors="ignore") as f:
                if FAIL_SIGNATURE in f.read():
                    failed = True
        else:
            missing_logs += 1

        r = results.setdefault(pair, {"jobs": 0, "failed_jobs": 0, "expected_points": 0})
        r["jobs"] += 1
        r["expected_points"] += expected
        if failed:
            r["failed_jobs"] += 1

    affected = {pair: r for pair, r in results.items() if r["failed_jobs"] > 0}

    print(f"{len(results)} pairs total, {len(affected)} hit the negative-yield signature\n")

    if affected:
        print(f"{'pair':30s} {'jobs':>6s} {'failed_jobs':>12s}")
        for pair, r in sorted(affected.items(), key=lambda kv: -kv[1]["failed_jobs"]):
            print(f"{pair:30s} {r['jobs']:6d} {r['failed_jobs']:12d}")
        print("\nComma-separated for --doOnly / readapt --pairs:")
        ops = sorted({op for pair in affected for op in pair.split("_", 1)})
        print("  affected pairs:", ",".join(sorted(affected)))
        print("  operators involved:", ",".join(ops))

    if missing_logs:
        print(f"\n[WARN] {missing_logs} jobs had no matching .out log (still running, or condor hasn't flushed logs yet)")


if __name__ == "__main__":
    main()
