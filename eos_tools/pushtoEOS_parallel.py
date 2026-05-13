#!/usr/bin/env python3
"""
pushtoEOS_parallel.py  <local_dir>  <eos_www_subpath>  [-j N]

Like pushtoEOS.py but uploads and index-plants in parallel using N worker
threads (default: 4).  Safe because every xrdcp call writes to a distinct
EOS path.

Requires: xrdcp and eos commands available (run cmsenv first).

EOS www base: /eos/user/a/aldufour/www/

Example:
    python pushtoEOS_parallel.py plots/syst_scan_v2  syst_scan_v2
    python pushtoEOS_parallel.py results/figures     DY2026/figures  -j 8
"""

import os
import sys
import subprocess
import threading
import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

EOS_WWW_BASE = "/eos/user/a/aldufour/www"
EOS_MGM_URL  = "root://eosuser.cern.ch"
INDEX_PHP    = Path(__file__).parent / "index.php"

_print_lock = threading.Lock()


def xrd(eos_abs_path: str) -> str:
    return f"{EOS_MGM_URL}/{eos_abs_path}"


def run(cmd: list[str], check: bool = True) -> subprocess.CompletedProcess:
    result = subprocess.run(cmd, check=False, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    with _print_lock:
        print("  $", " ".join(cmd))
        if result.stdout.strip():
            print(result.stdout, end="")
        if result.stderr.strip():
            print(result.stderr, end="", file=sys.stderr)
    if check and result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, cmd)
    return result


def upload_item(item: Path, eos_dest: str) -> tuple[Path, bool, str]:
    """Upload a single top-level item. Returns (item, success, error_msg)."""
    try:
        run(["xrdcp", "-r", "--silent", str(item), xrd(eos_dest + "/")])
        return item, True, ""
    except subprocess.CalledProcessError as e:
        return item, False, str(e)


def plant_index(eos_dir: str) -> tuple[str, bool]:
    """Plant index.php in one EOS directory. Returns (eos_dir, success)."""
    result = run(["xrdcp", "--silent", str(INDEX_PHP), xrd(f"{eos_dir}/index.php")],
                 check=False)
    return eos_dir, result.returncode == 0


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Push a local directory to EOS www with parallel uploads."
    )
    parser.add_argument("local_dir",    help="Local directory to upload")
    parser.add_argument("eos_subpath",  help="Destination under /eos/user/a/aldufour/www/")
    parser.add_argument("-j", "--workers", type=int, default=4,
                        help="Number of parallel workers (default: 4)")
    args = parser.parse_args()

    local_dir   = Path(args.local_dir).resolve()
    eos_subpath = args.eos_subpath.strip("/")
    eos_dest    = f"{EOS_WWW_BASE}/{eos_subpath}"
    workers     = args.workers

    if not local_dir.is_dir():
        sys.exit(f"Error: '{local_dir}' is not a directory")
    if not INDEX_PHP.is_file():
        sys.exit(f"Error: index.php not found at '{INDEX_PHP}'")

    os.environ["EOS_MGM_URL"] = EOS_MGM_URL

    # ------------------------------------------------------------------ #
    # 1. Check destination doesn't already exist, then create it
    # ------------------------------------------------------------------ #
    result = run(["eos", "ls", eos_dest], check=False)
    if result.returncode == 0:
        sys.exit(f"Error: '{eos_dest}' already exists on EOS — refusing to overwrite.")

    print(f"\n[1/3] Uploading '{local_dir}' -> {eos_dest}  (workers={workers})")
    run(["eos", "mkdir", "-p", eos_dest])

    items = sorted(local_dir.iterdir())
    failed_uploads = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(upload_item, item, eos_dest): item for item in items}
        for future in as_completed(futures):
            item, ok, err = future.result()
            if not ok:
                failed_uploads.append((item, err))

    if failed_uploads:
        print("\nERROR: some uploads failed:")
        for item, err in failed_uploads:
            print(f"  {item}: {err}")
        sys.exit(1)

    # ------------------------------------------------------------------ #
    # 2. Collect every directory created on EOS
    # ------------------------------------------------------------------ #
    print(f"\n[2/3] Collecting directories to index...")
    eos_dirs: list[str] = []
    for root, dirs, _ in os.walk(local_dir):
        dirs.sort()
        rel = Path(root).relative_to(local_dir)
        eos_dirs.append(eos_dest if str(rel) == "." else f"{eos_dest}/{rel}")
    print(f"  Found {len(eos_dirs)} directories")

    # ------------------------------------------------------------------ #
    # 3. Plant index.php in every directory (parallel)
    # ------------------------------------------------------------------ #
    print(f"\n[3/3] Planting index.php  (workers={workers})")
    failed_index = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(plant_index, d): d for d in eos_dirs}
        for future in as_completed(futures):
            eos_dir, ok = future.result()
            if not ok:
                failed_index.append(eos_dir)

    # ------------------------------------------------------------------ #
    # Summary
    # ------------------------------------------------------------------ #
    print()
    if failed_index:
        print(f"WARNING: index.php could not be planted in {len(failed_index)} directories:")
        for d in sorted(failed_index):
            print(f"  {d}")
    else:
        print(f"Done! index.php planted in {len(eos_dirs)} directories.")

    print(f"\nEOS path : {eos_dest}")
    print(f"Web URL  : https://aldufour.web.cern.ch/{eos_subpath}/")


if __name__ == "__main__":
    main()
