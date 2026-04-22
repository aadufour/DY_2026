#!/usr/bin/env python3
"""
pushtoEOS.py  <local_dir>  <eos_www_subpath>

Copies a local directory to the EOS www area and plants index.php in every
subdirectory so the CMS PHP plot browser works.

Requires: xrdcp and eos commands available (run cmsenv first).

EOS www base: /eos/user/a/aldufour/www/

Example:
    python pushtoEOS.py plots/syst_scan_v2  syst_scan_v2
    python pushtoEOS.py results/figures     DY2026/figures
"""

import os
import sys
import subprocess
from pathlib import Path

EOS_WWW_BASE = "/eos/user/a/aldufour/www"
EOS_MGM_URL  = "root://eosuser.cern.ch"
INDEX_PHP    = Path(__file__).parent / "index.php"


def xrd(eos_abs_path: str) -> str:
    return f"{EOS_MGM_URL}/{eos_abs_path}"


def run(cmd: list[str], check: bool = True) -> subprocess.CompletedProcess:
    print("  $", " ".join(cmd))
    result = subprocess.run(cmd, check=check, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.stdout.strip():
        print(result.stdout, end="")
    if result.stderr.strip():
        print(result.stderr, end="", file=sys.stderr)
    return result


def main() -> None:
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <local_dir> <eos_www_subpath>")
        print(f"  e.g. {sys.argv[0]} plots/syst_scan_v2  syst_scan_v2")
        sys.exit(1)

    local_dir   = Path(sys.argv[1]).resolve()
    eos_subpath = sys.argv[2].strip("/")
    eos_dest    = f"{EOS_WWW_BASE}/{eos_subpath}"

    if not local_dir.is_dir():
        sys.exit(f"Error: '{local_dir}' is not a directory")
    if not INDEX_PHP.is_file():
        sys.exit(f"Error: index.php not found at '{INDEX_PHP}'")

    # Make EOS_MGM_URL visible to child processes even if caller forgot to export it
    os.environ["EOS_MGM_URL"] = EOS_MGM_URL

    # ------------------------------------------------------------------ #
    # 1. Check destination doesn't already exist, then create it
    # ------------------------------------------------------------------ #
    result = run(["eos", "ls", eos_dest], check=False)
    if result.returncode == 0:
        sys.exit(f"Error: '{eos_dest}' already exists on EOS — refusing to overwrite.")

    print(f"\n[1/3] Uploading '{local_dir}' -> {eos_dest}")
    # xrdcp -r requires dest to exist, but if dest exists it nests src inside it.
    # Solution: create dest, then copy each item in local_dir into it individually.
    run(["eos", "mkdir", "-p", eos_dest])
    for item in sorted(local_dir.iterdir()):
        run(["xrdcp", "-r", "--silent", str(item), xrd(eos_dest + "/")])

    # ------------------------------------------------------------------ #
    # 2. Collect every directory that was just created on EOS by walking
    #    the local tree (structure is identical after the upload).
    # ------------------------------------------------------------------ #
    print(f"\n[2/3] Collecting directories to index...")
    eos_dirs: list[str] = []

    for root, dirs, _ in os.walk(local_dir):
        dirs.sort()  # deterministic order
        rel = Path(root).relative_to(local_dir)
        if str(rel) == ".":
            eos_dirs.append(eos_dest)
        else:
            eos_dirs.append(f"{eos_dest}/{rel}")

    print(f"  Found {len(eos_dirs)} directories")

    # ------------------------------------------------------------------ #
    # 3. Plant index.php in every directory
    # ------------------------------------------------------------------ #
    print(f"\n[3/3] Planting index.php ...")
    failed = []
    for eos_dir in eos_dirs:
        result = run(["xrdcp", "--silent", str(INDEX_PHP), xrd(f"{eos_dir}/index.php")],
                     check=False)
        if result.returncode != 0:
            failed.append(eos_dir)

    # ------------------------------------------------------------------ #
    # Summary
    # ------------------------------------------------------------------ #
    print()
    if failed:
        print(f"WARNING: index.php could not be planted in {len(failed)} directories:")
        for d in failed:
            print(f"  {d}")
    else:
        print(f"Done! index.php planted in {len(eos_dirs)} directories.")

    print(f"\nEOS path : {eos_dest}")
    print(f"Web URL  : https://aldufour.web.cern.ch/{eos_subpath}/")


if __name__ == "__main__":
    main()
