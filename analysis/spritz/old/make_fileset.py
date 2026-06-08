#!/usr/bin/env python3
"""
Build fileset.json for DY SMEFTsim LO NanoAOD samples on T2_FR_GRIF_LLR.

Run inside the spritz apptainer (xrootd + uproot available):
    python3 make_fileset.py

Output: fileset.json in the current directory.
"""

import json
import subprocess
import sys

XRD_SERVER = "eos.grif.fr:1094"
XRD_BASE   = "/eos/grif/cms/llr/store/user/aldufour/3DY_SMEFTsim_LO"
XRD_PREFIX = f"root://{XRD_SERVER}/"

BINS = {
    "DYSMEFTsim_LO_mll_50_120":   ("mll_50_120",   "260504_081708"),
    "DYSMEFTsim_LO_mll_120_200":  ("mll_120_200",  "260504_081714"),
    "DYSMEFTsim_LO_mll_200_400":  ("mll_200_400",  "260504_081725"),
    "DYSMEFTsim_LO_mll_400_600":  ("mll_400_600",  "260504_081732"),
    "DYSMEFTsim_LO_mll_600_800":  ("mll_600_800",  "260504_081739"),
    "DYSMEFTsim_LO_mll_800_1000": ("mll_800_1000", "260504_081745"),
    "DYSMEFTsim_LO_mll_1000_3000":("mll_1000_3000","260504_081752"),
}

CRAB_DATASET_NAME = "DYSMEFTMll-nanoaod18_SMEFTsim_mll_{bin}"


def xrdfs_ls(path):
    """Return list of entries under path on XRD_SERVER."""
    cmd = ["xrdfs", XRD_SERVER, "ls", path]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  WARNING xrdfs ls failed for {path}: {r.stderr.strip()}", file=sys.stderr)
        return []
    return [line.strip() for line in r.stdout.splitlines() if line.strip()]


def get_nevents(xrd_url):
    """Open a NanoAOD file via uproot and return number of events."""
    try:
        import uproot
        with uproot.open({xrd_url: "Events"}) as t:
            return t.num_entries
    except Exception as e:
        print(f"  WARNING could not get nevents for {xrd_url}: {e}", file=sys.stderr)
        return 500  # fallback: nominal events per job


def collect_files(sample_name, mll_bin, timestamp):
    crab_dir = CRAB_DATASET_NAME.format(bin=mll_bin.replace("mll_", ""))
    # CRAB output: XRD_BASE/crab_dir/crab_dir/timestamp/000N/file.root
    top = f"{XRD_BASE}/{crab_dir}/{crab_dir}/{timestamp}"
    subdirs = xrdfs_ls(top)  # e.g. /eos/.../0000, /eos/.../0001, ...
    if not subdirs:
        print(f"  No subdirs found under {top}", file=sys.stderr)
        return []

    files = []
    for subdir in sorted(subdirs):
        entries = xrdfs_ls(subdir)
        root_files = [e for e in entries if e.endswith(".root")]
        for rf in sorted(root_files):
            xrd_url = f"{XRD_PREFIX}{rf}"
            nevents = get_nevents(xrd_url)
            files.append({"path": [xrd_url], "nevents": nevents})
            print(f"    {rf.split('/')[-1]}  ({nevents} events)")

    return files


def main():
    fileset = {}
    for sample_name, (mll_bin, timestamp) in BINS.items():
        print(f"\n=== {sample_name} (bin={mll_bin}, ts={timestamp}) ===")
        files = collect_files(sample_name, mll_bin, timestamp)
        fileset[sample_name] = {
            "query": "",   # no DAS query for private samples
            "files": files,
        }
        print(f"  -> {len(files)} files")

    out = "fileset.json"
    with open(out, "w") as f:
        json.dump(fileset, f, indent=2)
    print(f"\nWritten {out}")


if __name__ == "__main__":
    main()
