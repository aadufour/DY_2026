#!/usr/bin/env python3
"""
merge_eft2d.py

Split 2D production: combine the merged results of the background run (propcorr_v3_pdfcut: data, backgrounds,
MiNNLO, runner_pdfcut.py) with the signal-only run (propcorr_v3_eft2d: EFT datasets with all 406 points,
runner_2d.py). EFT datasets are taken from the signal run only; everything else from the background run.
Writes <out>/condor/results_merged_new.pkl and <out>/config.py (signal config with SIGNAL_ONLY = False),
ready for spritz-postproc / spritz-cards in <out>.

Usage (spritz-env.sif, from spritz_fabian/configs):
    python3 merge_eft2d.py --bkg propcorr_v3_pdfcut --sig propcorr_v3_eft2d --out propcorr_v3_2d
"""
import argparse
import os

from spritz.framework.framework import read_chunks, write_chunks

EFT_PREFIX = "DYSMEFTsim_LO_propcorr"

p = argparse.ArgumentParser()
p.add_argument("--bkg", required=True)
p.add_argument("--sig", required=True)
p.add_argument("--out", required=True)
a = p.parse_args()

if os.path.exists(f"{a.out}/condor/results_merged_new.pkl"):
    raise SystemExit(f"{a.out}/condor/results_merged_new.pkl exists, not overwriting")

bkg = read_chunks(f"{a.bkg}/condor/results_merged_new.pkl")
sig = read_chunks(f"{a.sig}/condor/results_merged_new.pkl")
assert all(k.startswith(EFT_PREFIX) for k in sig), "signal run contains non-EFT datasets"
n_sm_bkg = len([k for k in bkg if k.startswith(EFT_PREFIX)])
merged = {k: v for k, v in bkg.items() if not k.startswith(EFT_PREFIX)}
merged.update(sig)
print(f"bkg: {len(bkg)} entries ({n_sm_bkg} EFT dropped), sig: {len(sig)} entries, merged: {len(merged)}")

os.makedirs(f"{a.out}/condor", exist_ok=True)
cfg = open(f"{a.sig}/config.py").read()
assert cfg.count("SIGNAL_ONLY = True") == 1
open(f"{a.out}/config.py", "w").write(cfg.replace("SIGNAL_ONLY = True", "SIGNAL_ONLY = False"))
if os.path.exists(f"{a.bkg}/data") and not os.path.exists(f"{a.out}/data"):
    os.symlink(os.path.abspath(f"{a.bkg}/data"), f"{a.out}/data")
write_chunks(merged, f"{a.out}/condor/results_merged_new.pkl")
print("written", f"{a.out}/condor/results_merged_new.pkl")
