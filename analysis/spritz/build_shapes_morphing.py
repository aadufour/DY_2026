#!/usr/bin/env python3
"""
build_shapes_morphing.py
========================
Converts histos.root (spritz-postproc output) to a shapes.root + datacard.txt
compatible with AnomalousCouplingMorphing.py (template_morphing branch).

Expected histogram names in histos.root (v7 config, morphing convention):
  histo_sm        = SM template  w(0)
  histo_w1_{op}   = c=+1 template  w(+1)
  histo_wm1_{op}  = c=-1 template  w(-1)

Output:
  shapes.root   — same histograms, plus histo_Data = real data (falls back to SM Asimov if absent)
  datacard.txt  — combine datacard with process/rate lines

Usage (from analysis_venv):
    python3 build_shapes_morphing.py
    python3 build_shapes_morphing.py --region inc_mm --variable mll
    python3 build_shapes_morphing.py --operators cHDD cHWB
"""

import argparse
import os

import uproot
import numpy as np

# --------------------------------------------------
# Operators (27, matching v7 config)
# --------------------------------------------------
OPERATORS = [
    "cHDD", "cHWB", "cbWRe", "cbBRe",
    "cHj1", "cHQ1", "cHj3",  "cHQ3",
    "cHu",  "cHd",  "cHbq",
    "cHl1", "cHl3", "cHe",
    "cll1",
    "clj1", "clj3",
    "cQl1", "cQl3",
    "ceu",  "ced",
    "cbe",  "cje",  "cQe",
    "clu",  "cld",  "cbl",
]

# Theory systematics written by post_process_eft.py
THEORY_SYSTS = ["QCDscale", "PDF"]  # → histo_{proc}_QCDscaleUp/Down, histo_{proc}_PDFUp/Down

# --------------------------------------------------
# Args
# --------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--input",    default="/grid_mnt/data__data.polcms/cms/adufour/spritz/configs/dy_smeftsim_v7/histos.root")
parser.add_argument("--outdir",   default="/grid_mnt/data__data.polcms/cms/adufour/spritz/configs/dy_smeftsim_v7/datacards_morphing")
parser.add_argument("--region",   default="inc_mm", choices=["inc_ee", "inc_mm", "inc_em"])
parser.add_argument("--variable", default="mll")
parser.add_argument("--operators", nargs="+", default=None,
                    help="Subset of operators. Default: all 27.")
parser.add_argument("--no-theory", action="store_true",
                    help="Skip theory systematics (QCDscale, PDF) even if present in histos.root.")
args = parser.parse_args()

ops      = args.operators if args.operators else OPERATORS
bin_name = f"{args.region}_{args.variable}"
outdir   = os.path.join(args.outdir, args.region, args.variable)
os.makedirs(outdir, exist_ok=True)

# --------------------------------------------------
# Read histos.root
# --------------------------------------------------
f         = uproot.open(args.input)
base_path = f"{args.region}/{args.variable}/nominal"

def get_vals(name, silent=False):
    key = f"{base_path}/histo_{name}"
    if key not in f:
        if not silent:
            print(f"  WARNING: {key} not found in ROOT file.")
        return None
    return f[key].values()

sm_vals   = get_vals("sm")
edges     = f[f"{base_path}/histo_sm"].axes[0].edges()
data_vals = get_vals("Data", silent=True)
if data_vals is None:
    print("  WARNING: histo_Data not found in histos.root — falling back to Asimov (SM)")
    data_vals = sm_vals

# --------------------------------------------------
# Build process list and histogram map
# --------------------------------------------------
# sm = index 0 (signal); EFT templates = negative indices
processes = ["sm"]
histo_map = {"sm": sm_vals}          # nominal only
syst_map  = {}                        # {(proc, syst, "Up"/"Down"): vals}

for op in ops:
    p1 = get_vals(f"w1_{op}", silent=True)
    m1 = get_vals(f"wm1_{op}", silent=True)
    if p1 is None or m1 is None:
        print(f"  WARNING: {op} not found in histos.root, skipping.")
        continue
    processes.append(f"w1_{op}")
    processes.append(f"wm1_{op}")
    histo_map[f"w1_{op}"]  = p1
    histo_map[f"wm1_{op}"] = m1

n_proc       = len(processes)
# AnomalousCouplingMorphing convention:
#   sm=1 (background/reference), w1_op1=0, wm1_op1=-1, w1_op2=-2, ...
# combine requires ≥1 positive index (background); sm fills that role.
proc_indices = [1] + list(range(0, -(n_proc - 1), -1))

# --------------------------------------------------
# Collect theory Up/Down histograms for each process
# --------------------------------------------------
active_systs = []   # systs present in histos.root for at least one process
if not args.no_theory:
    for syst in THEORY_SYSTS:
        found_any = False
        for proc in processes:
            up   = get_vals(f"{proc}_{syst}Up",   silent=True)
            down = get_vals(f"{proc}_{syst}Down", silent=True)
            if up is not None and down is not None:
                syst_map[(proc, syst, "Up")]   = up
                syst_map[(proc, syst, "Down")] = down
                found_any = True
            else:
                # Fall back to nominal (variation = 0) — combine treats it as unconstrained
                if found_any:
                    print(f"  WARNING: {syst} Up/Down missing for {proc}, using nominal as fallback.")
                syst_map[(proc, syst, "Up")]   = histo_map[proc].copy()
                syst_map[(proc, syst, "Down")] = histo_map[proc].copy()
        if found_any:
            active_systs.append(syst)
            print(f"  Theory syst '{syst}' found for at least one process.")

# --------------------------------------------------
# Write shapes.root
# --------------------------------------------------
shapes_path = os.path.join(outdir, "shapes.root")
with uproot.recreate(shapes_path) as out:
    # Nominal
    for name, vals in histo_map.items():
        out[f"histo_{name}"] = (vals, edges)
    out["histo_Data"] = (data_vals, edges)
    # Theory Up/Down
    for (proc, syst, ud), vals in syst_map.items():
        out[f"histo_{proc}_{syst}{ud}"] = (vals, edges)

print(f"Written : {shapes_path}")

# --------------------------------------------------
# Write datacard.txt
# --------------------------------------------------
sep   = "\t"
rates = [f"{histo_map[p].sum():.6f}" for p in processes]

lines = [
    "## EFT morphing datacard — AnomalousCouplingMorphing.py (template_morphing)",
    "imax 1 number of channels",
    "jmax * number of background",
    "kmax * number of nuisance parameters",
    "-" * 100,
    f"bin         {bin_name}",
    f"observation {data_vals.sum():.6f}",
    f"shapes  *         * shapes.root  histo_$PROCESS histo_$PROCESS_$SYSTEMATIC",
    f"shapes  data_obs  * shapes.root  histo_Data",
    "-" * 100,
    "bin"     + sep + sep.join([bin_name]  * n_proc),
    "process" + sep + sep.join(processes),
    "process" + sep + sep.join([str(i) for i in proc_indices]),
    "rate"    + sep + sep.join(rates),
    "-" * 100,
]

# Theory shape nuisances: one column per process, value = 1.0 if syst present else "-"
for syst in active_systs:
    cols = []
    for proc in processes:
        up   = get_vals(f"{proc}_{syst}Up",   silent=True)
        cols.append("1.0" if up is not None else "-")
    lines.append(f"{syst}{sep}shape{sep}" + sep.join(cols))

lines.append(f"{bin_name} autoMCStats 10 0 1")

datacard_path = os.path.join(outdir, "datacard.txt")
with open(datacard_path, "w") as dc:
    dc.write("\n".join(lines) + "\n")

print(f"Written : {datacard_path}")
print(f"Processes ({n_proc}): sm + {len(ops)} ops × 2 = {1 + len(ops)*2}")
if active_systs:
    print(f"Theory systs written: {active_systs}")
print()
print("Next steps (after dy_combine_morphing):")
print(f"  cd {outdir}")
print(f"  text2workspace.py datacard.txt \\")
print(f"    -P HiggsAnalysis.AnalyticAnomalousCoupling.AnomalousCouplingMorphing:analiticAnomalousCouplingMorphing \\")
print(f"    --PO eftOperators={','.join(ops)} \\")
print(f"    -o workspace.root")
print(f"  combine -M MultiDimFit workspace.root --algo=grid --points=100 \\")
print(f"    -P k_cHDD --setParameterRanges k_cHDD=-5,5 --redefineSignalPOIs k_cHDD")
