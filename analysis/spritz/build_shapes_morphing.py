#!/usr/bin/env python3
"""
build_shapes_morphing.py
========================
Converts v6 histos.root (DYSMEFTsim naming) to a shapes.root + datacard.txt
compatible with AnomalousCouplingMorphing.py (template_morphing branch).

Process naming convention expected by the physics model:
  histo_sm        = SM template  w(0)
  histo_w1_{op}   = c=+1 template  w(+1)  (SM + Lin + Quad)
  histo_wm1_{op}  = c=-1 template  w(-1)  (SM - Lin + Quad)
  histo_Data      = Asimov data_obs (= SM)

Usage (from analysis_venv):
    python3 build_shapes_morphing.py
    python3 build_shapes_morphing.py --region inc_ee --variable mll
    python3 build_shapes_morphing.py --operators cHDD cHWB
"""

import argparse
import os

import uproot
import numpy as np

# --------------------------------------------------
# Operators (27, matching v6 config)
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

# --------------------------------------------------
# Args
# --------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--input",    default="/grid_mnt/data__data.polcms/cms/adufour/spritz/configs/dy_smeftsim_v6/histos.root")
parser.add_argument("--outdir",   default="/grid_mnt/data__data.polcms/cms/adufour/spritz/configs/dy_smeftsim_v6/datacards_morphing")
parser.add_argument("--region",   default="inc_mm", choices=["inc_ee", "inc_mm", "inc_em"])
parser.add_argument("--variable", default="mll")
parser.add_argument("--operators", nargs="+", default=None,
                    help="Subset of operators. Default: all 27.")
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

def get_vals(sample):
    return f[f"{base_path}/histo_{sample}"].values()

sm_vals = get_vals("DYSMEFTsim_SM")
edges   = f[f"{base_path}/histo_DYSMEFTsim_SM"].axes[0].edges()

# --------------------------------------------------
# Build process list and histogram map
# --------------------------------------------------
# sm = index 0 (signal); EFT templates = negative indices
processes  = ["sm"]
histo_map  = {"sm": sm_vals}

for op in ops:
    try:
        p1 = get_vals(f"DYSMEFTsim_{op}")
        m1 = get_vals(f"DYSMEFTsim_{op}_m1")
    except Exception:
        print(f"  WARNING: {op} not found in histos.root, skipping.")
        continue
    processes.append(f"w1_{op}")
    processes.append(f"wm1_{op}")
    histo_map[f"w1_{op}"]  = p1
    histo_map[f"wm1_{op}"] = m1

n_proc       = len(processes)
# sm=0, then -1, -2, ... for EFT components
proc_indices = [0] + list(range(-1, -(n_proc), -1))

# --------------------------------------------------
# Write shapes.root
# --------------------------------------------------
shapes_path = os.path.join(outdir, "shapes.root")
with uproot.recreate(shapes_path) as out:
    for name, vals in histo_map.items():
        out[f"histo_{name}"] = (vals, edges)
    out["histo_Data"] = (sm_vals, edges)  # Asimov: data_obs = SM

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
    f"observation {sm_vals.sum():.6f}",
    f"shapes  *         * shapes.root  histo_$PROCESS histo_$PROCESS_$SYSTEMATIC",
    f"shapes  data_obs  * shapes.root  histo_Data",
    "-" * 100,
    "bin"     + sep + sep.join([bin_name]  * n_proc),
    "process" + sep + sep.join(processes),
    "process" + sep + sep.join([str(i) for i in proc_indices]),
    "rate"    + sep + sep.join(rates),
    "-" * 100,
    f"{bin_name} autoMCStats 10 0 1",
]

datacard_path = os.path.join(outdir, "datacard.txt")
with open(datacard_path, "w") as dc:
    dc.write("\n".join(lines) + "\n")

print(f"Written : {datacard_path}")
print(f"Processes ({n_proc}): sm + {len(ops)} ops × 2 = {1 + len(ops)*2}")
print()
print("Next steps (after dy_combine_morphing):")
print(f"  cd {outdir}")
print(f"  text2workspace.py datacard.txt \\")
print(f"    -P HiggsAnalysis.AnalyticAnomalousCoupling.AnomalousCouplingMorphing:analiticAnomalousCouplingMorphing \\")
print(f"    --PO eftOperators={','.join(ops)} \\")
print(f"    -o workspace.root")
print(f"  combine -M MultiDimFit workspace.root --algo=grid --points=100 \\")
print(f"    -P k_cHDD --setParameterRanges k_cHDD=-5,5 --redefineSignalPOIs k_cHDD")
