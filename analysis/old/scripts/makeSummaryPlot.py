#! /opt/anaconda3/envs/combine/bin/python3
"""
makeSummaryPlot.py

Read 1D scan ROOT files (higgsCombine.<op>.individual.MultiDimFit.mH125.root)
from the current directory, compute 68% and 95% CL intervals from the NLL
profile, and produce a horizontal summary plot sorted by interval width.

Usage:
    python3 makeSummaryPlot.py
    python3 makeSummaryPlot.py --input /path/to/scan/dir
    python3 makeSummaryPlot.py --output summary.png
    python3 makeSummaryPlot.py --doOnly ced,cHDD
    python3 makeSummaryPlot.py --sort descending
    python3 makeSummaryPlot.py --metadata other.json
"""

import argparse
import glob
import json
import os
import re

import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker

# ---- thresholds on 2*deltaNLL (same convention as mkEFTScan.py) ----
THRESH_1S = 1.0   # 68% CL
THRESH_2S = 4.0   # 95% CL

# ---- Args -------------------------------------------------------

parser = argparse.ArgumentParser(description="Summary plot of 1D EFT scan intervals.")
parser.add_argument("--input",    default=".",                help="Directory containing higgsCombine.*.root files (default: cwd)")
parser.add_argument("--output",   default="summary_plot.png", help="Output plot filename (default: summary_plot.png)")
parser.add_argument("--doOnly",   default="",                 help="Comma-separated list of operators to include (default: all found)")
parser.add_argument("--sort",     default="ascending", choices=["ascending", "descending", "none"],
                                                              help="Sort by 68%% CL interval width (default: ascending = tightest at top)")
parser.add_argument("--metadata", default="metadata.json",   help="metadata.json used to identify 1D operators (default: metadata.json)")
parser.add_argument("--maxNLL",   type=float, default=10.0,  help="Max 2*deltaNLL to read from scan (default: 10)")
args = parser.parse_args()

# ---- Helpers ----------------------------------------------------

def find_crossing(x, y2, threshold):
    """Linear interpolation: find where y2 crosses threshold left and right of the minimum.
    Returns (lo, best, hi); lo or hi may be None if not found within scan range.
    """
    i_min = int(np.argmin(y2))
    best  = float(x[i_min])
    lo = hi = None

    # left of minimum
    xl = x[:i_min+1][::-1]
    yl = y2[:i_min+1][::-1]
    for i in range(len(xl) - 1):
        if yl[i] <= threshold <= yl[i+1]:
            frac = (threshold - yl[i]) / (yl[i+1] - yl[i])
            lo = xl[i] + frac * (xl[i+1] - xl[i])
            break

    # right of minimum
    xr = x[i_min:]
    yr = y2[i_min:]
    for i in range(len(xr) - 1):
        if yr[i] <= threshold <= yr[i+1]:
            frac = (threshold - yr[i]) / (yr[i+1] - yr[i])
            hi = xr[i] + frac * (xr[i+1] - xr[i])
            break

    return lo, best, hi


def read_scan(filepath, poi):
    """Return (x, 2*deltaNLL) arrays from a combine TTree, NLL shifted so min=0."""
    with uproot.open(filepath) as f:
        tree = f["limit"]
        x    = tree[poi].array(library="np")
        dnll = tree["deltaNLL"].array(library="np")

    # entry 0 is the snapshot best-fit point — skip it
    x    = x[1:]
    dnll = dnll[1:]
    y2   = 2.0 * dnll

    # apply maxNLL filter, sort, deduplicate
    mask = y2 < args.maxNLL
    x, y2 = x[mask], y2[mask]
    idx   = np.argsort(x)
    x, y2 = x[idx], y2[idx]
    x, uid = np.unique(x, return_index=True)
    y2 = y2[uid]

    # shift minimum to 0
    y2 = y2 - np.min(y2)
    return x, y2


def process_file(filepath):
    """Return a result dict for one operator, or None on failure."""
    basename = os.path.basename(filepath)
    m = re.match(r"higgsCombine\.(.+)\.individual\.MultiDimFit\.mH125\.root$", basename)
    if not m:
        return None

    name = m.group(1)
    poi  = f"k_{name}"

    try:
        x, y2 = read_scan(filepath, poi)
    except Exception as e:
        print(f"  [SKIP] {basename}: {e}")
        return None

    lo1, best, hi1 = find_crossing(x, y2, THRESH_1S)
    lo2, _,    hi2 = find_crossing(x, y2, THRESH_2S)

    if lo1 is None or hi1 is None:
        print(f"  [WARN] {name}: 68% crossing not found (scan range may be too narrow)")
    if lo2 is None or hi2 is None:
        print(f"  [WARN] {name}: 95% crossing not found (scan range may be too narrow)")

    width1 = (hi1 - lo1) if (lo1 is not None and hi1 is not None) else None

    return dict(op=name, best=best,
                lo1=lo1, hi1=hi1,
                lo2=lo2, hi2=hi2,
                width1=width1)

# ---- Collect files ----------------------------------------------

pattern  = os.path.join(args.input, "higgsCombine.*.individual.MultiDimFit.mH125.root")
allfiles = sorted(glob.glob(pattern))

if not allfiles:
    raise SystemExit(f"No scan files found in '{args.input}' matching the expected pattern.")

# Load metadata to distinguish 1D operators from 2D pairs
meta_ops = set()
if os.path.isfile(args.metadata):
    with open(args.metadata) as fh:
        md = json.load(fh)
    meta_ops = set(md.get("operators", {}).keys())
    print(f"Loaded {len(meta_ops)} operators from {args.metadata}")
else:
    print(f"[WARN] {args.metadata} not found — will attempt all files (may include 2D pairs)")

whitelist = [op.strip() for op in args.doOnly.split(",") if op.strip()] if args.doOnly else []

results = []
for fp in allfiles:
    r = process_file(fp)
    if r is None:
        continue
    # skip multi-operator files (2D etc.) using metadata as reference
    if meta_ops and r["op"] not in meta_ops:
        continue
    if whitelist and r["op"] not in whitelist:
        continue
    results.append(r)

if not results:
    raise SystemExit("No valid 1D scan results found.")

# ---- Print table ------------------------------------------------

print(f"\nProcessed {len(results)} operators:")
print(f"  {'operator':<12}  {'best':>8}  {'68% CL interval':>22}  {'95% CL interval':>22}  {'width68':>8}")
for r in results:
    w   = f"{r['width1']:.5f}" if r['width1'] is not None else "N/A"
    lo1 = f"{r['lo1']:+.5f}" if r['lo1'] is not None else "   N/A"
    hi1 = f"{r['hi1']:+.5f}" if r['hi1'] is not None else "   N/A"
    lo2 = f"{r['lo2']:+.5f}" if r['lo2'] is not None else "   N/A"
    hi2 = f"{r['hi2']:+.5f}" if r['hi2'] is not None else "   N/A"
    print(f"  {r['op']:<12}  {r['best']:>+8.5f}  [{lo1}, {hi1}]  [{lo2}, {hi2}]  {w:>8}")

# ---- Sort -------------------------------------------------------

def sort_key(r):
    return r["width1"] if r["width1"] is not None else float("inf")

if args.sort == "ascending":
    results.sort(key=sort_key)        # tightest at top
elif args.sort == "descending":
    results.sort(key=sort_key, reverse=True)

# ---- Plot -------------------------------------------------------

n  = len(results)
ys = np.arange(n)

fig_h = max(4, 0.45 * n + 1.5)
fig, ax = plt.subplots(figsize=(8, fig_h))

COLOR_2S = "#FFCC00"   # yellow — 95% CL
COLOR_1S = "#00CC00"   # green  — 68% CL
COLOR_BF = "black"     # best-fit marker

for i, r in enumerate(results):
    y = ys[i]

    if r["lo2"] is not None and r["hi2"] is not None:
        ax.barh(y, r["hi2"] - r["lo2"], left=r["lo2"],
                height=0.6, color=COLOR_2S, zorder=2)

    if r["lo1"] is not None and r["hi1"] is not None:
        ax.barh(y, r["hi1"] - r["lo1"], left=r["lo1"],
                height=0.6, color=COLOR_1S, zorder=3)

    ax.plot(r["best"], y, marker="o", ms=5, color=COLOR_BF, zorder=4)

ax.axvline(0.0, color="gray", linestyle="--", linewidth=0.8, zorder=1)

ax.set_yticks(ys)
ax.set_yticklabels([r["op"] for r in results], fontsize=9)
ax.set_ylim(-0.8, n - 0.2)
ax.set_xlabel("k", fontsize=11)
ax.set_title("", fontsize=12)
# ax.set_title("EFT operator intervals (1D, Asimov)", fontsize=12)

patch_2s  = mpatches.Patch(color=COLOR_2S, label="95% CL")
patch_1s  = mpatches.Patch(color=COLOR_1S, label="68% CL")
marker_bf = plt.Line2D([0], [0], marker="o", color="w",
                        markerfacecolor=COLOR_BF, markersize=6, label="Best fit")
ax.legend(handles=[patch_2s, patch_1s, marker_bf],
          loc="lower right", fontsize=9, framealpha=0.8)

ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
ax.grid(axis="x", which="both", linestyle=":", linewidth=0.5, alpha=0.5)

plt.tight_layout()
plt.savefig(args.output, dpi=150)
print(f"\nSaved: {args.output}")

pdf_path = os.path.splitext(args.output)[0] + ".pdf"
plt.savefig(pdf_path)
print(f"Saved: {pdf_path}")
