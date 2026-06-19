#!/usr/bin/env python3
"""
makeSummaryThreeVar.py

Summary plot overlaying three observables (mll, rapll_abs, costhetastar)
on the same panel.  Operators are ordered by mll 95% CL sensitivity
(tightest first).  Only stat+syst scans are used.

Usage (run from the configs dir, e.g. eft_bkg_fullsyst_v2/datacards/inc_mm):

    makeSummaryThreeVar.py \
        --mll           mll/scan        \
        --rapll         rapll_abs/scan  \
        --costhetastar  costhetastar/scan
"""

import os
import glob
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep

from HiggsAnalysis.AnalyticAnomalousCoupling.utils.scan import scanEFT

# ============================================================
# PLOT CONFIGURATION
# ============================================================

FONT_SIZE       = 34
LABEL_SIZE      = 26
TICK_LABELSIZE  = 32
LEGEND_FONTSIZE = 20

FIG_HEIGHT      = 10      # inches
WIDTH_PER_OP    = 0.9     # inches per operator

# colours and display names for the three variables
VARS = {
    "mll":          {"color": "#2166ac", "label": r"$m_{\ell\ell}$"},
    "rapll":        {"color": "#d6604d", "label": r"$|y_{\ell\ell}|$"},
    "costhetastar": {"color": "#4dac26", "label": r"$\cos\theta^*$"},
}

# horizontal offset of each variable's markers within one operator slot
SHIFTS = {"mll": 0.0, "rapll": +0.22, "costhetastar": -0.22}

# height_ratios top (intervals) : bottom (Lambda)
HEIGHT_RATIOS = [2.5, 1.8]

# ============================================================

plt.style.use(hep.style.CMS)
plt.rcParams.update({
    "font.size":       FONT_SIZE,
    "axes.labelsize":  LABEL_SIZE,
    "xtick.labelsize": TICK_LABELSIZE,
    "ytick.labelsize": TICK_LABELSIZE,
    "legend.fontsize": LEGEND_FONTSIZE,
})

# -------------------------
# Helpers (same logic as makeSummary.py)
# -------------------------

def getBestFit(graphScan):
    x = list(graphScan.GetX())
    y = list(graphScan.GetY())
    x, y = zip(*sorted(zip(x, y)))
    y_b, x_b = zip(*sorted(zip(y, x)))
    return x_b[0]


def getLSintersections(graphScan, val, x_b=None):
    xings = []
    n = graphScan.GetN()
    x = list(graphScan.GetX())
    y = list(graphScan.GetY())
    x, y = zip(*sorted(zip(x, y)))

    for i in range(n):
        if y[i] == val:
            xings.append(x[i])
            continue
        if i > 0:
            if (y[i] - val) * (y[i-1] - val) < 0:
                xings.append(
                    x[i-1] + abs((y[i-1] - val) * (x[i] - x[i-1]) / (y[i] - y[i-1]))
                )

    if len(xings) < 2:
        return [min(x), max(x)]

    if x_b is not None:
        lo = [xi for xi in xings if xi <= x_b]
        hi = [xi for xi in xings if xi >= x_b]
        if lo and hi:
            return [max(lo), min(hi)]

    return xings[:2]


def extract_intervals(filepath, poi, maxNLL=10):
    scanUtil = scanEFT()
    scanUtil.setFile(filepath)
    scanUtil.setTree("limit")
    scanUtil.setPOI("k_" + poi)
    scanUtil.setupperNLLimit(maxNLL)
    scanUtil.setNuisanceStyle(False)

    gs  = scanUtil.getScan()
    x_b = getBestFit(gs)
    x1  = getLSintersections(gs, 1.0, x_b)
    x2  = getLSintersections(gs, 4.0, x_b)

    return {"best": x_b, "1sigma": x1, "2sigma": x2}


def discover_operators(scan_dir):
    pattern = os.path.join(scan_dir, "higgsCombine.*.individual.MultiDimFit.mH125.root")
    ops = []
    for f in glob.glob(pattern):
        m = re.search(r"higgsCombine\.(.+)\.individual", f)
        if m and "_stat" not in f:
            ops.append(m.group(1))
    return sorted(ops)


def load_results(scan_dir, operators):
    results = {}
    for op in operators:
        fp = os.path.join(scan_dir, f"higgsCombine.{op}.individual.MultiDimFit.mH125.root")
        if not os.path.exists(fp):
            print(f"  [skip] {op}: not found in {scan_dir}")
            continue
        try:
            results[op] = extract_intervals(fp, op)
        except Exception as e:
            print(f"  [skip] {op}: {e}")
    return results


# -------------------------
# Args
# -------------------------

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--mll",           required=True, help="Directory with mll scan ROOT files")
parser.add_argument("--rapll",         required=True, help="Directory with rapll_abs scan ROOT files")
parser.add_argument("--costhetastar",  required=True, help="Directory with costhetastar scan ROOT files")
parser.add_argument("--ops",  nargs="+", default=None, help="Subset of operators (default: all found in --mll dir)")
args = parser.parse_args()

# -------------------------
# Load
# -------------------------

dirs = {"mll": args.mll, "rapll": args.rapll, "costhetastar": args.costhetastar}

if args.ops:
    operators = sorted(args.ops)
else:
    operators = discover_operators(args.mll)

print(f"Operators found in mll dir: {len(operators)}")

all_results = {}
for var, d in dirs.items():
    all_results[var] = load_results(d, operators)

# keep only operators present in mll (the reference for ordering)
operators = [op for op in operators if op in all_results["mll"]]

# sort by mll 95% CL (tightest first)
operators.sort(key=lambda op: max(
    abs(all_results["mll"][op]["2sigma"][0]),
    abs(all_results["mll"][op]["2sigma"][1])
))

n_ops = len(operators)
pos   = np.arange(n_ops)

# -------------------------
# Plot
# -------------------------

fig_width = max(10, WIDTH_PER_OP * n_ops)
fig, (ax, ax2) = plt.subplots(
    nrows=2,
    figsize=(fig_width, FIG_HEIGHT),
    gridspec_kw={"height_ratios": HEIGHT_RATIOS},
    sharex=True,
)

for i, op in enumerate(operators):
    for var in ("mll", "rapll", "costhetastar"):
        if op not in all_results[var]:
            continue

        r     = all_results[var][op]
        color = VARS[var]["color"]
        p     = i + SHIFTS[var]

        x_b  = r["best"]
        x1   = r["1sigma"]
        x2s  = r["2sigma"]

        # ---- top panel: intervals ----
        ax.vlines(p, x2s[0], x2s[1], colors=color, linestyles="dashed", linewidth=1.5)
        ax.vlines(p, x1[0],  x1[1],  colors=color, linestyles="solid",  linewidth=3)
        ax.plot(p, x_b, "o", color=color, markersize=5)

        # ---- bottom panel: Lambda bars ----
        a = abs(x2s[0]) + abs(x2s[1])
        if a <= 0:
            continue
        lam1   = np.sqrt(1.0 / a)
        lam4pi = np.sqrt((4 * np.pi)**2 / a)

        ax2.bar(p, lam1,          width=0.18, color=color, alpha=0.9)
        ax2.bar(p, lam4pi - lam1, width=0.18, bottom=lam1, color=color, alpha=0.3)

# -------------------------
# Formatting top panel
# -------------------------

ax.set_xticks(pos)
ax.tick_params(axis="x", labelbottom=False)
ax.axhline(0, color="black", linestyle="--", linewidth=1)
ax.set_xlim(-1.0, n_ops + 0.5)
ax.set_ylabel("Wilson coefficient")

# -------------------------
# Formatting bottom panel
# -------------------------

ax2.set_xticks(pos)
ax2.set_xticklabels(operators, rotation=45, ha="right", rotation_mode="anchor")
ax2.set_ylabel(r"$\Lambda$ at 95% CL [TeV]")
ax2.set_yscale("log")

# -------------------------
# Legends
# -------------------------

from matplotlib.lines import Line2D
from matplotlib.patches import Patch

interval_handles = [
    Line2D([], [], color=VARS[v]["color"], lw=3,   label=VARS[v]["label"]) for v in VARS
] + [
    Line2D([], [], color="grey", lw=3,            label=r"$1\sigma$"),
    Line2D([], [], color="grey", lw=1.5, ls="--", label=r"$2\sigma$"),
]
ax.legend(handles=interval_handles, ncol=3, frameon=False,
          loc="upper right", columnspacing=2.0)

lambda_handles = [
    Patch(facecolor="grey", alpha=0.9, label=r"$c=1$"),
    Patch(facecolor="grey", alpha=0.3, label=r"$c=(4\pi)^2$"),
]
ax2.legend(handles=lambda_handles, ncol=1, frameon=False, loc="upper right")

# -------------------------
# CMS label + save
# -------------------------

hep.cms.label(ax=ax, data=True, label="Preliminary")

plt.tight_layout()
plt.savefig("eft_summary_three_var.pdf", bbox_inches="tight")
plt.savefig("eft_summary_three_var.png", dpi=150, bbox_inches="tight")
print("Saved: eft_summary_three_var.pdf / .png")
plt.show()
