#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep

from HiggsAnalysis.AnalyticAnomalousCoupling.utils.scan import scanEFT

plt.style.use(hep.style.CMS)
plt.rcParams.update({
    "font.size":        24,
    "axes.titlesize":   24,
    "axes.labelsize":   24,
    "xtick.labelsize":  22,
    "ytick.labelsize":  22,
    "legend.fontsize":  22,
})

# -------------------------
# Helper functions (yours)
# -------------------------

def getBestFit(graphScan):
    x = list(graphScan.GetX())
    y = list(graphScan.GetY())
    x, y = zip(*sorted(zip(x, y)))
    y_b, x_b = zip(*sorted(zip(y, x)))
    return x_b[0]

def getLSintersections(graphScan, val):
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
        xings = [min(x), max(x)]

    return xings[:2]


# -------------------------
# Extract intervals
# -------------------------

def extract_intervals(filepath, poi, maxNLL=10, isNuis=False):

    scanUtil = scanEFT()
    scanUtil.setFile(filepath)
    scanUtil.setTree("limit")
    scanUtil.setPOI("k_" + poi)
    scanUtil.setupperNLLimit(maxNLL)
    scanUtil.setNuisanceStyle(isNuis)

    gs = scanUtil.getScan()   # <-- correct TGraph

    x_b = getBestFit(gs)
    x1 = getLSintersections(gs, 1.0)
    x2 = getLSintersections(gs, 4.0)

    return {
        "best": x_b,
        "1sigma": x1,
        "2sigma": x2
    }


# -------------------------
# Scan directories
# -------------------------

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--indir", default=".")
parser.add_argument("--ops", nargs="+", default=None)
parser.add_argument("--horizontal", action="store_true",
                    help="Slide-friendly layout: operators on x-axis, panels stacked vertically")
args = parser.parse_args()

base_dir = args.indir

import glob, re
if args.ops:
    operators = sorted(args.ops)
else:
    # discover from full-scan files: higgsCombine.{op}.individual.MultiDimFit.mH125.root
    pattern = os.path.join(base_dir, "higgsCombine.*.individual.MultiDimFit.mH125.root")
    operators = sorted([
        re.search(r"higgsCombine\.(.+)\.individual", f).group(1)
        for f in glob.glob(pattern)
        if "_stat" not in f
    ])

results = {}

for op in operators:

    poi = op

    file_mc   = f"{base_dir}/higgsCombine.{op}.individual.MultiDimFit.mH125.root"
    file_nomc = f"{base_dir}/higgsCombine.{op}_stat.individual.MultiDimFit.mH125.root"

    if not (os.path.exists(file_mc) and os.path.exists(file_nomc)):
        print(f"Skipping {op} (missing files)")
        continue

    results[op] = {
        "MC": extract_intervals(file_mc, poi),
        "noMC": extract_intervals(file_nomc, poi)
    }


# -------------------------
# Plot
# -------------------------

n_ops = len(results)
pos = np.arange(n_ops)

if args.horizontal:
    # Slide-friendly: operators on x-axis, panels stacked top/bottom
    fig, (ax, ax2) = plt.subplots(
        nrows=2,
        figsize=(max(10, 0.9 * n_ops), 10),
        gridspec_kw={"height_ratios": [2.5, 1.8]},
        sharex=True
    )
else:
    fig, (ax, ax2) = plt.subplots(
        ncols=2,
        figsize=(12, 0.6 * n_ops),
        gridspec_kw={"width_ratios": [2.5, 1]},
        sharey=True
    )

for i, (op, res) in enumerate(results.items()):

    # =========================
    # FIRST PANEL (intervals)
    # =========================

    for shift, key, color in [(+0.15, "MC", "blue"), (-0.15, "noMC", "red")]:
        r = res[key]
        x_b = r["best"]
        x1   = r["1sigma"]
        x2s  = r["2sigma"]

        if args.horizontal:
            p = i + shift   # position along x (operator axis)
            ax.vlines(p, x2s[0], x2s[1], colors=color, linestyles='dashed', linewidth=1.5)
            ax.vlines(p, x1[0],  x1[1],  colors=color, linestyles='solid',  linewidth=3)
            ax.plot(p, x_b, 'o', color=color)
        else:
            p = i + shift   # position along y
            ax.hlines(p, x2s[0], x2s[1], colors=color, linestyles='dashed', linewidth=1.5)
            ax.hlines(p, x1[0],  x1[1],  colors=color, linestyles='solid',  linewidth=3)
            ax.plot(x_b, p, 'o', color=color)

    # =========================
    # SECOND PANEL (Lambda bars)
    # =========================

    for shift, key, color in [(+0.15, "MC", "blue"), (-0.15, "noMC", "red")]:
        r = res[key]
        p = i + shift

        a = abs(r["2sigma"][0]) + abs(r["2sigma"][1])
        if a <= 0:
            continue

        lam1   = np.sqrt(1.0 / a)
        lam4pi = np.sqrt((4*np.pi)**2 / a)

        if args.horizontal:
            ax2.bar(p, lam1,            width=0.25, color=color, alpha=0.8, label=r"$c=1$")
            ax2.bar(p, lam4pi - lam1,   width=0.25, bottom=lam1, color=color, alpha=0.3, label=r"$c=(4\pi)^2$")
        else:
            ax2.barh(p, lam1,           height=0.25, color=color, alpha=0.8, label=r"$c=1$")
            ax2.barh(p, lam4pi - lam1,  left=lam1, height=0.25, color=color, alpha=0.3, label=r"$c=(4\pi)^2$")


# =========================
# Formatting FIRST panel
# =========================

if args.horizontal:
    # ticks live on the shared axis — set once, show only on bottom panel
    ax.set_xticks(pos)
    ax.tick_params(axis='x', labelbottom=False)   # hide on top panel
    ax.axhline(0, color='black', linestyle='--', linewidth=1)
    ax.set_xlim(-1.0, n_ops + 0.5)
    ax.set_ylabel("Wilson coefficient")
else:
    ax.set_yticks(pos)
    ax.set_yticklabels(list(results.keys()))
    ax.axvline(0, color='black', linestyle='--', linewidth=1)
    ax.set_ylim(-1.0, n_ops + 2)
    ax.set_xlabel("Wilson coefficient")

# =========================
# Formatting SECOND panel
# =========================

if args.horizontal:
    ax2.set_ylabel(r"$\Lambda$ at 95% CL [TeV]")
    ax2.set_xticklabels(list(results.keys()), rotation=45, ha='left')
    ax2.set_yscale("log")
else:
    ax2.set_xlabel(r"$\Lambda$ at 95% CL [TeV]")
    ax2.tick_params(axis='y', left=False, labelleft=False)
    ax2.set_xscale("log")

# =========================
# Legend
# =========================

from matplotlib.lines import Line2D

handles = [
    Line2D([], [], color='blue', lw=3),
    Line2D([], [], color='blue', lw=1.5, linestyle='dashed'),
    Line2D([], [], color='red', lw=3),
    Line2D([], [], color='red', lw=1.5, linestyle='dashed'),
]

labels = [
    "Stat + Syst 1σ",
    "Stat + Syst 2σ",
    "Stat only 1σ",
    "Stat only 2σ",
]

ax.legend(
    handles, labels,
    ncol=2,
    columnspacing=3.0,
    loc='upper right' if args.horizontal else 'upper center',
    frameon=False
)

from matplotlib.patches import Patch

bar_handles = [
    Patch(facecolor='black', alpha=1.0, label=r'$c=1$'),
    Patch(facecolor='black', alpha=0.3, label=r'$c=(4\pi)^2$'),
]

ax2.legend(
    handles=bar_handles,
    ncol=1,
    loc='upper right',
    frameon=False
)

# =========================
# CMS label
# =========================

hep.cms.label(ax=ax, data=True, label="Preliminary")


# =========================
# Final
# =========================

suffix = "_horizontal" if args.horizontal else ""
plt.tight_layout()
plt.savefig(f"eft_summary_two_panel{suffix}.pdf", bbox_inches='tight')
plt.savefig(f"eft_summary_two_panel{suffix}.png", dpi=150, bbox_inches='tight')

plt.show()