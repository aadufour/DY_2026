#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep

from HiggsAnalysis.AnalyticAnomalousCoupling.utils.scan import scanEFT

plt.style.use(hep.style.CMS)

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
fig, (ax, ax2) = plt.subplots(
    ncols=2,
    figsize=(12, 0.6 * len(results)),
    gridspec_kw={"width_ratios": [2.5, 1]},
    sharey=True
)

y_positions = np.arange(len(results))

for i, (op, res) in enumerate(results.items()):

    # =========================
    # LEFT PANEL (intervals)
    # =========================

    # --- WITH MC (blue)
    r = res["MC"]
    y = i + 0.15

    x_b = r["best"]
    x1 = r["1sigma"]
    x2s = r["2sigma"]

    ax.hlines(y, x2s[0], x2s[1], colors='blue', linestyles='dashed', linewidth=1.5)
    ax.hlines(y, x1[0], x1[1], colors='blue', linestyles='solid', linewidth=3)
    ax.plot(x_b, y, 'o', color='blue')

    # --- WITHOUT MC (red)
    r = res["noMC"]
    y = i - 0.15

    x_b = r["best"]
    x1 = r["1sigma"]
    x2s = r["2sigma"]

    ax.hlines(y, x2s[0], x2s[1], colors='red', linestyles='dashed', linewidth=1.5)
    ax.hlines(y, x1[0], x1[1], colors='red', linestyles='solid', linewidth=3)
    ax.plot(x_b, y, 'o', color='red')


    # =========================
    # RIGHT PANEL (Lambda bars)
    # =========================

    for shift, key, color in [(+0.15, "MC", "blue"), (-0.15, "noMC", "red")]:

        r = res[key]
        y = i + shift

        # take upper 95% CL (>0)
        a = abs(r["2sigma"][0]) + abs(r["2sigma"][1])

        if a <= 0:
            continue

        # Lambda values
        lam1 = np.sqrt(1.0 / a)
        lam4pi = np.sqrt((4*np.pi)**2 / a)

        # stacked bar:
        ax2.barh(
            y,
            lam1,
            height=0.25,
            color=color,
            alpha=0.8,
            label=r"$c=1$"
        )

        ax2.barh(
            y,
            lam4pi - lam1,
            left=lam1,
            height=0.25,
            color=color,
            alpha=0.3,
            label=r"$c=(4\pi)^2$"
        )


# =========================
# Formatting LEFT
# =========================

n_ops = len(results)

ax.set_yticks(y_positions)
ax.set_yticklabels(list(results.keys()))

ax.axvline(0, color='black', linestyle='--', linewidth=1)
ax.set_ylim(-1.0, n_ops + 2)

ax.set_xlabel("Wilson coefficient")


# =========================
# Formatting RIGHT
# =========================

ax2.set_xlabel(r"$\Lambda$ at 95% CL [TeV]")

# remove duplicate y labels
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
    fontsize=17,
    columnspacing=3.0,
    loc='upper center',
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
    fontsize=17,
    loc='upper center',
    frameon=False
)

# =========================
# CMS label
# =========================

hep.cms.label(ax=ax, data=True, label="Preliminary")


# =========================
# Final
# =========================

plt.tight_layout()
plt.savefig("eft_summary_two_panel.pdf")
plt.savefig("eft_summary_two_panel.png", dpi=150)

plt.show()