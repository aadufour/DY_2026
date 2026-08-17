#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep

from HiggsAnalysis.AnalyticAnomalousCoupling.utils.scan import scanEFT

plt.style.use(hep.style.CMS)
plt.rcParams.update({
    "font.size":        34,
    "axes.titlesize":   26,
    "axes.labelsize":   26,
    "xtick.labelsize":  32,
    "ytick.labelsize":  32,
    "legend.fontsize":  20,
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
    import ROOT
    ROOT.gROOT.SetBatch(True)
    f = ROOT.TFile(filepath)
    t = f.Get("limit")
    xs, ys = [], []
    for ev in t:
        x = getattr(ev, "k_" + poi)
        y = 2 * ev.deltaNLL
        if y <= maxNLL:
            xs.append(x)
            ys.append(y)
    f.Close()

    xs, ys = zip(*sorted(zip(xs, ys)))
    # shift so minimum is 0
    y_min = min(ys)
    ys = [y - y_min for y in ys]

    x_b = xs[ys.index(min(ys))]
    x1 = getAllIntervals_xy(xs, ys, 1.0)
    x2 = getAllIntervals_xy(xs, ys, 4.0)

    return {
        "best": x_b,
        "1sigma": x1,
        "2sigma": x2
    }


def getAllIntervals_xy(xs, ys, val):
    """
    Return every disjoint (low, high) x-interval where the (linearly
    interpolated) scan sits at or below `val`, instead of collapsing to a
    single bracketing pair. A run that touches the first/last scan point
    stays below threshold there and is left open at the scan boundary
    (low=min(xs) / high=max(xs)), matching the old fallback behaviour.
    """
    n = len(xs)
    below = [y <= val for y in ys]

    intervals = []
    i = 0
    while i < n:
        if not below[i]:
            i += 1
            continue

        if i == 0:
            low = xs[0]
        else:
            low = xs[i-1] + (val - ys[i-1]) * (xs[i] - xs[i-1]) / (ys[i] - ys[i-1])

        j = i
        while j < n and below[j]:
            j += 1

        if j == n:
            high = xs[n-1]
        else:
            high = xs[j-1] + (val - ys[j-1]) * (xs[j] - xs[j-1]) / (ys[j] - ys[j-1])

        intervals.append((low, high))
        i = j

    if not intervals:
        intervals = [(min(xs), max(xs))]

    return intervals


def outerBounds(intervals):
    """Envelope (min low, max high) across a list of (low, high) intervals."""
    return min(iv[0] for iv in intervals), max(iv[1] for iv in intervals)


# -------------------------
# Scan directories
# -------------------------

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--indir", default=".")
parser.add_argument("--ops", nargs="+", default=None)
parser.add_argument("--horizontal", action="store_true",
                    help="Slide-friendly layout: operators on x-axis, panels stacked vertically")
parser.add_argument("--logscale",   action="store_true",    help="Symlog scale on the Wilson coefficient axis (handles negatives)")
parser.add_argument("--linthresh",  type=float, default=1e-2, help="Linear zone half-width for --logscale (default: 1e-2)")
parser.add_argument("--verbose",    action="store_true",    help="Print extracted best-fit and CL intervals for each operator")
args = parser.parse_args()

if args.logscale:
    linthresh = args.linthresh
    print(f"symlog linthresh = {linthresh:.2e}")

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

    res_mc   = extract_intervals(file_mc,   poi)
    res_nomc = extract_intervals(file_nomc, poi)
    results[op] = {"MC": res_mc, "noMC": res_nomc}
    if args.verbose:
        def fmt(intervals):
            return " U ".join(f"[{lo:+.4f}, {hi:+.4f}]" for lo, hi in intervals)
        print(f"{op}:")
        print(f"  MC   best={res_mc['best']:+.4f}  1s={fmt(res_mc['1sigma'])}  2s={fmt(res_mc['2sigma'])}")
        print(f"  noMC best={res_nomc['best']:+.4f}  1s={fmt(res_nomc['1sigma'])}  2s={fmt(res_nomc['2sigma'])}")


# sort by MC stat+syst 95% CL: smallest farthest bound from zero (across all
# disjoint intervals) = best sensitivity first
results = dict(sorted(
    results.items(),
    key=lambda kv: max(abs(v) for lo, hi in kv[1]["MC"]["2sigma"] for v in (lo, hi))
))

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
        intervals_1s = r["1sigma"]
        intervals_2s = r["2sigma"]

        if args.horizontal:
            p = i + shift   # position along x (operator axis)
            for lo, hi in intervals_2s:
                ax.vlines(p, lo, hi, colors=color, linestyles='dashed', linewidth=1.5)
            for lo, hi in intervals_1s:
                ax.vlines(p, lo, hi, colors=color, linestyles='solid',  linewidth=3)
            ax.plot(p, x_b, 'o', color=color)
        else:
            p = i + shift   # position along y
            for lo, hi in intervals_2s:
                ax.hlines(p, lo, hi, colors=color, linestyles='dashed', linewidth=1.5)
            for lo, hi in intervals_1s:
                ax.hlines(p, lo, hi, colors=color, linestyles='solid',  linewidth=3)
            ax.plot(x_b, p, 'o', color=color)

    # =========================
    # SECOND PANEL (Lambda bars)
    # =========================

    for shift, key, color in [(+0.15, "MC", "blue"), (-0.15, "noMC", "red")]:
        r = res[key]
        p = i + shift

        # envelope across all disjoint 2sigma intervals: the Lambda reach is
        # set by the furthest-out allowed |c|, regardless of any gap between
        # disjoint regions
        lo_env, hi_env = outerBounds(r["2sigma"])
        a = abs(lo_env) + abs(hi_env)
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

op_labels = list(results.keys())

if args.horizontal:
    ax.set_xticks(pos)
    ax.tick_params(axis='x', labelbottom=False)
    ax.axhline(0, color='black', linestyle='--', linewidth=1)
    ax.set_xlim(-1.0, n_ops + 0.5)
    ax.set_ylabel("Wilson coefficient")
    if args.logscale:
        ax.set_yscale("symlog", linthresh=linthresh)
else:
    ax.set_yticks(pos)
    ax.set_yticklabels(op_labels)
    ax.axvline(0, color='black', linestyle='--', linewidth=1)
    ax.set_ylim(-1.0, n_ops + 2)
    ax.set_xlabel("Wilson coefficient")
    if args.logscale:
        ax.set_xscale("symlog", linthresh=linthresh)

# =========================
# Formatting SECOND panel
# =========================

if args.horizontal:
    ax2.set_xticks(pos)
    ax2.set_xticklabels(op_labels, rotation=45, ha='right', rotation_mode='anchor')
    ax2.set_ylabel(r"$\Lambda$ at 95% CL [TeV]")
    ax2.set_yscale("log")
else:
    ax2.set_yticks(pos)
    ax2.set_yticklabels(op_labels)
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
plt.savefig(f"eft_summary_two_panel_multiInterval{suffix}.pdf", bbox_inches='tight')
plt.savefig(f"eft_summary_two_panel_multiInterval{suffix}.png", dpi=150, bbox_inches='tight')

plt.show()