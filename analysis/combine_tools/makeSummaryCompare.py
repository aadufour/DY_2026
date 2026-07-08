#!/usr/bin/env python3
"""
makeSummaryCompare.py

Same as makeSummary.py but compares two different scan folders instead of
MC vs noMC (stat-only) within a single folder. Stat-only files are ignored
in both folders; only the full stat+syst scan is used.

Usage:
    python3 makeSummaryCompare.py --indir1 scan_noprop --indir2 scan_prop \\
        --label1 "no prop" --label2 "prop"
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep

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
# Extract intervals
# -------------------------

def extract_intervals(filepath, poi, maxNLL=10):
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
    x1 = getLSintersections_xy(xs, ys, 1.0)
    x2 = getLSintersections_xy(xs, ys, 4.0)

    return {
        "best": x_b,
        "1sigma": x1,
        "2sigma": x2
    }


def getLSintersections_xy(xs, ys, val):
    xings = []
    for i in range(1, len(xs)):
        if ys[i] == val:
            xings.append(xs[i])
        elif (ys[i] - val) * (ys[i-1] - val) < 0:
            xings.append(
                xs[i-1] + (val - ys[i-1]) * (xs[i] - xs[i-1]) / (ys[i] - ys[i-1])
            )
    if len(xings) < 2:
        xings = [min(xs), max(xs)]
    return xings[:2]


# -------------------------
# Args
# -------------------------

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--indir1", required=True, help="First scan folder")
parser.add_argument("--indir2", required=True, help="Second scan folder")
parser.add_argument("--label1", default="no prop", help="Label for --indir1 (default: 'no prop')")
parser.add_argument("--label2", default="prop", help="Label for --indir2 (default: 'prop')")
parser.add_argument("--ops", nargs="+", default=None)
parser.add_argument("--horizontal", action="store_true",
                    help="Slide-friendly layout: operators on x-axis, panels stacked vertically")
parser.add_argument("--logscale",   action="store_true",    help="Symlog scale on the Wilson coefficient axis (handles negatives)")
parser.add_argument("--linthresh",  type=float, default=1e-2, help="Linear zone half-width for --logscale (default: 1e-2)")
parser.add_argument("--verbose",    action="store_true",    help="Print extracted best-fit and CL intervals for each operator")
parser.add_argument("-o", "--outname", default="eft_summary_compare", help="Output file base name")
args = parser.parse_args()

if args.logscale:
    linthresh = args.linthresh
    print(f"symlog linthresh = {linthresh:.2e}")

import glob, re

def discover_ops(base_dir):
    pattern = os.path.join(base_dir, "higgsCombine.*.individual.MultiDimFit.mH125.root")
    return set(
        re.search(r"higgsCombine\.(.+)\.individual", f).group(1)
        for f in glob.glob(pattern)
        if "_stat" not in f
    )

if args.ops:
    operators = sorted(args.ops)
else:
    operators = sorted(discover_ops(args.indir1) & discover_ops(args.indir2))

results = {}

for op in operators:

    poi = op

    file1 = f"{args.indir1}/higgsCombine.{op}.individual.MultiDimFit.mH125.root"
    file2 = f"{args.indir2}/higgsCombine.{op}.individual.MultiDimFit.mH125.root"

    if not (os.path.exists(file1) and os.path.exists(file2)):
        print(f"Skipping {op} (missing files)")
        continue

    res1 = extract_intervals(file1, poi)
    res2 = extract_intervals(file2, poi)
    results[op] = {"dir1": res1, "dir2": res2}
    if args.verbose:
        print(f"{op}:")
        print(f"  {args.label1:10s} best={res1['best']:+.4f}  1s=[{res1['1sigma'][0]:+.4f}, {res1['1sigma'][1]:+.4f}]  2s=[{res1['2sigma'][0]:+.4f}, {res1['2sigma'][1]:+.4f}]")
        print(f"  {args.label2:10s} best={res2['best']:+.4f}  1s=[{res2['1sigma'][0]:+.4f}, {res2['1sigma'][1]:+.4f}]  2s=[{res2['2sigma'][0]:+.4f}, {res2['2sigma'][1]:+.4f}]")


# sort by dir1 95% CL: smallest farthest bound from zero = best sensitivity first
results = dict(sorted(
    results.items(),
    key=lambda kv: max(abs(kv[1]["dir1"]["2sigma"][0]), abs(kv[1]["dir1"]["2sigma"][1]))
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

    for shift, key, color in [(+0.15, "dir1", "blue"), (-0.15, "dir2", "red")]:
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

    for shift, key, color in [(+0.15, "dir1", "blue"), (-0.15, "dir2", "red")]:
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
    f"{args.label1} 1σ",
    f"{args.label1} 2σ",
    f"{args.label2} 1σ",
    f"{args.label2} 2σ",
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
plt.savefig(f"{args.outname}{suffix}.pdf", bbox_inches='tight')
plt.savefig(f"{args.outname}{suffix}.png", dpi=150, bbox_inches='tight')

plt.show()
