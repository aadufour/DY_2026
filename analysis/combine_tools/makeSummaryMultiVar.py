#!/usr/bin/env python3
"""
makeSummaryMultiVar.py

Summary plot overlaying up to four observables (mll, rapll_abs, costhetastar,
triple_diff) on the same panel.  Operators are ordered by mll 95% CL sensitivity
(tightest first).  Only stat+syst scans are used.

Usage:
    makeSummaryMultiVar.py --mll mll/scan [--rapll rapll_abs/scan]
                           [--costhetastar costhetastar/scan]
                           [--triple_diff triple_diff/scan]
                           [--horizontal] [--logscale] [--verbose]
"""

import os
import glob
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import ROOT
ROOT.gROOT.SetBatch(True)

# ============================================================
# PLOT CONFIGURATION
# ============================================================

FONT_SIZE       = 34
LABEL_SIZE      = 24
TICK_LABELSIZE  = 25
LEGEND_FONTSIZE = 19

FIG_HEIGHT      = 10
WIDTH_PER_OP    = 0.9

VARS = {
    "mll":          {"color": "#2166ac", "label": r"$m_{\ell\ell}$"},
    "rapll":        {"color": "#d6604d", "label": r"$|y_{\ell\ell}|$"},
    "costhetastar": {"color": "#4dac26", "label": r"$\cos\theta^*$"},
    "triple_diff":  {"color": "#8B008B", "label": r"Triple-diff"},
}

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
# Interval extraction (reads ROOT tree directly)
# -------------------------

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
        return [min(xs), max(xs)]
    best = xs[ys.index(min(ys))]
    lo = [x for x in xings if x <= best]
    hi = [x for x in xings if x >= best]
    if lo and hi:
        return [max(lo), min(hi)]
    return xings[:2]


def extract_intervals(filepath, poi, maxNLL=10):
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
    y_min = min(ys)
    ys = [y - y_min for y in ys]
    x_b = xs[ys.index(min(ys))]
    x1  = getLSintersections_xy(xs, ys, 1.0)
    x2  = getLSintersections_xy(xs, ys, 4.0)
    return {"best": x_b, "1sigma": x1, "2sigma": x2}


def discover_operators(scan_dir):
    pattern = os.path.join(scan_dir, "higgsCombine.*.individual.MultiDimFit.mH125.root")
    ops = []
    for f in glob.glob(pattern):
        m = re.search(r"higgsCombine\.(.+)\.individual", f)
        if m and "_stat" not in f:
            ops.append(m.group(1))
    return sorted(ops)


def load_results(scan_dir, operators, verbose=False):
    results = {}
    for op in operators:
        fp = os.path.join(scan_dir, f"higgsCombine.{op}.individual.MultiDimFit.mH125.root")
        if not os.path.exists(fp):
            print(f"  [skip] {op}: not found in {scan_dir}")
            continue
        try:
            r = extract_intervals(fp, op)
            results[op] = r
            if verbose:
                print(f"  {op}: best={r['best']:+.4f}  "
                      f"1s=[{r['1sigma'][0]:+.4f}, {r['1sigma'][1]:+.4f}]  "
                      f"2s=[{r['2sigma'][0]:+.4f}, {r['2sigma'][1]:+.4f}]")
        except Exception as e:
            print(f"  [skip] {op}: {e}")
    return results


# -------------------------
# Args
# -------------------------

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--mll",           required=True, help="Directory with mll scan ROOT files")
parser.add_argument("--rapll",         default=None,  help="Directory with rapll_abs scan ROOT files")
parser.add_argument("--costhetastar",  default=None,  help="Directory with costhetastar scan ROOT files")
parser.add_argument("--triple_diff",   default=None,  help="Directory with triple_diff scan ROOT files")
parser.add_argument("--ops",        nargs="+", default=None)
parser.add_argument("--sort-by",    default=None,
                    help="Variable to use for operator ordering (default: mll). "
                         "Must be one of the active variables, e.g. --sort-by triple_diff")
parser.add_argument("--horizontal", action="store_true")
parser.add_argument("--logscale",   action="store_true")
parser.add_argument("--linthresh",  type=float, default=1e-2)
parser.add_argument("--verbose",    action="store_true")
args = parser.parse_args()

# -------------------------
# Load
# -------------------------

dirs = {"mll": args.mll}
if args.rapll:        dirs["rapll"]        = args.rapll
if args.costhetastar: dirs["costhetastar"] = args.costhetastar
if args.triple_diff:  dirs["triple_diff"]  = args.triple_diff

DISPLAY_ORDER = ["triple_diff", "mll", "costhetastar", "rapll"]
active_vars = [v for v in DISPLAY_ORDER if v in dirs]

n_vars  = len(active_vars)
offsets = np.linspace(-0.3, 0.3, n_vars) if n_vars > 1 else [0.0]
SHIFTS  = {v: offsets[i] for i, v in enumerate(active_vars)}

if args.ops:
    operators = sorted(args.ops)
else:
    operators = discover_operators(args.mll)

print(f"Operators: {len(operators)}   Variables: {active_vars}")

all_results = {}
for var, d in dirs.items():
    if args.verbose:
        print(f"\n--- {var} ---")
    all_results[var] = load_results(d, operators, verbose=args.verbose)

sort_var = args.sort_by if args.sort_by else "mll"
if sort_var not in all_results:
    raise ValueError(f"--sort-by '{sort_var}' not in active variables {list(all_results.keys())}")

operators = [op for op in operators if op in all_results[sort_var]]
operators.sort(key=lambda op: max(
    abs(all_results[sort_var][op]["2sigma"][0]),
    abs(all_results[sort_var][op]["2sigma"][1])
))
print(f"Ordering by: {sort_var}")

if args.logscale:
    linthresh = args.linthresh
    print(f"symlog linthresh = {linthresh:.2e}")

n_ops = len(operators)
pos   = np.arange(n_ops)

# -------------------------
# Plot
# -------------------------

if args.horizontal:
    fig_width = max(10, WIDTH_PER_OP * n_ops)
    fig, (ax, ax2) = plt.subplots(
        nrows=2,
        figsize=(fig_width, FIG_HEIGHT),
        gridspec_kw={"height_ratios": HEIGHT_RATIOS},
        sharex=True,
    )
else:
    fig, (ax, ax2) = plt.subplots(
        ncols=2,
        figsize=(12, max(6, 0.6 * n_ops)),
        gridspec_kw={"width_ratios": [2.5, 1]},
        sharey=True,
    )

for i, op in enumerate(operators):
    for var in active_vars:
        if op not in all_results[var]:
            continue

        r     = all_results[var][op]
        color = VARS[var]["color"]
        p     = i + SHIFTS[var]
        x_b   = r["best"]
        x1    = r["1sigma"]
        x2s   = r["2sigma"]

        if args.horizontal:
            ax.vlines(p, x2s[0], x2s[1], colors=color, linestyles="dashed", linewidth=1.5)
            ax.vlines(p, x1[0],  x1[1],  colors=color, linestyles="solid",  linewidth=3)
            ax.plot(p, x_b, "o", color=color, markersize=5)
        else:
            ax.hlines(p, x2s[0], x2s[1], colors=color, linestyles="dashed", linewidth=1.5)
            ax.hlines(p, x1[0],  x1[1],  colors=color, linestyles="solid",  linewidth=3)
            ax.plot(x_b, p, "o", color=color, markersize=5)

        a = abs(x2s[0]) + abs(x2s[1])
        if a <= 0:
            continue
        lam1   = np.sqrt(1.0 / a)
        lam4pi = np.sqrt((4 * np.pi)**2 / a)

        if args.horizontal:
            ax2.bar(p, lam1,          width=0.18, color=color, alpha=0.9)
            ax2.bar(p, lam4pi - lam1, width=0.18, bottom=lam1, color=color, alpha=0.3)
        else:
            ax2.barh(p, lam1,          height=0.18, color=color, alpha=0.9)
            ax2.barh(p, lam4pi - lam1, height=0.18, left=lam1,  color=color, alpha=0.3)

# -------------------------
# Formatting
# -------------------------

from matplotlib.lines import Line2D
from matplotlib.patches import Patch

if args.horizontal:
    ax.set_xticks(pos)
    ax.tick_params(axis="x", labelbottom=False)
    ax.axhline(0, color="black", linestyle="--", linewidth=1)
    ax.set_xlim(-1.0, n_ops + 0.5)
    ax.set_ylabel("Wilson coefficient")
    if args.logscale:
        ax.set_yscale("symlog", linthresh=linthresh)

    ax2.set_xticks(pos)
    ax2.set_xticklabels(operators, rotation=45, ha="right", rotation_mode="anchor")
    ax2.set_ylabel(r"$\Lambda$ at 95% CL [TeV]")
    ax2.set_yscale("log")
else:
    ax.set_yticks(pos)
    ax.set_yticklabels(operators)
    ax.axvline(0, color="black", linestyle="--", linewidth=1)
    ax.set_ylim(-1.0, n_ops + 0.5)
    ax.set_xlabel("Wilson coefficient", loc="left")
    if args.logscale:
        ax.set_xscale("symlog", linthresh=linthresh)

    ax2.set_yticks(pos)
    ax2.set_yticklabels(operators)
    ax2.tick_params(axis="y", left=False, labelleft=False)
    ax2.set_xlabel(r"$\Lambda$ at 95% CL [TeV]")
    ax2.set_xscale("log")

# -------------------------
# Legends
# -------------------------

from matplotlib.lines import Line2D
from matplotlib.patches import Patch

interval_handles = [
    Line2D([], [], color=VARS[v]["color"], lw=3, label=VARS[v]["label"])
    for v in active_vars
] + [
    Line2D([], [], color="grey", lw=3,            label="68% CL"),
    Line2D([], [], color="grey", lw=1.5, ls="--", label="95% CL"),
]

lambda_handles = [
    Patch(facecolor="grey", alpha=0.9, label=r"$c=1$"),
    Patch(facecolor="grey", alpha=0.3, label=r"$c=(4\pi)^2$"),
]
ax2.legend(handles=lambda_handles, ncol=1, frameon=False,
           loc="upper right" if args.horizontal else "upper center")

hep.cms.label(ax=ax, data=True, label="Preliminary")

# Reserve whitespace above the panels for the variable/CL legend so it sits
# in its own row above the CMS label instead of competing with the data.
plt.tight_layout()
top_margin = 0.84 if args.horizontal else 0.90
fig.subplots_adjust(top=top_margin)
legend_y = top_margin + 0.02
fig.legend(
    handles=interval_handles,
    loc="lower center",
    bbox_to_anchor=(0.5, legend_y),
    ncol=len(interval_handles),
    frameon=False,
    fontsize=LEGEND_FONTSIZE,
    columnspacing=1.4,
    handlelength=1.6,
    handletextpad=0.5,
)

suffix = "_horizontal" if args.horizontal else ""
plt.savefig(f"eft_summary_multivar{suffix}.pdf", bbox_inches="tight")
plt.savefig(f"eft_summary_multivar{suffix}.png", dpi=150, bbox_inches="tight")
print(f"Saved: eft_summary_multivar{suffix}.pdf / .png")
plt.show()
