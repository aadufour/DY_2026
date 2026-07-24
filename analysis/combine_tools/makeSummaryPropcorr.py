#!/usr/bin/env python3
"""
makeSummaryPropcorr.py

Same as makeSummaryMultiVar.py, but instead of overlaying all variables on
one panel, it produces one summary plot per variable (mll, rapll_abs,
costhetastar, triple_diff), each comparing two fits -- baseline and
propagator-corrected (propcorr) -- for every operator. Only stat+syst scans
are used.

Both configs are expected to have the same layout:
    <root>/<datacards-subpath>/<mll|rapll_abs|costhetastar|triple_diff>/
        higgsCombine.<op>.individual.MultiDimFit.mH125.root

Usage:
    makeSummaryPropcorr.py --bl-dir /path/to/eft_bkg_fullsyst_v9 \\
                            --pc-dir /path/to/propcorr_v1 \\
                            [--sort-by triple_diff] [--horizontal] \\
                            [--logscale] [--verbose]
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

BL_COLOR = "#2166ac"
PC_COLOR = "#d6604d"

VARS = {
    "mll":          {"dirname": "mll",          "label": r"$m_{\ell\ell}$"},
    "rapll":        {"dirname": "rapll_abs",     "label": r"$|y_{\ell\ell}|$"},
    "costhetastar": {"dirname": "costhetastar",  "label": r"$\cos\theta^*$"},
    "triple_diff":  {"dirname": "triple_diff",   "label": r"Triple-diff"},
}
DISPLAY_ORDER = ["triple_diff", "mll", "costhetastar", "rapll"]

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
parser.add_argument("--bl-dir", dest="bl_dir", required=True,
                    help="Root directory of the baseline (non-propcorr) config")
parser.add_argument("--pc-dir", dest="pc_dir", required=True,
                    help="Root directory of the propagator-corrected config")
parser.add_argument("--datacards-subpath", default="datacards/inc_mm",
                    help="Subpath from each root dir to the per-variable scan "
                         "folders (default: datacards/inc_mm)")
parser.add_argument("--label-bl", default="baseline", help="Legend label for --bl-dir")
parser.add_argument("--label-pc", default="propcorr", help="Legend label for --pc-dir")
parser.add_argument("--ops",        nargs="+", default=None)
parser.add_argument("--sort-by",    default="mll",
                    help="Variable used to order operators consistently across "
                         "all plots, ranked by baseline 95%% CL sensitivity "
                         "(default: mll)")
parser.add_argument("--horizontal", action="store_true")
parser.add_argument("--logscale",   action="store_true")
parser.add_argument("--linthresh",  type=float, default=1e-2)
parser.add_argument("--verbose",    action="store_true")
parser.add_argument("-o", "--outname", default="eft_summary_propcorr",
                    help="Output file base name prefix (variable name is appended)")
args = parser.parse_args()

if args.logscale:
    linthresh = args.linthresh
    print(f"symlog linthresh = {linthresh:.2e}")

# -------------------------
# Resolve per-variable directories, keep only variables present in both configs
# -------------------------

var_dirs = {}
for var in DISPLAY_ORDER:
    bl_path = os.path.join(args.bl_dir, args.datacards_subpath, VARS[var]["dirname"])
    pc_path = os.path.join(args.pc_dir, args.datacards_subpath, VARS[var]["dirname"])
    if not os.path.isdir(bl_path):
        print(f"  [skip var] {var}: baseline dir not found: {bl_path}")
        continue
    if not os.path.isdir(pc_path):
        print(f"  [skip var] {var}: propcorr dir not found: {pc_path}")
        continue
    var_dirs[var] = {"bl": bl_path, "pc": pc_path}

active_vars = [v for v in DISPLAY_ORDER if v in var_dirs]
if not active_vars:
    raise SystemExit("No variable found in both --bl-dir and --pc-dir")

if args.sort_by not in active_vars:
    raise ValueError(f"--sort-by '{args.sort_by}' not among available variables {active_vars}")

print(f"Variables: {active_vars}")

# -------------------------
# Load results per variable/dataset
# -------------------------

if args.ops:
    base_operators = sorted(args.ops)
else:
    base_operators = discover_operators(var_dirs[args.sort_by]["bl"])

all_results = {}
for var in active_vars:
    if args.verbose:
        print(f"\n--- {var} ({args.label_bl}) ---")
    bl_res = load_results(var_dirs[var]["bl"], base_operators, verbose=args.verbose)
    if args.verbose:
        print(f"\n--- {var} ({args.label_pc}) ---")
    pc_res = load_results(var_dirs[var]["pc"], base_operators, verbose=args.verbose)
    all_results[var] = {
        op: {"bl": bl_res[op], "pc": pc_res[op]}
        for op in base_operators if op in bl_res and op in pc_res
    }

# global operator order, ranked by baseline 95% CL sensitivity of --sort-by
sort_results = all_results[args.sort_by]
operators = [op for op in base_operators if op in sort_results]
operators.sort(key=lambda op: max(
    abs(sort_results[op]["bl"]["2sigma"][0]),
    abs(sort_results[op]["bl"]["2sigma"][1])
))
print(f"Ordering by: {args.sort_by} ({args.label_bl})")

n_ops = len(operators)
pos   = np.arange(n_ops)

# -------------------------
# Plot (one figure per variable)
# -------------------------

from matplotlib.lines import Line2D
from matplotlib.patches import Patch


def build_plot(var):
    results = all_results[var]
    ops     = [op for op in operators if op in results]
    n       = len(ops)
    p_all   = np.arange(n)

    if args.horizontal:
        fig_width = max(10, WIDTH_PER_OP * n)
        fig, (ax, ax2) = plt.subplots(
            nrows=2,
            figsize=(fig_width, FIG_HEIGHT),
            gridspec_kw={"height_ratios": HEIGHT_RATIOS},
            sharex=True,
        )
    else:
        fig, (ax, ax2) = plt.subplots(
            ncols=2,
            figsize=(12, max(6, 0.6 * n)),
            gridspec_kw={"width_ratios": [2.5, 1]},
            sharey=True,
        )

    for i, op in enumerate(ops):
        for shift, key, color in [(+0.18, "bl", BL_COLOR), (-0.18, "pc", PC_COLOR)]:
            r   = results[op][key]
            p   = i + shift
            x_b = r["best"]
            x1  = r["1sigma"]
            x2s = r["2sigma"]

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

    if args.horizontal:
        ax.set_xticks(p_all)
        ax.tick_params(axis="x", labelbottom=False)
        ax.axhline(0, color="black", linestyle="--", linewidth=1)
        ax.set_xlim(-1.0, n + 0.5)
        ax.set_ylabel("Wilson coefficient")
        if args.logscale:
            ax.set_yscale("symlog", linthresh=linthresh)

        ax2.set_xticks(p_all)
        ax2.set_xticklabels(ops, rotation=45, ha="right", rotation_mode="anchor")
        ax2.set_ylabel(r"$\Lambda$ at 95% CL [TeV]")
        ax2.set_yscale("log")
    else:
        ax.set_yticks(p_all)
        ax.set_yticklabels(ops)
        ax.axvline(0, color="black", linestyle="--", linewidth=1)
        ax.set_ylim(-1.0, n + 0.5)
        ax.set_xlabel("Wilson coefficient")
        if args.logscale:
            ax.set_xscale("symlog", linthresh=linthresh)

        ax2.set_yticks(p_all)
        ax2.set_yticklabels(ops)
        ax2.tick_params(axis="y", left=False, labelleft=False)
        ax2.set_xlabel(r"$\Lambda$ at 95% CL" + "\n[TeV]")
        ax2.set_xscale("log")

    # -------------------------
    # Legends
    # -------------------------

    interval_handles = [
        Line2D([], [], color=BL_COLOR, lw=3, label=args.label_bl),
        Line2D([], [], color=PC_COLOR, lw=3, label=args.label_pc),
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
    fig.suptitle(VARS[var]["label"], y=top_margin + 0.09)

    suffix = "_horizontal" if args.horizontal else ""
    outname = f"{args.outname}_{var}{suffix}"
    plt.savefig(f"{outname}.pdf", bbox_inches="tight")
    plt.savefig(f"{outname}.png", dpi=150, bbox_inches="tight")
    print(f"Saved: {outname}.pdf / .png")


for var in active_vars:
    build_plot(var)

plt.show()
