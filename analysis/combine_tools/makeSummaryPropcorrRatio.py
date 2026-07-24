#!/usr/bin/env python3
"""
makeSummaryPropcorrRatio.py

Companion to makeSummaryPropcorr.py: for each variable (mll, rapll_abs,
costhetastar, triple_diff) and each operator, computes the ratio of the 68%
CL interval width between the propagator-corrected (propcorr) and baseline
fits. Only stat+syst scans are used.

    ratio = width_pc / width_bl

ratio == 1: no change. ratio < 1: propcorr narrows the interval.
ratio > 1: propcorr widens it.

Both configs are expected to have the same layout:
    <root>/<datacards-subpath>/<mll|rapll_abs|costhetastar|triple_diff>/
        higgsCombine.<op>.individual.MultiDimFit.mH125.root

Usage:
    makeSummaryPropcorrRatio.py --bl-dir /path/to/eft_bkg_fullsyst_v9 \\
                                 --pc-dir /path/to/propcorr_v1 \\
                                 [--horizontal] [--verbose]
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

BAR_COLOR = "#2166ac"

VARS = {
    "mll":          {"dirname": "mll",          "label": r"$m_{\ell\ell}$"},
    "rapll":        {"dirname": "rapll_abs",     "label": r"$|y_{\ell\ell}|$"},
    "costhetastar": {"dirname": "costhetastar",  "label": r"$\cos\theta^*$"},
    "triple_diff":  {"dirname": "triple_diff",   "label": r"Triple-diff"},
}
DISPLAY_ORDER = ["triple_diff", "mll", "costhetastar", "rapll"]

# Operators that actually enter the propagator correction for this process,
# per notes/propagator_correction.md section 9/10 (code-derived via
# auto_detect_operators_propcorr.py, cross-checked against the SMEFTsim
# practical guide). Everything else is unaffected by construction and would
# just show ~0% change, so it's excluded from the default operator set.
PROPCORR_OPS = [
    "cHDD", "cHWB",
    "cHj1", "cHj3", "cHQ1", "cHQ3", "cHu", "cHd", "cHbq", "cbWRe", "cbBRe",
    "cHl1", "cHl3", "cHe",
    "cll1",
]

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
parser.add_argument("--ops",        nargs="+", default=None,
                    help="Explicit operator list, overrides the default "
                         "propcorr-relevant filter")
parser.add_argument("--all-ops",    action="store_true",
                    help="Disable the default propcorr-relevant filter and "
                         "use every operator discovered in the scan dir")
parser.add_argument("--sort-by",    default="mll",
                    help="Variable used only to discover the master operator "
                         "list when --ops is not given (default: mll).")
parser.add_argument("--horizontal", action="store_true")
parser.add_argument("--verbose",    action="store_true")
parser.add_argument("-o", "--outname", default="eft_widthratio_propcorr",
                    help="Output file base name prefix (variable name is appended)")
args = parser.parse_args()

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
# Load results and compute 68% CL width ratios per variable
# -------------------------

if args.ops:
    base_operators = sorted(args.ops)
else:
    discovered = discover_operators(var_dirs[args.sort_by]["bl"])
    if args.all_ops:
        base_operators = discovered
    else:
        base_operators = [op for op in PROPCORR_OPS if op in discovered]
        missing = [op for op in PROPCORR_OPS if op not in discovered]
        if missing:
            print(f"  [note] propcorr-relevant ops not found in discovery dir: {missing}")
print(f"Operators: {base_operators}")

all_ratios = {}
for var in active_vars:
    if args.verbose:
        print(f"\n--- {var} (baseline) ---")
    bl_res = load_results(var_dirs[var]["bl"], base_operators, verbose=args.verbose)
    if args.verbose:
        print(f"\n--- {var} (propcorr) ---")
    pc_res = load_results(var_dirs[var]["pc"], base_operators, verbose=args.verbose)

    ratios = {}
    for op in base_operators:
        if op not in bl_res or op not in pc_res:
            continue
        w_bl = bl_res[op]["1sigma"][1] - bl_res[op]["1sigma"][0]
        w_pc = pc_res[op]["1sigma"][1] - pc_res[op]["1sigma"][0]
        if w_bl <= 0:
            continue
        ratio = w_pc / w_bl
        ratios[op] = ratio
        if args.verbose:
            print(f"  {op}: width_bl={w_bl:.4f}  width_pc={w_pc:.4f}  ratio={ratio:.3f}")
    all_ratios[var] = ratios

# -------------------------
# Plot (one figure per variable)
# -------------------------

def build_plot(var):
    ratios = all_ratios[var]
    ops = sorted(ratios.keys(), key=lambda op: abs(ratios[op] - 1.0), reverse=True)
    n   = len(ops)
    pos = np.arange(n)
    vals = [ratios[op] for op in ops]

    if args.horizontal:
        fig_width = max(10, WIDTH_PER_OP * n)
        fig, ax = plt.subplots(figsize=(fig_width, FIG_HEIGHT))
        ax.bar(pos, vals, width=0.6, color=BAR_COLOR)
        ax.axhline(1, color="black", linestyle="--", linewidth=1)
        ax.set_xticks(pos)
        ax.set_xticklabels(ops, rotation=45, ha="right", rotation_mode="anchor")
        ax.set_ylabel("68% CL width ratio\n(propcorr / baseline)")
    else:
        fig, ax = plt.subplots(figsize=(12, max(6, 0.6 * n)))
        ax.barh(pos, vals, height=0.6, color=BAR_COLOR)
        ax.axvline(1, color="black", linestyle="--", linewidth=1)
        ax.set_yticks(pos)
        ax.set_yticklabels(ops)
        ax.set_xlabel("68% CL width ratio (propcorr / baseline)")

    hep.cms.label(ax=ax, data=True, label="Preliminary")

    plt.tight_layout()
    top_margin = 0.84 if args.horizontal else 0.90
    fig.subplots_adjust(top=top_margin)
    fig.suptitle(VARS[var]["label"], x=0.25, ha="left", y=top_margin + 0.05)

    suffix = "_horizontal" if args.horizontal else ""
    outname = f"{args.outname}_{var}{suffix}"
    plt.savefig(f"{outname}.pdf", bbox_inches="tight")
    plt.savefig(f"{outname}.png", dpi=150, bbox_inches="tight")
    print(f"Saved: {outname}.pdf / .png")


for var in active_vars:
    build_plot(var)

plt.show()
