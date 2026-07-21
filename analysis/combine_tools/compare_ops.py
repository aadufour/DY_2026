#!/usr/bin/env python3
"""
compare_ops.py
==============
General-purpose propcorr-vs-baseline comparison for LHE reweight caches.
Weight decomposition matches notes/lhe.md section 3.1 exactly:
    w_lin   = 0.5*(w_p1 - w_m1)
    w_quad  = 0.5*(w_p1 + w_m1) - w_SM
    w_slq   = w_SM + C*w_lin + C^2*w_quad          ("complete", single op)
    w_inter = w_pp[op1,op2] - w_p1[op1] - w_p1[op2] + w_SM   ("mixed", two ops)

Usage:
    python3 compare_ops.py --op1 cHDD --component lin
    python3 compare_ops.py --op1 cHDD --component quad
    python3 compare_ops.py --op1 cHDD --component complete
    python3 compare_ops.py --op1 cHDD --op2 cHWB --component lin      # both, side by side
    python3 compare_ops.py --op1 cHDD --op2 cHWB --component mixed    # cross term
    python3 compare_ops.py --op1 cHDD --lo 50 --hi 3000 --binwidth 10 --component complete
"""

import argparse
import os
import pickle
from copy import deepcopy

import numpy as np
import matplotlib
matplotlib.use('Agg')   # no display on the cluster -- write straight to file
import matplotlib.pyplot as plt
import mplhep as hep

PROP_CACHE = '/grid_mnt/data__data.polcms/cms/adufour/LHE/propcorr/CACHE/lhe_cache_propcorr_parallel.pkl'
BASE_CACHE = '/grid_mnt/data__data.polcms/cms/adufour/LHE/SYST_slc7/CACHE/lhe_cache_parallel.pkl'

FONT_SIZE       = 22
LABEL_SIZE      = 20
TICK_LABELSIZE  = 18
LEGEND_FONTSIZE = 16
FIG_DPI         = 200

style = deepcopy(hep.style.CMS)
style["font.size"]       = FONT_SIZE
style["axes.labelsize"]  = LABEL_SIZE
style["xtick.labelsize"] = TICK_LABELSIZE
style["ytick.labelsize"] = TICK_LABELSIZE
style["legend.fontsize"] = LEGEND_FONTSIZE
plt.style.use(style)

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--op1', required=True)
parser.add_argument('--op2', default=None, help='second operator; without --component mixed, just runs op1 and op2 independently side by side')
parser.add_argument('--component', choices=['lin', 'quad', 'complete', 'mixed'], default='complete')
parser.add_argument('--lo', type=float, default=80)
parser.add_argument('--hi', type=float, default=102)
parser.add_argument('--binwidth', type=float, default=0.5)
parser.add_argument('--C', type=float, default=1.0, help='Wilson coefficient value for --component complete')
parser.add_argument('--prop-cache', default=PROP_CACHE)
parser.add_argument('--base-cache', default=BASE_CACHE)
parser.add_argument('--outdir', default='./compare_ops_plots', help='where to save the .pdf/.png')
parser.add_argument('--no-normalize', action='store_true',
                     help='plot raw weighted bin sums instead of normalizing each histogram to sum=1')
parser.add_argument('--print-stats', action='store_true',
                     help='print the per-bin propcorr/baseline/pull table and chi2/dof summary (off by default)')
parser.add_argument('--lims', type=float, default=None,
                     help='fix the ratio-panel y-axis to (1-LIMS, 1+LIMS) for every plot, instead of '
                          'auto-scaling per operator -- use to get a common scale across a batch run')
args = parser.parse_args()
os.makedirs(args.outdir, exist_ok=True)

with open(args.prop_cache, 'rb') as f:
    prop = pickle.load(f)
with open(args.base_cache, 'rb') as f:
    base = pickle.load(f)

# N_gen normalization fix, matching build_datacard_new.py exactly:
# cache weights sum to sigma_sample * N_GEN_PER_SAMPLE, not sigma_sample.
# Dividing by N_GEN_PER_SAMPLE puts everything in genuine pb-scale cross
# section units -- required for --no-normalize to be physically meaningful
# and for propcorr/baseline to be directly comparable regardless of any
# difference in how many events survived cuts in each sample.
N_GEN_PER_SAMPLE = 100_000

def _apply_ngen_fix(cache):
    cache['w_SM'] = cache['w_SM'] / N_GEN_PER_SAMPLE
    cache['w_p1'] = {op: w / N_GEN_PER_SAMPLE for op, w in cache['w_p1'].items()}
    cache['w_m1'] = {op: w / N_GEN_PER_SAMPLE for op, w in cache['w_m1'].items()}
    cache['w_pp'] = {k: w / N_GEN_PER_SAMPLE for k, w in cache.get('w_pp', {}).items()}
    return cache

prop = _apply_ngen_fix(prop)
base = _apply_ngen_fix(base)

edges = np.arange(args.lo, args.hi + args.binwidth, args.binwidth)


def get_weight(cache, component, op1, op2, C):
    w_SM = cache['w_SM']
    if component == 'mixed':
        if op2 is None:
            raise ValueError('--op2 required for --component mixed')
        w_pp_all = cache.get('w_pp', {})
        pair = (op1, op2) if (op1, op2) in w_pp_all else (op2, op1)
        if pair not in w_pp_all:
            raise KeyError(f"No cross-weight for ({op1},{op2}) in cache "
                            f"-- was it built with --nodoubles?")
        return w_pp_all[pair] - cache['w_p1'][op1] - cache['w_p1'][op2] + w_SM

    wp1 = cache['w_p1'][op1]
    wm1 = cache['w_m1'][op1]
    w_lin  = 0.5 * (wp1 - wm1)
    w_quad = 0.5 * (wp1 + wm1) - w_SM
    if component == 'lin':
        return w_lin
    if component == 'quad':
        return w_quad
    if component == 'complete':
        return w_SM + C * w_lin + C**2 * w_quad
    raise ValueError(component)


def shape_and_sigma(cache, component, op1, op2, C):
    """Bin content and its exact statistical uncertainty, sigma = sqrt(sum(w_i^2))
    per bin (ROOT Sumw2 / boost_histogram Weight() convention). If normalizing,
    both content and sigma are divided by the same total -- an approximation
    that ignores the correlation between a bin and the total it's part of,
    but with --no-normalize (the recommended mode) this is the exact result."""
    w = get_weight(cache, component, op1, op2, C)
    mll = cache['mll']
    h,  _ = np.histogram(mll, bins=edges, weights=w)
    h2, _ = np.histogram(mll, bins=edges, weights=w**2)
    sigma = np.sqrt(h2)
    if args.no_normalize:
        return h, sigma
    total = h.sum()
    if total == 0:
        return h, sigma
    return h / total, sigma / total


def run_comparison(component, op1, op2):
    sp, sig_p = shape_and_sigma(prop, component, op1, op2, args.C)
    sb, sig_b = shape_and_sigma(base, component, op1, op2, args.C)

    sigma = np.sqrt(sig_p**2 + sig_b**2)
    sigma = np.where(sigma == 0, np.inf, sigma)
    pull = (sp - sb) / sigma
    chi2 = np.nansum(pull**2)
    ndof = np.sum(np.isfinite(pull))

    label = op1 if component != 'mixed' else f'{op1}_{op2}'
    norm_label = 'normalized (sum=1)' if not args.no_normalize else 'cross section [pb], N_gen-corrected'
    if args.print_stats:
        print(f"Component: {component}   Operator(s): {label}   [{norm_label}]")
        print(f"{'bin':>16} {'propcorr':>12} {'+-sig_p':>10} {'baseline':>12} {'+-sig_b':>10} {'ratio':>8} {'pull':>7}")
        for i in range(len(edges) - 1):
            lo, hi = edges[i], edges[i + 1]
            r = sp[i] / sb[i] if sb[i] != 0 else float('nan')
            print(f"{lo:6.1f}-{hi:<6.1f}   {sp[i]:12.5g} {sig_p[i]:10.3g} {sb[i]:12.5g} {sig_b[i]:10.3g} {r:8.4f} {pull[i]:7.2f}")
        print(f"\nchi2/dof = {chi2/ndof:.2f} over {ndof} bins\n")

    # --- plot: CMS/mplhep style, shape overlay on top, ratio panel below ---
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, sharex=True, dpi=FIG_DPI,
        gridspec_kw={"height_ratios": [3, 1]})
    fig.tight_layout(pad=-0.5)
    hep.cms.label("Preliminary", data=False, ax=ax_top)

    edges_rep = np.repeat(edges, 2)[1:-1]   # doubled-edges trick for step-shaped fill_between

    ax_top.fill_between(
        edges_rep, np.repeat(sp - sig_p, 2), np.repeat(sp + sig_p, 2),
        step="pre", alpha=0.3, color="crimson", linewidth=0, zorder=1)
    ax_top.fill_between(
        edges_rep, np.repeat(sb - sig_b, 2), np.repeat(sb + sig_b, 2),
        step="pre", alpha=0.3, color="steelblue", linewidth=0, zorder=0)
    ax_top.stairs(sp, edges=edges, color="crimson",   linewidth=2.0, label="propcorr", fill=False, zorder=3)
    ax_top.stairs(sb, edges=edges, color="steelblue", linewidth=2.0, label="baseline", fill=False, zorder=2)
    ax_top.set_ylabel("Normalized shape" if not args.no_normalize else "Cross section [pb / bin]")
    ax_top.tick_params(labelbottom=False)
    ax_top.legend(loc="upper right", framealpha=0.8)
    ax_top.text(
        0.03, 0.97, f"{component}   {label}",
        transform=ax_top.transAxes, va="top", ha="left", fontsize=LEGEND_FONTSIZE,
    )

    ratio = np.divide(sp, sb, out=np.full_like(sp, np.nan), where=sb != 0)
    # propagate independent Poisson uncertainties through the division
    ratio_sigma = np.abs(ratio) * np.sqrt(
        np.divide(sig_p, sp, out=np.zeros_like(sp), where=sp != 0)**2
        + np.divide(sig_b, sb, out=np.zeros_like(sb), where=sb != 0)**2
    )
    ax_bot.fill_between(
        edges_rep, np.repeat(ratio - ratio_sigma, 2), np.repeat(ratio + ratio_sigma, 2),
        step="pre", alpha=0.3, color="gray", linewidth=0, zorder=0)
    ax_bot.stairs(ratio, edges=edges, color="black", linewidth=1.2)
    ax_bot.axhline(1.0, color="gray", linestyle="dashed", linewidth=1)
    if args.lims is not None:
        half = args.lims
    else:
        finite = ratio[np.isfinite(ratio)]
        half = max(np.max(np.abs(finite - 1.0)) * 1.2, 0.05) if finite.size else 0.3
    ax_bot.set_ylim(1.0 - half, 1.0 + half)
    ax_bot.set_ylabel("propcorr / baseline")
    ax_bot.set_xlabel(r"$m_{\ell\ell}$ [GeV]")
    ax_bot.set_xlim(edges[0], edges[-1])

    suffix = "" if not args.no_normalize else "_raw"
    base_name = os.path.join(args.outdir, f"compare_{label}_{component}{suffix}")
    for ext in ("png", "pdf"):
        fig.savefig(f"{base_name}.{ext}", facecolor="white", pad_inches=0.1, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved plot: {base_name}.png / .pdf\n")


if args.component == 'mixed':
    run_comparison('mixed', args.op1, args.op2)
elif args.op2 is not None:
    run_comparison(args.component, args.op1, None)
    run_comparison(args.component, args.op2, None)
else:
    run_comparison(args.component, args.op1, None)
