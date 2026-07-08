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

PROP_CACHE = '/grid_mnt/data__data.polcms/cms/adufour/LHE/propcorr/CACHE/lhe_cache_propcorr.pkl'
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
args = parser.parse_args()
os.makedirs(args.outdir, exist_ok=True)

with open(args.prop_cache, 'rb') as f:
    prop = pickle.load(f)
with open(args.base_cache, 'rb') as f:
    base = pickle.load(f)

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


def shape_and_n(cache, component, op1, op2, C):
    w = get_weight(cache, component, op1, op2, C)
    mll = cache['mll']
    h, _ = np.histogram(mll, bins=edges, weights=w)
    n, _ = np.histogram(mll, bins=edges)
    total = h.sum()
    return (h / total if total != 0 else h), n


def run_comparison(component, op1, op2):
    sp, n_prop = shape_and_n(prop, component, op1, op2, args.C)
    sb, n_base = shape_and_n(base, component, op1, op2, args.C)

    sig_p = sp / np.sqrt(np.maximum(n_prop, 1))
    sig_b = sb / np.sqrt(np.maximum(n_base, 1))
    sigma = np.sqrt(sig_p**2 + sig_b**2)
    sigma = np.where(sigma == 0, np.inf, sigma)
    pull = (sp - sb) / sigma
    chi2 = np.nansum(pull**2)
    ndof = np.sum(np.isfinite(pull))

    label = op1 if component != 'mixed' else f'{op1}_{op2}'
    print(f"Component: {component}   Operator(s): {label}")
    print(f"{'bin':>16} {'propcorr':>10} {'baseline':>10} {'ratio':>8} {'pull':>7} {'N_prop':>8} {'N_base':>8}")
    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        r = sp[i] / sb[i] if sb[i] != 0 else float('nan')
        print(f"{lo:6.1f}-{hi:<6.1f}   {sp[i]:10.5f} {sb[i]:10.5f} {r:8.4f} {pull[i]:7.2f} {n_prop[i]:8d} {n_base[i]:8d}")
    print(f"\nchi2/dof = {chi2/ndof:.2f} over {ndof} bins\n")

    # --- plot: CMS/mplhep style, shape overlay on top, ratio panel below ---
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, sharex=True, dpi=FIG_DPI,
        gridspec_kw={"height_ratios": [3, 1]})
    fig.tight_layout(pad=-0.5)
    hep.cms.label("Simulation", data=False, ax=ax_top)

    ax_top.stairs(sp, edges=edges, color="crimson",   linewidth=2.0, label="propcorr", fill=False, zorder=3)
    ax_top.stairs(sb, edges=edges, color="steelblue", linewidth=2.0, label="baseline", fill=False, zorder=2)
    ax_top.set_ylabel("Normalized shape")
    ax_top.tick_params(labelbottom=False)
    ax_top.legend(loc="upper right", framealpha=0.8)
    ax_top.text(
        0.03, 0.97, f"{component}   {label}",
        transform=ax_top.transAxes, va="top", ha="left", fontsize=LEGEND_FONTSIZE,
    )

    ratio = np.divide(sp, sb, out=np.full_like(sp, np.nan), where=sb != 0)
    ax_bot.stairs(ratio, edges=edges, color="black", linewidth=1.2)
    ax_bot.axhline(1.0, color="gray", linestyle="dashed", linewidth=1)
    finite = ratio[np.isfinite(ratio)]
    half = max(np.max(np.abs(finite - 1.0)) * 1.2, 0.05) if finite.size else 0.3
    ax_bot.set_ylim(1.0 - half, 1.0 + half)
    ax_bot.set_ylabel("propcorr / baseline")
    ax_bot.set_xlabel(r"$m_{\ell\ell}$ [GeV]")
    ax_bot.set_xlim(edges[0], edges[-1])

    base_name = os.path.join(args.outdir, f"compare_{label}_{component}")
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
