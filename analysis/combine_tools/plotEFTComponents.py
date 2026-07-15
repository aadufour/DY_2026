#!/usr/bin/env python3
"""
plotEFTComponents.py
====================
Plot the EFT decomposition (sm, lin, quad) nominal histograms for each operator.
No nuisance variations — clean "components" plots suitable for presentations.

Three separate figures are produced per operator:
  sm_full_{op}  —  SM + full prediction at c=1 (and extra c values if given)
  lin_{op}      —  linear term at c=1  (never log-scaled; can be negative)
  quad_{op}     —  quadratic term at c=1

RECO morphing mode (flat keys, "histo_" prefix):
    sm   = histo_sm
    lin  = 0.5 * (histo_w1_op  - histo_wm1_op)
    quad = 0.5 * (histo_w1_op  + histo_wm1_op) - histo_sm

LHE mode (nested under a channel directory):
    sm   = {channel}/sm
    quad = {channel}/quad_op
    lin  = {channel}/sm_lin_quad_op - sm - quad

The pipeline is auto-detected from the file structure.

Usage:
    plotEFTComponents.py --shapes shapes.root --outdir plots/eft_components
    plotEFTComponents.py --shapes shapes.root --operators cHDD cHWB
    plotEFTComponents.py --shapes shapes.root --c-values 0.5 1.0 2.0 --logy
"""

import argparse
import os
from multiprocessing import Pool

import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import uproot

hep.style.use("CMS")

OPERATORS = [
    "cHDD", "cHWB", "cbWRe", "cbBRe", "cHj1", "cHQ1", "cHj3", "cHQ3",
    "cHu", "cHd", "cHbq", "cHl1", "cHl3", "cHe", "cll1", "clj1", "clj3",
    "cQl1", "cQl3", "ceu", "ced", "cbe", "cje", "cQe", "clu", "cld", "cbl",
]

SM_COLOR   = "#5790fc"
LIN_COLOR  = "#f89c20"
QUAD_COLOR = "#e42536"
# colors for additional c values beyond c=1
EXTRA_COLORS = ["#9467bd", "#8c564b", "#17becf"]

VAR_XLABELS = {
    "mll":          r"$m_{\ell\ell}$ (GeV)",
    "costhetastar": r"$\cos\theta^*$",
    "rapll_abs":    r"$|y_{\ell\ell}|$",
    "triple_diff":  "Unrolled bin",
}

FIG_STYLE = {"figsize": (10, 7)}


# ---------------------------------------------------------------------------
# helpers shared with drawNuisancesEFT
# ---------------------------------------------------------------------------

def detect_mode(f):
    keys = [k.split(";")[0] for k in f.keys()]
    if any(k == "histo_sm" for k in keys):
        return "morphing", None
    for k in keys:
        if k.endswith("/sm"):
            return "lhe", k.rsplit("/", 1)[0]
    raise RuntimeError(
        "Could not detect pipeline: no 'histo_sm' (morphing) or '<channel>/sm' "
        "(LHE) key found in file."
    )


def get_vals(f, key):
    return f[key].values().copy()


def get_edges(f, key):
    return f[key].axes[0].edges()


def decompose(f, op, mode, channel):
    if mode == "morphing":
        sm   = get_vals(f, "histo_sm")
        w1   = get_vals(f, f"histo_w1_{op}")
        wm1  = get_vals(f, f"histo_wm1_{op}")
        lin  = 0.5 * (w1 - wm1)
        quad = 0.5 * (w1 + wm1) - sm
    else:
        sm   = get_vals(f, f"{channel}/sm")
        quad = get_vals(f, f"{channel}/quad_{op}")
        slq  = get_vals(f, f"{channel}/sm_lin_quad_{op}")
        lin  = slq - sm - quad
    return sm, lin, quad


def autodetect_variable(f, mode, channel):
    sm_key = "histo_sm" if mode == "morphing" else f"{channel}/sm"
    edges  = get_edges(f, sm_key)
    n = len(edges) - 1
    if n == 200:
        return "triple_diff"
    elif n == 50 and edges[-1] > 2.0:
        return "rapll_abs"
    elif n == 50:
        return "costhetastar"
    return "mll"


def _save(fig, stem):
    fig.savefig(f"{stem}.png", bbox_inches="tight", facecolor="white")
    fig.savefig(f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _stairs(ax, vals, edges, widths, color, label, ls="-", lw=2.0):
    ax.stairs(vals / widths, edges=edges, color=color,
              linewidth=lw, linestyle=ls, label=label, fill=False)


def _decorate(ax, xlabel, title, logy=False):
    ax.set_ylabel("Events / unit")
    ax.set_xlabel(xlabel)
    ax.text(0.97, 0.97, title, transform=ax.transAxes,
            ha="right", va="top", fontsize=20, fontweight="bold")
    ax.legend(loc="upper left", fontsize=16)
    if logy:
        ax.set_yscale("log")
    hep.cms.label(loc=0, label="Preliminary", data=False, ax=ax)


# ---------------------------------------------------------------------------
# per-operator plotting
# ---------------------------------------------------------------------------

def plot_task(d):
    shapes_file = d["shapes"]
    op          = d["op"]
    outdir      = d["outdir"]
    logy        = d["logy"]
    mode        = d["mode"]
    channel     = d["channel"]
    xlabel      = d["xlabel"]
    c_values    = d["c_values"]   # list of floats; c=1 is always included first

    try:
        f      = uproot.open(shapes_file)
        sm_key = "histo_sm" if mode == "morphing" else f"{channel}/sm"
        edges  = get_edges(f, sm_key)
        widths = np.diff(edges)
        sm, lin, quad = decompose(f, op, mode, channel)
        f.close()
    except Exception as e:
        print(f"  [skip] {op}: {e}")
        return

    # ---- figure 1: SM + full prediction(s) --------------------------------
    fig1, ax1 = plt.subplots(**FIG_STYLE)
    _stairs(ax1, sm, edges, widths, SM_COLOR, "SM", lw=2.5)
    for cv, col in zip(c_values, [LIN_COLOR] + EXTRA_COLORS):
        full = sm + cv * lin + cv**2 * quad
        label = fr"SM + EFT  ($c={cv}$)" if len(c_values) > 1 else r"SM + EFT  ($c=1$)"
        _stairs(ax1, full, edges, widths, col, label, ls="--", lw=2.0)
    _decorate(ax1, xlabel, op, logy=logy)
    _save(fig1, os.path.join(outdir, f"sm_full_{op}"))

    # ---- figure 2: linear term at c=1 (no log scale — can be negative) ----
    fig2, ax2 = plt.subplots(**FIG_STYLE)
    _stairs(ax2, lin, edges, widths, LIN_COLOR, r"linear term  ($c=1$)", lw=2.5)
    ax2.axhline(0, color="black", linewidth=0.8, linestyle="dashed")
    _decorate(ax2, xlabel, op, logy=False)
    _save(fig2, os.path.join(outdir, f"lin_{op}"))

    # ---- figure 3: quadratic term at c=1 ----------------------------------
    fig3, ax3 = plt.subplots(**FIG_STYLE)
    _stairs(ax3, quad, edges, widths, QUAD_COLOR, r"quadratic term  ($c=1$)", lw=2.5)
    _decorate(ax3, xlabel, op, logy=logy)
    _save(fig3, os.path.join(outdir, f"quad_{op}"))

    print(f"  {op:12s}  ->  sm_full / lin / quad  [{outdir}]")


# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--shapes",    required=True,
                        help="Path to shapes.root or histograms.root")
    parser.add_argument("--outdir",    default="plots/eft_components")
    parser.add_argument("--variable",  default=None, choices=list(VAR_XLABELS.keys()),
                        help="x-axis label variable (auto-detected if not given)")
    parser.add_argument("--operators", nargs="+", default=OPERATORS)
    parser.add_argument("--logy",      action="store_true",
                        help="Log y-scale on sm_full and quad panels (never applied to lin)")
    parser.add_argument("--c-values",  nargs="+", type=float, default=[1.0],
                        help="Wilson coefficient values for the full-prediction overlay. "
                             "Default: 1.0")
    parser.add_argument("--ncores",    type=int, default=4)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    f = uproot.open(args.shapes)
    mode, channel = detect_mode(f)
    variable = args.variable or autodetect_variable(f, mode, channel)
    xlabel   = VAR_XLABELS[variable]
    f.close()

    print(f"Shapes file : {args.shapes}")
    print(f"Pipeline    : {mode}" + (f"  (channel: {channel})" if mode == "lhe" else ""))
    print(f"Variable    : {variable}  →  {xlabel}")
    print(f"Operators   : {args.operators}")
    print(f"c values    : {args.c_values}")
    print(f"Output      : {args.outdir}\n")

    tasks = [
        {
            "shapes":   args.shapes,
            "op":       op,
            "outdir":   args.outdir,
            "logy":     args.logy,
            "mode":     mode,
            "channel":  channel,
            "xlabel":   xlabel,
            "c_values": args.c_values,
        }
        for op in args.operators
    ]

    with Pool(processes=args.ncores) as pool:
        pool.map(plot_task, tasks)

    print("\nDone.")


if __name__ == "__main__":
    main()
