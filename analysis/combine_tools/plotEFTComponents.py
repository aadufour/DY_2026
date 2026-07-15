#!/usr/bin/env python3
"""
plotEFTComponents.py
====================
Plot the EFT decomposition (sm, lin, quad) nominal histograms for each operator.
No nuisance variations — clean "components" plots suitable for presentations.

RECO morphing mode  (flat keys, "histo_" prefix):
    sm   = histo_sm
    lin  = 0.5 * (histo_w1_op  - histo_wm1_op)
    quad = 0.5 * (histo_w1_op  + histo_wm1_op) - histo_sm

LHE mode  (nested under a channel directory):
    sm   = {channel}/sm
    quad = {channel}/quad_op
    lin  = {channel}/sm_lin_quad_op - sm - quad

The pipeline is auto-detected from the file structure.

Usage:
    plotEFTComponents.py --shapes shapes.root --outdir plots/eft_components
    plotEFTComponents.py --shapes shapes.root --operators cHDD cHWB --logy
    plotEFTComponents.py --shapes shapes.root --c-values 0.5 1.0 2.0 --full
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
FULL_COLORS = ["#34a853", "#9467bd", "#8c564b"]  # for c-values overlay

VAR_XLABELS = {
    "mll":          r"$m_{\ell\ell}$ (GeV)",
    "costhetastar": r"$\cos\theta^*$",
    "rapll_abs":    r"$|y_{\ell\ell}|$",
    "triple_diff":  "Unrolled bin",
}


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
    edges = get_edges(f, sm_key)
    n = len(edges) - 1
    if n == 200:
        return "triple_diff"
    elif n == 50 and edges[-1] > 2.0:
        return "rapll_abs"
    elif n == 50:
        return "costhetastar"
    return "mll"


def plot_task(d):
    shapes_file = d["shapes"]
    op          = d["op"]
    outdir      = d["outdir"]
    logy        = d["logy"]
    mode        = d["mode"]
    channel     = d["channel"]
    xlabel      = d["xlabel"]
    c_values    = d["c_values"]
    show_full   = d["show_full"]

    try:
        f = uproot.open(shapes_file)
        sm_key = "histo_sm" if mode == "morphing" else f"{channel}/sm"
        edges  = get_edges(f, sm_key)
        widths = np.diff(edges)

        sm, lin, quad = decompose(f, op, mode, channel)
        f.close()
    except Exception as e:
        print(f"  [skip] {op}: {e}")
        return

    centers = 0.5 * (edges[:-1] + edges[1:])

    if show_full:
        fig, (ax, rax) = plt.subplots(
            2, 1, sharex=True,
            figsize=(10, 10),
            gridspec_kw={"height_ratios": (3, 1)},
        )
        fig.subplots_adjust(hspace=0.07)
    else:
        fig, ax = plt.subplots(figsize=(10, 7))
        rax = None

    def stairs(ax_, vals, color, label, ls="-", lw=2.0):
        ax_.stairs(vals / widths, edges=edges, color=color,
                   linewidth=lw, linestyle=ls, label=label, fill=False)

    stairs(ax, sm,   SM_COLOR,   "SM",   lw=2.5)
    stairs(ax, lin,  LIN_COLOR,  "lin",  lw=2.0, ls="--")
    stairs(ax, quad, QUAD_COLOR, "quad", lw=2.0, ls=":")

    if show_full:
        for cv, col in zip(c_values, FULL_COLORS):
            full = sm + cv * lin + cv**2 * quad
            stairs(ax, full, col, fr"full  $c={cv}$", lw=1.5, ls="-.")

    ax.set_ylabel("Events / unit")
    ax.text(0.97, 0.97, op, transform=ax.transAxes,
            ha="right", va="top", fontsize=22, fontweight="bold")
    ax.legend(loc="upper left", fontsize=16)
    if logy:
        ax.set_yscale("log")

    hep.cms.label(loc=0, label="Preliminary", data=False, ax=ax)

    if show_full and rax is not None:
        safe = np.where(np.abs(sm) > 0, sm, np.nan)
        for cv, col in zip(c_values, FULL_COLORS):
            full = sm + cv * lin + cv**2 * quad
            rax.stairs(full / safe, edges=edges, color=col,
                       linewidth=1.5, linestyle="-.", label=fr"$c={cv}$")
        rax.axhline(1.0, color="black", linewidth=0.8, linestyle="dashed")
        rax.set_ylim(0.5, 1.5)
        rax.set_ylabel("Full / SM")
        rax.set_xlabel(xlabel)
        rax.legend(loc="upper left", fontsize=14)
        rax.autoscale(axis="x", tight=True)
    elif rax is None:
        ax.set_xlabel(xlabel)

    stem = os.path.join(outdir, f"components_{op}")
    fig.savefig(f"{stem}.png", bbox_inches="tight", facecolor="white")
    fig.savefig(f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  {op:12s}  ->  {stem}.png")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--shapes",    required=True, help="Path to shapes.root or histograms.root")
    parser.add_argument("--outdir",    default="plots/eft_components")
    parser.add_argument("--variable",  default=None, choices=list(VAR_XLABELS.keys()),
                        help="x-axis label variable (auto-detected if not given)")
    parser.add_argument("--operators", nargs="+", default=OPERATORS)
    parser.add_argument("--logy",      action="store_true")
    parser.add_argument("--full",      action="store_true",
                        help="Also overlay sm+c*lin+c²*quad for each --c-values")
    parser.add_argument("--c-values",  nargs="+", type=float, default=[1.0],
                        help="Wilson coefficient values to use for the full prediction overlay "
                             "(only relevant with --full). Default: 1.0")
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
    print(f"Output      : {args.outdir}\n")

    tasks = [
        {
            "shapes":    args.shapes,
            "op":        op,
            "outdir":    args.outdir,
            "logy":      args.logy,
            "mode":      mode,
            "channel":   channel,
            "xlabel":    xlabel,
            "c_values":  args.c_values,
            "show_full": args.full,
        }
        for op in args.operators
    ]

    with Pool(processes=args.ncores) as pool:
        pool.map(plot_task, tasks)

    print("\nDone.")


if __name__ == "__main__":
    main()
