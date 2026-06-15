#!/usr/bin/env python3
"""
drawNuisancesEFT.py
===================
Plot EFT decomposition (sm, lin, quad) Up/Down nuisance variations from shapes.root.

For each operator × nuisance, computes:
    sm   = histo_sm
    lin  = 0.5 * (histo_w1_op  - histo_wm1_op)
    quad = 0.5 * (histo_w1_op  + histo_wm1_op) - histo_sm

And for Up/Down variations:
    sm_Up   = histo_sm_{nuis}Up
    lin_Up  = 0.5 * (histo_w1_op_{nuis}Up - histo_wm1_op_{nuis}Up)
    quad_Up = 0.5 * (histo_w1_op_{nuis}Up + histo_wm1_op_{nuis}Up) - histo_sm_{nuis}Up

Usage:
    drawNuisancesEFT.py --shapes shapes.root --outdir plots/eft_nuisances
    drawNuisancesEFT.py --shapes shapes.root --operators cHDD cHWB --nuisances QCDScale PDFweight
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

NUISANCES = [
    "QCDScale", "PDFweight", "alphaS", "PSWeight",
    "mu_reco", "mu_idiso", "mu_trig", "PU", "prefireWeight",
    "rochester_stat", "rochester_syst",
]

SM_COLOR   = "#5790fc"
LIN_COLOR  = "#f89c20"
QUAD_COLOR = "#e42536"

ratio_fig_style = {
    "figsize": (10, 10),
    "gridspec_kw": {"height_ratios": (3, 1)},
}


def get_vals(f, key):
    return f[key].values().copy()


def get_variances(f, key):
    return f[key].variances().copy()


def get_edges(f, key):
    return f[key].axes[0].edges()


def decompose(f, op, suffix=""):
    """Return (sm, lin, quad) values for a given operator and histogram suffix (e.g. '_QCDScaleUp')."""
    sm   = get_vals(f, f"histo_sm{suffix}")
    w1   = get_vals(f, f"histo_w1_{op}{suffix}")
    wm1  = get_vals(f, f"histo_wm1_{op}{suffix}")
    lin  = 0.5 * (w1 - wm1)
    quad = 0.5 * (w1 + wm1) - sm
    return sm, lin, quad


def plot_task(d):
    shapes_file = d["shapes"]
    op          = d["op"]
    nuis        = d["nuis"]
    outdir      = d["outdir"]
    logy        = d["logy"]

    try:
        f = uproot.open(shapes_file)

        edges  = get_edges(f, "histo_sm")
        widths = np.diff(edges)
        centers = 0.5 * (edges[:-1] + edges[1:])

        # nominal decomposition
        sm_nom,  lin_nom,  quad_nom  = decompose(f, op)

        # up/down decomposition
        sm_up,   lin_up,   quad_up   = decompose(f, op, f"_{nuis}Up")
        sm_down, lin_down, quad_down = decompose(f, op, f"_{nuis}Down")

    except Exception as e:
        print(f"  [skip] {op} / {nuis}: {e}")
        return

    _UP_COLOR   = "#f89c20"
    _DOWN_COLOR = "#e42536"

    def make_axes(vals_nom, vals_up, vals_down, label, color):
        """Return (fig, ax_top, ax_bot) for one component."""
        fig, (ax, rax) = plt.subplots(2, 1, sharex=True, **ratio_fig_style)
        fig.subplots_adjust(hspace=0.07)

        # top panel — events/GeV
        def stairs(ax_, vals, ls="-", lw=1.5, col=color, lab=None):
            ax_.stairs(vals / widths, edges=edges, color=col, linewidth=lw,
                       linestyle=ls, label=lab, fill=False)

        stairs(ax, vals_nom,  ls="-",  lw=2.0, lab="Nominal")
        stairs(ax, vals_up,   ls="--", lw=1.5, col=_UP_COLOR,   lab=f"{nuis} Up")
        stairs(ax, vals_down, ls=":",  lw=1.5, col=_DOWN_COLOR,  lab=f"{nuis} Down")

        ax.set_ylabel("Events / GeV")
        ax.set_title(f"{label}  [{op}]  —  {nuis}", fontsize=13)
        ax.legend(loc="best")
        # lin can be negative → always linear; sm/quad use --logy if requested
        if logy and label != "lin":
            ax.set_yscale("log")

        hep.cms.label(loc=0, label="Preliminary", data=False, ax=ax)

        # ratio panel
        safe = np.where(np.abs(vals_nom) > 0, vals_nom, np.nan)
        rax.stairs(vals_up   / safe, edges=edges, color=_UP_COLOR,   linewidth=1.5, linestyle="--", label=f"{nuis} Up")
        rax.stairs(vals_down / safe, edges=edges, color=_DOWN_COLOR, linewidth=1.5, linestyle=":",  label=f"{nuis} Down")
        rax.axhline(1.0, color="black", linewidth=0.8, linestyle="dashed")
        rax.set_ylim(0.7, 1.3)
        rax.set_ylabel("Var / Nom.")
        rax.set_xlabel(r"$m_{\ell\ell}$ (GeV)")
        rax.autoscale(axis="x", tight=True)

        return fig

    for component, vals_nom, vals_up, vals_down, color in [
        ("sm",   sm_nom,   sm_up,   sm_down,   SM_COLOR),
        ("lin",  lin_nom,  lin_up,  lin_down,  LIN_COLOR),
        ("quad", quad_nom, quad_up, quad_down, QUAD_COLOR),
    ]:
        fig = make_axes(vals_nom, vals_up, vals_down, component, color)
        stem = os.path.join(outdir, f"{component}_{op}_{nuis}")
        fig.savefig(f"{stem}.png", bbox_inches="tight", facecolor="white")
        fig.savefig(f"{stem}.pdf", bbox_inches="tight", facecolor="white")
        plt.close(fig)
        print(f"  {component:4s}  {op:10s}  {nuis:20s}  ->  {stem}.png")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--shapes",     required=True,        help="Path to shapes.root")
    parser.add_argument("--outdir",     default="plots/eft_nuisances")
    parser.add_argument("--operators",  nargs="+", default=OPERATORS)
    parser.add_argument("--nuisances",  nargs="+", default=NUISANCES)
    parser.add_argument("--logy",       action="store_true")
    parser.add_argument("--ncores",     type=int, default=4)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # check which nuisances actually exist in shapes.root
    f = uproot.open(args.shapes)
    available_keys = [k.split(";")[0] for k in f.keys()]
    available_nuis = [n for n in args.nuisances
                      if any(f"_sm_{n}Up" in k or f"histo_sm_{n}Up" in k for k in available_keys)]
    if not available_nuis:
        # fallback: check without histo_ prefix
        available_nuis = [n for n in args.nuisances
                          if any(f"{n}Up" in k for k in available_keys)]
    print(f"Shapes file : {args.shapes}")
    print(f"Operators   : {args.operators}")
    print(f"Nuisances   : {available_nuis}")
    print(f"Output      : {args.outdir}\n")

    tasks = [
        {"shapes": args.shapes, "op": op, "nuis": nuis,
         "outdir": args.outdir, "logy": args.logy}
        for op in args.operators
        for nuis in available_nuis
    ]

    with Pool(processes=args.ncores) as pool:
        pool.map(plot_task, tasks)

    print("\nDone.")


if __name__ == "__main__":
    main()
