#!/usr/bin/env python3
"""
drawNuisancesEFT.py
===================
Plot EFT decomposition (sm, lin, quad) Up/Down nuisance variations from
shapes.root (RECO morphing pipeline) or histograms.root (LHE pipeline).

RECO morphing mode (flat keys, "histo_" prefix):
    sm   = histo_sm
    lin  = 0.5 * (histo_w1_op  - histo_wm1_op)
    quad = 0.5 * (histo_w1_op  + histo_wm1_op) - histo_sm
    nuisance suffix: histo_sm_{nuis}Up / histo_w1_op_{nuis}Up / ...

LHE mode (keys nested under a channel directory, e.g. "triple_DY/...",
processes already decomposed by build_datacard_new.py):
    sm   = {channel}/sm
    quad = {channel}/quad_op
    lin  = {channel}/sm_lin_quad_op - sm - quad   (C=1 reference point)
    nuisance suffix: {channel}/sm_{nuis}Up / {channel}/quad_op_{nuis}Up / ...
    (nuisance names are "muf_scale"/"pdf", not "QCDScale"/"PDFweight" —
    common aliases are accepted, see NUISANCE_ALIASES)

The pipeline is auto-detected from the file structure; no flag needed.

Usage:
    drawNuisancesEFT.py --shapes shapes.root      --outdir plots/eft_nuisances
    drawNuisancesEFT.py --shapes histograms.root  --outdir plots/eft_nuisances_lhe
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

NUISANCES_MORPHING = [
    "QCDScale", "PDFweight", "alphaS", "PSWeight",
    "mu_reco", "mu_idiso", "mu_trig", "PU", "prefireWeight",
    "rochester_stat", "rochester_syst",
]

NUISANCES_LHE = ["muf_scale", "pdf"]

# accepted aliases when running against the LHE pipeline, so the same
# --nuisances spelling used for the morphing pipeline still works
NUISANCE_ALIASES_LHE = {
    "QCDScale": "muf_scale", "QCDscale": "muf_scale", "qcdscale": "muf_scale",
    "PDFweight": "pdf", "PDF": "pdf", "pdfweight": "pdf",
}


def detect_mode(f):
    """Return ('morphing', None) or ('lhe', channel_name) based on file layout."""
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

SM_COLOR   = "#5790fc"
LIN_COLOR  = "#f89c20"
QUAD_COLOR = "#e42536"

VAR_XLABELS = {
    "mll":          r"$m_{\ell\ell}$ (GeV)",
    "costhetastar": r"$\cos\theta^*$",
    "rapll_abs":    r"$|y_{\ell\ell}|$",
    "triple_diff":  "Unrolled bin",
}

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


def decompose(f, op, mode, channel, suffix=""):
    """Return (sm, lin, quad) values for a given operator and histogram suffix
    (e.g. '_QCDScaleUp' for morphing, '_muf_scaleUp' for LHE)."""
    if mode == "morphing":
        sm   = get_vals(f, f"histo_sm{suffix}")
        w1   = get_vals(f, f"histo_w1_{op}{suffix}")
        wm1  = get_vals(f, f"histo_wm1_{op}{suffix}")
        lin  = 0.5 * (w1 - wm1)
        quad = 0.5 * (w1 + wm1) - sm
    else:  # lhe
        sm   = get_vals(f, f"{channel}/sm{suffix}")
        quad = get_vals(f, f"{channel}/quad_{op}{suffix}")
        slq  = get_vals(f, f"{channel}/sm_lin_quad_{op}{suffix}")
        lin  = slq - sm - quad   # sm_lin_quad = sm + C*lin + C^2*quad, C=1 reference
    return sm, lin, quad


def render_panel(edges, widths, vals_nom, vals_up, vals_down, nuis, logy, title, color, xlabel):
    """Build and return a (Nominal/Up/Down + ratio) figure for one histogram component."""
    _UP_COLOR   = "#f89c20"
    _DOWN_COLOR = "#e42536"

    fig, (ax, rax) = plt.subplots(2, 1, sharex=True, **ratio_fig_style)
    fig.subplots_adjust(hspace=0.07)

    def stairs(ax_, vals, ls="-", lw=1.5, col=color, lab=None):
        ax_.stairs(vals / widths, edges=edges, color=col, linewidth=lw,
                   linestyle=ls, label=lab, fill=False)

    stairs(ax, vals_nom,  ls="-",  lw=2.0, lab="Nominal")
    stairs(ax, vals_up,   ls="--", lw=1.5, col=_UP_COLOR,   lab=f"{nuis} Up")
    stairs(ax, vals_down, ls=":",  lw=1.5, col=_DOWN_COLOR,  lab=f"{nuis} Down")

    ax.set_ylabel("Events / GeV")
    # title in top-right corner inside the axes — avoids overlap with CMS label (top-left)
    ax.text(0.97, 0.97, title, transform=ax.transAxes,
            ha="right", va="top", fontsize=18)
    ax.legend(loc="upper left")
    # lin can be negative → always linear; sm/quad use --logy if requested
    if logy and not title.startswith("lin"):
        ax.set_yscale("log")

    hep.cms.label(loc=0, label="Preliminary", data=False, ax=ax)

    safe = np.where(np.abs(vals_nom) > 0, vals_nom, np.nan)
    rax.stairs(vals_up   / safe, edges=edges, color=_UP_COLOR,   linewidth=1.5, linestyle="--", label=f"{nuis} Up")
    rax.stairs(vals_down / safe, edges=edges, color=_DOWN_COLOR, linewidth=1.5, linestyle=":",  label=f"{nuis} Down")
    rax.axhline(1.0, color="black", linewidth=0.8, linestyle="dashed")
    rax.set_ylim(0.7, 1.3)
    rax.set_ylabel("Var / Nom.")
    rax.set_xlabel(xlabel)
    rax.autoscale(axis="x", tight=True)

    return fig


def plot_sm_task(d):
    """sm doesn't depend on the operator (C=0 reference) — plot it once per nuisance."""
    shapes_file = d["shapes"]
    nuis        = d["nuis"]
    outdir      = d["outdir"]
    logy        = d["logy"]
    mode        = d["mode"]
    channel     = d["channel"]
    xlabel      = d["xlabel"]

    try:
        f = uproot.open(shapes_file)
        sm_key = "histo_sm" if mode == "morphing" else f"{channel}/sm"
        edges  = get_edges(f, sm_key)
        widths = np.diff(edges)

        sm_nom  = get_vals(f, sm_key)
        sm_up   = get_vals(f, f"{sm_key}_{nuis}Up")
        sm_down = get_vals(f, f"{sm_key}_{nuis}Down")
    except Exception as e:
        print(f"  [skip] sm / {nuis}: {e}")
        return

    fig = render_panel(edges, widths, sm_nom, sm_up, sm_down, nuis, logy,
                        title=f"sm  —  {nuis}", color=SM_COLOR, xlabel=xlabel)
    stem = os.path.join(outdir, f"sm_{nuis}")
    fig.savefig(f"{stem}.png", bbox_inches="tight", facecolor="white")
    fig.savefig(f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  {'sm':4s}  {'(shared)':10s}  {nuis:20s}  ->  {stem}.png")


def plot_task(d):
    """lin/quad depend on the operator — one plot per (op, nuisance)."""
    shapes_file = d["shapes"]
    op          = d["op"]
    nuis        = d["nuis"]
    outdir      = d["outdir"]
    logy        = d["logy"]
    mode        = d["mode"]
    channel     = d["channel"]
    xlabel      = d["xlabel"]

    try:
        f = uproot.open(shapes_file)

        sm_key = "histo_sm" if mode == "morphing" else f"{channel}/sm"
        edges  = get_edges(f, sm_key)
        widths = np.diff(edges)

        # nominal decomposition
        sm_nom,  lin_nom,  quad_nom  = decompose(f, op, mode, channel)

        # up/down decomposition
        sm_up,   lin_up,   quad_up   = decompose(f, op, mode, channel, f"_{nuis}Up")
        sm_down, lin_down, quad_down = decompose(f, op, mode, channel, f"_{nuis}Down")

    except Exception as e:
        print(f"  [skip] {op} / {nuis}: {e}")
        return

    for component, vals_nom, vals_up, vals_down, color in [
        ("lin",  lin_nom,  lin_up,  lin_down,  LIN_COLOR),
        ("quad", quad_nom, quad_up, quad_down, QUAD_COLOR),
    ]:
        title = f"{component}  [{op}]  —  {nuis}"
        fig = render_panel(edges, widths, vals_nom, vals_up, vals_down, nuis, logy,
                            title=title, color=color, xlabel=xlabel)
        stem = os.path.join(outdir, f"{component}_{op}_{nuis}")
        fig.savefig(f"{stem}.png", bbox_inches="tight", facecolor="white")
        fig.savefig(f"{stem}.pdf", bbox_inches="tight", facecolor="white")
        plt.close(fig)
        print(f"  {component:4s}  {op:10s}  {nuis:20s}  ->  {stem}.png")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--shapes",     required=True,        help="Path to shapes.root or histograms.root")
    parser.add_argument("--outdir",     default="plots/eft_nuisances")
    parser.add_argument("--variable",   default=None,
                        choices=list(VAR_XLABELS.keys()),
                        help="Observable plotted (sets x-axis label). "
                             "Auto-detected from the shapes path if not given.")
    parser.add_argument("--operators",  nargs="+", default=OPERATORS)
    parser.add_argument("--nuisances",  nargs="+", default=None,
                        help="Default: all morphing nuisances, or "
                             "['muf_scale', 'pdf'] in LHE mode. "
                             "QCDScale/PDFweight aliases accepted in LHE mode.")
    parser.add_argument("--logy",       action="store_true")
    parser.add_argument("--ncores",     type=int, default=4)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    f = uproot.open(args.shapes)
    mode, channel = detect_mode(f)
    available_keys = [k.split(";")[0] for k in f.keys()]

    if mode == "lhe":
        # strip channel prefix so suffix-matching below works the same as morphing mode
        available_keys = [k.split("/", 1)[1] if k.startswith(f"{channel}/") else k
                          for k in available_keys]
        requested = args.nuisances if args.nuisances is not None else NUISANCES_LHE
        requested = [NUISANCE_ALIASES_LHE.get(n, n) for n in requested]
    else:
        requested = args.nuisances if args.nuisances is not None else NUISANCES_MORPHING

    # check which nuisances actually exist in the file
    available_nuis = [n for n in requested
                      if any(f"_sm_{n}Up" in k or f"sm_{n}Up" in k for k in available_keys)]
    if not available_nuis:
        # fallback: looser match, no "sm_" anchor
        available_nuis = [n for n in requested
                          if any(f"{n}Up" in k for k in available_keys)]

    print(f"Shapes file : {args.shapes}")
    print(f"Pipeline    : {mode}" + (f"  (channel: {channel})" if mode == "lhe" else ""))
    print(f"Operators   : {args.operators}")
    print(f"Nuisances   : {available_nuis}")
    print(f"Output      : {args.outdir}\n")

    # auto-detect variable from histogram shape (works even when running from the variable dir)
    if args.variable:
        variable = args.variable
    else:
        f_tmp = uproot.open(args.shapes)
        sm_key_tmp = "histo_sm" if mode == "morphing" else f"{channel}/sm"
        edges_tmp  = get_edges(f_tmp, sm_key_tmp)
        n = len(edges_tmp) - 1
        if n == 200:
            variable = "triple_diff"
        elif n == 50 and edges_tmp[-1] > 2.0:
            variable = "rapll_abs"
        elif n == 50:
            variable = "costhetastar"
        else:
            variable = "mll"
        f_tmp.close()

    xlabel = VAR_XLABELS[variable]
    print(f"Variable    : {variable}  →  x-label: {xlabel}")

    sm_tasks = [
        {"shapes": args.shapes, "nuis": nuis,
         "outdir": args.outdir, "logy": args.logy,
         "mode": mode, "channel": channel, "xlabel": xlabel}
        for nuis in available_nuis
    ]
    eft_tasks = [
        {"shapes": args.shapes, "op": op, "nuis": nuis,
         "outdir": args.outdir, "logy": args.logy,
         "mode": mode, "channel": channel, "xlabel": xlabel}
        for op in args.operators
        for nuis in available_nuis
    ]

    with Pool(processes=args.ncores) as pool:
        pool.map(plot_sm_task, sm_tasks)
        pool.map(plot_task, eft_tasks)

    print("\nDone.")


if __name__ == "__main__":
    main()
