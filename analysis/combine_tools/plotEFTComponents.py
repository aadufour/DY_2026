#!/usr/bin/env python3
"""
plotEFTComponents.py
====================
Plot the EFT decomposition (sm, lin, quad) nominal histograms for each operator.
MC statistical uncertainty band shown as a shaded region.
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
    plotEFTComponents.py --shapes shapes.root --c-values 0.5 1.0 2.0
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
EXTRA_COLORS = ["#9467bd", "#8c564b", "#17becf"]

VAR_XLABELS = {
    "mll":          r"$m_{\ell\ell}$ (GeV)",
    "costhetastar": r"$\cos\theta^*$",
    "rapll_abs":    r"$|y_{\ell\ell}|$",
    "triple_diff":  "Unrolled bin",
}

VAR_LOGX = {
    "mll":          True,
    "costhetastar": False,
    "rapll_abs":    False,
    "triple_diff":  False,
}

FIG_STYLE = {
    "figsize": (10, 10),
    "gridspec_kw": {"height_ratios": (3, 1)},
}


# ---------------------------------------------------------------------------
# helpers
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


def get_variances(f, key):
    return f[key].variances().copy()


def get_edges(f, key):
    return f[key].axes[0].edges()


def decompose(f, op, mode, channel):
    """Return (sm, lin, quad) values and their MC stat variances."""
    if mode == "morphing":
        sm_key  = "histo_sm"
        w1_key  = f"histo_w1_{op}"
        wm1_key = f"histo_wm1_{op}"
        sm   = get_vals(f, sm_key);   var_sm  = get_variances(f, sm_key)
        w1   = get_vals(f, w1_key);   var_w1  = get_variances(f, w1_key)
        wm1  = get_vals(f, wm1_key);  var_wm1 = get_variances(f, wm1_key)
        lin  = 0.5 * (w1 - wm1)
        quad = 0.5 * (w1 + wm1) - sm
        # error propagation (independent MC samples)
        var_lin  = 0.25 * (var_w1 + var_wm1)
        var_quad = 0.25 * (var_w1 + var_wm1) + var_sm
    else:
        sm_key   = f"{channel}/sm"
        quad_key = f"{channel}/quad_{op}"
        slq_key  = f"{channel}/sm_lin_quad_{op}"
        sm   = get_vals(f, sm_key);   var_sm   = get_variances(f, sm_key)
        quad = get_vals(f, quad_key); var_quad = get_variances(f, quad_key)
        slq  = get_vals(f, slq_key);  var_slq  = get_variances(f, slq_key)
        lin  = slq - sm - quad
        var_lin = var_slq + var_sm + var_quad
    return sm, lin, quad, var_sm, var_lin, var_quad


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


def _stairs(ax, vals, edges, color, label, ls="-", lw=2.0):
    ax.stairs(vals, edges=edges, color=color,
              linewidth=lw, linestyle=ls, label=label, fill=False)


def _band(ax, vals, variances, edges, color):
    """Draw a MC stat uncertainty band (±1σ) as a shaded step region."""
    sigma = np.sqrt(np.abs(variances))
    lo = vals - sigma
    hi = vals + sigma
    x    = np.repeat(edges, 2)[1:-1]
    ax.fill_between(x, np.repeat(lo, 2), np.repeat(hi, 2),
                    color=color, alpha=0.25, linewidth=0, label="MC stat. unc.")


def _ratio_band(rax, vals, variances, edges, color):
    """Bottom panel: draw 1 ± sigma/|nominal| band."""
    sigma = np.sqrt(np.abs(variances))
    safe  = np.where(np.abs(vals) > 0, np.abs(vals), np.nan)
    rel   = sigma / safe
    x     = np.repeat(edges, 2)[1:-1]
    rax.fill_between(x, np.repeat(1 - rel, 2), np.repeat(1 + rel, 2),
                     color=color, alpha=0.35, linewidth=0)
    rax.axhline(1.0, color="black", linewidth=0.8, linestyle="dashed")
    rax.set_ylabel("Stat. unc.")
    rax.set_ylim(0.5, 1.5)


def _make_fig(logx):
    fig, (ax, rax) = plt.subplots(2, 1, sharex=True, **FIG_STYLE)
    fig.subplots_adjust(hspace=0.07)
    if logx:
        ax.set_xscale("log")
        rax.set_xscale("log")
    return fig, ax, rax


def _decorate(ax, rax, xlabel, title, logy=False):
    ax.set_ylabel("Events")
    ax.text(0.97, 0.97, title, transform=ax.transAxes,
            ha="right", va="top", fontsize=20, fontweight="bold")
    ax.legend(loc="upper left", fontsize=16)
    if logy:
        ax.set_yscale("log")
    rax.set_xlabel(xlabel)
    rax.autoscale(axis="x", tight=True)
    hep.cms.label(loc=0, label="Preliminary", data=False, ax=ax)


# ---------------------------------------------------------------------------
# per-operator plotting
# ---------------------------------------------------------------------------

def plot_task(d):
    shapes_file = d["shapes"]
    op          = d["op"]
    outdir      = d["outdir"]
    mode        = d["mode"]
    channel     = d["channel"]
    xlabel      = d["xlabel"]
    logx        = d["logx"]
    c_values    = d["c_values"]

    try:
        f      = uproot.open(shapes_file)
        sm_key = "histo_sm" if mode == "morphing" else f"{channel}/sm"
        edges  = get_edges(f, sm_key)
        sm, lin, quad, var_sm, var_lin, var_quad = decompose(f, op, mode, channel)
        f.close()
    except Exception as e:
        print(f"  [skip] {op}: {e}")
        return

    # ---- figure 1: SM + full prediction(s) --------------------------------
    fig1, ax1, rax1 = _make_fig(logx)
    _stairs(ax1, sm, edges, SM_COLOR, "SM", lw=2.5)
    _band(ax1, sm, var_sm, edges, SM_COLOR)
    _ratio_band(rax1, sm, var_sm, edges, SM_COLOR)
    for cv, col in zip(c_values, [LIN_COLOR] + EXTRA_COLORS):
        full = sm + cv * lin + cv**2 * quad
        var_full = var_sm * (1 - cv**2)**2 + var_lin * (2*cv)**2 + var_quad * (2*cv**2)**2
        label = fr"SM + EFT  ($c={cv}$)" if len(c_values) > 1 else r"SM + EFT  ($c=1$)"
        _stairs(ax1, full, edges, col, label, ls="--", lw=2.0)
        _band(ax1, full, var_full, edges, col)
        _ratio_band(rax1, full, var_full, edges, col)
    _decorate(ax1, rax1, xlabel, op, logy=True)
    _save(fig1, os.path.join(outdir, f"sm_full_{op}"))

    # ---- figure 2: linear term (always linear scale — can be negative) -----
    fig2, ax2, rax2 = _make_fig(logx)
    _stairs(ax2, lin, edges, LIN_COLOR, r"linear term  ($c=1$)", lw=2.5)
    _band(ax2, lin, var_lin, edges, LIN_COLOR)
    _ratio_band(rax2, lin, var_lin, edges, LIN_COLOR)
    ax2.axhline(0, color="black", linewidth=0.8, linestyle="dashed")
    _decorate(ax2, rax2, xlabel, op, logy=False)
    _save(fig2, os.path.join(outdir, f"lin_{op}"))

    # ---- figure 3: quadratic term ------------------------------------------
    fig3, ax3, rax3 = _make_fig(logx)
    _stairs(ax3, quad, edges, QUAD_COLOR, r"quadratic term  ($c=1$)", lw=2.5)
    _band(ax3, quad, var_quad, edges, QUAD_COLOR)
    _ratio_band(rax3, quad, var_quad, edges, QUAD_COLOR)
    _decorate(ax3, rax3, xlabel, op, logy=True)
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
    logx     = VAR_LOGX[variable]
    f.close()

    print(f"Shapes file : {args.shapes}")
    print(f"Pipeline    : {mode}" + (f"  (channel: {channel})" if mode == "lhe" else ""))
    print(f"Variable    : {variable}  →  {xlabel}  (log x: {logx})")
    print(f"Operators   : {args.operators}")
    print(f"c values    : {args.c_values}")
    print(f"Output      : {args.outdir}\n")

    tasks = [
        {
            "shapes":   args.shapes,
            "op":       op,
            "outdir":   args.outdir,
            "mode":     mode,
            "channel":  channel,
            "xlabel":   xlabel,
            "logx":     logx,
            "c_values": args.c_values,
        }
        for op in args.operators
    ]

    with Pool(processes=args.ncores) as pool:
        pool.map(plot_task, tasks)

    print("\nDone.")


if __name__ == "__main__":
    main()
