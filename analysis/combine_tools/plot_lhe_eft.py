#!/usr/bin/env python3
"""
plot_lhe_eft.py
===============
Plot EFT decomposition (SM, lin, quad) from the LHE cache pkl.

Cache format (from build_cache_propcorr.py):
    cache['mll']       — array of mll values (GeV)
    cache['rap']       — array of |y_ll|
    cache['cstar']     — array of cos(theta*)_CS
    cache['w_SM']      — array of SM weights
    cache['w_p1'][op]  — array of weights at c_{op} = +1
    cache['w_m1'][op]  — array of weights at c_{op} = -1

Morphing:
    lin  = 0.5 * (w_p1 - w_m1)
    quad = 0.5 * (w_p1 + w_m1) - w_SM

Three figures per operator per variable:
    {var}_sm_full_{op}  — SM + SM+EFT at c=1 (log y)
    {var}_lin_{op}      — linear term at c=1  (linear y, can be negative)
    {var}_quad_{op}     — quadratic term at c=1 (log y)

Plus for triple_diff: an additional reconstructed-canvas figure.

Usage:
    python plot_lhe_eft.py --cache lhe_cache_propcorr_parallel.pkl --outdir plots/lhe_eft
    python plot_lhe_eft.py --cache /path/to/cache.pkl --operators cHDD cHWB --c-values 0.5 1.0
"""

import argparse
import os
import pickle
import sys

import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np

hep.style.use("CMS")

# ---------------------------------------------------------------------------
# Binning
# 1D variables: fine binning from config.py (good for thesis plots)
# triple_diff:  coarse binning from config_v9.py (14×5×4 = 280 bins)
# ---------------------------------------------------------------------------

# fine 1D binning
MLL_EDGES   = np.array([
    *range(50, 76, 5),    # 50–75:  5 GeV steps
    *range(76, 106, 2),   # 76–105: 2 GeV steps (Z peak)
    *range(106, 120, 5),  # 106–119: 5 GeV steps
    120, 150, 200, 250, 300, 400, 600, 800, 1000, 1500, 3000,
], dtype=float)
COSTH_EDGES = np.linspace(-1.0, 1.0, 51)   # 50 uniform bins
RAP_EDGES   = np.linspace(0.0,  2.5, 51)   # 50 uniform bins

# coarse triple_diff binning
TD_MLL_EDGES   = np.array([40, 60, 80, 100, 120, 140, 180, 220, 270, 350, 500, 700, 1000, 1500, 3000], dtype=float)
TD_COSTH_EDGES = np.array([-1.0, -0.6, -0.2, 0.2, 0.6, 1.0])
TD_RAP_EDGES   = np.array([0.0, 0.48, 0.96, 1.44, 2.4])

N_TD_MLL   = len(TD_MLL_EDGES)   - 1  # 14
N_TD_COSTH = len(TD_COSTH_EDGES) - 1  # 5
N_TD_RAP   = len(TD_RAP_EDGES)   - 1  # 4
N_TD       = N_TD_MLL * N_TD_COSTH * N_TD_RAP  # 280

OPERATORS = [
    "cHDD", "cHWB", "cbWRe", "cbBRe", "cHj1", "cHQ1", "cHj3", "cHQ3",
    "cHu", "cHd", "cHbq", "cHl1", "cHl3", "cHe", "cll1", "clj1", "clj3",
    "cQl1", "cQl3", "ceu", "ced", "cbe", "cje", "cQe", "clu", "cld", "cbl",
]

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------
SM_COLOR     = "#5790fc"
LIN_COLOR    = "#f89c20"
QUAD_COLOR   = "#e42536"
EXTRA_COLORS = ["#9467bd", "#8c564b", "#17becf"]

VAR_META = {
    "mll":          {"xlabel": r"$m_{\ell\ell}$ [GeV]",  "logx": True,  "unit": True},
    "costhetastar": {"xlabel": r"$\cos\theta^*$",         "logx": False, "unit": False},
    "rapll_abs":    {"xlabel": r"$|y_{\ell\ell}|$",       "logx": False, "unit": False},
    "triple_diff":  {"xlabel": "Unrolled bin",            "logx": False, "unit": False},
}

FIG_STYLE = {"figsize": (10, 10), "gridspec_kw": {"height_ratios": (3, 1)}}


# ---------------------------------------------------------------------------
# Histogram filling from cache arrays
# ---------------------------------------------------------------------------

def fill_histograms(cache, op):
    """Return dict of {var: (sm, lin, quad, sm_v, lin_v, quad_v, edges)}."""
    mll   = cache["mll"]
    rap   = cache["rap"]
    costh = cache["cstar"]
    w_sm  = cache["w_SM"]
    w_p1  = cache["w_p1"][op]
    w_m1  = cache["w_m1"][op]

    w_lin  = 0.5 * (w_p1 - w_m1)
    w_quad = 0.5 * (w_p1 + w_m1) - w_sm

    def fill1d(vals, edges):
        sm_h,   _ = np.histogram(vals, bins=edges, weights=w_sm)
        lin_h,  _ = np.histogram(vals, bins=edges, weights=w_lin)
        quad_h, _ = np.histogram(vals, bins=edges, weights=w_quad)
        sm_v,   _ = np.histogram(vals, bins=edges, weights=w_sm**2)
        lin_v,  _ = np.histogram(vals, bins=edges, weights=w_lin**2)
        quad_v, _ = np.histogram(vals, bins=edges, weights=w_quad**2)
        return (sm_h, lin_h, quad_h, sm_v, lin_v, quad_v)

    out = {}
    out["mll"]          = fill1d(mll,   MLL_EDGES)   + (MLL_EDGES,)
    out["costhetastar"] = fill1d(costh, COSTH_EDGES) + (COSTH_EDGES,)
    out["rapll_abs"]    = fill1d(rap,   RAP_EDGES)   + (RAP_EDGES,)

    # triple_diff uses its own coarser binning (14×5×4 = 280 bins)
    irap  = np.clip(np.searchsorted(TD_RAP_EDGES,   rap,   side="right") - 1, 0, N_TD_RAP   - 1)
    icos  = np.clip(np.searchsorted(TD_COSTH_EDGES, costh, side="right") - 1, 0, N_TD_COSTH - 1)
    imll  = np.clip(np.searchsorted(TD_MLL_EDGES,   mll,   side="right") - 1, 0, N_TD_MLL   - 1)
    flat  = (irap * N_TD_COSTH * N_TD_MLL + icos * N_TD_MLL + imll).astype(int)

    def fill_td(w):
        h = np.zeros(N_TD)
        np.add.at(h, flat, w)
        return h

    td_edges = np.arange(N_TD + 1, dtype=float)
    out["triple_diff"] = (
        fill_td(w_sm),   fill_td(w_lin),   fill_td(w_quad),
        fill_td(w_sm**2), fill_td(w_lin**2), fill_td(w_quad**2),
        td_edges,
    )
    return out


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def _make_fig(logx):
    fig, (ax, rax) = plt.subplots(2, 1, sharex=True, **FIG_STYLE)
    fig.subplots_adjust(hspace=0.07)
    if logx:
        ax.set_xscale("log"); rax.set_xscale("log")
    return fig, ax, rax


def _stairs(ax, vals, edges, color, label, ls="-", lw=2.0):
    ax.stairs(vals, edges=edges, color=color, linewidth=lw, linestyle=ls,
              label=label, fill=False)


def _band(ax, vals, variances, edges, color):
    sigma = np.sqrt(np.abs(variances))
    x = np.repeat(edges, 2)[1:-1]
    ax.fill_between(x, np.repeat(vals - sigma, 2), np.repeat(vals + sigma, 2),
                    color=color, alpha=0.25, linewidth=0)


def _ratio_band(rax, vals, variances, edges, color):
    """MC stat uncertainty band around 1 (for lin/quad panels)."""
    sigma = np.sqrt(np.abs(variances))
    safe  = np.where(np.abs(vals) > 0, np.abs(vals), np.nan)
    rel   = sigma / safe
    x     = np.repeat(edges, 2)[1:-1]
    rax.fill_between(x, np.repeat(1 - rel, 2), np.repeat(1 + rel, 2),
                     color=color, alpha=0.35, linewidth=0)
    rax.axhline(1.0, color="black", linewidth=0.8, linestyle="dashed")
    rax.set_ylabel("MC stat.")
    rax.set_ylim(0.5, 1.5)


def _eft_ratio_panel(rax, sm, sm_v, full, full_v, edges, c_values, colors):
    """
    Bottom panel for sm_full figure: EFT/SM ratio with MC stat unc bands.
    One curve per c value. SM unc band shown in grey around 1.
    """
    x    = np.repeat(edges, 2)[1:-1]
    safe = np.where(sm > 0, sm, np.nan)

    # SM uncertainty band (grey, around 1)
    sm_rel = np.sqrt(np.abs(sm_v)) / safe
    rax.fill_between(x, np.repeat(1 - sm_rel, 2), np.repeat(1 + sm_rel, 2),
                     color="grey", alpha=0.3, linewidth=0, label="SM MC stat.")

    # EFT/SM ratio curve + its uncertainty band per c value
    for (cv, col, full_cv, full_v_cv) in zip(c_values, colors, full, full_v):
        ratio     = full_cv / safe
        ratio_err = np.sqrt(np.abs(full_v_cv)) / safe
        rax.stairs(ratio, edges=edges, color=col, linewidth=1.8, linestyle="--",
                   label=fr"SM+EFT ($c={cv}$)" if len(c_values) > 1 else r"SM+EFT ($c=1$)")
        rax.fill_between(x,
                         np.repeat(ratio - ratio_err, 2),
                         np.repeat(ratio + ratio_err, 2),
                         color=col, alpha=0.2, linewidth=0)

    rax.axhline(1.0, color="black", linewidth=0.8, linestyle="dashed")
    rax.set_ylabel("EFT / SM", fontsize=12)
    rax.legend(loc="upper left", fontsize=9)
    # auto-range with a sensible cap
    rax.autoscale(axis="y")
    lo, hi = rax.get_ylim()
    rax.set_ylim(max(lo, 0.5), min(hi, 3.0))


def _decorate(ax, rax, ylabel, xlabel, op, logy=False):
    ax.set_ylabel(ylabel, fontsize=14)
    ax.text(0.97, 0.97, op, transform=ax.transAxes,
            ha="right", va="top", fontsize=20, fontweight="bold")
    ax.legend(loc="upper left", fontsize=12)
    if logy:
        ax.set_yscale("log")
    rax.set_xlabel(xlabel)
    rax.autoscale(axis="x", tight=True)
    hep.cms.label(loc=0, label="Preliminary", data=False, ax=ax)


def _save(fig, stem):
    for ext in ("png", "pdf"):
        fig.savefig(f"{stem}.{ext}", bbox_inches="tight", facecolor="white")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Per-variable 1D plot
# ---------------------------------------------------------------------------

def plot_variable(var, meta, histos, op, c_values, outdir):
    sm, lin, quad, sm_v, lin_v, quad_v, edges = histos[var]

    # divide by bin width for Events/GeV axes
    if meta["unit"]:
        w = np.diff(edges)
        sm, lin, quad    = sm/w,   lin/w,   quad/w
        sm_v, lin_v, quad_v = sm_v/w**2, lin_v/w**2, quad_v/w**2

    logx  = meta["logx"]
    ylab  = "Events / GeV" if meta["unit"] else "Events"
    xlab  = meta["xlabel"]
    stem  = os.path.join(outdir, f"{var}_sm_full_{op}")

    # SM + full (bottom panel: EFT/SM ratio)
    fig, ax, rax = _make_fig(logx)
    _stairs(ax, sm, edges, SM_COLOR, "SM", lw=2.5)
    _band(ax, sm, sm_v, edges, SM_COLOR)
    eft_colors  = [LIN_COLOR] + EXTRA_COLORS
    full_list   = []
    full_v_list = []
    for cv, col in zip(c_values, eft_colors):
        full_cv   = sm + cv * lin + cv**2 * quad
        full_v_cv = sm_v + cv**2 * lin_v + cv**4 * quad_v
        lbl       = fr"SM+EFT ($c={cv}$)" if len(c_values) > 1 else r"SM+EFT ($c=1$)"
        _stairs(ax, full_cv, edges, col, lbl, ls="--", lw=2.0)
        _band(ax, full_cv, full_v_cv, edges, col)
        full_list.append(full_cv)
        full_v_list.append(full_v_cv)
    _eft_ratio_panel(rax, sm, sm_v, full_list, full_v_list, edges,
                     c_values=c_values, colors=eft_colors[:len(c_values)])
    _decorate(ax, rax, ylab, xlab, op, logy=True)
    _save(fig, stem)

    # linear term
    fig, ax, rax = _make_fig(logx)
    _stairs(ax, lin, edges, LIN_COLOR, r"linear ($c=1$)", lw=2.5)
    _band(ax, lin, lin_v, edges, LIN_COLOR)
    _ratio_band(rax, lin, lin_v, edges, LIN_COLOR)
    ax.axhline(0, color="black", linewidth=0.8, linestyle="dashed")
    _decorate(ax, rax, ylab, xlab, op, logy=False)
    _save(fig, os.path.join(outdir, f"{var}_lin_{op}"))

    # quadratic term
    fig, ax, rax = _make_fig(logx)
    _stairs(ax, quad, edges, QUAD_COLOR, r"quadratic ($c=1$)", lw=2.5)
    _band(ax, quad, quad_v, edges, QUAD_COLOR)
    _ratio_band(rax, quad, quad_v, edges, QUAD_COLOR)
    _decorate(ax, rax, ylab, xlab, op, logy=True)
    _save(fig, os.path.join(outdir, f"{var}_quad_{op}"))


# ---------------------------------------------------------------------------
# Triple-diff reconstructed canvas
# Grid: N_RAP columns x N_COSTH rows, mll on x-axis (log scale, Events/GeV)
# ---------------------------------------------------------------------------

def plot_triple_diff_2d(histos, op, c_values, outdir):
    sm, lin, quad, sm_v, lin_v, quad_v, _ = histos["triple_diff"]
    widths = np.diff(TD_MLL_EDGES)

    fig, axes = plt.subplots(
        N_TD_COSTH, N_TD_RAP,
        figsize=(4.5 * N_TD_RAP, 3.5 * N_TD_COSTH),
        sharex=True, sharey=False,
    )
    fig.subplots_adjust(hspace=0.08, wspace=0.35)

    handles, labels = [], []

    for irap in range(N_TD_RAP):
        for icos in range(N_TD_COSTH):
            ax = axes[N_TD_COSTH - 1 - icos][irap]  # flip: icos=0 at bottom

            sl = slice(
                irap * N_TD_COSTH * N_TD_MLL + icos * N_TD_MLL,
                irap * N_TD_COSTH * N_TD_MLL + icos * N_TD_MLL + N_TD_MLL,
            )
            sm_s    = sm[sl]   / widths
            lin_s   = lin[sl]  / widths
            quad_s  = quad[sl] / widths
            sm_v_s  = sm_v[sl] / widths**2

            is_legend_cell = (irap == 0 and icos == N_TD_COSTH - 1)

            h1, = ax.step(np.append(TD_MLL_EDGES[:-1], TD_MLL_EDGES[-1]),
                          np.append(sm_s, sm_s[-1]),
                          where="post", color=SM_COLOR, linewidth=1.8,
                          label="SM")
            _band(ax, sm_s, sm_v_s, TD_MLL_EDGES, SM_COLOR)

            for cv, col in zip(c_values, [LIN_COLOR] + EXTRA_COLORS):
                full   = sm_s + cv * lin_s + cv**2 * quad_s
                full_v = sm_v_s + cv**2 * (quad_v[sl] / widths**2)
                lbl    = fr"SM+EFT ($c={cv}$)" if len(c_values) > 1 else r"SM+EFT ($c=1$)"
                h2, = ax.step(np.append(TD_MLL_EDGES[:-1], TD_MLL_EDGES[-1]),
                              np.append(full, full[-1]),
                              where="post", color=col, linewidth=1.5,
                              linestyle="--", label=lbl)
                if is_legend_cell:
                    handles.append(h2); labels.append(lbl)

            if is_legend_cell:
                handles.insert(0, h1); labels.insert(0, "SM")

            ax.set_xscale("log")
            ax.set_yscale("log")

            rap_lo, rap_hi   = TD_RAP_EDGES[irap],   TD_RAP_EDGES[irap + 1]
            cos_lo, cos_hi   = TD_COSTH_EDGES[icos], TD_COSTH_EDGES[icos + 1]
            lbl_text = (fr"$|y|$: [{rap_lo},{rap_hi}]" + "\n" +
                        fr"$\cos\theta^*$: [{cos_lo},{cos_hi}]")
            ax.text(0.97, 0.97, lbl_text, transform=ax.transAxes,
                    ha="right", va="top", fontsize=7)

            if icos == 0:
                ax.set_xlabel(r"$m_{\ell\ell}$ [GeV]", fontsize=9)
            if irap == 0:
                ax.set_ylabel("Events / GeV", fontsize=9)

    fig.legend(handles, labels, loc="upper center", ncol=len(c_values) + 1,
               fontsize=10, bbox_to_anchor=(0.5, 1.02))
    hep.cms.label(loc=0, label="Preliminary", data=False, ax=axes[0][0])
    fig.suptitle(f"Triple-diff — {op}", y=1.05, fontsize=13)

    _save(fig, os.path.join(outdir, f"triple_diff_2d_{op}"))


# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cache",     required=True,
                        help="Path to lhe_cache_propcorr_parallel.pkl")
    parser.add_argument("--outdir",    default="plots/lhe_eft")
    parser.add_argument("--operators", nargs="+", default=OPERATORS)
    parser.add_argument("--c-values",  nargs="+", type=float, default=[1.0])
    parser.add_argument("--variables", nargs="+",
                        default=["mll", "costhetastar", "rapll_abs",
                                 "triple_diff", "triple_diff_2d"],
                        choices=["mll", "costhetastar", "rapll_abs",
                                 "triple_diff", "triple_diff_2d"])
    args = parser.parse_args()

    print(f"Loading cache: {args.cache}")
    with open(args.cache, "rb") as f:
        cache = pickle.load(f)
    n = len(cache["mll"])
    print(f"Events in cache : {n:,}")
    print(f"Operators       : {args.operators}")
    print(f"c values        : {args.c_values}")
    print(f"Variables       : {args.variables}\n")

    os.makedirs(args.outdir, exist_ok=True)

    for op in args.operators:
        if op not in cache["w_p1"]:
            print(f"  [skip] {op}: not in cache")
            continue

        histos = fill_histograms(cache, op)

        for var in args.variables:
            if var == "triple_diff_2d":
                plot_triple_diff_2d(histos, op, args.c_values, args.outdir)
            else:
                plot_variable(var, VAR_META[var], histos, op, args.c_values, args.outdir)

        print(f"  {op:12s}  done")

    print("\nDone.")


if __name__ == "__main__":
    main()
