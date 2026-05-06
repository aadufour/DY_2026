#!/usr/bin/env python3
"""
operator_sensitivity.py

Per-operator EFT sensitivity study from histograms.root.

For each operator produces a 3-panel figure:
  Top    : absolute mll shapes (SM + EFT) with QCD scale and PDF syst bands
  Middle : fractional deviation (EFT/SM - 1) with syst and stat (sqrt-N) bands
  Bottom : S/sqrt(B) per bin, where S = |sm_lin_quad_{op} - sm|, B = sm

Summary outputs:
  summary_sensitivity_heatmap.pdf  — S/sqrt(B) for all ops x bins
  summary_ratio_heatmap.pdf        — signed EFT/SM-1 for all ops x bins
  summary_operator_ranking.pdf     — ranking bar chart, colour = tail/bulk fraction

Usage:
  python operator_sensitivity.py --input histograms.root [--outdir plots/sensitivity/]
                                  [--operators cHDD cHWB ...]
                                  [--channel triple_DY]
                                  [--lumi 59740]
                                  [--tail-threshold 200]
                                  [--no-per-op]
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import uproot

# -- operator list (Warsaw basis, SMEFTsim) ---------------------------------
ALL_OPS = [
    "cHDD", "cHWB", "cbWRe", "cbBRe",
    "cHj1", "cHQ1", "cHj3", "cHQ3",
    "cHu", "cHd", "cHbq",
    "cHl1", "cHl3", "cHe",
    "cll1", "clj1", "clj3",
    "cQl1", "cQl3",
    "ceu", "ced", "cbe", "cje", "cQe",
    "clu", "cld", "cbl",
]

# -- colour palette ---------------------------------------------------------
C_SM   = "black"
C_EFT  = "crimson"
C_QCD  = "steelblue"
C_PDF  = "darkorange"
C_STAT = "forestgreen"


# -- helpers ----------------------------------------------------------------

def read_hist(f, channel, name):
    """Return (values, edges) from a ROOT histogram, or (None, None) if missing."""
    key = f"{channel}/{name}"
    if key not in f:
        return None, None
    vals, edges = f[key].to_numpy()
    return vals.copy(), edges.copy()


def get_vals(f, channel, name):
    v, _ = read_hist(f, channel, name)
    return v


def make_bin_labels(edges):
    lo = edges[:-1].astype(int)
    hi = edges[1:].astype(int)
    return [f"[{l},{h}]" for l, h in zip(lo, hi)]


def fractional(num, denom, fill=0.0):
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(denom > 0, num / denom, fill)


def fill_step(ax, edges, lo, hi, **kw):
    """Step-function fill_between using histogram bin edges."""
    ax.fill_between(
        edges,
        np.append(lo, lo[-1]),
        np.append(hi, hi[-1]),
        step="post", **kw,
    )


# -- per-operator plot ------------------------------------------------------

def make_operator_plot(f, op, channel, outdir, lumi_fb, tail_thr):
    """Three-panel sensitivity figure for one operator.
    Returns a summary dict, or None if histograms are missing.
    """

    sm, edges = read_hist(f, channel, "sm")
    eft = get_vals(f, channel, f"sm_lin_quad_{op}")

    if sm is None or eft is None:
        print(f"    [skip] missing sm or sm_lin_quad_{op}")
        return None

    # -- systematics (fall back to SM nominal if absent) ---------------
    def _syst(name): v = get_vals(f, channel, name); return v if v is not None else sm.copy()
    qcd_up = _syst("sm_qcd_scaleUp")
    qcd_dn = _syst("sm_qcd_scaleDown")
    pdf_up = _syst("sm_pdfUp")
    pdf_dn = _syst("sm_pdfDown")

    # -- derived quantities --------------------------------------------
    signal  = eft - sm
    stat    = np.sqrt(np.maximum(sm, 0))

    ratio    = fractional(signal, sm)
    # asymmetric syst bands: measure actual up/down deviation from SM nominal
    qcd_up_rel = fractional(qcd_up - sm, sm)
    qcd_dn_rel = fractional(qcd_dn - sm, sm)
    pdf_up_rel = fractional(pdf_up - sm, sm)
    pdf_dn_rel = fractional(pdf_dn - sm, sm)
    stat_rel   = fractional(stat, sm)
    # total syst: quadrature of the larger deviation per source per bin
    tot_syst_rel = np.sqrt(
        np.maximum(np.abs(qcd_up_rel), np.abs(qcd_dn_rel))**2 +
        np.maximum(np.abs(pdf_up_rel), np.abs(pdf_dn_rel))**2
    )

    pull = fractional(np.abs(signal), stat)   # |EFT - SM| / sqrt(SM)
    pull_pos = np.where(signal >= 0, pull, 0.)
    pull_neg = np.where(signal <  0, pull, 0.)

    # -- figure: 3 panels sharing the mll x-axis -----------------------
    fig, (ax1, ax2, ax3) = plt.subplots(
        3, 1, figsize=(9, 11),
        gridspec_kw={"height_ratios": [4, 2, 2], "hspace": 0.05},
        sharex=True,
    )

    # -- Panel 1: absolute shapes --------------------------------------
    fill_step(ax1, edges, qcd_dn, qcd_up, alpha=0.25, color=C_QCD, label="QCD scale")
    fill_step(ax1, edges, pdf_dn, pdf_up, alpha=0.25, color=C_PDF, label="PDF")
    hep.histplot(sm,  bins=edges, ax=ax1, color=C_SM,  linewidth=1.8,
                 histtype="step", label="SM")
    hep.histplot(eft, bins=edges, ax=ax1, color=C_EFT, linewidth=1.8,
                 histtype="step", linestyle="--",
                 label=rf"SM+lin+quad  ({op}=1.0)")

    ax1.semilogy()
    ax1.set_xscale("log")
    ax1.set_xlim(edges[0], edges[-1])
    ax1.set_ylabel("Events / bin")
    ax1.legend(frameon=False, fontsize=9, ncol=2, loc="upper right")
    hep.cms.label(f"{op} = 1.0", ax=ax1, data=True, lumi=lumi_fb, loc=0)

    # -- Panel 2: fractional deviation (asymmetric syst bands) --------
    ax2.axhline(0, color="black", linewidth=0.9)
    fill_step(ax2, edges, qcd_dn_rel, qcd_up_rel, alpha=0.25, color=C_QCD, label="QCD scale")
    fill_step(ax2, edges, pdf_dn_rel, pdf_up_rel, alpha=0.25, color=C_PDF,  label="PDF")
    fill_step(ax2, edges, -stat_rel,  +stat_rel,  alpha=0.15, color=C_STAT, label=r"Stat $\sqrt{N}$")
    hep.histplot(ratio, bins=edges, ax=ax2, color=C_EFT, linewidth=1.8,
                 histtype="step", label="EFT / SM - 1")

    ax2.set_ylabel("EFT / SM - 1")
    ax2.legend(frameon=False, fontsize=8, ncol=4, loc="best")

    # -- Panel 3: |EFT - SM| / sqrt(SM) -------------------------------
    ax3.axhline(1, color="gray", linestyle="--", linewidth=0.8, alpha=0.7)
    hep.histplot(pull_pos, bins=edges, ax=ax3, color=C_EFT, alpha=0.75,
                 histtype="fill", label="EFT > SM")
    hep.histplot(pull_neg, bins=edges, ax=ax3, color="navy", alpha=0.75,
                 histtype="fill", label="EFT < SM")

    ax3.set_ylabel(r"$|\mathrm{EFT}-\mathrm{SM}|\,/\,\sqrt{\mathrm{SM}}$")
    ax3.set_xlabel(r"$m_{\ell\ell}$ [GeV]")
    ax3.set_xticks(edges)
    ax3.set_xticklabels([str(int(e)) for e in edges], rotation=45, ha="right", fontsize=7)
    ax3.legend(frameon=False, fontsize=8, loc="upper left")

    # bin-edge separators on all panels
    for e in edges[1:-1]:
        for ax in (ax1, ax2, ax3):
            ax.axvline(e, color="gray", linewidth=0.4, linestyle=":", alpha=0.6)

    plt.tight_layout()
    outpath = os.path.join(outdir, f"sensitivity_{op}.pdf")
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"    saved: {outpath}")

    # -- summary statistics ---------------------------------------------
    x_centers = 0.5 * (edges[:-1] + edges[1:])
    tail_mask = x_centers >= tail_thr
    tail_frac = pull[tail_mask].sum() / (pull.sum() + 1e-12)

    return {
        "op":       op,
        "max_pull": float(pull.max()),
        "peak_bin": int(pull.argmax()),
        "peak_mll": float(x_centers[pull.argmax()]),
        "tail_frac": float(tail_frac),
        "pull":     pull.tolist(),
        "ratio":    ratio.tolist(),
    }


# -- summary plots ----------------------------------------------------------

def make_heatmap(results, edges, outdir):
    """Heatmap: operators (rows) x mll bins (cols), cells = |EFT-SM|/sqrt(SM)."""
    ops    = [r["op"] for r in results]
    mat    = np.array([r["pull"] for r in results])
    labels = make_bin_labels(edges)

    fig, ax = plt.subplots(figsize=(max(8, 1.1 * len(labels)), max(5, 0.35 * len(ops))))
    vmax = np.percentile(mat[mat > 0], 95) if mat.max() > 0 else 1.0
    im = ax.imshow(mat, aspect="auto", cmap="YlOrRd", vmin=0, vmax=vmax)
    cbar = plt.colorbar(im, ax=ax)
    cbar.ax.tick_params(labelsize=7)
    cbar.set_label(r"$|\mathrm{EFT}-\mathrm{SM}|\,/\,\sqrt{\mathrm{SM}}$", fontsize=9)

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    ax.set_yticks(range(len(ops)))
    ax.set_yticklabels(ops, fontsize=8)
    ax.set_xlabel(r"$m_{\ell\ell}$ bin [GeV]", fontsize=9)

    for i, row in enumerate(mat):
        for j, val in enumerate(row):
            text_col = "white" if val > 0.65 * vmax else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                    fontsize=6, color=text_col)

    plt.tight_layout()
    path = os.path.join(outdir, "summary_sensitivity_heatmap.pdf")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"    saved: {path}")


def make_ratio_heatmap(results, edges, outdir):
    """Heatmap of EFT/SM - 1 (signed) — shows where operators deviate and in which direction."""
    ops    = [r["op"] for r in results]
    mat    = np.array([r["ratio"] for r in results])
    labels = make_bin_labels(edges)

    vext = max(np.percentile(np.abs(mat), 98), 1e-3)
    fig, ax = plt.subplots(figsize=(max(8, 1.1 * len(labels)), max(5, 0.35 * len(ops))))
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-vext, vmax=vext)
    plt.colorbar(im, ax=ax, label="EFT / SM - 1")

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    ax.set_yticks(range(len(ops)))
    ax.set_yticklabels(ops, fontsize=8)
    ax.set_xlabel(r"$m_{\ell\ell}$ bin [GeV]")
    ax.set_title(r"EFT / SM - 1 per operator per $m_{\ell\ell}$ bin  ($C=1$)", fontsize=12)

    for i, row in enumerate(mat):
        for j, val in enumerate(row):
            ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                    fontsize=5.5, color="black")

    plt.tight_layout()
    path = os.path.join(outdir, "summary_ratio_heatmap.pdf")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"    saved: {path}")


def make_ranking(results, outdir):
    """Horizontal bar chart ranking operators by max |EFT-SM|/sqrt(SM).
    Colour encodes tail_frac: green = tail-dominated, red = bulk-dominated.
    """
    results_s = sorted(results, key=lambda r: r["max_pull"], reverse=True)
    ops   = [r["op"] for r in results_s]
    vals  = [r["max_pull"] for r in results_s]
    fracs = [r["tail_frac"] for r in results_s]

    cmap   = plt.cm.RdYlGn
    colors = [cmap(f) for f in fracs]

    fig, ax = plt.subplots(figsize=(9, max(4, 0.38 * len(ops))))
    ax.barh(range(len(ops)), vals, color=colors, alpha=0.85, edgecolor="white", linewidth=0.5)
    ax.set_yticks(range(len(ops)))
    ax.set_yticklabels(ops, fontsize=9)
    ax.set_xscale("log")
    ax.set_xlabel(r"max $|\mathrm{EFT}-\mathrm{SM}|\,/\,\sqrt{\mathrm{SM}}$ across $m_{\ell\ell}$ bins")
    ax.set_title(
        r"Operator ranking by EFT sensitivity  ($C=1$)"
        "\n"
        r"colour: green = tail-dominated ($m_{\ell\ell} \gg M_Z$), "
        r"red = Z-peak dominated",
        fontsize=10,
    )
    ax.axvline(1, color="gray", linestyle="--", linewidth=0.8)
    ax.invert_yaxis()
    ax.grid(axis="x", linestyle=":", alpha=0.4)

    sm_p = mpatches.Patch(color=cmap(0.05), label=r"bulk ($m_{\ell\ell} \sim M_Z$)")
    tl_p = mpatches.Patch(color=cmap(0.95), label=r"tail ($m_{\ell\ell} \gg M_Z$)")
    ax.legend(handles=[sm_p, tl_p], frameon=False, fontsize=9, loc="lower right")

    plt.tight_layout()
    path = os.path.join(outdir, "summary_operator_ranking.pdf")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"    saved: {path}")


# -- main -------------------------------------------------------------------

def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input",    required=True,
                    help="Path to histograms.root (built with --all_op --lumi 59740 --pdf-flavour 5)")
    ap.add_argument("--outdir",   default="plots/sensitivity/",
                    help="Output directory for PDF figures")
    ap.add_argument("--channel",  default="triple_DY")
    ap.add_argument("--operators", nargs="+", default=None,
                    help="Subset of operators to process (default: all 27)")
    ap.add_argument("--lumi",     type=float, default=59740.,
                    help="Luminosity in pb^-1 used when building the histograms")
    ap.add_argument("--tail-threshold", type=float, default=200.,
                    help="mll [GeV] above which a bin is considered 'tail' (default: 200)")
    ap.add_argument("--no-per-op", action="store_true",
                    help="Skip per-operator plots; only produce summary figures")
    return ap.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    hep.style.use("CMS")

    lumi_fb = args.lumi / 1000.
    ops     = args.operators or ALL_OPS

    print(f"Input   : {args.input}")
    print(f"Channel : {args.channel}")
    print(f"Lumi    : {lumi_fb:.2f} fb⁻¹")
    print(f"Ops     : {len(ops)}")
    print(f"Outdir  : {args.outdir}\n")

    with uproot.open(args.input) as f:

        # read bin edges once from SM histogram
        _, edges = read_hist(f, args.channel, "sm")
        if edges is None:
            raise KeyError(f"Histogram {args.channel}/sm not found in {args.input}")
        print(f"Bins    : {len(edges)-1}  edges = {edges.tolist()}\n")

        results = []
        for op in ops:
            print(f"  {op}")
            if args.no_per_op:
                sm  = get_vals(f, args.channel, "sm")
                eft = get_vals(f, args.channel, f"sm_lin_quad_{op}")
                if sm is None or eft is None:
                    print("    [skip]"); continue
                signal  = eft - sm
                stat    = np.sqrt(np.maximum(sm, 0))
                pull    = fractional(np.abs(signal), stat)
                ratio   = fractional(signal, sm)
                centers = 0.5 * (edges[:-1] + edges[1:])
                tail_frac = pull[centers >= args.tail_threshold].sum() / (pull.sum() + 1e-12)
                results.append({
                    "op": op, "max_pull": float(pull.max()),
                    "peak_bin": int(pull.argmax()),
                    "peak_mll": float(centers[pull.argmax()]),
                    "tail_frac": float(tail_frac),
                    "pull": pull.tolist(),
                    "ratio": ratio.tolist(),
                })
            else:
                res = make_operator_plot(
                    f, op, args.channel, args.outdir, lumi_fb, args.tail_threshold,
                )
                if res is not None:
                    results.append(res)

    if not results:
        print("No results — check that histograms.root contains sm_lin_quad_* entries.")
        return

    print(f"\nBuilding summary plots ({len(results)} operators) ...")
    make_heatmap(results, edges, args.outdir)
    make_ratio_heatmap(results, edges, args.outdir)
    make_ranking(results, args.outdir)

    # -- text summary -------------------------------------------------
    print("\n{:<12}  {:>12}  {:>10}  {:>12}  {:>10}".format(
        "Operator", "max |Δ|/√σSM", "peak mll", "tail_frac", "dominance"))
    print("-" * 62)
    for r in sorted(results, key=lambda x: x["max_pull"], reverse=True):
        dom = "tail" if r["tail_frac"] > 0.5 else "bulk"
        print("{:<12}  {:>12.3f}  {:>9.0f}  {:>11.1%}  {:>10}".format(
            r["op"], r["max_pull"], r["peak_mll"], r["tail_frac"], dom))

    print("\nDone.")


if __name__ == "__main__":
    main()
