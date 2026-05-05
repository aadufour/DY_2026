#!/usr/bin/env python3
"""
operator_sensitivity.py

Per-operator EFT sensitivity study from histograms.root.

For each operator produces a 3-panel figure:
  Top    : absolute mll shapes (SM + EFT) with QCD scale and PDF syst bands
  Middle : fractional deviation (EFT/SM − 1) with syst and stat (sqrt-N) bands
  Bottom : S/sqrt(B) per bin, where S = |sm_lin_quad_{op} − sm|

Summary outputs:
  summary_sensitivity_heatmap.pdf  — S/sqrt(B) for all ops × bins
  summary_operator_ranking.pdf     — ranking bar chart, colour = tail/bulk fraction

Usage:
  python operator_sensitivity.py --input histograms.root [--outdir plots/sensitivity/]
                                  [--operators cHDD cHWB ...]
                                  [--channel triple_DY]
                                  [--lumi 59740]
                                  [--tail-threshold 200]
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import uproot

# ── operator list (Warsaw basis, SMEFTsim) ─────────────────────────────────
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

# ── colour palette ─────────────────────────────────────────────────────────
C_SM    = "black"
C_EFT   = "crimson"
C_QCD   = "steelblue"
C_PDF   = "darkorange"
C_STAT  = "forestgreen"


# ══ helpers ════════════════════════════════════════════════════════════════

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


def step_xy(edges, values):
    """Build x/y arrays suitable for ax.step(..., where='post').
    Returns (x, y) where x has n+1 points and y has n+1 points (last value repeated).
    """
    return edges, np.append(values, values[-1])


def make_bin_labels(edges):
    lo = edges[:-1].astype(int)
    hi = edges[1:].astype(int)
    return [f"[{l},{h}]" for l, h in zip(lo, hi)]


def fractional(num, denom, fill=0.0):
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(denom > 0, num / denom, fill)


# ══ per-operator plot ══════════════════════════════════════════════════════

def make_operator_plot(f, op, channel, outdir, lumi_fb, tail_thr):
    """Three-panel sensitivity figure for one operator.
    Returns a summary dict, or None if histograms are missing.
    """

    sm, edges = read_hist(f, channel, "sm")
    eft = get_vals(f, channel, f"sm_lin_quad_{op}")

    if sm is None or eft is None:
        print(f"    [skip] missing sm or sm_lin_quad_{op}")
        return None

    n_bins = len(sm)

    # ── systematics (fall back to SM nominal if absent) ───────────────
    qcd_up = get_vals(f, channel, "sm_qcd_scaleUp")   or sm.copy()
    qcd_dn = get_vals(f, channel, "sm_qcd_scaleDown") or sm.copy()
    pdf_up = get_vals(f, channel, "sm_pdfUp")         or sm.copy()
    pdf_dn = get_vals(f, channel, "sm_pdfDown")       or sm.copy()

    # ── derived quantities ────────────────────────────────────────────
    signal   = eft - sm
    abs_sig  = np.abs(signal)
    stat     = np.sqrt(np.maximum(sm, 0))

    ratio    = fractional(signal, sm)
    qcd_rel  = fractional(qcd_up - qcd_dn, 2 * sm)   # half-band / SM
    pdf_rel  = fractional(pdf_up - pdf_dn, 2 * sm)
    stat_rel = fractional(stat, sm)

    s_sqrtB  = fractional(abs_sig, stat)              # S / sqrt(B)

    # ── figure ────────────────────────────────────────────────────────
    fig, (ax1, ax2, ax3) = plt.subplots(
        3, 1, figsize=(9, 11),
        gridspec_kw={"height_ratios": [3, 2, 2], "hspace": 0.04},
        sharex=False,
    )

    # shared log x-axis (mll)
    x_centers = 0.5 * (edges[:-1] + edges[1:])

    def _step_axes(ax, fill=True):
        """Configure mll x-axis with bin-edge lines."""
        ax.set_xscale("log")
        ax.set_xlim(edges[0], edges[-1])
        for e in edges[1:-1]:
            ax.axvline(e, color="gray", lw=0.4, ls=":")
        ax.set_xticks(x_centers)
        ax.set_xticklabels([])
        if fill:
            ax.minorticks_off()

    def _fill_step(ax, lo, hi, color, alpha, **kw):
        ax.fill_between(edges,
                        np.append(lo, lo[-1]),
                        np.append(hi, hi[-1]),
                        step="post", alpha=alpha, color=color, **kw)

    def _line_step(ax, vals, **kw):
        ax.step(edges, np.append(vals, vals[-1]), where="post", **kw)

    # ── Panel 1: absolute shapes ──────────────────────────────────────
    _step_axes(ax1)
    _fill_step(ax1, qcd_dn, qcd_up, C_QCD, 0.25, label="QCD scale env.")
    _fill_step(ax1, pdf_dn, pdf_up, C_PDF, 0.25, label="PDF band")
    _line_step(ax1, sm,  color=C_SM,  lw=1.8, label="SM")
    _line_step(ax1, eft, color=C_EFT, lw=1.8, ls="--",
               label=rf"SM+lin+quad  ($C_{{{op}}}=1$)")

    ax1.set_yscale("log")
    ax1.set_ylabel("Events / bin", fontsize=11)
    ax1.set_title(
        rf"Operator: {op}   |   $\mathcal{{L}}$ = {lumi_fb:.0f} fb$^{{-1}}$",
        fontsize=12, pad=6,
    )
    ax1.legend(fontsize=9, ncol=2, loc="upper right")
    ax1.grid(axis="y", ls=":", alpha=0.35)
    plt.setp(ax1.get_xticklabels(), visible=False)

    # ── Panel 2: fractional deviation ────────────────────────────────
    _step_axes(ax2)
    ax2.axhline(0, color="black", lw=0.9)
    _fill_step(ax2, -stat_rel,  +stat_rel,  C_STAT, 0.15, label=r"Stat $\sqrt{N}$")
    _fill_step(ax2, -pdf_rel,   +pdf_rel,   C_PDF,  0.25, label="PDF")
    _fill_step(ax2, -qcd_rel,   +qcd_rel,   C_QCD,  0.25, label="QCD scale")
    _line_step(ax2, ratio, color=C_EFT, lw=1.8, label="EFT / SM − 1")

    ax2.set_ylabel("EFT / SM − 1", fontsize=10)
    ax2.legend(fontsize=8, ncol=4, loc="best")
    ax2.grid(axis="y", ls=":", alpha=0.35)
    plt.setp(ax2.get_xticklabels(), visible=False)

    # add zero-crossings annotation if relevant
    crosses = np.where(np.diff(np.sign(ratio)) != 0)[0]
    for c in crosses:
        mid = np.sqrt(x_centers[c] * x_centers[c + 1]) if c + 1 < n_bins else x_centers[c]
        ax2.axvline(mid, color=C_EFT, lw=0.7, ls="--", alpha=0.5)

    # ── Panel 3: S/sqrt(B) ────────────────────────────────────────────
    ax3.set_xscale("log")
    ax3.set_xlim(edges[0], edges[-1])
    for e in edges[1:-1]:
        ax3.axvline(e, color="gray", lw=0.4, ls=":")
    ax3.axhline(1, color="gray", ls="--", lw=0.8, alpha=0.7)

    widths = edges[1:] - edges[:-1]
    bar_colors = [C_EFT if s >= 0 else "navy" for s in signal]
    ax3.bar(x_centers, s_sqrtB, width=0.6 * widths,
            color=bar_colors, alpha=0.75, label="S / √B")

    # overlay total syst band for reference
    tot_syst_rel = np.sqrt(qcd_rel**2 + pdf_rel**2)
    ax3.step(edges, np.append(tot_syst_rel, tot_syst_rel[-1]),
             where="post", color="gray", lw=1.0, ls=":", label="Total syst / SM")

    ax3.set_ylabel("S / √B", fontsize=10)
    ax3.set_xlabel(r"$m_{\ell\ell}$ [GeV]", fontsize=11)
    ax3.set_xticks(x_centers)
    ax3.set_xticklabels(make_bin_labels(edges), rotation=30, ha="right", fontsize=7)
    ax3.legend(fontsize=8, loc="upper left")
    ax3.grid(axis="y", ls=":", alpha=0.35)

    plt.tight_layout()
    outpath = os.path.join(outdir, f"sensitivity_{op}.pdf")
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print(f"    saved: {outpath}")

    # ── summary statistics ─────────────────────────────────────────────
    tail_mask = x_centers >= tail_thr
    tail_sum  = s_sqrtB[tail_mask].sum()
    total_sum = s_sqrtB.sum()
    tail_frac = tail_sum / (total_sum + 1e-12)

    return {
        "op":          op,
        "max_s_sqrtB": float(s_sqrtB.max()),
        "peak_bin":    int(s_sqrtB.argmax()),
        "peak_mll":    float(x_centers[s_sqrtB.argmax()]),
        "tail_frac":   float(tail_frac),
        "s_sqrtB":     s_sqrtB.tolist(),
        "ratio":       ratio.tolist(),
    }


# ══ summary plots ══════════════════════════════════════════════════════════

def make_heatmap(results, edges, outdir):
    """Heatmap: operators (rows) × mll bins (cols), cells = S/√B."""
    ops   = [r["op"] for r in results]
    mat   = np.array([r["s_sqrtB"] for r in results])
    labels = make_bin_labels(edges)

    fig, ax = plt.subplots(figsize=(max(8, 1.1 * len(labels)), max(5, 0.35 * len(ops))))
    vmax = np.percentile(mat[mat > 0], 95) if mat.max() > 0 else 1.0
    im = ax.imshow(mat, aspect="auto", cmap="YlOrRd", vmin=0, vmax=vmax)
    plt.colorbar(im, ax=ax, label="S / √B")

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    ax.set_yticks(range(len(ops)))
    ax.set_yticklabels(ops, fontsize=8)
    ax.set_xlabel(r"$m_{\ell\ell}$ bin [GeV]", fontsize=11)
    ax.set_title(r"S / $\sqrt{B}$ per operator per $m_{\ell\ell}$ bin  ($C=1$)", fontsize=12)

    for i, row in enumerate(mat):
        for j, val in enumerate(row):
            text_col = "white" if val > 0.65 * vmax else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                    fontsize=6, color=text_col)

    plt.tight_layout()
    path = os.path.join(outdir, "summary_sensitivity_heatmap.pdf")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"    saved: {path}")


def make_ranking(results, outdir):
    """Horizontal bar chart ranking operators by max S/√B.
    Colour encodes tail_frac: green=tail-dominated, red=bulk-dominated.
    """
    results_s = sorted(results, key=lambda r: r["max_s_sqrtB"], reverse=True)
    ops   = [r["op"] for r in results_s]
    vals  = [r["max_s_sqrtB"] for r in results_s]
    fracs = [r["tail_frac"] for r in results_s]

    cmap   = plt.cm.RdYlGn
    colors = [cmap(f) for f in fracs]

    fig, ax = plt.subplots(figsize=(9, max(4, 0.38 * len(ops))))
    ax.barh(range(len(ops)), vals, color=colors, alpha=0.85, edgecolor="white", lw=0.5)
    ax.set_yticks(range(len(ops)))
    ax.set_yticklabels(ops, fontsize=9)
    ax.set_xlabel(r"max  S / $\sqrt{B}$  across $m_{\ell\ell}$ bins", fontsize=11)
    ax.set_title(
        "Operator ranking by EFT sensitivity  ($C=1$)\n"
        r"colour: green = tail-dominated ($m_{\ell\ell} \gg M_Z$), "
        r"red = Z-peak dominated",
        fontsize=10,
    )
    ax.axvline(1, color="gray", ls="--", lw=0.8)
    ax.invert_yaxis()

    sm_p = mpatches.Patch(color=cmap(0.05), label=r"bulk ($m_{\ell\ell} \sim M_Z$)")
    tl_p = mpatches.Patch(color=cmap(0.95), label=r"tail ($m_{\ell\ell} \gg M_Z$)")
    ax.legend(handles=[sm_p, tl_p], fontsize=9, loc="lower right")
    ax.grid(axis="x", ls=":", alpha=0.4)

    plt.tight_layout()
    path = os.path.join(outdir, "summary_operator_ranking.pdf")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"    saved: {path}")


def make_ratio_heatmap(results, edges, outdir):
    """Heatmap of EFT/SM − 1 (signed) to show where operators deviate."""
    ops    = [r["op"] for r in results]
    mat    = np.array([r["ratio"] for r in results])
    labels = make_bin_labels(edges)

    vext = np.percentile(np.abs(mat), 98)
    fig, ax = plt.subplots(figsize=(max(8, 1.1 * len(labels)), max(5, 0.35 * len(ops))))
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-vext, vmax=vext)
    plt.colorbar(im, ax=ax, label="EFT / SM − 1")

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    ax.set_yticks(range(len(ops)))
    ax.set_yticklabels(ops, fontsize=8)
    ax.set_xlabel(r"$m_{\ell\ell}$ bin [GeV]", fontsize=11)
    ax.set_title(r"EFT / SM − 1 per operator per $m_{\ell\ell}$ bin  ($C=1$)", fontsize=12)

    for i, row in enumerate(mat):
        for j, val in enumerate(row):
            ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                    fontsize=5.5, color="black")

    plt.tight_layout()
    path = os.path.join(outdir, "summary_ratio_heatmap.pdf")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"    saved: {path}")


# ══ main ═══════════════════════════════════════════════════════════════════

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
                # lightweight path: only collect summary stats
                sm  = get_vals(f, args.channel, "sm")
                eft = get_vals(f, args.channel, f"sm_lin_quad_{op}")
                if sm is None or eft is None:
                    print("    [skip]"); continue
                signal  = eft - sm
                stat    = np.sqrt(np.maximum(sm, 0))
                s_sqrtB = fractional(np.abs(signal), stat)
                ratio   = fractional(signal, sm)
                centers = 0.5 * (edges[:-1] + edges[1:])
                tail_mask = centers >= args.tail_threshold
                tail_frac = s_sqrtB[tail_mask].sum() / (s_sqrtB.sum() + 1e-12)
                results.append({
                    "op": op, "max_s_sqrtB": float(s_sqrtB.max()),
                    "peak_bin": int(s_sqrtB.argmax()),
                    "peak_mll": float(centers[s_sqrtB.argmax()]),
                    "tail_frac": float(tail_frac),
                    "s_sqrtB": s_sqrtB.tolist(),
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

    # ── text summary ─────────────────────────────────────────────────
    print("\n{:<12}  {:>10}  {:>10}  {:>12}  {:>10}".format(
        "Operator", "max S/√B", "peak mll", "tail_frac", "dominance"))
    print("-" * 60)
    for r in sorted(results, key=lambda x: x["max_s_sqrtB"], reverse=True):
        dom = "tail" if r["tail_frac"] > 0.5 else "bulk"
        print("{:<12}  {:>10.3f}  {:>9.0f}  {:>11.1%}  {:>10}".format(
            r["op"], r["max_s_sqrtB"], r["peak_mll"], r["tail_frac"], dom))

    print("\nDone.")


if __name__ == "__main__":
    main()
