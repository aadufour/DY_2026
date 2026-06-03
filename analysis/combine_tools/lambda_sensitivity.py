#!/usr/bin/env python3
"""
lambda_sensitivity.py

Per-operator EFT sensitivity expressed as a reach in the EFT cutoff scale Λ.

Physics:  σ = σ_SM + (C/Λ²) σ_lin + (C/Λ²)² σ_quad
Templates are generated at C=1, Λ=Λ_ref (default 1 TeV).
In the linear regime the signal scales as 1/Λ², so the pull at arbitrary Λ is
    pull(Λ) = pull(Λ_ref) × (Λ_ref/Λ)²
Setting pull(Λ*) = 1 gives the operator's Λ reach per bin:
    Λ* [TeV] = Λ_ref × sqrt( |EFT - SM| / sqrt(SM) )

Using the linear template (sm_lin_{op}) is preferred because this 1/Λ²
rescaling is exact; using sm_lin_quad_{op} mixes in a 1/Λ⁴ quadratic
contribution that breaks the scaling (use --quad to opt in if lin is absent).

Outputs (per operator):
  lambda_reach_{op}.pdf   — Λ* per mll bin (bar chart) with syst bands
Summary outputs:
  lambda_heatmap.pdf      — Λ* for all operators × bins
  lambda_ranking.pdf      — ranking bar chart, colour = tail/bulk fraction

Usage:
  python lambda_sensitivity.py --input histograms.root [--outdir plots/lambda/]
                                [--operators cHDD cHWB ...]
                                [--channel triple_DY]
                                [--lumi 59740]
                                [--lambda-ref 1.0]
                                [--tail-threshold 200]
                                [--quad]
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

C_SM   = "black"
C_EFT  = "crimson"
C_QCD  = "steelblue"
C_PDF  = "darkorange"
C_STAT = "forestgreen"


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def read_hist(f, channel, name):
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


def save_both(fig, path):
    fig.savefig(path, dpi=150, bbox_inches="tight")
    fig.savefig(path.replace(".pdf", ".png"), dpi=150, bbox_inches="tight")


def fill_step(ax, edges, lo, hi, **kw):
    ax.fill_between(
        edges,
        np.append(lo, lo[-1]),
        np.append(hi, hi[-1]),
        step="post", **kw,
    )


def lambda_reach(eft, sm, lambda_ref):
    """Per-bin Λ* [TeV] = λ_ref × sqrt( |EFT−SM| / sqrt(SM) ).
    Bins with SM≤0 or |EFT−SM|=0 return 0.
    """
    with np.errstate(invalid="ignore", divide="ignore"):
        pull = np.where(sm > 0, np.abs(eft - sm) / np.sqrt(sm), 0.0)
    return lambda_ref * np.sqrt(np.maximum(pull, 0.0))


# ---------------------------------------------------------------------------
# per-operator plot
# ---------------------------------------------------------------------------

def make_operator_plot(f, op, channel, outdir, lumi_fb, tail_thr, lambda_ref, use_quad):
    sm, edges = read_hist(f, channel, "sm")

    template_name = f"sm_lin_quad_{op}" if use_quad else f"sm_lin_{op}"
    eft = get_vals(f, channel, template_name)

    # fall back to lin+quad if linear template is missing
    if eft is None and not use_quad:
        eft = get_vals(f, channel, f"sm_lin_quad_{op}")
        if eft is not None:
            template_name = f"sm_lin_quad_{op}"
            print(f"    [warn] sm_lin_{op} absent, using sm_lin_quad_{op} (Λ scaling approximate)")

    if sm is None or eft is None:
        print(f"    [skip] missing sm or {template_name}")
        return None

    def _syst(name): v = get_vals(f, channel, name); return v if v is not None else sm.copy()
    qcd_up = _syst("sm_qcd_scaleUp")
    qcd_dn = _syst("sm_qcd_scaleDown")
    pdf_up = _syst("sm_pdfUp")
    pdf_dn = _syst("sm_pdfDown")

    lam     = lambda_reach(eft, sm, lambda_ref)
    lam_qup = lambda_reach(qcd_up, sm, lambda_ref)   # syst shifts Λ* by changing SM prediction
    lam_qdn = lambda_reach(qcd_dn, sm, lambda_ref)
    lam_pup = lambda_reach(pdf_up, sm, lambda_ref)
    lam_pdn = lambda_reach(pdf_dn, sm, lambda_ref)

    centers  = 0.5 * (edges[:-1] + edges[1:])
    x_pos    = np.arange(len(centers))
    bin_labs = make_bin_labels(edges)

    # -- figure ----------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(max(9, 0.7 * len(centers)), 5), layout="constrained")

    bar_w = 0.6
    ax.bar(x_pos, lam, width=bar_w, color=C_EFT, alpha=0.75, label=r"$\Lambda^*$")

    # syst uncertainty on Λ* as error bars (envelope of QCD + PDF shifts)
    lam_err_up = np.maximum(lam_qup - lam, 0) + np.maximum(lam_pup - lam, 0)
    lam_err_dn = np.maximum(lam - lam_qdn, 0) + np.maximum(lam - lam_pdn, 0)
    ax.errorbar(x_pos, lam, yerr=[lam_err_dn, lam_err_up],
                fmt="none", color="gray", linewidth=1.2, capsize=3, label="Syst (QCD+PDF)")

    ax.axhline(lambda_ref, color="black", linestyle="--", linewidth=1.0,
               label=rf"$\Lambda_\mathrm{{ref}} = {lambda_ref:.0f}$ TeV")

    ax.set_xticks(x_pos)
    ax.set_xticklabels(bin_labs, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(r"$\Lambda^*$ [TeV]", fontsize=12)
    ax.set_xlabel(r"$m_{\ell\ell}$ bin [GeV]", fontsize=11)
    ax.legend(frameon=False, fontsize=9)
    ax.grid(axis="y", linestyle=":", alpha=0.4)

    mode_str = "lin+quad (approx)" if "quad" in template_name else "lin"
    hep.cms.label(
        rf"{op} ($C=1$, {mode_str})",
        ax=ax, data=True, lumi=lumi_fb, loc=0,
    )

    outpath = os.path.join(outdir, f"lambda_reach_{op}.pdf")
    save_both(fig, outpath)
    plt.close(fig)
    print(f"    saved: {outpath}")

    tail_mask = centers >= tail_thr
    tail_frac = lam[tail_mask].sum() / (lam.sum() + 1e-12)

    return {
        "op":         op,
        "max_lambda": float(lam.max()),
        "peak_bin":   int(lam.argmax()),
        "peak_mll":   float(centers[lam.argmax()]),
        "tail_frac":  float(tail_frac),
        "lambda":     lam.tolist(),
    }


# ---------------------------------------------------------------------------
# summary plots
# ---------------------------------------------------------------------------

def make_lambda_heatmap(results, edges, outdir, lambda_ref):
    ops    = [r["op"] for r in results]
    mat    = np.array([r["lambda"] for r in results])
    labels = make_bin_labels(edges)

    fig, ax = plt.subplots(figsize=(max(8, 1.1 * len(labels)), max(5, 0.35 * len(ops))))
    vmax = np.percentile(mat[mat > 0], 95) if mat.max() > 0 else 1.0
    im = ax.imshow(mat, aspect="auto", cmap="YlOrRd", vmin=0, vmax=vmax)
    cbar = plt.colorbar(im, ax=ax)
    cbar.ax.tick_params(labelsize=7)
    cbar.set_label(r"$\Lambda^*$ [TeV]", fontsize=9)

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    ax.set_yticks(range(len(ops)))
    ax.set_yticklabels(ops, fontsize=8)
    ax.set_xlabel(r"$m_{\ell\ell}$ bin [GeV]", fontsize=9)
    ax.set_title(
        rf"$\Lambda^*$ reach per operator per $m_{{\ell\ell}}$ bin  "
        rf"($C=1$, $\Lambda_\mathrm{{ref}}={lambda_ref:.0f}$ TeV)",
        fontsize=10,
    )

    for i, row in enumerate(mat):
        for j, val in enumerate(row):
            text_col = "white" if val > 0.65 * vmax else "black"
            ax.text(j, i, f"{val:.1f}", ha="center", va="center",
                    fontsize=6, color=text_col)

    path = os.path.join(outdir, "lambda_heatmap.pdf")
    save_both(fig, path)
    plt.close(fig)
    print(f"    saved: {path}")


def make_lambda_ranking(results, outdir, lambda_ref):
    results_s = sorted(results, key=lambda r: r["max_lambda"], reverse=True)
    ops   = [r["op"] for r in results_s]
    vals  = [r["max_lambda"] for r in results_s]
    fracs = [r["tail_frac"] for r in results_s]

    cmap   = plt.cm.RdYlGn
    colors = [cmap(f) for f in fracs]

    fig, ax = plt.subplots(figsize=(9, max(4, 0.38 * len(ops))))
    ax.barh(range(len(ops)), vals, color=colors, alpha=0.85, edgecolor="white", linewidth=0.5)
    ax.set_yticks(range(len(ops)))
    ax.set_yticklabels(ops, fontsize=9)
    ax.set_xlabel(r"max $\Lambda^*$ [TeV]  (best bin)")
    ax.set_title(
        rf"Operator $\Lambda$ reach ranking  ($C=1$, $\Lambda_\mathrm{{ref}}={lambda_ref:.0f}$ TeV)"
        "\n"
        r"colour: green = tail-dominated ($m_{\ell\ell} \gg M_Z$), "
        r"red = Z-peak dominated",
        fontsize=10,
    )
    ax.axvline(lambda_ref, color="gray", linestyle="--", linewidth=0.8,
               label=rf"$\Lambda_\mathrm{{ref}}={lambda_ref:.0f}$ TeV")
    ax.invert_yaxis()
    ax.grid(axis="x", linestyle=":", alpha=0.4)
    ax.legend(frameon=False, fontsize=9)

    sm_p = mpatches.Patch(color=cmap(0.05), label=r"bulk ($m_{\ell\ell} \sim M_Z$)")
    tl_p = mpatches.Patch(color=cmap(0.95), label=r"tail ($m_{\ell\ell} \gg M_Z$)")
    ax.legend(handles=[sm_p, tl_p], frameon=False, fontsize=9, loc="lower right")

    path = os.path.join(outdir, "lambda_ranking.pdf")
    save_both(fig, path)
    plt.close(fig)
    print(f"    saved: {path}")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input",     required=True)
    ap.add_argument("--outdir",    default="plots/lambda/")
    ap.add_argument("--channel",   default="triple_DY")
    ap.add_argument("--operators", nargs="+", default=None)
    ap.add_argument("--lumi",      type=float, default=59740.)
    ap.add_argument("--lambda-ref", type=float, default=1.0,
                    help="Reference Λ in TeV used when generating templates (default: 1.0)")
    ap.add_argument("--tail-threshold", type=float, default=200.)
    ap.add_argument("--quad", action="store_true",
                    help="Force use of sm_lin_quad template (Λ scaling only approximate)")
    ap.add_argument("--no-per-op", action="store_true")
    return ap.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    hep.style.use("CMS")

    lumi_fb    = args.lumi / 1000.
    ops        = args.operators or ALL_OPS
    lambda_ref = args.lambda_ref

    print(f"Input      : {args.input}")
    print(f"Channel    : {args.channel}")
    print(f"Lumi       : {lumi_fb:.2f} fb⁻¹")
    print(f"Λ_ref      : {lambda_ref:.1f} TeV")
    print(f"Template   : {'sm_lin_quad (forced)' if args.quad else 'sm_lin (preferred)'}")
    print(f"Ops        : {len(ops)}")
    print(f"Outdir     : {args.outdir}\n")

    with uproot.open(args.input) as f:
        _, edges = read_hist(f, args.channel, "sm")
        if edges is None:
            raise KeyError(f"{args.channel}/sm not found in {args.input}")
        print(f"Bins : {len(edges)-1}  edges = {edges.tolist()}\n")

        results = []
        for op in ops:
            print(f"  {op}")
            if args.no_per_op:
                sm = get_vals(f, args.channel, "sm")
                template_name = f"sm_lin_quad_{op}" if args.quad else f"sm_lin_{op}"
                eft = get_vals(f, args.channel, template_name)
                if eft is None and not args.quad:
                    eft = get_vals(f, args.channel, f"sm_lin_quad_{op}")
                    if eft is not None:
                        print(f"    [warn] using sm_lin_quad_{op}")
                if sm is None or eft is None:
                    print("    [skip]"); continue
                lam       = lambda_reach(eft, sm, lambda_ref)
                centers   = 0.5 * (edges[:-1] + edges[1:])
                tail_frac = lam[centers >= args.tail_threshold].sum() / (lam.sum() + 1e-12)
                results.append({
                    "op": op, "max_lambda": float(lam.max()),
                    "peak_bin": int(lam.argmax()),
                    "peak_mll": float(centers[lam.argmax()]),
                    "tail_frac": float(tail_frac),
                    "lambda": lam.tolist(),
                })
            else:
                res = make_operator_plot(
                    f, op, args.channel, args.outdir, lumi_fb,
                    args.tail_threshold, lambda_ref, args.quad,
                )
                if res is not None:
                    results.append(res)

    if not results:
        print("No results.")
        return

    print(f"\nBuilding summary plots ({len(results)} operators) ...")
    make_lambda_heatmap(results, edges, args.outdir, lambda_ref)
    make_lambda_ranking(results, args.outdir, lambda_ref)

    print("\n{:<12}  {:>12}  {:>10}  {:>12}  {:>10}".format(
        "Operator", "max Λ* [TeV]", "peak mll", "tail_frac", "dominance"))
    print("-" * 62)
    for r in sorted(results, key=lambda x: x["max_lambda"], reverse=True):
        dom = "tail" if r["tail_frac"] > 0.5 else "bulk"
        print("{:<12}  {:>12.2f}  {:>9.0f}  {:>11.1%}  {:>10}".format(
            r["op"], r["max_lambda"], r["peak_mll"], r["tail_frac"], dom))

    print("\nDone.")


if __name__ == "__main__":
    main()
