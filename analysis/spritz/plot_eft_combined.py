#!/usr/bin/env python3
"""
EFT operator plots with k-factor rescaling.

For each of the 27 SMEFT operators, produces a stacked-background plot with the
EFT signal at c=1 overlaid (replacing DYll / MiNNLO), plus a ratio panel (EFT/SM).

Run from the config dir (where histos.root lives), inside apptainer + analysis_venv:
    python3 /grid_mnt/.../DY_2026/analysis/spritz/plot_eft_combined.py \
        [--input histos.root] [--region inc_mm] [--variable mll] \
        [--outdir plots/eft_operators]

Or via the shell alias:
    spritz-plot-eft [--region inc_mm] [--variable mll] [--outdir plots/eft_operators]
"""

import argparse
import os
from copy import deepcopy

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import uproot

# -- operator list (27 operators) ------------------------------
OPERATORS = [
    "cHDD", "cHWB", "cbWRe", "cbBRe", "cHj1", "cHQ1", "cHj3", "cHQ3",
    "cHu", "cHd", "cHbq", "cHl1", "cHl3", "cHe", "cll1", "clj1", "clj3",
    "cQl1", "cQl3", "ceu", "ced", "cbe", "cje", "cQe", "clu", "cld", "cbl",
]

# background stack order (bottom -> top)
BKG_STACK = ["DYtt", "GGToLL", "Single Top", "ZZ", "WZ", "WW", "TT"]

# fallback colours (Petroff palette) if config.py doesn't provide them
DEFAULT_COLORS = {
    "DYll":        "#5790fc",
    "TT":          "#f89c20",
    "WW":          "#e42536",
    "WZ":          "#964a8b",
    "ZZ":          "#9c9ca1",
    "Single Top":  "#7a21dd",
    "GGToLL":      "#92dadd",
    "DYtt":        "#2ca02c",
    "Data":        "black",
}


def _darker(color):
    r, g, b, _ = mpl.colors.to_rgba(color)
    f = 0.75
    return (r * f, g * f, b * f)


def get_vals(directory, name):
    vals, _ = directory[f"histo_{name}"].to_numpy()
    return vals.copy()


def get_variances(directory, name):
    return directory[f"histo_{name}"].variances().copy()


def get_edges(directory, name):
    _, edges = directory[f"histo_{name}"].to_numpy()
    return edges


def main():
    parser = argparse.ArgumentParser(description="EFT operator plots with k-factor")
    parser.add_argument("--input",    default="histos.root",        help="Path to histos.root")
    parser.add_argument("--region",   default="inc_mm",             help="Region key in histos.root")
    parser.add_argument("--variable", default="mll",                help="Variable key in histos.root")
    parser.add_argument("--outdir",   default="plots/eft_operators", help="Output directory for PNGs")
    parser.add_argument("--blind-above", type=float, default=500.,
                        help="Blind data above this mll value in GeV (default: 500)")
    parser.add_argument("--no-blind", action="store_true",
                        help="Disable blinding and show all data")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # -- try to get colours / lumi from config.py in CWD ---------------------
    colors = dict(DEFAULT_COLORS)
    lumi = 59.74
    year_label = "2018"
    try:
        from spritz.framework.framework import get_analysis_dict
        ad = get_analysis_dict()
        colors.update(ad.get("colors", {}))
        lumi = ad.get("lumi", lumi)
        year_label = ad.get("year_label", year_label)
    except Exception:
        pass  # analysis_venv may not have spritz_fabian in PYTHONPATH; defaults are fine

    # -- open histos.root -----------------------------------------------------
    f = uproot.open(args.input)
    directory = f[f"{args.region}/{args.variable}"]

    edges   = get_edges(directory, "sm")
    widths  = edges[1:] - edges[:-1]


    # -- read backgrounds -----------------------------------------------------
    bkg_vals = {}
    for s in BKG_STACK:
        try:
            bkg_vals[s] = get_vals(directory, s)
        except Exception:
            pass  # sample absent in this histos.root; silently skip

    dyll      = get_vals(directory, "DYll")
    sm        = get_vals(directory, "sm")
    data      = get_vals(directory, "Data")
    data_var  = get_variances(directory, "Data")
    centers   = 0.5 * (edges[:-1] + edges[1:])

    # -- bin-by-bin k-factor --------------------------------------------------
    k = np.where(sm > 0, dyll / sm, 1.0)

    # -- background cumulative stack (shared across all operators) -------------
    present = [s for s in BKG_STACK if s in bkg_vals]
    stack   = np.array([bkg_vals[s] for s in present])   # shape (n_samples, n_bins)
    cumsum  = np.cumsum(stack, axis=0)                    # shape (n_samples, n_bins)
    bkg_total = cumsum[-1]

    # -- matplotlib style -----------------------------------------------------
    style = deepcopy(hep.style.CMS)
    style["font.size"] = 12
    style["figure.figsize"] = (6, 6)
    plt.style.use(style)

    # -- one plot per operator ------------------------------------------------
    for op in OPERATORS:
        try:
            w1  = get_vals(directory, f"w1_{op}")
            wm1 = get_vals(directory, f"wm1_{op}")
        except Exception:
            print(f"  [skip] {op}: histo_w1_{op} / histo_wm1_{op} not found in {args.input}")
            continue

        # k-factor rescaling
        sm_k  = sm  * k             #becomes MiNNLO
        w1_k  = w1  * k             #becomes consistent with MiNNLO
        wm1_k = wm1 * k 

        # EFT decomposition
        lin  = 0.5 * (w1_k - wm1_k)
        quad = 0.5 * (w1_k + wm1_k) - sm_k

        # EFT at c=1: sm_k + lin + quad = w1_k  (sanity: c=0 -> sm_k = dyll )
        eft = sm_k + lin + quad

        sm_total  = bkg_total + dyll
        eft_total = bkg_total + eft

        # -- figure -----------------------------------------------------------
        fig, (ax_top, ax_bot) = plt.subplots(
            2, 1, sharex=True,
            gridspec_kw={"height_ratios": [3, 1]},
            dpi=200,
        )
        fig.tight_layout(pad=-0.5)
        hep.cms.label("", data=False, lumi=round(lumi, 2), ax=ax_top, year=year_label)

        # draw stacked backgrounds (bottom -> top)
        for i, name in enumerate(present):
            base = cumsum[i - 1] if i > 0 else np.zeros_like(bkg_total)
            ax_top.stairs(
                cumsum[i] / widths,
                edges=edges,
                baseline=base / widths,
                fill=True,
                color=colors.get(name, "grey"),
                edgecolor=_darker(colors.get(name, "grey")),
                linewidth=0.5,
                label=name,
                zorder=1.0 - i * 0.01,
            )

        # SM total (bkg + MiNNLO DYll) as dashed reference line
        ax_top.stairs(
            sm_total / widths, edges=edges,
            color=colors.get("DYll", DEFAULT_COLORS["DYll"]),
            linewidth=1.5, linestyle="dashed",
            label="SM (MiNNLO)", fill=False, zorder=2,
        )

        # EFT total (bkg + EFT signal at c=1)
        ax_top.stairs(
            eft_total / widths, edges=edges,
            color="crimson", linewidth=1.5,
            label=f"EFT {op} (c=1)", fill=False, zorder=3,
        )

        # Data as black dots with Poisson error bars (blind above threshold)
        data_unc = np.sqrt(np.abs(data_var))
        if not args.no_blind:
            blind_mask = centers > args.blind_above
            data_plot     = np.where(blind_mask, np.nan, data)
            data_unc_plot = np.where(blind_mask, np.nan, data_unc)
        else:
            data_plot     = data
            data_unc_plot = data_unc
        blind_label = f" [blind > {int(args.blind_above)} GeV]" if not args.no_blind else ""
        ax_top.errorbar(
            centers,
            data_plot / widths,
            yerr=data_unc_plot / widths,
            fmt="o", markersize=4, color="black",
            label=f"Data [{int(round(data.sum()))}]{blind_label}",
            zorder=4,
        )

        # y-axis range
        ymax = max(np.max(sm_total / widths), np.max(eft_total / widths), np.max(data / widths))
        pos_vals = np.concatenate([v[v > 0] / widths[v > 0] for v in stack if np.any(v > 0)])
        ymin = max(1e-4, 0.3 * np.min(pos_vals)) if pos_vals.size else 1e-4

        ax_top.set_yscale("log")
        ax_top.set_ylim(ymin, ymax * 3e3)
        ax_top.set_ylabel("Events / GeV")
        ax_top.tick_params(labelbottom=False)
        ax_top.legend(loc="upper right", fontsize=7, ncols=2, framealpha=0.8)

        # ratio panel: EFT/SM and Data/SM
        denom = np.where(sm_total > 0, sm_total, 1e-30)
        ratio_eft  = eft_total / denom
        ratio_data = data / denom

        ax_bot.stairs(ratio_eft, edges=edges, color="crimson", linewidth=1.2, label="EFT/SM")
        ax_bot.errorbar(
            centers,
            np.where(blind_mask if not args.no_blind else False, np.nan, ratio_data),
            yerr=np.where(blind_mask if not args.no_blind else False, np.nan, data_unc / denom),
            fmt="o", markersize=3, color="black", label="Data/SM",
            zorder=4,
        )
        ax_bot.axhline(1.0, color="black", linewidth=0.8, linestyle="dashed")
        ax_bot.set_ylabel("/ SM")
        # auto-range: pad 20% around the actual spread, minimum ±0.05 around 1
        finite = ratio_eft[np.isfinite(ratio_eft)]
        rmin = min(finite.min(), 1.0) - 0.05
        rmax = max(finite.max(), 1.0) + 0.05
        pad  = 0.2 * (rmax - rmin)
        ax_bot.set_ylim(rmin - pad, rmax + pad)
        ax_bot.set_xlabel(r"$m_{\ell\ell}$ (GeV)")
        ax_bot.set_xscale("log")
        ax_bot.set_xlim(edges[0], edges[-1])
        ax_bot.legend(loc="upper right", fontsize=7, framealpha=0.8)

        stem = os.path.join(args.outdir, f"eft_{op}")
        for ext in ("png", "pdf"):
            fig.savefig(f"{stem}.{ext}", facecolor="white", pad_inches=0.1, bbox_inches="tight")
        plt.close()
        print(f"  {op:10s}  ->  {stem}.png / .pdf")

    print("Done.")


if __name__ == "__main__":
    main()
