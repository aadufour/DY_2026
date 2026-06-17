#!/usr/bin/env python3
"""
EFT operator plots with k-factor rescaling.

For each of the 27 SMEFT operators, produces a stacked-background plot with the
EFT signal at c=±1 overlaid (replacing DYll / MiNNLO), plus a ratio panel (EFT/SM).

Run from the config dir (where histos.root lives), inside apptainer + analysis_venv:
    python3 .../plot_eft.py [--input histos.root] [--region inc_mm]
                            [--variable mll|costhetastar|rapll_abs|all]
                            [--outdir plots/eft_operators]

--variable all  loops over every variable found in histos.root for the given region,
                writing each to a separate subdirectory under --outdir.
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

# ============================================================
# PLOT CONFIGURATION
# ============================================================

FONT_SIZE        = 22
LEGEND_FONTSIZE  = 16
TICK_LABELSIZE   = 18
LABEL_SIZE       = 20

FIG_SIZE         = (8, 8)
FIG_DPI          = 200

LEGEND_NCOLS_TOP = 2
LEGEND_LOC_TOP   = "upper right"
LEGEND_LOC_BOT   = "upper right"

MARKER_SIZE_DATA = 5
MARKER_SIZE_RATIO= 4

LINE_WIDTH_EFT   = 2.0
LINE_WIDTH_SM    = 1.5
LINE_WIDTH_BKG   = 0.8

# ============================================================

# Per-variable display metadata: label, whether to use log x-scale,
# whether to blind data (only meaningful for mll).
# These are defaults; get_analysis_dict() can override labels if spritz is available.
VAR_META_DEFAULT = {
    "mll":          {"label": r"$m_{\ell\ell}$ (GeV)", "log_x": True,  "blind": True},
    "costhetastar": {"label": r"$\cos\theta^*$",        "log_x": False, "blind": False},
    "rapll_abs":    {"label": r"$|y_{\ell\ell}|$",      "log_x": False, "blind": False},
}

OPERATORS = [
    "cHDD", "cHWB", "cbWRe", "cbBRe", "cHj1", "cHQ1", "cHj3", "cHQ3",
    "cHu", "cHd", "cHbq", "cHl1", "cHl3", "cHe", "cll1", "clj1", "clj3",
    "cQl1", "cQl3", "ceu", "ced", "cbe", "cje", "cQe", "clu", "cld", "cbl",
]

BKG_STACK = ["DYtt", "GGToLL", "Single Top", "ZZ", "WZ", "WW", "TT"]

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


def plot_one_variable(
    f, region, variable, var_meta, outdir,
    colors, lumi, year_label,
    shapes_path, blind_above, no_blind,
):
    """Produce one PNG+PDF per operator for a single (region, variable)."""
    directory = f[f"{region}/{variable}"]

    edges   = get_edges(directory, "sm")
    widths  = edges[1:] - edges[:-1]
    centers = 0.5 * (edges[:-1] + edges[1:])

    xlabel  = var_meta["label"]
    log_x   = var_meta["log_x"]
    do_blind = var_meta["blind"] and not no_blind

    # backgrounds
    bkg_vals = {}
    for s in BKG_STACK:
        try:
            bkg_vals[s] = get_vals(directory, s)
        except Exception:
            pass

    dyll     = get_vals(directory, "DYll")
    sm       = get_vals(directory, "sm")
    data     = get_vals(directory, "Data")
    data_var = get_variances(directory, "Data")

    k = np.where(sm > 0, dyll / sm, 1.0)

    present = [s for s in BKG_STACK if s in bkg_vals]
    stack   = np.array([bkg_vals[s] for s in present])
    cumsum  = np.cumsum(stack, axis=0)
    bkg_total = cumsum[-1]

    # systematic band from shapes.root
    syst_up_sm   = np.zeros_like(bkg_total)
    syst_down_sm = np.zeros_like(bkg_total)
    if shapes_path is not None:
        try:
            fs = uproot.open(shapes_path)
            bkg_nuisances = [
                ("QCDScale",      ["Single_Top", "TT", "WW", "DYtt", "DYll"]),
                ("PDFweight",     ["TT", "WW", "DYtt", "DYll"]),
                ("alphaS",        ["DYtt", "DYll"]),
                ("PSWeight",      ["Single_Top", "TT", "WW", "WZ", "ZZ", "DYtt", "DYll"]),
                ("mu_reco",       ["GGToLL", "Single_Top", "TT", "WW", "WZ", "ZZ", "DYtt"]),
                ("mu_idiso",      ["GGToLL", "Single_Top", "TT", "WW", "WZ", "ZZ", "DYtt"]),
                ("mu_trig",       ["GGToLL", "Single_Top", "TT", "WW", "WZ", "ZZ", "DYtt"]),
                ("PU",            ["GGToLL", "Single_Top", "TT", "WW", "WZ", "ZZ", "DYtt"]),
                ("prefireWeight", ["GGToLL", "Single_Top", "TT", "WW", "WZ", "ZZ", "DYtt"]),
                ("tt_ptrw",       ["TT"]),
                ("rochester_stat",["GGToLL", "Single_Top", "TT", "WW", "WZ", "ZZ", "DYtt"]),
                ("rochester_syst",["GGToLL", "Single_Top", "TT", "WW", "WZ", "ZZ", "DYtt"]),
                ("lumi",          ["Single_Top", "TT", "WW", "WZ", "ZZ", "DYtt", "DYll", "GGToLL"]),
            ]
            sm_samples_shapes = list(bkg_vals.keys()) + ["DYll"]
            nom_total = np.zeros_like(bkg_total)
            for s in sm_samples_shapes:
                try:
                    nom_total += fs[f"histo_{s.replace(' ', '_')}"].values().copy()
                except Exception:
                    pass

            for nuis_name, affected in bkg_nuisances:
                tot_up   = nom_total.copy()
                tot_down = nom_total.copy()
                for s in affected:
                    skey = s.replace(" ", "_")
                    try:
                        nom_s = fs[f"histo_{skey}"].values().copy()
                        if nuis_name == "lumi":
                            tot_up   += (1.0084 - 1.0) * nom_s
                            tot_down += (1.0 / 1.0084 - 1.0) * nom_s
                        else:
                            up_s = fs[f"histo_{skey}_{nuis_name}Up"].values().copy()
                            do_s = fs[f"histo_{skey}_{nuis_name}Down"].values().copy()
                            tot_up   += (up_s - nom_s)
                            tot_down += (do_s - nom_s)
                    except Exception:
                        pass
                syst_up_sm   += np.square(tot_up   - nom_total)
                syst_down_sm += np.square(tot_down - nom_total)
            syst_up_sm   = np.sqrt(syst_up_sm)
            syst_down_sm = np.sqrt(syst_down_sm)
            print(f"  Syst band loaded from {shapes_path}")
        except Exception as e:
            print(f"  [warn] Could not compute syst band from {shapes_path}: {e}")

    # matplotlib style
    style = deepcopy(hep.style.CMS)
    style["font.size"]        = FONT_SIZE
    style["axes.labelsize"]   = LABEL_SIZE
    style["xtick.labelsize"]  = TICK_LABELSIZE
    style["ytick.labelsize"]  = TICK_LABELSIZE
    style["legend.fontsize"]  = LEGEND_FONTSIZE
    style["figure.figsize"]   = FIG_SIZE
    plt.style.use(style)

    # blinding mask — only applied when variable supports it
    if do_blind:
        blind_mask = centers > blind_above
    else:
        blind_mask = np.zeros(len(centers), dtype=bool)

    data_unc = np.sqrt(np.abs(data_var))

    for op in OPERATORS:
        try:
            w1  = get_vals(directory, f"w1_{op}")
            wm1 = get_vals(directory, f"wm1_{op}")
        except Exception:
            print(f"  [skip] {op}: histo_w1_{op} / histo_wm1_{op} not found")
            continue

        sm_k  = sm  * k
        w1_k  = w1  * k
        wm1_k = wm1 * k

        lin  = 0.5 * (w1_k - wm1_k)
        quad = 0.5 * (w1_k + wm1_k) - sm_k

        eft_p1 = sm_k + lin + quad   # = w1_k
        eft_m1 = sm_k - lin + quad   # = wm1_k

        sm_total   = bkg_total + dyll
        eft_total  = bkg_total + eft_p1
        eftm_total = bkg_total + eft_m1

        fig, (ax_top, ax_bot) = plt.subplots(
            2, 1, sharex=True,
            gridspec_kw={"height_ratios": [3, 1]},
            dpi=FIG_DPI,
        )
        fig.tight_layout(pad=-0.5)
        hep.cms.label("", data=False, lumi=round(lumi, 2), ax=ax_top, year=year_label)

        for i, name in enumerate(present):
            base = cumsum[i - 1] if i > 0 else np.zeros_like(bkg_total)
            ax_top.stairs(
                cumsum[i] / widths, edges=edges,
                baseline=base / widths,
                fill=True,
                color=colors.get(name, "grey"),
                edgecolor=_darker(colors.get(name, "grey")),
                linewidth=LINE_WIDTH_BKG,
                label=name, zorder=1.0 - i * 0.01,
            )

        ax_top.stairs(
            sm_total / widths, edges=edges,
            color=colors.get("DYll", DEFAULT_COLORS["DYll"]),
            linewidth=LINE_WIDTH_SM, linestyle="dashed",
            label="SM (MiNNLO)", fill=False, zorder=2,
        )

        if shapes_path is not None and np.any(syst_up_sm > 0):
            ax_top.fill_between(
                np.repeat(edges, 2)[1:-1],
                np.repeat((sm_total - syst_down_sm) / widths, 2),
                np.repeat((sm_total + syst_up_sm)   / widths, 2),
                step="pre", alpha=0.30, color="grey",
                hatch="///", linewidth=0, label="Syst. unc.", zorder=2,
            )

        ax_top.stairs(
            eft_total / widths, edges=edges,
            color="crimson", linewidth=LINE_WIDTH_EFT,
            label=f"EFT {op} (c=+1)", fill=False, zorder=3,
        )
        ax_top.stairs(
            eftm_total / widths, edges=edges,
            color="steelblue", linewidth=LINE_WIDTH_EFT,
            label=f"EFT {op} (c=-1)", fill=False, zorder=3,
        )

        data_plot     = np.where(blind_mask, np.nan, data)
        data_unc_plot = np.where(blind_mask, np.nan, data_unc)
        blind_label = f" [blind > {int(blind_above)} GeV]" if do_blind else ""
        ax_top.errorbar(
            centers, data_plot / widths,
            yerr=data_unc_plot / widths,
            fmt="o", markersize=MARKER_SIZE_DATA, color="black",
            label=f"Data [{int(round(data.sum()))}]{blind_label}",
            zorder=4,
        )

        ymax = max(
            np.nanmax(sm_total / widths),
            np.nanmax(eft_total / widths),
            np.nanmax(eftm_total / widths),
            np.nanmax(data / widths),
        )
        pos_vals_list = [v[v > 0] / widths[v > 0] for v in stack if np.any(v > 0)]
        pos_vals = np.concatenate(pos_vals_list) if pos_vals_list else np.array([1e-4])
        ymin = max(1e-4, 0.3 * np.min(pos_vals))

        ax_top.set_yscale("log")
        ax_top.set_ylim(ymin, ymax * 3e3)
        ax_top.set_ylabel("Events / GeV")
        ax_top.tick_params(labelbottom=False)
        ax_top.legend(loc=LEGEND_LOC_TOP, ncols=LEGEND_NCOLS_TOP, framealpha=0.8)

        denom      = np.where(sm_total > 0, sm_total, 1e-30)
        ratio_eft  = eft_total  / denom
        ratio_eftm = eftm_total / denom
        ratio_data = data / denom

        ax_bot.stairs(ratio_eft,  edges=edges, color="crimson",   linewidth=1.2, label="c=+1 / SM")
        ax_bot.stairs(ratio_eftm, edges=edges, color="steelblue", linewidth=1.2, label="c=-1 / SM")
        ax_bot.errorbar(
            centers,
            np.where(blind_mask, np.nan, ratio_data),
            yerr=np.where(blind_mask, np.nan, data_unc / denom),
            fmt="o", markersize=MARKER_SIZE_RATIO, color="black", label="Data/SM",
            zorder=4,
        )
        ax_bot.axhline(1.0, color="black", linewidth=0.8, linestyle="dashed")
        if shapes_path is not None and np.any(syst_up_sm > 0):
            ax_bot.fill_between(
                np.repeat(edges, 2)[1:-1],
                np.repeat((sm_total - syst_down_sm) / denom, 2),
                np.repeat((sm_total + syst_up_sm)   / denom, 2),
                step="pre", alpha=0.30, color="grey",
                hatch="///", linewidth=0, zorder=0,
            )

        ax_bot.set_ylabel("Ratio")
        finite = np.concatenate([ratio_eft[np.isfinite(ratio_eft)], ratio_eftm[np.isfinite(ratio_eftm)]])
        half = max(np.max(np.abs(finite - 1.0)) * 1.2, 0.05) if finite.size else 0.3
        ax_bot.set_ylim(1.0 - half, 1.0 + half)
        ax_bot.set_xlabel(xlabel)
        if log_x:
            ax_bot.set_xscale("log")
        ax_bot.set_xlim(edges[0], edges[-1])
        ax_bot.legend(loc=LEGEND_LOC_BOT, framealpha=0.8)

        stem = os.path.join(outdir, f"eft_{op}")
        for ext in ("png", "pdf"):
            fig.savefig(f"{stem}.{ext}", facecolor="white", pad_inches=0.1, bbox_inches="tight")
        plt.close()
        print(f"  {op:12s}  ->  {stem}.png / .pdf")


def main():
    parser = argparse.ArgumentParser(description="EFT operator plots with k-factor")
    parser.add_argument("--input",    default="histos.root",         help="Path to histos.root")
    parser.add_argument("--shapes",   default=None,                  help="Path to shapes.root (for syst band)")
    parser.add_argument("--region",   default="inc_mm",              help="Region key in histos.root")
    parser.add_argument("--variable", default="mll",
                        help="Variable key in histos.root, or 'all' to loop over every variable in the region")
    parser.add_argument("--outdir",   default="plots/eft_operators",  help="Output directory for PNGs")
    parser.add_argument("--blind-above", type=float, default=500.,
                        help="Blind data above this value (only for variables with blinding enabled, default: 500 GeV for mll)")
    parser.add_argument("--no-blind", action="store_true",
                        help="Disable blinding for all variables")
    args = parser.parse_args()

    # -- try to get colours / lumi / variable labels from config.py ------------
    colors = dict(DEFAULT_COLORS)
    lumi = 59.74
    year_label = "2018"
    var_meta_config = {}
    try:
        from spritz.framework.framework import get_analysis_dict
        ad = get_analysis_dict()
        colors.update(ad.get("colors", {}))
        lumi = ad.get("lumi", lumi)
        year_label = ad.get("year_label", year_label)
        for vname, vdict in ad.get("variables", {}).items():
            label = vdict.get("label", vname)
            unit  = vdict.get("unit", "")
            if unit:
                label = f"{label} ({unit})"
            var_meta_config[vname] = {"label": label}
    except Exception:
        pass

    # -- resolve which variables to plot ---------------------------------------
    f = uproot.open(args.input)
    if args.variable == "all":
        region_dir = f[args.region]
        variables = [
            k.split(";")[0] for k, cls in region_dir.classnames().items()
            if "TDirectory" in cls and "/" not in k
        ]
    else:
        variables = [args.variable]

    # -- loop over variables ---------------------------------------------------
    for variable in variables:
        # build var_meta: start from hardcoded defaults, overlay config labels
        meta = dict(VAR_META_DEFAULT.get(variable, {"label": variable, "log_x": False, "blind": False}))
        if variable in var_meta_config:
            meta["label"] = var_meta_config[variable]["label"]

        if args.variable == "all":
            outdir = os.path.join(args.outdir, variable)
        else:
            outdir = args.outdir
        os.makedirs(outdir, exist_ok=True)

        print(f"\n=== {args.region} / {variable} ===")
        plot_one_variable(
            f, args.region, variable, meta, outdir,
            colors, lumi, year_label,
            args.shapes, args.blind_above, args.no_blind,
        )

    print("\nDone.")


if __name__ == "__main__":
    main()
