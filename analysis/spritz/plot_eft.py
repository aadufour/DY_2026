#!/usr/bin/env python3
"""
EFT operator plots — read-only, no data modification.

For each of the 27 SMEFT operators, produces a stacked-background plot with the
EFT signal at c=±1 overlaid, plus a ratio panel (EFT/SM).

The k-factor (MiNNLO / SMEFTsim LO SM) is applied upstream in spritz-postproc-eft
(post_process.py). Templates read from shapes.root are already k-scaled — this
script does NOT compute or apply any k-factor.

Sources:
  --input   histos.root   background stack + EFT templates (from spritz-postproc-eft)
  --shapes  shapes.root   systematic band on SM prediction (from spritz-cards, optional)

Run from the config dir, with analysis_venv active:
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
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
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
# Physical bin edges for the triple_diff unrolled histogram.
# Unrolling order: rapll_abs slowest, costhetastar middle, mll fastest (.T.flatten())
# Flat index for (irapll, icos): irapll * N_COSTH * N_MLL + icos * N_MLL
MLL_EDGES   = np.array([40, 60, 80, 100, 120, 140, 180, 220, 270, 350, 500, 700, 1000, 1500, 3000], dtype=float)
COSTH_EDGES = np.array([-1, -0.6, -0.2, 0.2, 0.6, 1], dtype=float)
RAPLL_EDGES = np.array([0, 0.48, 0.96, 1.44, 2.4], dtype=float)
N_MLL   = len(MLL_EDGES)   - 1   # 14
N_COSTH = len(COSTH_EDGES) - 1   # 5
N_RAPLL = len(RAPLL_EDGES) - 1   # 4


def _td_slice(irapll, icos):
    """Return the slice into the 280-bin triple_diff array for panel (irapll, icos)."""
    start = irapll * N_COSTH * N_MLL + icos * N_MLL
    return slice(start, start + N_MLL)


VAR_META_DEFAULT = {
    # blind_all=True: hide all data (variable is mll-integrated, can't threshold-cut)
    # blind_all=False: hide data above blind_above (only meaningful for mll)
    # range_max: ignore bins above this value when auto-ranging the ratio y-axis
    "mll":          {"label": r"$m_{\ell\ell}$ (GeV)", "log_x": True,  "blind": True,  "blind_all": False, "range_max": None},
    "costhetastar": {"label": r"$\cos\theta^*$",        "log_x": False, "blind": False, "blind_all": False, "range_max": None},
    "rapll_abs":    {"label": r"$|y_{\ell\ell}|$",      "log_x": False, "blind": False, "blind_all": False, "range_max": 2.4},
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
    do_blind  = var_meta["blind"] and not no_blind
    blind_all = var_meta.get("blind_all", False)

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

    # blinding mask
    if not do_blind:
        blind_mask = np.zeros(len(centers), dtype=bool)
    elif blind_all:
        # variable is mll-integrated: can't threshold-cut, hide all data
        blind_mask = np.ones(len(centers), dtype=bool)
    else:
        blind_mask = centers > blind_above

    data_unc = np.sqrt(np.abs(data_var))

    for op in OPERATORS:
        try:
            w1  = get_vals(directory, f"w1_{op}")
            wm1 = get_vals(directory, f"wm1_{op}")
        except Exception:
            print(f"  [skip] {op}: histo_w1_{op} / histo_wm1_{op} not found")
            continue

        lin  = 0.5 * (w1 - wm1)
        quad = 0.5 * (w1 + wm1) - sm

        eft_p1 = sm + lin + quad   # = w1
        eft_m1 = sm - lin + quad   # = wm1

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
        if do_blind:
            blind_label = f" [blind > {int(blind_above)} GeV]"
        elif variable in ("costhetastar", "rapll_abs"):
            blind_label = " [mll < 500 GeV]"
        else:
            blind_label = ""
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
        ax_top.set_ylabel("Events / GeV" if var_meta.get("unit") else "Events")
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
            fmt="o", markersize=MARKER_SIZE_RATIO, color="black", label="_nolegend_",
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
        ratio_data_visible = np.where(blind_mask, np.nan, ratio_data)
        candidates = [ratio_eft, ratio_eftm, ratio_data_visible]
        if shapes_path is not None and np.any(syst_up_sm > 0):
            candidates += [(sm_total + syst_up_sm) / denom, (sm_total - syst_down_sm) / denom]
        # mask bins outside selection acceptance (e.g. |y_ll| > 2.4 always empty)
        range_mask = centers <= (var_meta.get("range_max") or np.inf)
        all_finite = np.concatenate([a[range_mask][np.isfinite(a[range_mask])] for a in candidates])
        half = max(np.max(np.abs(all_finite - 1.0)) * 1.2, 0.05) if all_finite.size else 0.3
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


def plot_triple_diff(f, region, outdir, colors, lumi, year_label, shapes_path):
    """Multi-panel EFT plot for the triple_diff unrolled histogram.

    Creates one figure per operator: a 4-row (rapll_abs) × 5-col (costhetastar)
    grid.  Each panel shows the 10 mll bins for that (rapll_abs, costhetastar) slice,
    with a stacked background + SM + c=±1 EFT overlay and a ratio panel below.
    """
    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

    directory = f[f"{region}/triple_diff"]
    mll_widths = MLL_EDGES[1:] - MLL_EDGES[:-1]
    mll_centers = 0.5 * (MLL_EDGES[:-1] + MLL_EDGES[1:])

    def _get(name):
        return directory[f"histo_{name}"].values().copy()

    def _getvar(name):
        return directory[f"histo_{name}"].variances().copy()

    try:
        sm   = _get("sm")
        dyll = _get("DYll")
        data     = _get("Data")
        data_var = _getvar("Data")
    except Exception as e:
        print(f"  [skip] triple_diff: cannot read sm/DYll/Data: {e}")
        return

    bkg_vals = {}
    for s in BKG_STACK:
        try:
            bkg_vals[s] = _get(s)
        except Exception:
            pass
    present  = [s for s in BKG_STACK if s in bkg_vals]
    n_td = N_MLL * N_COSTH * N_RAPLL
    stack    = np.array([bkg_vals[s] for s in present])   # (n_bkg, n_td)
    cumsum   = np.cumsum(stack, axis=0)                    # (n_bkg, n_td)
    bkg_total = cumsum[-1]                                 # (n_td,)

    # --- systematic band from shapes.root ---
    syst_up   = np.zeros(n_td)
    syst_down = np.zeros(n_td)
    if shapes_path is not None:
        try:
            fs = uproot.open(shapes_path)

            def _sget(name):
                return fs[f"histo_{name}"].values().copy()

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
            nom_total = np.zeros(n_td)
            for s in sm_samples_shapes:
                try:
                    nom_total += _sget(s.replace(" ", "_"))
                except Exception:
                    pass

            for nuis_name, affected in bkg_nuisances:
                tot_up   = nom_total.copy()
                tot_down = nom_total.copy()
                for s in affected:
                    skey = s.replace(" ", "_")
                    try:
                        nom_s = _sget(skey)
                        if nuis_name == "lumi":
                            tot_up   += (1.0084 - 1.0) * nom_s
                            tot_down += (1.0 / 1.0084 - 1.0) * nom_s
                        else:
                            up_s = _sget(f"{skey}_{nuis_name}Up")
                            do_s = _sget(f"{skey}_{nuis_name}Down")
                            tot_up   += (up_s - nom_s)
                            tot_down += (do_s - nom_s)
                    except Exception:
                        pass
                syst_up   += np.square(tot_up   - nom_total)
                syst_down += np.square(tot_down - nom_total)
            syst_up   = np.sqrt(syst_up)
            syst_down = np.sqrt(syst_down)
            print(f"  Syst band loaded from {shapes_path}")
        except Exception as e:
            print(f"  [warn] Could not compute syst band: {e}")

    style = deepcopy(hep.style.CMS)
    style["font.size"]       = 7
    style["axes.labelsize"]  = 6
    style["xtick.labelsize"] = 5
    style["ytick.labelsize"] = 5
    style["legend.fontsize"] = 5
    plt.style.use(style)

    PANEL_W, PANEL_H = 3.5, 3.5
    fig_w = N_COSTH * PANEL_W
    fig_h = N_RAPLL * PANEL_H

    for op in OPERATORS:
        try:
            w1  = _get(f"w1_{op}")
            wm1 = _get(f"wm1_{op}")
        except Exception:
            print(f"  [skip] triple_diff {op}: w1/wm1 not found")
            continue

        fig = plt.figure(figsize=(fig_w, fig_h))
        outer = GridSpec(N_RAPLL, N_COSTH, figure=fig, hspace=0.38, wspace=0.22)

        for irapll in range(N_RAPLL):
            for icos in range(N_COSTH):
                sl = _td_slice(irapll, icos)

                inner = GridSpecFromSubplotSpec(
                    2, 1, subplot_spec=outer[irapll, icos],
                    height_ratios=[3, 1], hspace=0.06,
                )
                ax_top = fig.add_subplot(inner[0])
                ax_bot = fig.add_subplot(inner[1], sharex=ax_top)

                bkg_sl      = cumsum[:, sl]
                bkg_tot_sl  = bkg_total[sl]
                sm_sl       = sm[sl]
                dyll_sl     = dyll[sl]
                eft_p_sl    = w1[sl]
                eft_m_sl    = wm1[sl]
                data_sl     = data[sl]
                data_unc_sl = np.sqrt(np.abs(data_var[sl]))
                syst_up_sl  = syst_up[sl]
                syst_dn_sl  = syst_down[sl]

                sm_total_sl  = bkg_tot_sl + dyll_sl
                eft_tot_sl   = bkg_tot_sl + eft_p_sl
                eftm_tot_sl  = bkg_tot_sl + eft_m_sl

                is_first = (irapll == 0 and icos == 0)

                for i, name in enumerate(present):
                    base = bkg_sl[i - 1] if i > 0 else np.zeros(N_MLL)
                    ax_top.stairs(
                        bkg_sl[i] / mll_widths, edges=MLL_EDGES,
                        baseline=base / mll_widths,
                        fill=True,
                        color=colors.get(name, "grey"),
                        edgecolor=_darker(colors.get(name, "grey")),
                        linewidth=0.4,
                        label=name if is_first else "_nolegend_",
                        zorder=1.0 - i * 0.01,
                    )
                ax_top.stairs(
                    sm_total_sl / mll_widths, edges=MLL_EDGES,
                    color=colors.get("DYll", DEFAULT_COLORS["DYll"]),
                    linewidth=0.8, linestyle="dashed",
                    label="SM (MiNNLO)" if is_first else "_nolegend_",
                    fill=False, zorder=2,
                )

                has_syst = np.any(syst_up_sl > 0)
                if has_syst:
                    ax_top.fill_between(
                        np.repeat(MLL_EDGES, 2)[1:-1],
                        np.repeat((sm_total_sl - syst_dn_sl) / mll_widths, 2),
                        np.repeat((sm_total_sl + syst_up_sl) / mll_widths, 2),
                        step="pre", alpha=0.30, color="grey",
                        hatch="///", linewidth=0,
                        label="Syst. unc." if is_first else "_nolegend_", zorder=2,
                    )

                ax_top.stairs(
                    eft_tot_sl / mll_widths, edges=MLL_EDGES,
                    color="crimson", linewidth=0.8,
                    label=f"{op} c=+1" if is_first else "_nolegend_",
                    fill=False, zorder=3,
                )
                ax_top.stairs(
                    eftm_tot_sl / mll_widths, edges=MLL_EDGES,
                    color="steelblue", linewidth=0.8,
                    label=f"{op} c=-1" if is_first else "_nolegend_",
                    fill=False, zorder=3,
                )

                _data_blind = mll_centers > 500
                ax_top.errorbar(
                    mll_centers,
                    np.where(_data_blind, np.nan, data_sl / mll_widths),
                    yerr=np.where(_data_blind, np.nan, data_unc_sl / mll_widths),
                    fmt="o", markersize=2, color="black", linewidth=0.6,
                    label=f"Data [{int(data_sl.sum())}] [mll < 500 GeV]" if is_first else "_nolegend_",
                    zorder=4,
                )

                ax_top.set_yscale("log")
                ax_top.set_xscale("log")
                ax_top.tick_params(labelbottom=False)

                # Per-panel bin range label — top right to avoid covering the peak
                cos_lo = COSTH_EDGES[icos]
                cos_hi = COSTH_EDGES[icos + 1]
                rap_lo = RAPLL_EDGES[irapll]
                rap_hi = RAPLL_EDGES[irapll + 1]
                ax_top.text(
                    0.97, 0.97,
                    (f"$\\cos\\theta^* \\in [{cos_lo:.1f},{cos_hi:.1f}]$\n"
                     f"$|y_{{\\ell\\ell}}| \\in [{rap_lo:.2f},{rap_hi:.2f}]$"),
                    transform=ax_top.transAxes,
                    va="top", ha="right", fontsize=5, linespacing=1.4,
                )

                if is_first:
                    hep.cms.label(
                        "", data=True, lumi=round(lumi, 2),
                        ax=ax_top, year=year_label, fontsize=6,
                    )

                if icos == 0:
                    ax_top.set_ylabel("Events/GeV", fontsize=5)

                # Ratio panel
                denom      = np.where(sm_total_sl > 0, sm_total_sl, 1e-30)
                ratio_p    = eft_tot_sl  / denom
                ratio_m    = eftm_tot_sl / denom
                ratio_data = data_sl / denom
                data_blind_mask = mll_centers > 500
                ratio_data_plot = np.where(data_blind_mask, np.nan, ratio_data)
                data_unc_plot   = np.where(data_blind_mask, np.nan, data_unc_sl / denom)

                ax_bot.stairs(ratio_p, edges=MLL_EDGES, color="crimson",   linewidth=0.8)
                ax_bot.stairs(ratio_m, edges=MLL_EDGES, color="steelblue", linewidth=0.8)
                ax_bot.errorbar(
                    mll_centers, ratio_data_plot,
                    yerr=data_unc_plot,
                    fmt="o", markersize=2, color="black", linewidth=0.6, zorder=4,
                )
                ax_bot.axhline(1.0, color="black", linewidth=0.6, linestyle="dashed")

                if has_syst:
                    ax_bot.fill_between(
                        np.repeat(MLL_EDGES, 2)[1:-1],
                        np.repeat((sm_total_sl - syst_dn_sl) / denom, 2),
                        np.repeat((sm_total_sl + syst_up_sl) / denom, 2),
                        step="pre", alpha=0.30, color="grey",
                        hatch="///", linewidth=0, zorder=0,
                    )

                candidates = [ratio_p, ratio_m, ratio_data_plot]
                if has_syst:
                    candidates += [
                        (sm_total_sl + syst_up_sl) / denom,
                        (sm_total_sl - syst_dn_sl) / denom,
                    ]
                finite_r = np.concatenate([a[np.isfinite(a)] for a in candidates])
                half = max(np.max(np.abs(finite_r - 1.0)) * 1.2, 0.05) if finite_r.size else 0.3
                ax_bot.set_ylim(1.0 - half, 1.0 + half)
                ax_bot.set_xlim(MLL_EDGES[0], MLL_EDGES[-1])

                if irapll == N_RAPLL - 1:
                    ax_bot.set_xlabel(r"$m_{\ell\ell}$ [GeV]", fontsize=5)
                if icos == 0:
                    ax_bot.set_ylabel("Ratio", fontsize=5)

        fig.suptitle(f"EFT {op} (c=±1)  —  triple-diff", fontsize=9, y=0.995)

        # Common legend across all panels
        legend_handles = (
            [Patch(color=colors.get(s, "grey"), label=s) for s in present]
            + [
                Line2D([0], [0], color=colors.get("DYll", DEFAULT_COLORS["DYll"]), linestyle="dashed", linewidth=0.8, label="SM (MiNNLO)"),
                Line2D([0], [0], color="crimson",   linewidth=0.8, label=f"{op} c=+1"),
                Line2D([0], [0], color="steelblue", linewidth=0.8, label=f"{op} c=-1"),
                Line2D([0], [0], color="black", marker="o", markersize=2, linewidth=0.6, label=f"Data [mll < 500 GeV]"),
            ]
        )
        fig.legend(
            handles=legend_handles,
            fontsize=5, ncols=4, framealpha=0.7,
            loc="lower center", bbox_to_anchor=(0.5, -0.02),
            handlelength=1.2, handletextpad=0.4, columnspacing=0.8,
        )

        stem = os.path.join(outdir, f"eft_triple_diff_{op}")
        for ext in ("png", "pdf"):
            fig.savefig(
                f"{stem}.{ext}", facecolor="white",
                pad_inches=0.05, bbox_inches="tight", dpi=FIG_DPI,
            )
        plt.close()
        print(f"  {op:12s}  ->  {stem}.png / .pdf")


def main():
    parser = argparse.ArgumentParser(description="EFT operator plots with k-factor")
    parser.add_argument("--input",    default="histos.root",         help="Path to histos.root")
    parser.add_argument("--shapes",   default=None,
                        help="Path to shapes.root (for syst band). When --variable all is used "
                             "and this is omitted, auto-detected from --datacards-dir/<variable>/shapes.root")
    parser.add_argument("--datacards-dir", default=None,
                        help="Directory containing per-variable datacard subdirs "
                             "(e.g. datacards/inc_mm). Used to auto-find shapes.root per variable.")
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
            var_meta_config[vname] = {"label": label, "unit": unit}
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
        meta = dict(VAR_META_DEFAULT.get(variable, {"label": variable, "log_x": False, "blind": False, "blind_all": False}))
        if variable in var_meta_config:
            meta["label"] = var_meta_config[variable]["label"]
            meta["unit"]  = var_meta_config[variable].get("unit", "")

        if args.variable == "all":
            outdir = os.path.join(args.outdir, variable)
        else:
            outdir = args.outdir
        os.makedirs(outdir, exist_ok=True)

        # auto-detect shapes.root per variable when --datacards-dir is given
        shapes_path = args.shapes
        if shapes_path is None and args.datacards_dir is not None:
            candidate = os.path.join(args.datacards_dir, variable, "shapes.root")
            if os.path.isfile(candidate):
                shapes_path = candidate

        print(f"\n=== {args.region} / {variable} ===")
        if variable == "triple_diff":
            plot_triple_diff(f, args.region, outdir, colors, lumi, year_label, shapes_path)
        else:
            plot_one_variable(
                f, args.region, variable, meta, outdir,
                colors, lumi, year_label,
                shapes_path, args.blind_above, args.no_blind,
            )

    print("\nDone.")


if __name__ == "__main__":
    main()
