#!/usr/bin/env python3
"""
plot_eft_propcorr.py

Same as plot_eft.py, but overlays the propagator-corrected (propcorr) EFT
signal on top of the baseline one, so both c=+-1 curves are visible together
for every operator.

Backgrounds, data, and the SM (MiNNLO) reference are read only from the
baseline histos.root: config_v9.py and config_propcorr_v1.py are identical
in backgrounds/data/nuisances/binning (verified by diff) -- the *only*
difference between the two configs is which DY-SMEFTsim production feeds
the sm/w1_op/wm1_op templates. So there is exactly one background stack, one
SM line, one systematic band, and one Data overlay; only the EFT c=+-1
curves are drawn twice (solid = baseline, dashed = propcorr).

Sources:
  --bl-input   histos.root from the baseline config (also supplies backgrounds/data)
  --pc-input   histos.root from the propcorr config (only sm/w1_op/wm1_op are used)
  --shapes     shapes.root (systematic band on SM prediction, from spritz-cards, optional)

Run from a config dir, with analysis_venv active:
    python3 .../plot_eft_propcorr.py --bl-input /path/to/eft_bkg_fullsyst_v9/histos.root \\
                                      --pc-input /path/to/propcorr_v1/histos.root \\
                                      [--region inc_mm] [--variable mll|costhetastar|rapll_abs|triple_diff|all]
                                      [--outdir plots/eft_operators_propcorr]

--variable all  loops over every variable found in the baseline histos.root for the given region.
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
LEGEND_FONTSIZE  = 12
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

COLOR_P1 = "crimson"
COLOR_M1 = "steelblue"

# ============================================================

MLL_EDGES   = np.array([40, 60, 80, 100, 120, 140, 180, 220, 270, 350, 500, 700, 1000, 1500, 3000], dtype=float)
COSTH_EDGES = np.array([-1, -0.6, -0.2, 0.2, 0.6, 1], dtype=float)
RAPLL_EDGES = np.array([0, 0.48, 0.96, 1.44, 2.4], dtype=float)
N_MLL   = len(MLL_EDGES)   - 1
N_COSTH = len(COSTH_EDGES) - 1
N_RAPLL = len(RAPLL_EDGES) - 1


def _td_slice(irapll, icos):
    start = irapll * N_COSTH * N_MLL + icos * N_MLL
    return slice(start, start + N_MLL)


VAR_META_DEFAULT = {
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

BKG_NUISANCES = [
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


def compute_syst_band(shapes_path, bkg_vals_keys, n_bins, get_shape_fn):
    """Quadrature-sum up/down background systematic band from shapes.root."""
    syst_up   = np.zeros(n_bins)
    syst_down = np.zeros(n_bins)
    if shapes_path is None:
        return syst_up, syst_down
    try:
        sm_samples_shapes = list(bkg_vals_keys) + ["DYll"]
        nom_total = np.zeros(n_bins)
        for s in sm_samples_shapes:
            try:
                nom_total += get_shape_fn(s.replace(" ", "_"))
            except Exception:
                pass

        for nuis_name, affected in BKG_NUISANCES:
            tot_up   = nom_total.copy()
            tot_down = nom_total.copy()
            for s in affected:
                skey = s.replace(" ", "_")
                try:
                    nom_s = get_shape_fn(skey)
                    if nuis_name == "lumi":
                        tot_up   += (1.0084 - 1.0) * nom_s
                        tot_down += (1.0 / 1.0084 - 1.0) * nom_s
                    else:
                        up_s = get_shape_fn(f"{skey}_{nuis_name}Up")
                        do_s = get_shape_fn(f"{skey}_{nuis_name}Down")
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
        print(f"  [warn] Could not compute syst band from {shapes_path}: {e}")
    return syst_up, syst_down


def plot_one_variable(
    f_bl, f_pc, region, variable, var_meta, outdir,
    colors, lumi, year_label,
    shapes_path, blind_above, no_blind,
):
    """Produce one PNG+PDF per operator, overlaying baseline and propcorr EFT curves."""
    directory_bl = f_bl[f"{region}/{variable}"]
    directory_pc = f_pc[f"{region}/{variable}"]

    edges   = get_edges(directory_bl, "sm")
    widths  = edges[1:] - edges[:-1]
    centers = 0.5 * (edges[:-1] + edges[1:])

    xlabel  = var_meta["label"]
    log_x   = var_meta["log_x"]
    do_blind  = var_meta["blind"] and not no_blind
    blind_all = var_meta.get("blind_all", False)

    bkg_vals = {}
    for s in BKG_STACK:
        try:
            bkg_vals[s] = get_vals(directory_bl, s)
        except Exception:
            pass

    dyll     = get_vals(directory_bl, "DYll")
    data     = get_vals(directory_bl, "Data")
    data_var = get_variances(directory_bl, "Data")

    present = [s for s in BKG_STACK if s in bkg_vals]
    stack   = np.array([bkg_vals[s] for s in present])
    cumsum  = np.cumsum(stack, axis=0)
    bkg_total = cumsum[-1]
    sm_total = bkg_total + dyll

    fs = uproot.open(shapes_path) if shapes_path is not None else None
    syst_up_sm, syst_down_sm = compute_syst_band(
        shapes_path, bkg_vals.keys(), len(bkg_total),
        lambda name: fs[f"histo_{name}"].values().copy(),
    ) if fs is not None else (np.zeros_like(bkg_total), np.zeros_like(bkg_total))

    style = deepcopy(hep.style.CMS)
    style["font.size"]        = FONT_SIZE
    style["axes.labelsize"]   = LABEL_SIZE
    style["xtick.labelsize"]  = TICK_LABELSIZE
    style["ytick.labelsize"]  = TICK_LABELSIZE
    style["legend.fontsize"]  = LEGEND_FONTSIZE
    style["figure.figsize"]   = FIG_SIZE
    plt.style.use(style)

    if not do_blind:
        blind_mask = np.zeros(len(centers), dtype=bool)
    elif blind_all:
        blind_mask = np.ones(len(centers), dtype=bool)
    else:
        blind_mask = centers > blind_above

    data_unc = np.sqrt(np.abs(data_var))

    for op in OPERATORS:
        try:
            w1_bl  = get_vals(directory_bl, f"w1_{op}")
            wm1_bl = get_vals(directory_bl, f"wm1_{op}")
        except Exception:
            print(f"  [skip] {op}: baseline histo_w1_{op} / histo_wm1_{op} not found")
            continue
        try:
            w1_pc  = get_vals(directory_pc, f"w1_{op}")
            wm1_pc = get_vals(directory_pc, f"wm1_{op}")
        except Exception:
            print(f"  [skip] {op}: propcorr histo_w1_{op} / histo_wm1_{op} not found")
            continue

        eft_p1_bl = bkg_total + w1_bl
        eft_m1_bl = bkg_total + wm1_bl
        eft_p1_pc = bkg_total + w1_pc
        eft_m1_pc = bkg_total + wm1_pc

        fig, (ax_top, ax_bot) = plt.subplots(
            2, 1, sharex=True,
            gridspec_kw={"height_ratios": [3, 1]},
            dpi=FIG_DPI,
        )
        fig.tight_layout(pad=-0.5)
        hep.cms.label("Preliminary", data=False, lumi=round(lumi, 2), ax=ax_top, year=year_label)

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

        ax_top.stairs(eft_p1_bl / widths, edges=edges, color=COLOR_P1, linewidth=LINE_WIDTH_EFT,
                      linestyle="solid", label=f"{op} baseline (c=+1)", fill=False, zorder=3)
        ax_top.stairs(eft_m1_bl / widths, edges=edges, color=COLOR_M1, linewidth=LINE_WIDTH_EFT,
                      linestyle="solid", label=f"{op} baseline (c=-1)", fill=False, zorder=3)
        ax_top.stairs(eft_p1_pc / widths, edges=edges, color=COLOR_P1, linewidth=LINE_WIDTH_EFT,
                      linestyle="dashed", label=f"{op} propcorr (c=+1)", fill=False, zorder=3)
        ax_top.stairs(eft_m1_pc / widths, edges=edges, color=COLOR_M1, linewidth=LINE_WIDTH_EFT,
                      linestyle="dashed", label=f"{op} propcorr (c=-1)", fill=False, zorder=3)

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
            np.nanmax(eft_p1_bl / widths), np.nanmax(eft_m1_bl / widths),
            np.nanmax(eft_p1_pc / widths), np.nanmax(eft_m1_pc / widths),
        )
        pos_vals_list = [v[v > 0] / widths[v > 0] for v in stack if np.any(v > 0)]
        pos_vals = np.concatenate(pos_vals_list) if pos_vals_list else np.array([1e-4])
        ymin = max(1e-4, 0.3 * np.min(pos_vals))

        ax_top.set_yscale("log")
        ax_top.set_ylim(ymin, ymax * 3e3)
        ax_top.set_ylabel("Events / GeV" if var_meta.get("unit") else "Events")
        ax_top.tick_params(labelbottom=False)
        ax_top.legend(loc=LEGEND_LOC_TOP, ncols=LEGEND_NCOLS_TOP, framealpha=0.8)

        denom = np.where(sm_total > 0, sm_total, 1e-30)
        ratio_p1_bl = eft_p1_bl / denom
        ratio_m1_bl = eft_m1_bl / denom
        ratio_p1_pc = eft_p1_pc / denom
        ratio_m1_pc = eft_m1_pc / denom
        ratio_data  = data / denom

        ax_bot.stairs(ratio_p1_bl, edges=edges, color=COLOR_P1, linewidth=1.2, linestyle="solid",  label="baseline c=+1 / SM")
        ax_bot.stairs(ratio_m1_bl, edges=edges, color=COLOR_M1, linewidth=1.2, linestyle="solid",  label="baseline c=-1 / SM")
        ax_bot.stairs(ratio_p1_pc, edges=edges, color=COLOR_P1, linewidth=1.2, linestyle="dashed", label="propcorr c=+1 / SM")
        ax_bot.stairs(ratio_m1_pc, edges=edges, color=COLOR_M1, linewidth=1.2, linestyle="dashed", label="propcorr c=-1 / SM")
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
        candidates = [ratio_p1_bl, ratio_m1_bl, ratio_p1_pc, ratio_m1_pc]
        range_mask = centers <= (var_meta.get("range_max") or np.inf)
        all_finite = np.concatenate([a[range_mask][np.isfinite(a[range_mask])] for a in candidates])
        half = max(np.max(np.abs(all_finite - 1.0)) * 1.2, 0.05) if all_finite.size else 0.3
        ax_bot.set_ylim(1.0 - half, 1.0 + half)
        ax_bot.set_xlabel(xlabel)
        if log_x:
            ax_bot.set_xscale("log")
        ax_bot.set_xlim(edges[0], edges[-1])
        ax_bot.legend(loc=LEGEND_LOC_BOT, framealpha=0.8, ncols=2, fontsize=LEGEND_FONTSIZE - 2)

        stem = os.path.join(outdir, f"eft_propcorr_{op}")
        for ext in ("png", "pdf"):
            fig.savefig(f"{stem}.{ext}", facecolor="white", pad_inches=0.1, bbox_inches="tight")
        plt.close()
        print(f"  {op:12s}  ->  {stem}.png / .pdf")


def plot_triple_diff(f_bl, f_pc, region, outdir, colors, lumi, year_label, shapes_path):
    """Multi-panel EFT plot for the triple_diff unrolled histogram, baseline vs propcorr."""
    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

    directory_bl = f_bl[f"{region}/triple_diff"]
    directory_pc = f_pc[f"{region}/triple_diff"]
    mll_widths = MLL_EDGES[1:] - MLL_EDGES[:-1]
    mll_centers = 0.5 * (MLL_EDGES[:-1] + MLL_EDGES[1:])

    def _get(directory, name):
        return directory[f"histo_{name}"].values().copy()

    def _getvar(directory, name):
        return directory[f"histo_{name}"].variances().copy()

    try:
        dyll     = _get(directory_bl, "DYll")
        data     = _get(directory_bl, "Data")
        data_var = _getvar(directory_bl, "Data")
    except Exception as e:
        print(f"  [skip] triple_diff: cannot read DYll/Data: {e}")
        return

    bkg_vals = {}
    for s in BKG_STACK:
        try:
            bkg_vals[s] = _get(directory_bl, s)
        except Exception:
            pass
    present  = [s for s in BKG_STACK if s in bkg_vals]
    n_td = N_MLL * N_COSTH * N_RAPLL
    stack    = np.array([bkg_vals[s] for s in present])
    cumsum   = np.cumsum(stack, axis=0)
    bkg_total = cumsum[-1]
    sm_total  = bkg_total + dyll

    fs = uproot.open(shapes_path) if shapes_path is not None else None
    syst_up, syst_down = compute_syst_band(
        shapes_path, bkg_vals.keys(), n_td,
        lambda name: fs[f"histo_{name}"].values().copy(),
    ) if fs is not None else (np.zeros(n_td), np.zeros(n_td))

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
            w1_bl  = _get(directory_bl, f"w1_{op}")
            wm1_bl = _get(directory_bl, f"wm1_{op}")
        except Exception:
            print(f"  [skip] triple_diff {op}: baseline w1/wm1 not found")
            continue
        try:
            w1_pc  = _get(directory_pc, f"w1_{op}")
            wm1_pc = _get(directory_pc, f"wm1_{op}")
        except Exception:
            print(f"  [skip] triple_diff {op}: propcorr w1/wm1 not found")
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
                sm_total_sl = sm_total[sl]
                data_sl     = data[sl]
                data_unc_sl = np.sqrt(np.abs(data_var[sl]))
                syst_up_sl  = syst_up[sl]
                syst_dn_sl  = syst_down[sl]

                eft_p1_bl_sl = bkg_tot_sl + w1_bl[sl]
                eft_m1_bl_sl = bkg_tot_sl + wm1_bl[sl]
                eft_p1_pc_sl = bkg_tot_sl + w1_pc[sl]
                eft_m1_pc_sl = bkg_tot_sl + wm1_pc[sl]

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

                ax_top.stairs(eft_p1_bl_sl / mll_widths, edges=MLL_EDGES, color=COLOR_P1, linewidth=0.8,
                              linestyle="solid", label=f"{op} bl c=+1" if is_first else "_nolegend_", fill=False, zorder=3)
                ax_top.stairs(eft_m1_bl_sl / mll_widths, edges=MLL_EDGES, color=COLOR_M1, linewidth=0.8,
                              linestyle="solid", label=f"{op} bl c=-1" if is_first else "_nolegend_", fill=False, zorder=3)
                ax_top.stairs(eft_p1_pc_sl / mll_widths, edges=MLL_EDGES, color=COLOR_P1, linewidth=0.8,
                              linestyle="dashed", label=f"{op} pc c=+1" if is_first else "_nolegend_", fill=False, zorder=3)
                ax_top.stairs(eft_m1_pc_sl / mll_widths, edges=MLL_EDGES, color=COLOR_M1, linewidth=0.8,
                              linestyle="dashed", label=f"{op} pc c=-1" if is_first else "_nolegend_", fill=False, zorder=3)

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
                        "Preliminary", data=False, lumi=round(lumi, 2),
                        ax=ax_top, year=year_label, fontsize=6,
                    )

                if icos == 0:
                    ax_top.set_ylabel("Events/GeV", fontsize=5)

                denom      = np.where(sm_total_sl > 0, sm_total_sl, 1e-30)
                ratio_p1_bl = eft_p1_bl_sl / denom
                ratio_m1_bl = eft_m1_bl_sl / denom
                ratio_p1_pc = eft_p1_pc_sl / denom
                ratio_m1_pc = eft_m1_pc_sl / denom
                ratio_data  = data_sl / denom
                data_blind_mask = mll_centers > 500
                ratio_data_plot = np.where(data_blind_mask, np.nan, ratio_data)
                data_unc_plot   = np.where(data_blind_mask, np.nan, data_unc_sl / denom)

                ax_bot.stairs(ratio_p1_bl, edges=MLL_EDGES, color=COLOR_P1, linewidth=0.8, linestyle="solid")
                ax_bot.stairs(ratio_m1_bl, edges=MLL_EDGES, color=COLOR_M1, linewidth=0.8, linestyle="solid")
                ax_bot.stairs(ratio_p1_pc, edges=MLL_EDGES, color=COLOR_P1, linewidth=0.8, linestyle="dashed")
                ax_bot.stairs(ratio_m1_pc, edges=MLL_EDGES, color=COLOR_M1, linewidth=0.8, linestyle="dashed")
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

                candidates = [ratio_p1_bl, ratio_m1_bl, ratio_p1_pc, ratio_m1_pc]
                finite_r = np.concatenate([a[np.isfinite(a)] for a in candidates])
                half = max(np.max(np.abs(finite_r - 1.0)) * 1.2, 0.05) if finite_r.size else 0.3
                ax_bot.set_ylim(1.0 - half, 1.0 + half)
                ax_bot.set_xlim(MLL_EDGES[0], MLL_EDGES[-1])

                if irapll == N_RAPLL - 1:
                    ax_bot.set_xlabel(r"$m_{\ell\ell}$ [GeV]", fontsize=5)
                if icos == 0:
                    ax_bot.set_ylabel("Ratio", fontsize=5)

        fig.suptitle(f"EFT {op} (c=±1): baseline (solid) vs propcorr (dashed)", fontsize=11, y=0.93)

        legend_handles = (
            [Patch(color=colors.get(s, "grey"), label=s) for s in present]
            + [
                Line2D([0], [0], color=colors.get("DYll", DEFAULT_COLORS["DYll"]), linestyle="dashed", linewidth=0.8, label="SM (MiNNLO)"),
                Line2D([0], [0], color=COLOR_P1, linestyle="solid",  linewidth=0.8, label=f"{op} baseline c=+1"),
                Line2D([0], [0], color=COLOR_M1, linestyle="solid",  linewidth=0.8, label=f"{op} baseline c=-1"),
                Line2D([0], [0], color=COLOR_P1, linestyle="dashed", linewidth=0.8, label=f"{op} propcorr c=+1"),
                Line2D([0], [0], color=COLOR_M1, linestyle="dashed", linewidth=0.8, label=f"{op} propcorr c=-1"),
                Line2D([0], [0], color="black", marker="o", markersize=2, linewidth=0.6, label="Data [mll < 500 GeV]"),
            ]
        )
        fig.legend(
            handles=legend_handles,
            fontsize=7, ncols=4, framealpha=0.7,
            loc="lower center", bbox_to_anchor=(0.7, 0.9),
            handlelength=1.2, handletextpad=0.4, columnspacing=0.8,
        )

        stem = os.path.join(outdir, f"eft_triple_diff_propcorr_{op}")
        for ext in ("png", "pdf"):
            fig.savefig(
                f"{stem}.{ext}", facecolor="white",
                pad_inches=0.05, bbox_inches="tight", dpi=FIG_DPI,
            )
        plt.close()
        print(f"  {op:12s}  ->  {stem}.png / .pdf")


def main():
    parser = argparse.ArgumentParser(description="EFT operator plots comparing baseline vs propcorr")
    parser.add_argument("--bl-input", dest="bl_input", required=True,
                        help="Path to baseline histos.root (also supplies backgrounds/data)")
    parser.add_argument("--pc-input", dest="pc_input", required=True,
                        help="Path to propcorr histos.root (only sm/w1_op/wm1_op are used)")
    parser.add_argument("--shapes",   default=None,
                        help="Path to shapes.root (for syst band). When --variable all is used "
                             "and this is omitted, auto-detected from --datacards-dir/<variable>/shapes.root")
    parser.add_argument("--datacards-dir", default=None,
                        help="Directory containing per-variable datacard subdirs "
                             "(e.g. datacards/inc_mm), from the baseline config. Used to "
                             "auto-find shapes.root per variable.")
    parser.add_argument("--region",   default="inc_mm",              help="Region key in histos.root")
    parser.add_argument("--variable", default="mll",
                        help="Variable key in histos.root, or 'all' to loop over every variable in the region")
    parser.add_argument("--outdir",   default="plots/eft_operators_propcorr",  help="Output directory for PNGs")
    parser.add_argument("--blind-above", type=float, default=500.,
                        help="Blind data above this value (only for variables with blinding enabled, default: 500 GeV for mll)")
    parser.add_argument("--no-blind", action="store_true",
                        help="Disable blinding for all variables")
    args = parser.parse_args()

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

    f_bl = uproot.open(args.bl_input)
    f_pc = uproot.open(args.pc_input)

    if args.variable == "all":
        region_dir = f_bl[args.region]
        variables = [
            k.split(";")[0] for k, cls in region_dir.classnames().items()
            if "TDirectory" in cls and "/" not in k
        ]
    else:
        variables = [args.variable]

    for variable in variables:
        meta = dict(VAR_META_DEFAULT.get(variable, {"label": variable, "log_x": False, "blind": False, "blind_all": False}))
        if variable in var_meta_config:
            meta["label"] = var_meta_config[variable]["label"]
            meta["unit"]  = var_meta_config[variable].get("unit", "")

        if args.variable == "all":
            outdir = os.path.join(args.outdir, variable)
        else:
            outdir = args.outdir
        os.makedirs(outdir, exist_ok=True)

        shapes_path = args.shapes
        if shapes_path is None and args.datacards_dir is not None:
            candidate = os.path.join(args.datacards_dir, variable, "shapes.root")
            if os.path.isfile(candidate):
                shapes_path = candidate

        print(f"\n=== {args.region} / {variable} ===")
        if variable == "triple_diff":
            plot_triple_diff(f_bl, f_pc, args.region, outdir, colors, lumi, year_label, shapes_path)
        else:
            plot_one_variable(
                f_bl, f_pc, args.region, variable, meta, outdir,
                colors, lumi, year_label,
                shapes_path, args.blind_above, args.no_blind,
            )

    print("\nDone.")


if __name__ == "__main__":
    main()
