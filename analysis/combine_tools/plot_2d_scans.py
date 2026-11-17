#!/usr/bin/env python3
"""
Publication-style plots of 2D EFT grid scans (combine MultiDimFit --algo grid).

One input  -> 2*deltaNLL colour map + 68% / 95% CL contours.
Several    -> overlay of the 68% (filled + solid) and 95% (dashed) CL contours,
              one colour per variable (e.g. mll vs triple_diff).

Contours: 2*deltaNLL = 2.30 (68%) and 5.99 (95%), from the full scan grid
(no re-binning; missing points stay empty). The grid is smoothly upsampled
before contouring when it has no holes.

Usage:
  plot_2d_scans.py --pairs cHj3_cHu,clj3_cQe --inputs triple_diff=/path/td -o out/
  plot_2d_scans.py --pairs cHj3_cHu --inputs mll=/p/mll triple_diff=/p/td -o out/ --prefix overlay_
  # axis range: --frame metadata_pairs.json (use that file's per-pair box), default = grid range
"""
import argparse, json, os, re
import numpy as np
import ROOT
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap
import mplhep as hep
from scipy.ndimage import zoom, distance_transform_edt

ROOT.gROOT.SetBatch(True)
plt.style.use(hep.style.CMS)

LEVELS = (2.30, 5.99)
VAR_LABEL = {"mll": r"$m_{\ell\ell}$", "triple_diff": r"$m_{\ell\ell} \times |y_{\ell\ell}| \times \cos\theta^*$",
             "rapll_abs": r"$|y_{\ell\ell}|$", "costhetastar": r"$\cos\theta^*$"}
VAR_COLOR = {"mll": "#3f90da", "triple_diff": "#e42536", "rapll_abs": "#f89c20", "costhetastar": "#5a9a3a"}


# mkEFTScan look (HiggsAnalysis/AnalyticAnomalousCoupling utils/scan.py): white -> orange (2.30) -> blue (5.99) -> dark blue,
# 80% opacity on white; red contours (solid 1 s.d., long-dashed 2 s.d.); red cross "Expected"; axis titles k_<op>
def mkeft_cmap(zmax):
    stops = [0.0, min(2.30 / zmax, 1), min(5.99 / zmax, 1), 1.0]
    cols = [(1.0, 1.0, 1.0), (1.0, 0.549, 0.0), (0.145, 0.459, 1.0), (0.0, 0.0, 1.0)]
    cols = [tuple(0.8 * c + 0.2 for c in col) for col in cols]
    return LinearSegmentedColormap.from_list("mkeft", list(zip(stops, cols)), N=200)


LONGDASH = (0, (8, 4))


def load(fn, o1, o2):
    d = ROOT.RDataFrame("limit", fn).Filter("quantileExpected > -0.5").AsNumpy([f"k_{o1}", f"k_{o2}", "deltaNLL"])
    x, y, z = d[f"k_{o1}"], d[f"k_{o2}"], 2 * d["deltaNLL"]
    xs, ys = np.unique(np.round(x, 10)), np.unique(np.round(y, 10))
    Z = np.full((len(ys), len(xs)), np.nan)
    Z[np.searchsorted(ys, np.round(y, 10)), np.searchsorted(xs, np.round(x, 10))] = z
    Z -= np.nanmin(Z)
    # failed / missing grid points: fill with the nearest valid point (display only)
    holes = np.isnan(Z)
    if holes.any():
        idx = distance_transform_edt(holes, return_distances=False, return_indices=True)
        Z = Z[tuple(idx)]
        print(f"  {os.path.basename(fn)}: filled {int(holes.sum())} missing grid points (nearest neighbour)")
    return xs, ys, Z


def smooth(xs, ys, Z, f=4):
    if np.isnan(Z).any() or len(xs) < 4:
        return xs, ys, Z
    return (np.linspace(xs[0], xs[-1], len(xs) * f), np.linspace(ys[0], ys[-1], len(ys) * f),
            np.clip(zoom(Z, f, order=3), 0, None))


def decorate(ax, o1, o2, lumi, mcolor):
    hep.cms.label("Preliminary", data=True, lumi=lumi, com=13, ax=ax, fontsize=20)
    ax.set_xlabel(f"k_{o1}"); ax.set_ylabel(f"k_{o2}")
    ax.plot(0, 0, marker="P", color=mcolor, ms=14, ls="none", zorder=5)


p = argparse.ArgumentParser()
p.add_argument("--pairs", required=True)
p.add_argument("--inputs", nargs="+", required=True, help="label=dir")
p.add_argument("-o", "--outdir", default="plots_2d")
p.add_argument("--prefix", default="")
p.add_argument("--frame", default="", help="metadata json whose 'pairs' boxes set the axis range")
p.add_argument("--zmax", type=float, default=10.0)
p.add_argument("--lumi", type=float, default=59.74)
a = p.parse_args()
os.makedirs(a.outdir, exist_ok=True)
inputs = [x.split("=", 1) for x in a.inputs]
frames = json.load(open(a.frame)).get("pairs", {}) if a.frame else {}

for pair in a.pairs.split(","):
    o1, o2 = pair.split("_")
    scans = []
    for lab, d in inputs:
        fn = os.path.join(d, f"higgsCombine.{pair}.individual.MultiDimFit.mH125.root")
        if os.path.exists(fn):
            scans.append((lab, *load(fn, o1, o2)))
    if not scans:
        print("skip", pair); continue

    fig, ax = plt.subplots(figsize=(10, 9) if len(scans) == 1 else (10, 10))
    if len(scans) == 1:
        lab, xs, ys, Z = scans[0]
        Zc = np.clip(Z, 0, a.zmax)
        pc = ax.pcolormesh(xs, ys, Zc, cmap=mkeft_cmap(a.zmax), vmin=0, vmax=a.zmax, shading="gouraud", rasterized=True)
        cb = fig.colorbar(pc, ax=ax, pad=0.02); cb.set_label(r"$-2\ \Delta$ lnL")
        sx, sy, sZ = smooth(xs, ys, Z)
        ax.contour(sx, sy, sZ, levels=[LEVELS[0]], colors="red", linewidths=2, linestyles="solid")
        ax.contour(sx, sy, sZ, levels=[LEVELS[1]], colors="red", linewidths=2, linestyles=[LONGDASH])
        handles = [Line2D([], [], color="red", marker="P", ls="none", ms=14, label="Expected"),
                   Line2D([], [], color="red", lw=2, ls="solid", label=r"$\pm$ 1 s.d."),
                   Line2D([], [], color="red", lw=2, ls=LONGDASH, label=r"$\pm$ 2 s.d.")]
        ax.legend(handles=handles, loc="upper right", frameon=True, facecolor="white", edgecolor="none",
                  framealpha=1.0, fontsize=18)
        ax.set_xlim(xs[0], xs[-1]); ax.set_ylim(ys[0], ys[-1])
    else:
        handles = []
        for lab, xs, ys, Z in scans:
            c = VAR_COLOR.get(lab, "gray")
            sx, sy, sZ = smooth(xs, ys, Z)
            ax.contourf(sx, sy, sZ, levels=[-1, LEVELS[0]], colors=[c], alpha=0.25)
            ax.contour(sx, sy, sZ, levels=[LEVELS[0]], colors=[c], linewidths=2.5, linestyles="solid")
            ax.contour(sx, sy, sZ, levels=[LEVELS[1]], colors=[c], linewidths=2.5, linestyles=[LONGDASH])
            handles.append(Patch(facecolor=c, alpha=0.5, edgecolor=c, lw=2, label=VAR_LABEL.get(lab, lab)))
        handles += [Line2D([], [], color="black", marker="P", ls="none", ms=14, label="Expected"),
                    Line2D([], [], color="black", lw=2.5, ls="solid", label=r"$\pm$ 1 s.d."),
                    Line2D([], [], color="black", lw=2.5, ls=LONGDASH, label=r"$\pm$ 2 s.d.")]
        ax.legend(handles=handles, loc="upper right", frameon=True, facecolor="white", edgecolor="none",
                  framealpha=0.9, fontsize=16)
        xs0 = np.concatenate([s[1] for s in scans]); ys0 = np.concatenate([s[2] for s in scans])
        ax.set_xlim(xs0.min(), xs0.max()); ax.set_ylim(ys0.min(), ys0.max())
    if pair in frames:
        ax.set_xlim(*frames[pair][o1]); ax.set_ylim(*frames[pair][o2])
    decorate(ax, o1, o2, a.lumi, "red" if len(scans) == 1 else "black")
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(a.outdir, f"{a.prefix}{pair}.{ext}"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", pair)
