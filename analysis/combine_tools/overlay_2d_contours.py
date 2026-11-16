#!/usr/bin/env python3
"""
Overlay 68% / 95% CL contours of 2D grid scans of the same operator pair from
several variables (e.g. mll vs triple_diff).

Contours at 2*deltaNLL = 2.30 (68%) and 5.99 (95%), drawn from each scan's own
grid. Axis range: --zoom times the 95% extent of the tightest input (clipped to
the scan box), so a closed contour stays visible next to a flat band.

Usage:
  overlay_2d_contours.py --pairs cHj3_cHu,cbBRe_cHbq \
      --inputs mll=/path/mll triple_diff=/path/triple_diff -o outdir
"""
import argparse, os
import numpy as np
import ROOT
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT.gROOT.SetBatch(True)
LABELS = {"mll": r"$m_{\ell\ell}$", "triple_diff": "Triple-diff",
          "rapll_abs": r"$|y_{\ell\ell}|$", "costhetastar": r"$\cos\theta^*$"}
COLORS = ["#2166ac", "#b2182b", "#1b7837", "#762a83"]

p = argparse.ArgumentParser()
p.add_argument("--pairs", required=True)
p.add_argument("--inputs", nargs="+", required=True, help="label=dir")
p.add_argument("-o", "--outdir", default="overlay_2d")
p.add_argument("--zoom", type=float, default=3.0)
p.add_argument("--lumi", default="59.74")
a = p.parse_args()
os.makedirs(a.outdir, exist_ok=True)
inputs = [x.split("=", 1) for x in a.inputs]


def load(fn, o1, o2):
    d = ROOT.RDataFrame("limit", fn).Filter("quantileExpected > -0.5 && deltaNLL < 1e5").AsNumpy([f"k_{o1}", f"k_{o2}", "deltaNLL"])
    x, y, z = d[f"k_{o1}"], d[f"k_{o2}"], 2 * d["deltaNLL"]
    xs, ys = np.unique(np.round(x, 9)), np.unique(np.round(y, 9))
    Z = np.full((len(ys), len(xs)), np.nan)
    Z[np.searchsorted(ys, np.round(y, 9)), np.searchsorted(xs, np.round(x, 9))] = z
    Z -= np.nanmin(Z)
    return xs, ys, Z


for pair in a.pairs.split(","):
    o1, o2 = pair.split("_")
    scans = []
    for lab, d in inputs:
        fn = os.path.join(d, f"higgsCombine.{pair}.individual.MultiDimFit.mH125.root")
        if os.path.exists(fn):
            scans.append((lab, *load(fn, o1, o2)))
    if len(scans) < 2:
        print("skip", pair, "(need >=2 inputs)")
        continue

    # range: zoom x the smallest 95% extent among inputs, clipped to the box
    best = None
    for lab, xs, ys, Z in scans:
        Y, X = np.nonzero(Z <= 5.99)
        ext = (np.ptp(xs[X]) if len(X) else np.inf) * (np.ptp(ys[Y]) if len(Y) else np.inf)
        if best is None or ext < best[0]:
            best = (ext, xs[X], ys[Y], xs, ys)
    _, bx, by, xs0, ys0 = best
    cx, cy = (bx.max() + bx.min()) / 2, (by.max() + by.min()) / 2
    hx, hy = a.zoom * np.ptp(bx) / 2, a.zoom * np.ptp(by) / 2
    xlim = (max(cx - hx, min(s[1].min() for s in scans)), min(cx + hx, max(s[1].max() for s in scans)))
    ylim = (max(cy - hy, min(s[2].min() for s in scans)), min(cy + hy, max(s[2].max() for s in scans)))

    fig, ax = plt.subplots(figsize=(7, 6.5))
    handles = []
    for i, (lab, xs, ys, Z) in enumerate(scans):
        c = COLORS[i % len(COLORS)]
        ax.contour(xs, ys, Z, levels=[2.30], colors=[c], linestyles="solid", linewidths=2)
        ax.contour(xs, ys, Z, levels=[5.99], colors=[c], linestyles="dashed", linewidths=2)
        handles.append(Line2D([], [], color=c, lw=2, label=LABELS.get(lab, lab)))
    handles += [Line2D([], [], color="gray", lw=2, ls="solid", label="68% CL"),
                Line2D([], [], color="gray", lw=2, ls="dashed", label="95% CL"),
                Line2D([], [], color="black", marker="+", ls="none", ms=10, label="SM")]
    ax.plot(0, 0, "k+", ms=12, mew=2)
    ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    ax.set_xlabel(o1, fontsize=16); ax.set_ylabel(o2, fontsize=16)
    ax.legend(handles=handles, loc="upper right", fontsize=11, frameon=True)
    ax.set_title(r"$\bf{CMS}$ $\it{Preliminary}$", loc="left", fontsize=15)
    ax.set_title(rf"{a.lumi} fb$^{{-1}}$ (13 TeV)", loc="right", fontsize=13)
    ax.grid(alpha=0.3)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(a.outdir, f"overlay_{pair}.{ext}"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", pair)
