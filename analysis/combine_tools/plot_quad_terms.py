#!/usr/bin/env python3
"""
Plot the per-bin linear and quadratic EFT terms relative to SM, for a few operators,
comparing datacard shapes (e.g. propcorr vs no-propcorr):
  lin/sm = (w1-wm1)/2/sm,  quad/sm = ((w1+wm1)/2 - sm)/sm
Negative quad bins are marked.

Usage: plot_quad_terms.py --ops cHj3,cll1 --shapes propcorr=a/shapes.root v9=b/shapes.root -o out [--xlabel m_ll]
"""
import argparse
import numpy as np
import ROOT
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep

ROOT.gROOT.SetBatch(True)
p = argparse.ArgumentParser()
p.add_argument("--ops", required=True)
p.add_argument("--shapes", nargs="+", required=True, help="label=path/shapes.root")
p.add_argument("--xlabel", default=r"$m_{\ell\ell}$ [GeV]")
p.add_argument("-o", "--out", required=True)
a = p.parse_args()
hep.style.use("CMS")

ops = a.ops.split(",")
fig, axes = plt.subplots(2, len(ops), figsize=(7 * len(ops), 9), sharex=True, squeeze=False)
for lab_path, col in zip(a.shapes, ["#3f90da", "#e42536", "#94a4a2"]):
    lab, path = lab_path.split("=", 1)
    f = ROOT.TFile.Open(path)
    h = f.Get("histo_sm")
    n = h.GetNbinsX()
    edges = np.array([h.GetBinLowEdge(i) for i in range(1, n + 2)])
    A = lambda name: np.array([f.Get(name).GetBinContent(i) for i in range(1, n + 1)])
    sm = A("histo_sm")
    for j, o in enumerate(ops):
        w1, wm = A(f"histo_w1_{o}"), A(f"histo_wm1_{o}")
        lin, quad = (w1 - wm) / 2 / sm, ((w1 + wm) / 2 - sm) / sm
        axes[0, j].stairs(lin, edges, color=col, lw=2, label=lab)
        axes[1, j].stairs(quad, edges, color=col, lw=2, label=lab)
        neg = quad < 0
        if neg.any():
            ctr = 0.5 * (edges[1:] + edges[:-1])
            axes[1, j].scatter(ctr[neg], quad[neg], color=col, marker="v", s=80, zorder=5,
                               label=f"{lab}: quad<0 in {neg.sum()} bin(s)")
for j, o in enumerate(ops):
    axes[0, j].set_title(f"$c_{{{o[1:]}}}$ = 1", fontsize=20)
    axes[0, j].set_ylabel("lin / SM")
    axes[1, j].set_ylabel("quad / SM")
    axes[1, j].set_xlabel(a.xlabel)
    axes[1, j].set_xscale("log")
    axes[1, j].set_yscale("symlog", linthresh=1e-4)
    for ax in axes[:, j]:
        ax.axhline(0, color="k", lw=0.8)
        ax.legend(fontsize=13)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{a.out}.{ext}")
print("wrote", a.out)
