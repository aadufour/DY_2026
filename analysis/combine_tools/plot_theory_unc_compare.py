#!/usr/bin/env python3
"""
plot_theory_unc_compare.py

Compare the theory uncertainties on the EFT templates between two datacard sets,
e.g. MG5 LO variations (old apply_kfactor) vs MiNNLO variations (post_process_minnlo).
For each nuisance: relative Up/Down variation (%) of one process per bin.

Usage (dy_analysis env):
    python3 plot_theory_unc_compare.py \
        --dir1 .../propcorr_v1/datacards_single --label1 "MG5 LO" \
        --dir2 .../propcorr_v1_minnlo/datacards  --label2 "MiNNLO" \
        --var mll [--proc sm] [--outdir .]
"""
import argparse
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep

plt.style.use(hep.style.CMS)

p = argparse.ArgumentParser()
p.add_argument("--dir1", required=True)
p.add_argument("--dir2", required=True)
p.add_argument("--label1", default="MG5 LO")
p.add_argument("--label2", default="MiNNLO")
p.add_argument("--var", default="mll")
p.add_argument("--region", default="inc_mm")
p.add_argument("--proc", default="sm")
p.add_argument("--nuis", nargs="+", default=["QCDScale", "PDFweight", "alphaS"])
p.add_argument("--outdir", default=".")
a = p.parse_args()

XLABEL = {"mll": r"$m_{\ell\ell}$ [GeV]", "rapll_abs": r"$|y_{\ell\ell}|$",
          "costhetastar": r"$\cos\theta^*$", "triple_diff": "bin index"}
NAME = {"QCDScale": r"QCD scale ($\mu_R$, $\mu_F$)", "PDFweight": "PDF", "alphaS": r"$\alpha_S$"}


def rel(path, nuis):
    f = uproot.open(f"{path}/{a.region}/{a.var}/shapes.root")
    h = f[f"histo_{a.proc}"]
    n = h.values()
    safe = np.where(n != 0, n, 1.0)
    up = 100 * (f[f"histo_{a.proc}_{nuis}Up"].values() / safe - 1)
    dn = 100 * (f[f"histo_{a.proc}_{nuis}Down"].values() / safe - 1)
    return h.axis().edges(), up, dn


fig, axes = plt.subplots(len(a.nuis), 1, figsize=(12, 4.5 * len(a.nuis)), sharex=True)
axes = np.atleast_1d(axes)
for ax, nuis in zip(axes, a.nuis):
    for d, lab, col in [(a.dir1, a.label1, "tab:blue"), (a.dir2, a.label2, "tab:red")]:
        edges, up, dn = rel(d, nuis)
        if a.var == "triple_diff":
            edges = np.arange(len(up) + 1)
        hep.histplot(up, edges, ax=ax, color=col, label=f"{lab} up")
        hep.histplot(dn, edges, ax=ax, color=col, ls="--", label=f"{lab} down")
    ax.axhline(0, color="gray", lw=0.8)
    ax.set_ylabel("variation [%]")
    ax.set_title(NAME.get(nuis, nuis), fontsize=20, loc="left")
    ax.legend(fontsize=14, ncol=2)
    if a.var == "mll":
        ax.set_xscale("log")
axes[-1].set_xlabel(XLABEL.get(a.var, a.var))
hep.cms.label("Simulation Preliminary", data=False, ax=axes[0], fontsize=20)
fig.tight_layout()
out = f"{a.outdir}/theory_unc_compare_{a.proc}_{a.var}"
for ext in ("pdf", "png"):
    fig.savefig(f"{out}.{ext}")
print("saved", out)
