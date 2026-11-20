#!/usr/bin/env python3
"""
LHE-level check of candidate binnings for the same-binning study (mll, costhetastar, rapll_abs, triple_diff with
the same number of bins N each), compared with the current v3 binning.

Input: the propcorr LHE cache (build_cache_propcorr.py; fiducial |eta|<2.4, pT>25/10, 50<mll<3000), weights / N_gen.
Per binning and observable, one plot (unrolled for triple_diff):
  1. SM events per bin at 59.74 fb^-1 (parton level, no reco efficiency)
  2. LHE MC stat: sqrt(sum w^2)/sum w of the SM per bin, and the raw number of LHE events
  3. shape of the linear EFT term, lin/SM per bin, scaled to max |lin/SM| = 1, for a few operators
Plus a table: stat-only expected 1 sigma per operator at linear order, 1/sqrt(sum_i lin_i^2 / SM_i), per binning
(lin_i^2 corrected for the LHE MC variance of lin_i).
It ignores reco efficiency, backgrounds and systematics: use it to compare binnings, not as a limit.

Binnings are defined in BINNINGS below (v3 and the N=40 candidate, two triple_diff splits A and B).
Usage: lhe_binning_check.py [--cache ...] [-o outdir]
"""
import argparse, json, os, pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep

p = argparse.ArgumentParser()
p.add_argument("--cache", default="/grid_mnt/data__data.polcms/cms/adufour/LHE/propcorr/CACHE/lhe_cache_propcorr_parallel.pkl")
p.add_argument("-o", "--outdir", default="lhe_binning")
p.add_argument("--lumi", type=float, default=59740.0, help="pb^-1")
a = p.parse_args()
os.makedirs(a.outdir, exist_ok=True)
plt.style.use(hep.style.CMS)

NGEN = 100_000
c = pickle.load(open(a.cache, "rb"))
mll, y, cs = c["mll"], np.abs(c["rap"]), c["cstar"]
wsm = c["w_SM"] / NGEN * a.lumi
OPS = list(c["w_p1"])
lin = {op: 0.5 * (c["w_p1"][op] - c["w_m1"][op]) / NGEN * a.lumi for op in OPS}
del c
SHOW = ["cHj3", "cHe", "clj1", "ced", "clu"]

_o = np.argsort(y)
y_med = float(np.round(y[_o][np.searchsorted(np.cumsum(wsm[_o]), 0.5 * wsm.sum())] / 0.05) * 0.05)  # SM-weighted median
V3_MLL = [50, 60, 80, 100, 120, 140, 180, 220, 270, 350, 500, 700, 1000, 1500, 3000]
# 40 bins in the spirit of eft_bkg_fullsyst_v5 (5 GeV below the peak, 2 GeV over the Z peak, wide at high mll)
MLL40 = [50, 55, 60, 65, 70, 75, 78, *range(80, 107, 2), 110, 115, 120, 130, 140, 150, 175, 200, 250, 300, 350, 400,
         500, 600, 700, 800, 1000, 1200, 1500, 3000]
BINNINGS = {
    "v3": {"mll": V3_MLL, "costhetastar": [-1.0, -0.6, -0.2, 0.2, 0.6, 1.0], "rapll_abs": [0.0, 0.48, 0.96, 1.44, 2.4],
           "triple_diff": (V3_MLL, [-1.0, -0.6, -0.2, 0.2, 0.6, 1.0], [0.0, 0.48, 0.96, 1.44, 2.4])},
    "N40": {"mll": MLL40,
            "costhetastar": list(np.round(np.linspace(-1, 1, 41), 4)),
            "rapll_abs": list(np.round(np.linspace(0, 2.4, 41), 4)),
            # A: 10 mll x 2 cos (backward/forward) x 2 |y| (split at the SM median)
            "triple_diff": ([50, 70, 80, 86, 90, 94, 100, 120, 200, 500, 3000], [-1.0, 0.0, 1.0], [0.0, y_med, 2.4]),
            # B: 5 mll x 4 cos x 2 |y|
            "triple_diff_B": ([50, 80, 100, 200, 500, 3000], [-1.0, -0.5, 0.0, 0.5, 1.0], [0.0, y_med, 2.4]),
            # C: like A, fewer peak bins, more at high mll
            "triple_diff_C": ([50, 76, 86, 96, 106, 150, 250, 400, 600, 1000, 3000], [-1.0, 0.0, 1.0], [0.0, y_med, 2.4])},
}
for b in BINNINGS.values():
    for v, e in b.items():
        n = np.prod([len(x) - 1 for x in e]) if v.startswith("triple_diff") else len(e) - 1
        assert b is BINNINGS["v3"] or n == 40, (v, n)
VAL = {"mll": mll, "costhetastar": cs, "rapll_abs": y}
XLAB = {"mll": r"$m_{\ell\ell}$ [GeV]", "costhetastar": r"$\cos\theta^*$", "rapll_abs": r"$|y_{\ell\ell}|$",
        "triple_diff": "unrolled bin ($m_{\\ell\\ell}$ inner, then $\\cos\\theta^*$, $|y_{\\ell\\ell}|$ outer; as hist_unroll)"}


def binidx(var, edges):
    if not var.startswith("triple_diff"):
        e = np.asarray(edges, float)
        i = np.digitize(VAL[var], e) - 1
        return i, len(e) - 1, e
    em, ec, ey = (np.asarray(x, float) for x in edges)
    im, ic, iy = np.digitize(mll, em) - 1, np.digitize(cs, ec) - 1, np.digitize(y, ey) - 1
    ok = (im >= 0) & (im < len(em) - 1) & (ic >= 0) & (ic < len(ec) - 1) & (iy >= 0) & (iy < len(ey) - 1)
    # same order as hist_unroll in post_process_minnlo.py (view().T.flatten()): mll fastest, then cos, |y| outermost
    i = np.where(ok, (iy * (len(ec) - 1) + ic) * (len(em) - 1) + im, -1)
    return i, (len(em) - 1) * (len(ec) - 1) * (len(ey) - 1), None


def hsum(i, n, w):
    ok = (i >= 0) & (i < n)
    return np.bincount(i[ok], weights=w[ok], minlength=n)


summary = {}
for bname, bdef in BINNINGS.items():
    for var, edges in bdef.items():
        i, n, e = binidx(var, edges)
        sm = hsum(i, n, wsm)
        sw2 = hsum(i, n, wsm ** 2)
        nraw = hsum(i, n, np.ones_like(wsm))
        mcrel = np.sqrt(sw2) / np.where(sm > 0, sm, 1)
        L = {op: hsum(i, n, lin[op]) for op in OPS}
        # MC fluctuations of lin add to lin^2 on average: subtract the MC variance sum(w_lin^2) per bin (unbiased)
        LV = {op: hsum(i, n, lin[op] ** 2) for op in OPS}
        sig = {op: 1 / np.sqrt(max(np.sum(np.clip(L[op] ** 2 - LV[op], 0, None) / np.where(sm > 0, sm, np.inf)), 1e-300))
               for op in OPS}
        sig_raw = {op: 1 / np.sqrt(np.sum(L[op] ** 2 / np.where(sm > 0, sm, np.inf))) for op in OPS}
        summary[(bname, var)] = dict(nbins=n, sig=sig, sig_raw=sig_raw, mc_max=float(mcrel.max()), nraw_min=int(nraw.min()),
                                     sm_min=float(sm.min()))

        fig, ax = plt.subplots(3, 1, figsize=(13, 13), sharex=True, gridspec_kw=dict(height_ratios=[2, 1.2, 1.6], hspace=0.06))
        x = np.arange(n + 1) if e is None else e
        hep.histplot(sm, x, ax=ax[0], color="k", lw=2, label="SM (LHE, fiducial)")
        ax[0].set_yscale("log"); ax[0].set_ylabel("Events / bin")
        ax[0].legend(loc="lower left" if var.startswith("triple_diff") else "upper right", fontsize=16)
        ax[0].set_title(f"{var}: {bname} binning, {n} bins", fontsize=20, loc="left")
        hep.histplot(100 * mcrel, x, ax=ax[1], color="tab:blue", lw=2, label="LHE MC stat of SM [%]")
        ax[1].axhline(1, color="gray", ls=":")
        ax[1].set_ylabel("MC stat [%]"); ax[1].set_ylim(0, max(2, 1.2 * 100 * mcrel.max()))
        a2 = ax[1].twinx(); a2.step(x, np.append(nraw, nraw[-1]), where="post", color="tab:orange", lw=1.2)
        a2.set_yscale("log"); a2.set_ylabel("LHE events", color="tab:orange", fontsize=16); a2.tick_params(axis="y", colors="tab:orange", labelsize=12)
        for op in SHOW:
            r = L[op] / np.where(sm > 0, sm, 1)
            hep.histplot(r / np.abs(r).max(), x, ax=ax[2], lw=2, label=f"{op} (max |lin/SM| = {np.abs(r).max():.2g})")
        ax[2].axhline(0, color="gray", lw=0.8)
        ax[2].set_ylabel("lin/SM, scaled"); ax[2].set_ylim(-1.15, 1.6)
        ax[2].legend(ncol=2, fontsize=12, loc="upper left")
        ax[2].set_xlabel(XLAB["triple_diff" if var.startswith("triple_diff") else var])
        if var == "mll":
            for aa in ax: aa.set_xscale("log")
        if var.startswith("triple_diff"):
            nmb, nc, ny = (len(x) - 1 for x in bdef[var])
            ec, ey = bdef[var][1], bdef[var][2]
            for k in range(1, nc * ny):
                for aa in ax: aa.axvline(k * nmb, color="gray", lw=0.8 if k % nc else 1.6, ls="--" if k % nc else "-")
            for k in range(nc * ny):
                ic, iy = k % nc, k // nc
                ax[0].text((k + 0.5) * nmb, 0.97, f"cos [{ec[ic]:g},{ec[ic+1]:g}]\n|y| [{ey[iy]:g},{ey[iy+1]:g}]",
                           transform=ax[0].get_xaxis_transform(), ha="center", va="top", fontsize=9 if nc * ny <= 8 else 6)
        fig.savefig(f"{a.outdir}/{bname}_{var}.png", dpi=100, bbox_inches="tight")
        plt.close(fig)
        print("wrote", bname, var, n, "bins; MC stat max", f"{100*mcrel.max():.2f}%", "min LHE events", int(nraw.min()))

# ---- table: expected stat-only 1 sigma (linear) per operator, per binning
bnames = list(BINNINGS); vars_ = ["mll", "costhetastar", "rapll_abs", "triple_diff", "triple_diff_B", "triple_diff_C"]
lines = ["| op | " + " | ".join(f"{v} {b} ({summary[(b, v if (b, v) in summary else 'triple_diff')]['nbins']})" for v in vars_ for b in bnames) + " |",
         "|---|" + "---|" * (len(vars_) * len(bnames))]
for op in OPS:
    lines.append(f"| {op} | " + " | ".join(f"{summary[(b, v if (b, v) in summary else 'triple_diff')]['sig'][op]:.3g}" for v in vars_ for b in bnames) + " |")
lines.append("")
lines.append("ratio N/v3 of the expected 1 sigma (>1 = the new binning loses information):")
lines.append("| op | " + " | ".join(vars_) + " |"); lines.append("|---|" + "---|" * len(vars_))
for op in OPS:
    lines.append(f"| {op} | " + " | ".join(f"{summary[(bnames[1], v)]['sig'][op] / summary[(bnames[0], v if (bnames[0], v) in summary else 'triple_diff')]['sig'][op]:.2f}" for v in vars_) + " |")
lines.append("")
lines.append("MC-noise bias: ratio (raw, no variance subtraction) / (corrected) of the 1 sigma, N binning:")
lines.append("| op | " + " | ".join(vars_) + " |"); lines.append("|---|" + "---|" * len(vars_))
for op in OPS:
    lines.append(f"| {op} | " + " | ".join(f"{summary[(bnames[1], v)]['sig_raw'][op] / summary[(bnames[1], v)]['sig'][op]:.2f}" for v in vars_) + " |")
open(f"{a.outdir}/sigma_table.md", "w").write("\n".join(lines) + "\n")
json.dump({b: {v: (bdef[v] if v != "triple_diff" else list(bdef[v])) for v in bdef} for b, bdef in BINNINGS.items()},
          open(f"{a.outdir}/binnings.json", "w"), indent=1, default=float)
print("\n".join(lines))
