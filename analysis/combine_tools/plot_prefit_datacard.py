#!/usr/bin/env python3
"""
Prefit plot of exactly what combine fits: the datacard (datacard.txt + shapes.root), at the SM point
(EFT processes w1/wm1/w11 have zero coefficient at c=0, so the prediction is sm + backgrounds).

Top    : stacked backgrounds + SM DY (sm), Asimov (= total prediction, what the -t -1 fits use)
Middle : ratio to the total with the prefit uncertainty: systematics (shape + lnN, each nuisance
         = max(|up|,|down|) shift of the total, summed in quadrature) and MC stat (autoMCStats input)
Bottom : relative shift of the total for +1 sigma of each nuisance (sign kept). With combine's
         vertical interpolation a bin with relative shift d goes to zero near theta ~ -1/d (linear
         extrapolation beyond 1 sigma), which is where 'Number of events is negative' / FASTEXIT appears.

Usage (dy_combine_morphing env, from the datacard dir):
  plot_prefit_datacard.py [--card datacard.txt] [--var mll] [-o prefit/prefit_mll] [--per-gev]
"""
import argparse
from optparse import OptionParser
import numpy as np
import ROOT
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
from HiggsAnalysis.CombinedLimit.DatacardParser import parseCard, addDatacardParserOptions

ROOT.gROOT.SetBatch(True)
p = argparse.ArgumentParser()
p.add_argument("--card", default="datacard.txt")
p.add_argument("--var", default="mll")
p.add_argument("-o", "--out", default=None)
p.add_argument("--per-gev", action="store_true", help="divide by bin width (mll)")
p.add_argument("--nnuis", type=int, default=8, help="nuisances shown in the bottom panel")
a = p.parse_args()
out = a.out or f"prefit_{a.var}"

op = OptionParser(); addDatacardParserOptions(op); o, _ = op.parse_args([])
DC = parseCard(open(a.card), o)
b = DC.bins[0]
EFT = ("w1_", "wm1_", "w11_")
procs = [q for q in DC.exp[b] if not q.startswith(EFT)]
fname, pattern, spattern = DC.shapeMap["*"]["*"]
f = ROOT.TFile.Open(fname)


def H(name):
    h = f.Get(name)
    if not h:
        return None
    n = h.GetNbinsX()
    return (np.array([h.GetBinContent(i) for i in range(1, n + 1)]),
            np.array([h.GetBinError(i) for i in range(1, n + 1)]),
            np.array([h.GetBinLowEdge(i) for i in range(1, n + 2)]))


nom, err = {}, {}
for q in procs:
    v, e, edges = H(pattern.replace("$PROCESS", q))
    nom[q], err[q] = v, e
tot = sum(nom.values())
mcstat = np.sqrt(sum(e**2 for e in err.values()))

# per-nuisance relative shift of the total (up, down)
shifts = {}
for name, nofloat, pdf, args, errline in DC.systs:
    up = np.zeros_like(tot); dn = np.zeros_like(tot)
    for q in procs:
        val = errline[b].get(q, 0)
        if not val:
            continue
        if pdf == "lnN":
            if isinstance(val, (list, tuple)):
                kd, ku = val
            else:
                ku, kd = val, 1 / val
            up += nom[q] * (ku - 1); dn += nom[q] * (kd - 1)
        elif pdf.startswith("shape"):
            hu = H(spattern.replace("$PROCESS", q).replace("$SYSTEMATIC", name + "Up"))
            hd = H(spattern.replace("$PROCESS", q).replace("$SYSTEMATIC", name + "Down"))
            up += (hu[0] - nom[q]) * val; dn += (hd[0] - nom[q]) * val
    shifts[name] = (up / tot, dn / tot)
syst = np.sqrt(sum(np.maximum(np.abs(u), np.abs(d))**2 for u, d in shifts.values()))
mcrel = mcstat / tot

nb = len(tot)
unrolled = a.var == "triple_diff" or edges[0] == 0 and np.allclose(np.diff(edges), 1)
x = np.arange(nb + 1) if unrolled else edges
w = np.diff(x) if (a.per_gev and not unrolled) else np.ones(nb)

hep.style.use("CMS")
fig, (ax, rx, nx) = plt.subplots(3, 1, figsize=(12 if unrolled else 10, 13), sharex=True,
                                 gridspec_kw=dict(height_ratios=[3, 1.2, 1.6], hspace=0.05))
order = sorted([q for q in procs if q != "sm"], key=lambda q: nom[q].sum())
cols = dict(zip(order, ["#94a4a2", "#832db6", "#a96b59", "#e76300", "#b9ac70", "#717581", "#92dadd", "#3f90da"]))
labels = {"sm": "DY (SM, MiNNLO k)", "DYtt": "DY#tau#tau".replace("#tau", r"$\tau$"), "Single_Top": "Single top"}
hep.histplot([nom[q] / w for q in order] + [nom["sm"] / w], x, stack=True, histtype="fill", ax=ax,
             color=[cols[q] for q in order] + ["#ffa90e"], label=[labels.get(q, q) for q in order] + [labels["sm"]])
ctr = 0.5 * (x[1:] + x[:-1])
ax.errorbar(ctr, tot / w, yerr=np.sqrt(tot) / w, fmt="ko", ms=4, label="Asimov (fitted, = total)")
ax.stairs((tot + tot * syst) / w, x, baseline=(tot - tot * syst) / w, fill=True, color="none", hatch="////",
          edgecolor="gray", lw=0, label="Syst. (prefit)")
ax.set_yscale("log")
ax.set_ylabel("Events / GeV" if (a.per_gev and not unrolled) else "Events / bin")
ax.legend(ncol=2, fontsize=14, loc="upper right")
hep.cms.label("Preliminary", data=False, lumi=59.74, year=2018, ax=ax, fontsize=18)
ax.set_title(f"prefit, {a.card.split('/')[-1]} [{b}]", fontsize=14, loc="center", pad=30)

rx.stairs(1 + syst, x, baseline=1 - syst, fill=True, color="#c0c0c0", alpha=0.8, label="syst")
rx.stairs(1 + mcrel, x, baseline=1 - mcrel, fill=True, color="#3f90da", alpha=0.6, label="MC stat")
rx.errorbar(ctr, np.ones(nb), yerr=1 / np.sqrt(tot), fmt="k.", ms=0, label="data stat (Asimov)")
rx.axhline(1, color="k", lw=0.8)
rx.set_ylabel("Unc. / total")
lim = min(max(0.02, 1.3 * np.max(syst)), 0.5)
rx.set_ylim(1 - lim, 1 + lim)
rx.legend(ncol=3, fontsize=12, loc="upper left")

big = sorted(shifts, key=lambda n: -np.max(np.abs(shifts[n][0])))[: a.nnuis]
pal = plt.cm.tab10(np.arange(len(big)))
for n, c in zip(big, pal):
    nx.stairs(100 * shifts[n][0], x, color=c, lw=1.8, label=n)
nx.axhline(0, color="k", lw=0.8)
nx.set_ylabel(r"+1$\sigma$ shift [%]")
nx.legend(ncol=4, fontsize=11, loc="upper left")
if unrolled:
    nx.set_xlabel(f"{a.var} bin (unrolled)")
else:
    nx.set_xlabel({"mll": r"$m_{\ell\ell}$ [GeV]", "rapll_abs": r"$|y_{\ell\ell}|$", "costhetastar": r"$\cos\theta^*$"}.get(a.var, a.var))
    if a.var == "mll":
        nx.set_xscale("log")
for axx in (ax, rx, nx):
    axx.set_xlim(x[0], x[-1])
fig.savefig(out + ".png", dpi=110, bbox_inches="tight")
fig.savefig(out + ".pdf", bbox_inches="tight")

print(f"{a.var}: {nb} bins, total {tot.sum():.4g} events, processes {procs}")
print(f"  {'bin':>3s} {'edges':>18s} {'total':>10s} {'syst%':>6s} {'MCstat%':>7s}  largest +1sigma shifts (%), theta where total -> 0 (linear)")
for i in range(nb):
    top = sorted(shifts, key=lambda n: -abs(shifts[n][0][i]))[:3]
    s = "  ".join(f"{n} {100 * shifts[n][0][i]:+.2f} (theta~{-1 / shifts[n][0][i]:+.0f})" for n in top if shifts[n][0][i])
    if i < 20 or not unrolled:
        print(f"  {i:3d} {edges[i]:8.4g}-{edges[i + 1]:<8.4g} {tot[i]:10.4g} {100 * syst[i]:6.2f} {100 * mcrel[i]:7.3f}  {s}")
print("wrote", out + ".png/.pdf")
