#!/usr/bin/env python3
"""
Linearised Asimov Fisher check of a 1D shapes.root (spritz-cards-eft layout): sigma(k_op) stat-only and profiled over the
13 nuisances (unit Gaussian priors, lumi lnN), reproducing combine's 1D fits to ~20-40%.
Used to show (2026-11-21, notes/binning_study.md) that the N80 stat+syst gains come mostly from MC noise of the EFT
linear templates (the fit reads bin-to-bin noise as shape once the nuisances have absorbed the smooth part).

Usage: fisher_template_check.py <shapes.root> <op> <lumi lnN, e.g. 1.0084> [group G] [smoothing none|eft|all]
  G        : merge G adjacent bins (80 -> 80/G)
  smoothing: none; eft = EFT linear template replaced by a deg-6 polynomial fit of lin/nominal (at full binning,
             before merging); all = also every nuisance variation
"""
import sys, ROOT, numpy as np
ROOT.gROOT.SetBatch(True)
# linearised Asimov Fisher: parameters (k_op, 13 nuisances with unit Gaussian priors); bins: Poisson ~ Gaussian with var = nominal
f = ROOT.TFile.Open(sys.argv[1]); op = sys.argv[2]; lumi = float(sys.argv[3]); G = int(sys.argv[4]) if len(sys.argv) > 4 else 1; SM = sys.argv[5] if len(sys.argv) > 5 else "none"
PROC = ["sm", "GGToLL", "Single_Top", "TT", "WW", "WZ", "ZZ", "DYtt"]
NU = ["QCDScale", "PDFweight", "alphaS", "PSWeight", "mu_reco", "mu_idiso", "mu_trig", "PU", "prefireWeight", "tt_ptrw", "rochester_stat", "rochester_syst"]
def h80(n):
    o = f.Get(n)
    return None if not o else np.array([o.GetBinContent(i) for i in range(1, o.GetNbinsX() + 1)])
def h(n):
    v = h80(n)
    return None if v is None else v.reshape(-1, G).sum(1)
nom80 = sum(h80(f"histo_{p}") for p in PROC); x80 = np.linspace(-1, 1, len(nom80))
def smooth(d80):
    r = d80 / nom80; c = np.polyfit(x80, r, 6, w=np.sqrt(nom80)); return np.polyval(c, x80) * nom80
def grp(v): return v.reshape(-1, G).sum(1)
nom = sum(h(f"histo_{p}") for p in PROC)
L80 = 0.5 * (h80(f"histo_w1_{op}") - h80(f"histo_wm1_{op}"))
J = {f"k_{op}": grp(smooth(L80)) if SM in ("eft", "all") else grp(L80), "lumi": nom * np.log(lumi)}
for nu in NU:
    up = sum((h(f"histo_{p}_{nu}Up") if h(f"histo_{p}_{nu}Up") is not None else h(f"histo_{p}")) for p in PROC)
    dn = sum((h(f"histo_{p}_{nu}Down") if h(f"histo_{p}_{nu}Down") is not None else h(f"histo_{p}")) for p in PROC)
    if SM == "all":
        up80 = sum((h80(f"histo_{p}_{nu}Up") if h80(f"histo_{p}_{nu}Up") is not None else h80(f"histo_{p}")) for p in PROC)
        dn80 = sum((h80(f"histo_{p}_{nu}Down") if h80(f"histo_{p}_{nu}Down") is not None else h80(f"histo_{p}")) for p in PROC)
        J[nu] = grp(smooth(0.5 * (up80 - dn80)))
    else:
        J[nu] = 0.5 * (up - dn)
names = list(J); A = np.array([J[n] for n in names])
F = (A / nom) @ A.T
prior = np.diag([0] + [1] * (len(names) - 1))
C = np.linalg.inv(F + prior)
s_stat = 1 / np.sqrt(F[0, 0]); s_syst = np.sqrt(C[0, 0])
print(f"[{SM:4s}] {len(nom)} bins: sigma(k_{op}) stat-only {s_stat:.4g}   stat+syst {s_syst:.4g}   ratio {s_syst/s_stat:.1f}")
print("   post-fit nuisance sigma: " + " ".join(f"{n}={np.sqrt(C[i,i]):.2f}" for i, n in enumerate(names) if i))
# which nuisances matter: sigma(k) when removing (fixing) one nuisance at a time
out = []
for j in range(1, len(names)):
    keep = [i for i in range(len(names)) if i != j]
    Cj = np.linalg.inv((F + prior)[np.ix_(keep, keep)])
    out.append((np.sqrt(Cj[0, 0]), names[j]))
print("   sigma(k) with one nuisance frozen: " + " ".join(f"{n}:{s:.3g}" for s, n in sorted(out)[:6]))
