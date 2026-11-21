#!/usr/bin/env python3
"""
Which operators can be profiled together? Linear Fisher-matrix estimate from the datacard templates.

  lin_op[b] = (w1 - wm1)/2,  sigma_b^2 = sm + bkg (stat only)
  F_ij = sum_b lin_i[b] lin_j[b] / sigma_b^2
Optionally a free overall normalisation (--free-norm: worst case for lumi / xsec nuisances) is profiled too.
--systs: the datacard nuisances are added as parameters with a unit Gaussian constraint and profiled
  (Schur complement): shape nuisances d mu/d theta = sum over sm+bkg of (Up-Down)/2, lnN (uniform value
  assumed) ln(k)*mu, and the MC stat error of sm+bkg added to the per-bin variance (autoMCStats-like).
  Then sigma_individual = only that operator floating (nuisances profiled), sigma_profiled = all selected
  operators floating (nuisances profiled).

Per operator: sigma_individual = 1/sqrt(F_ii), sigma_profiled = sqrt((F^-1)_ii), ratio = degradation from profiling.
Eigen-decomposition of the correlation-normalised F gives the weakest directions (near-flat combinations).
--greedy: drop, one at a time, the operator with the worst profiled/individual ratio until all ratios < --max-ratio.

Linear only: a guide to choose the subset, not a result (quadratic terms can close flat directions far out).
Usage: fisher_profiling.py --shapes triple_diff/shapes.root [--ops a,b,...] [--skip cHQ3,cQl3] [--free-norm]
       [--systs [--datacard datacard.txt]] [--greedy --max-ratio 3]
"""
import argparse, os, re
import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
p = argparse.ArgumentParser()
p.add_argument("--shapes", required=True)
p.add_argument("--ops", default="")
p.add_argument("--skip", default="cHQ3,cQl3")
p.add_argument("--free-norm", action="store_true")
p.add_argument("--systs", action="store_true", help="profile the datacard nuisances (shape, lnN, MC stat)")
p.add_argument("--datacard", default="", help="default: datacard.txt next to --shapes")
p.add_argument("--greedy", action="store_true")
p.add_argument("--max-ratio", type=float, default=3.0)
p.add_argument("--nweak", type=int, default=4)
p.add_argument("--full", action="store_true", help="print all coefficients of the weakest directions, also in operator units")
a = p.parse_args()

f = ROOT.TFile.Open(a.shapes)
names = {k.GetName() for k in f.GetListOfKeys()}
A = lambda n: np.array([f.Get(n).GetBinContent(i) for i in range(1, f.Get(n).GetNbinsX() + 1)])
sm = A("histo_sm")
bkg = sum(A(n) for n in names if n.startswith("histo_") and not any(t in n for t in ("Up", "Down", "Data", "sm", "w1_", "wm1_", "w11_")))
var = sm + bkg
procs = ["sm"] + [n[6:] for n in names if n.startswith("histo_") and not any(t in n for t in ("Up", "Down", "Data", "sm", "w1_", "wm1_", "w11_"))]
NU = []  # (name, d mu / d theta, prior)
if a.systs:
    E = lambda n: np.array([f.Get(n).GetBinError(i) for i in range(1, f.Get(n).GetNbinsX() + 1)])
    var = var + sum(E(f"histo_{pr}") ** 2 for pr in procs)
    for l in open(a.datacard or os.path.join(os.path.dirname(a.shapes), "datacard.txt")):
        w = l.split()
        if len(w) > 2 and w[1] == "shape":
            d = sum((A(f"histo_{pr}_{w[0]}Up") - A(f"histo_{pr}_{w[0]}Down")) / 2 for pr in procs if f"histo_{pr}_{w[0]}Up" in names)
            NU.append((w[0], d, 1.0))
        elif len(w) > 2 and w[1] == "lnN":
            NU.append((w[0], np.log(float(w[2])) * (sm + bkg), 1.0))
sig = np.sqrt(np.maximum(var, 1e-9))
allops = sorted({m.group(1) for n in names for m in [re.fullmatch(r"histo_w1_([A-Za-z0-9]+)", n)] if m})
ops = [o for o in (a.ops.split(",") if a.ops else allops) if o not in a.skip.split(",")]
L = {o: (A(f"histo_w1_{o}") - A(f"histo_wm1_{o}")) / 2 / sig for o in ops}
norm = (sm + bkg) / sig  # derivative of the yield w.r.t. a relative normalisation shift
N = np.array([d / sig for _, d, _ in NU] + ([norm] if a.free_norm else []))
prior = np.array([pr for _, _, pr in NU] + ([0.0] if a.free_norm else []))


def analyse(sel, verbose):
    M = np.array([L[o] for o in sel])
    F = M @ M.T
    if len(N):  # profile nuisances / free norm: Schur complement
        Fon = M @ N.T
        F = F - Fon @ np.linalg.solve(N @ N.T + np.diag(prior), Fon.T)
    C = np.linalg.pinv(F)
    n = len(sel)
    ind = 1 / np.sqrt(np.diag(F)[:n])
    prof = np.sqrt(np.clip(np.diag(C)[:n], 0, None))
    if verbose:
        D = np.diag(1 / np.sqrt(np.diag(F)))
        ev, V = np.linalg.eigh(D @ F @ D)
        print(f"  condition number (corr.-normalised F): {ev[-1] / max(ev[0], 1e-300):.3g}")
        print(f"  {'op':6s} {'sigma_ind':>10s} {'sigma_prof':>10s} {'ratio':>7s}")
        for i in np.argsort(-prof / ind):
            print(f"  {sel[i]:6s} {ind[i]:10.4g} {prof[i]:10.4g} {prof[i] / ind[i]:7.1f}")
        labels = sel
        print(f"  weakest directions (eigenvalue; sigma along it = 1/sqrt(eigenvalue) in units of individual sigmas):")
        for k in range(min(a.nweak, len(ev))):
            v = V[:, k]
            top = np.argsort(-np.abs(v))[:5]
            print(f"    lambda={ev[k]:.2e} (x{1 / np.sqrt(max(ev[k], 1e-300)):.0f}): " + " ".join(f"{v[j]:+.2f}*{labels[j]}" for j in top))
            if a.full:
                # v is in units of individual sigmas; in operator units the flat direction is c = D v
                c = np.diag(D) * v
                c = c / np.abs(c).max()
                order = np.argsort(-np.abs(v))
                print("      " + f"{'op':6s} {'sigma units':>12s} {'operator units':>15s}")
                for j in order:
                    if abs(v[j]) > 0.01:
                        print(f"      {labels[j]:6s} {v[j]:+12.3f} {c[j]:+15.4g}")
                print(f"      residual |sum_i c_i lin_i| / |lin| = {np.linalg.norm(c @ M) / np.linalg.norm(M, axis=1).max():.1e}")
    return dict(zip(sel, prof / ind))


print(f"{a.shapes}: {len(sm)} bins, {'stat + ' + str(len(NU)) + ' nuisances + MC stat' if a.systs else 'stat only'}{' + free normalisation' if a.free_norm else ''}, linear terms")
print(f"== all {len(ops)} operators")
r = analyse(ops, True)
if a.greedy:
    sel = list(ops)
    dropped = []
    while max(r.values()) > a.max_ratio and len(sel) > 1:
        worst = max(r, key=r.get)
        dropped.append(f"{worst}({r[worst]:.1f})")
        sel.remove(worst)
        r = analyse(sel, False)
    print(f"\n== greedy: largest set with all profiled/individual ratios < {a.max_ratio}: {len(sel)} operators")
    print("  dropped (ratio when dropped): " + ", ".join(dropped))
    analyse(sel, True)
    print("  --doOnly " + ",".join(sel))
