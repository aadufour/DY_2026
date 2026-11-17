#!/usr/bin/env python3
"""
Which operators can be profiled together? Linear Fisher-matrix estimate from the datacard templates.

  lin_op[b] = (w1 - wm1)/2,  sigma_b^2 = sm + bkg (stat only)
  F_ij = sum_b lin_i[b] lin_j[b] / sigma_b^2
Optionally a free overall normalisation (--free-norm: worst case for lumi / xsec nuisances) is profiled too.

Per operator: sigma_individual = 1/sqrt(F_ii), sigma_profiled = sqrt((F^-1)_ii), ratio = degradation from profiling.
Eigen-decomposition of the correlation-normalised F gives the weakest directions (near-flat combinations).
--greedy: drop, one at a time, the operator with the worst profiled/individual ratio until all ratios < --max-ratio.

Linear only and stat only (+ optional free norm): a guide to choose the subset, not a result.
Usage: fisher_profiling.py --shapes triple_diff/shapes.root [--ops a,b,...] [--skip cHQ3,cQl3] [--free-norm] [--greedy --max-ratio 3]
"""
import argparse, re
import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
p = argparse.ArgumentParser()
p.add_argument("--shapes", required=True)
p.add_argument("--ops", default="")
p.add_argument("--skip", default="cHQ3,cQl3")
p.add_argument("--free-norm", action="store_true")
p.add_argument("--greedy", action="store_true")
p.add_argument("--max-ratio", type=float, default=3.0)
p.add_argument("--nweak", type=int, default=4)
a = p.parse_args()

f = ROOT.TFile.Open(a.shapes)
names = {k.GetName() for k in f.GetListOfKeys()}
A = lambda n: np.array([f.Get(n).GetBinContent(i) for i in range(1, f.Get(n).GetNbinsX() + 1)])
sm = A("histo_sm")
bkg = sum(A(n) for n in names if n.startswith("histo_") and not any(t in n for t in ("Up", "Down", "Data", "sm", "w1_", "wm1_", "w11_")))
sig = np.sqrt(np.maximum(sm + bkg, 1e-9))
allops = sorted({m.group(1) for n in names for m in [re.fullmatch(r"histo_w1_([A-Za-z0-9]+)", n)] if m})
ops = [o for o in (a.ops.split(",") if a.ops else allops) if o not in a.skip.split(",")]
L = {o: (A(f"histo_w1_{o}") - A(f"histo_wm1_{o}")) / 2 / sig for o in ops}
norm = (sm + bkg) / sig  # derivative of the yield w.r.t. a relative normalisation shift


def analyse(sel, verbose):
    M = np.array([L[o] for o in sel] + ([norm] if a.free_norm else []))
    F = M @ M.T
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
        labels = sel + (["norm"] if a.free_norm else [])
        print(f"  weakest directions (eigenvalue; sigma along it = 1/sqrt(eigenvalue) in units of individual sigmas):")
        for k in range(min(a.nweak, len(ev))):
            v = V[:, k]
            top = np.argsort(-np.abs(v))[:5]
            print(f"    lambda={ev[k]:.2e} (x{1 / np.sqrt(max(ev[k], 1e-300)):.0f}): " + " ".join(f"{v[j]:+.2f}*{labels[j]}" for j in top))
    return dict(zip(sel, prof / ind))


print(f"{a.shapes}: {len(sm)} bins, stat only{' + free normalisation' if a.free_norm else ''}, linear terms")
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
