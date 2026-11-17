#!/usr/bin/env python3
"""
2D positivity of the pure-EFT (quadratic) part, per bin, for every pair with a w11 template.

For a pair (a,b) the quadratic part of the yield is the quadratic form
  Q(ca,cb) = qa*ca^2 + qb*cb^2 + m*ca*cb,   m = mixed term
It is >= 0 for all directions only if the 2x2 matrix [[qa, m/2],[m/2, qb]] is positive semidefinite.
If not, the yield in that bin decreases without bound along some direction (non-convex NLL,
negative yields far from the origin, several minima).

Prints, per pair: bins with a negative eigenvalue, and the worst min-eigenvalue relative to sm in that bin.
Usage: check_quad_psd.py --shapes shapes.root [--pairs a_b,...] [--top 40]
"""
import argparse, re
import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
p = argparse.ArgumentParser()
p.add_argument("--shapes", required=True)
p.add_argument("--pairs", default="")
p.add_argument("--top", type=int, default=400)
a = p.parse_args()

f = ROOT.TFile.Open(a.shapes)
names = {k.GetName() for k in f.GetListOfKeys()}
A = lambda n: np.array([f.Get(n).GetBinContent(i) for i in range(1, f.Get(n).GetNbinsX() + 1)])
sm = A("histo_sm")
ops = sorted({m.group(1) for n in names for m in [re.fullmatch(r"histo_w1_(\w+)", n)] if m})
lin, quad = {}, {}
for o in ops:
    w1, wm = A(f"histo_w1_{o}"), A(f"histo_wm1_{o}")
    lin[o], quad[o] = (w1 - wm) / 2, (w1 + wm) / 2 - sm

pairs = a.pairs.split(",") if a.pairs else [n[len("histo_w11_"):] for n in names if re.fullmatch(r"histo_w11_[A-Za-z0-9]+_[A-Za-z0-9]+", n)]
rows = []
for pair in pairs:
    o1, o2 = pair.split("_")
    w11 = f"histo_w11_{o1}_{o2}" if f"histo_w11_{o1}_{o2}" in names else f"histo_w11_{o2}_{o1}"
    if w11 not in names:
        continue
    qa, qb = quad[o1], quad[o2]
    m = A(w11) - sm - lin[o1] - lin[o2] - qa - qb
    # min eigenvalue of [[qa, m/2],[m/2, qb]]
    emin = (qa + qb) / 2 - np.sqrt(((qa - qb) / 2) ** 2 + (m / 2) ** 2)
    scale = np.maximum(qa, 0) + np.maximum(qb, 0) + np.abs(m) / 2
    rel = np.where(scale > 0, emin / np.where(scale > 0, scale, 1), 0)
    bad = rel < -1e-6
    rows.append((rel.min(), pair, int(bad.sum()), len(sm), int((qa < 0).sum()), int((qb < 0).sum()), int(np.argmin(rel))))

rows.sort()
print(f"{'pair':16s} {'nonPSD bins':>11s} {'qa<0':>5s} {'qb<0':>5s} {'worst emin/scale':>16s} {'bin':>4s}")
for r in rows[: a.top]:
    print(f"{r[1]:16s} {r[2]:5d}/{r[3]:<5d} {r[4]:5d} {r[5]:5d} {r[0]:16.3g} {r[6]:4d}")
print(f"... {len(rows)} pairs; with non-PSD bins: {sum(r[2] > 0 for r in rows)}")
