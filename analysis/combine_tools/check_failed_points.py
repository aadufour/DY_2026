#!/usr/bin/env python3
"""
Why do grid points of a 2D scan fail (are missing from the limit tree)?

For each pair, rebuild the full grid from the scan(s), and at every point compute the
predicted signal yield per bin from the morphing templates at nominal nuisances:
  N(c1,c2) = sm + c1*lin1 + c2*lin2 + c1^2*quad1 + c2^2*quad2 + c1*c2*mixed
  lin = (w1-wm1)/2, quad = (w1+wm1)/2 - sm, mixed = w11_12 - sm - lin1 - lin2 - quad1 - quad2
Then compare failed vs successful points: fraction with min-bin yield < 0 (+ backgrounds).

Usage: check_failed_points.py --shapes shapes.root --scan higgsCombine.<pair>...root [--grid-from other.root] --pairs a_b,...
       (--scan-dir dir: uses higgsCombine.<pair>.individual.MultiDimFit.mH125.root inside)
"""
import argparse, os
import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
p = argparse.ArgumentParser()
p.add_argument("--shapes", required=True)
p.add_argument("--scan-dir", required=True)
p.add_argument("--grid-dir", default="", help="another scan of the same pair/box, to complete the grid")
p.add_argument("--pairs", required=True)
a = p.parse_args()

f = ROOT.TFile.Open(a.shapes)
names = {k.GetName() for k in f.GetListOfKeys()}
A = lambda n: np.array([f.Get(n).GetBinContent(i) for i in range(1, f.Get(n).GetNbinsX() + 1)])
sm = A("histo_sm")
bkg = sum(A(n) for n in names if n.startswith("histo_") and not any(t in n for t in ("Up", "Down", "Data", "sm", "w1_", "wm1_", "w11_")))


def pts(fn, o1, o2):
    d = ROOT.RDataFrame("limit", fn).Filter("quantileExpected>-0.5").AsNumpy([f"k_{o1}", f"k_{o2}"])
    return {(round(float(x), 8), round(float(y), 8)) for x, y in zip(d[f"k_{o1}"], d[f"k_{o2}"])}


for pair in a.pairs.split(","):
    o1, o2 = pair.split("_")
    fn = os.path.join(a.scan_dir, f"higgsCombine.{pair}.individual.MultiDimFit.mH125.root")
    ok = pts(fn, o1, o2)
    allp = set(ok)
    if a.grid_dir:
        allp |= pts(os.path.join(a.grid_dir, f"higgsCombine.{pair}.individual.MultiDimFit.mH125.root"), o1, o2)
    xs = sorted({x for x, _ in allp}); ys = sorted({y for _, y in allp})
    grid = [(x, y) for x in xs for y in ys]
    w1a, wma, w1b, wmb = A(f"histo_w1_{o1}"), A(f"histo_wm1_{o1}"), A(f"histo_w1_{o2}"), A(f"histo_wm1_{o2}")
    la, qa, lb, qb = (w1a - wma) / 2, (w1a + wma) / 2 - sm, (w1b - wmb) / 2, (w1b + wmb) / 2 - sm
    w11 = f"histo_w11_{o1}_{o2}" if f"histo_w11_{o1}_{o2}" in names else f"histo_w11_{o2}_{o1}"
    mix = A(w11) - sm - la - lb - qa - qb if w11 in names else 0 * sm

    def minyield(c1, c2):
        return np.min(sm + c1 * la + c2 * lb + c1 ** 2 * qa + c2 ** 2 * qb + c1 * c2 * mix)

    def mintot(c1, c2):
        return np.min(sm + bkg + c1 * la + c2 * lb + c1 ** 2 * qa + c2 ** 2 * qb + c1 * c2 * mix)

    failed = [g for g in grid if g not in ok]
    good = [g for g in grid if g in ok]
    fneg = np.mean([minyield(*g) < 0 for g in failed]) if failed else 0
    gneg = np.mean([minyield(*g) < 0 for g in good])
    fnegt = np.mean([mintot(*g) < 0 for g in failed]) if failed else 0
    print(f"{pair:14s} grid {len(grid)} failed {len(failed):5d} | signal yield<0 in some bin: failed {fneg:6.1%}  ok {gneg:6.1%}"
          f" | signal+bkg<0: failed {fnegt:6.1%} | mixed term {'yes' if w11 in names else 'MISSING'}")
