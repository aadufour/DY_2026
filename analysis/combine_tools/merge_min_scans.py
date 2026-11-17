#!/usr/bin/env python3
"""
Combine several grid scans of the same pair/box (e.g. loose and strict minimizer settings)
by taking, at every grid point, the lowest 2*deltaNLL found (the best minimum = best estimate
of the profile likelihood; deltaNLL is relative to the same initial fit in every scan).
A point is missing only if it failed in all inputs.

Writes <outdir>/higgsCombine.<pair>.individual.MultiDimFit.mH125.root with a "limit" tree
(k_<op1>, k_<op2>, deltaNLL, quantileExpected), readable by plot_2d_scans.py etc.

Usage: merge_min_scans.py --pairs a_b,c_d --inputs dir1 dir2 [...] -o outdir
"""
import argparse, os
import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
p = argparse.ArgumentParser()
p.add_argument("--pairs", required=True)
p.add_argument("--inputs", nargs="+", required=True)
p.add_argument("-o", "--outdir", required=True)
a = p.parse_args()
os.makedirs(a.outdir, exist_ok=True)

for pair in a.pairs.split(","):
    o1, o2 = pair.split("_")
    best, src = {}, {}
    for i, d in enumerate(a.inputs):
        fn = os.path.join(d, f"higgsCombine.{pair}.individual.MultiDimFit.mH125.root")
        if not os.path.exists(fn):
            continue
        r = ROOT.RDataFrame("limit", fn).Filter("quantileExpected > -0.5").AsNumpy([f"k_{o1}", f"k_{o2}", "deltaNLL"])
        for x, y, z in zip(r[f"k_{o1}"], r[f"k_{o2}"], r["deltaNLL"]):
            k = (round(float(x), 9), round(float(y), 9))
            if np.isfinite(z) and (k not in best or z < best[k]):
                best[k], src[k] = float(z), i
    keys = sorted(best)
    cols = {f"k_{o1}": np.array([k[0] for k in keys], dtype=np.float32),
            f"k_{o2}": np.array([k[1] for k in keys], dtype=np.float32),
            "deltaNLL": np.array([best[k] for k in keys], dtype=np.float64),
            "quantileExpected": np.zeros(len(keys), dtype=np.float32)}
    out = os.path.join(a.outdir, f"higgsCombine.{pair}.individual.MultiDimFit.mH125.root")
    ROOT.RDF.FromNumpy(cols).Snapshot("limit", out)
    n_from = np.bincount([src[k] for k in keys], minlength=len(a.inputs))
    print(f"{pair:14s} {len(keys)} points; best value taken from: " + ", ".join(f"{os.path.basename(d.rstrip('/'))}={n}" for d, n in zip(a.inputs, n_from)))
