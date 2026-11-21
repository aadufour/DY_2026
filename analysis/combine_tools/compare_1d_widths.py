#!/usr/bin/env python3
"""
Width of the 1D 95% intervals (2dNLL = 4 crossings), production B vs production A, for every variable and operator,
stat+syst and stat-only. Flags scans whose crossing is not clean:
  edge  : fewer than 2 crossings (interval = scan range edge, i.e. the scan range is too small or the curve is broken)
  multi : more than 2 crossings (non-parabolic / several minima)

Usage: compare_1d_widths.py <config_dir_A> <config_dir_B> [--labels v3 N80] [--vars mll rapll_abs ...] [-o out.md]
"""
import argparse, glob, os
import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
p = argparse.ArgumentParser()
p.add_argument("A"); p.add_argument("B")
p.add_argument("--labels", nargs=2, default=["A", "B"])
p.add_argument("--vars", nargs="+", default=["mll", "rapll_abs", "costhetastar", "triple_diff"])
p.add_argument("-o", "--out", default="")
a = p.parse_args()


def interval(f, op):
    d = ROOT.RDataFrame("limit", f).AsNumpy([f"k_{op}", "deltaNLL"])
    x, y = d[f"k_{op}"], 2 * d["deltaNLL"]
    o = np.argsort(x); x, y = x[o], y[o] - y.min()
    s = np.sign(y - 4)
    idx = np.where(s[1:] * s[:-1] < 0)[0]
    xs = [x[i] + (4 - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i]) for i in idx]
    flag = "edge" if len(xs) < 2 else ("multi" if len(xs) > 2 else "")
    lo, hi = (min(xs), max(xs)) if len(xs) >= 2 else (x.min() if not xs or xs[0] > 0 else xs[0], x.max() if not xs or xs[0] < 0 else xs[0])
    return lo, hi, flag


lines = []
for v in a.vars:
    da, db = [f"{c}/datacards/inc_mm/{v}" for c in (a.A, a.B)]
    ops = sorted({os.path.basename(f).split(".")[1] for f in glob.glob(f"{db}/higgsCombine.*.individual.MultiDimFit.mH125.root")
                  if "_stat" not in f})
    lines += [f"\n### {v}: 95% width ratio {a.labels[1]}/{a.labels[0]} (< 1 = {a.labels[1]} tighter)", "",
              f"| op | {a.labels[0]} 95% (stat+syst) | {a.labels[1]} 95% (stat+syst) | ratio | ratio stat-only | flags |",
              "|---|---|---|---|---|---|"]
    for op in ops:
        row, rat, fl = [], [], []
        for suf in ("", "_stat"):
            r = []
            for lab, d in zip(a.labels, (da, db)):
                f = f"{d}/higgsCombine.{op}{suf}.individual.MultiDimFit.mH125.root"
                if not os.path.exists(f):
                    r.append(None); continue
                lo, hi, flag = interval(f, op)
                r.append((lo, hi))
                if flag: fl.append(f"{lab}{suf or ''}:{flag}")
            if suf == "": row = r
            rat.append((r[1][1] - r[1][0]) / (r[0][1] - r[0][0]) if None not in r else np.nan)
        fmt = lambda t: f"[{t[0]:+.3g}, {t[1]:+.3g}]" if t else "-"
        lines.append(f"| {op} | {fmt(row[0])} | {fmt(row[1])} | {rat[0]:.2f} | {rat[1]:.2f} | {' '.join(fl)} |")
txt = "\n".join(lines) + "\n"
print(txt)
if a.out:
    open(a.out, "w").write(txt)
