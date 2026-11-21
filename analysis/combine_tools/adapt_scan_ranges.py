#!/usr/bin/env python3
"""
Per-operator 1D scan ranges sized to the 95% interval, from an existing (possibly too coarse) scan.

Problem it fixes: with a fixed range much wider than the interval, the grid points next to the minimum already sit at
2dNLL >> 4 and the 2dNLL = 4 "crossing" is a linear interpolation between 0 and e.g. 800, i.e. an interval that scales
like sigma^2 instead of sigma (2026-11-21, N80 rapll_abs cHe stat-only: points at +-0.05 with 2dNLL 800 -> +-0.0002,
true +-0.0035).

Per side of the minimum x0:
  resolved   (>= 3 scan points with 2dNLL < 4): w = distance of the 2dNLL = 4 crossing
  unresolved (< 3 points)                    : parabola through the nearest point (d, y): w = 2 d / sqrt(y)
  no crossing (2dNLL < 4 up to the range edge): keep the original edge
New range = [x0 - F wL, x0 + F wR] (F = --factor, default 3.5 -> ~28 of 100 points inside the 95% interval).

Writes <vardir>/<outdir>/metadata.json (from the stat+syst scans) and metadata_stat.json (from the _stat scans) and
symlinks the initialFit snapshots there, so that `runScans.py 1 scan [--stat] --metadata ...` can run in <outdir>.

Usage: adapt_scan_ranges.py <vardir> [--outdir rescan] [--factor 3.5] [--src <dir with the scans to size from>]
"""
import argparse, glob, json, os
import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
p = argparse.ArgumentParser()
p.add_argument("vardir")
p.add_argument("--outdir", default="rescan")
p.add_argument("--factor", type=float, default=3.5)
p.add_argument("--src", default="", help="folder with the scans to size the ranges from (default: vardir)")
a = p.parse_args()
vd = os.path.abspath(a.vardir); src = os.path.abspath(a.src or a.vardir); od = os.path.join(vd, a.outdir)
os.makedirs(od, exist_ok=True)
meta = json.load(open(os.path.join(vd, "metadata.json")))


def side_width(x, y, x0, sign, edge):
    m = (x - x0) * sign > 0
    xs, ys = np.abs(x[m] - x0), y[m]
    o = np.argsort(xs); xs, ys = xs[o], ys[o]
    if len(xs) == 0 or ys.max() < 4:
        return None, "edge"
    k = np.argmax(ys >= 4)  # first point beyond the crossing
    if k >= 3:
        x1, y1 = (xs[k - 1], ys[k - 1])
        return x1 + (4 - y1) * (xs[k] - x1) / (ys[k] - y1), "resolved"
    j = np.argmax(ys > 0)
    return 2 * xs[j] / np.sqrt(ys[j]), "parabola"


for suf, fname in (("", "metadata.json"), ("_stat", "metadata_stat.json")):
    out = dict(meta); out["operators"] = {}
    for op, (lo0, hi0) in meta["operators"].items():
        f = f"{src}/higgsCombine.{op}{suf}.individual.MultiDimFit.mH125.root"
        if not os.path.exists(f):
            out["operators"][op] = [lo0, hi0]; continue
        d = ROOT.RDataFrame("limit", f).AsNumpy([f"k_{op}", "deltaNLL"])
        x, y = d[f"k_{op}"].astype(float), 2 * d["deltaNLL"].astype(float)
        y = y - y.min(); x0 = x[np.argmin(y)]
        wl, tl = side_width(x, y, x0, -1, lo0); wr, tr = side_width(x, y, x0, +1, hi0)
        lo = lo0 if wl is None else x0 - a.factor * wl
        hi = hi0 if wr is None else x0 + a.factor * wr
        out["operators"][op] = [float(f"{lo:.4g}"), float(f"{hi:.4g}")]
        print(f"{op}{suf:5s}: old [{lo0:.3g},{hi0:.3g}] -> [{lo:.3g},{hi:.3g}]   ({tl}/{tr})")
    json.dump(out, open(os.path.join(od, fname), "w"), indent=2)

for f in glob.glob(f"{vd}/higgsCombine.initialFit_*.MultiDimFit.mH125.root") + [f"{vd}/{n}" for n in ("datacard.txt", "jsonComb.json", "shapes.root")]:
    t = os.path.join(od, os.path.basename(f))
    if os.path.exists(f) and not os.path.lexists(t):
        os.symlink(f, t)
print("wrote", od)
