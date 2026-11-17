#!/usr/bin/env python3
"""
Rank 2D grid scans by how "sketchy" they look, to pick the ones to inspect by hand.

Per pair (merged higgsCombine.<pair>.individual.MultiDimFit.mH125.root in --scan-dir):
  missing : grid points absent from the tree (failed fits / negative yields)
  spikes  : points inside the interesting region (2dNLL < 15) that differ from the median of their
            3x3 neighbourhood by > --spike (default 1.0 in 2dNLL) -> fits that jumped to another minimum
  tiny    : 68%/95% contour pieces enclosing < --tiny-cells grid cells (noise islands, not real lobes)
  edge    : 95% region touches the box edge (open contour)
  in95    : missing points inside the 95% region (nearest-neighbour filled)
Score = spikes + 2*tiny + in95 + 20*edge. Pairs listed worst first.

Usage: rank_2d_scan_quality.py --scan-dir full/final --skip cHQ1_cHQ3,cQl1_cQl3 [--dedup cHQ3:cHQ1,cQl3:cQl1] [--top 40]
"""
import argparse, glob, os, re
import numpy as np
import ROOT
from scipy.ndimage import median_filter, label, distance_transform_edt

ROOT.gROOT.SetBatch(True)
p = argparse.ArgumentParser()
p.add_argument("--scan-dir", required=True)
p.add_argument("--skip", default="")
p.add_argument("--dedup", default="", help="op aliases: pairs containing the first op are dropped (duplicates)")
p.add_argument("--spike", type=float, default=1.0)
p.add_argument("--tiny-cells", type=int, default=6)
p.add_argument("--top", type=int, default=40)
p.add_argument("--csv", default="")
a = p.parse_args()
skip = set(filter(None, a.skip.split(",")))
dup = [x.split(":")[0] for x in filter(None, a.dedup.split(","))]

rows = []
for fn in sorted(glob.glob(os.path.join(a.scan_dir, "higgsCombine.*.individual.MultiDimFit.mH125.root"))):
    pair = re.search(r"higgsCombine\.(\w+?)\.individual", os.path.basename(fn)).group(1)
    o1, o2 = pair.split("_")
    if pair in skip or o1 in dup or o2 in dup:
        continue
    d = ROOT.RDataFrame("limit", fn).Filter("quantileExpected > -0.5").AsNumpy([f"k_{o1}", f"k_{o2}", "deltaNLL"])
    x, y, z = d[f"k_{o1}"], d[f"k_{o2}"], 2 * d["deltaNLL"]
    xs, ys = np.unique(np.round(x, 10)), np.unique(np.round(y, 10))
    Z = np.full((len(ys), len(xs)), np.nan)
    Z[np.searchsorted(ys, np.round(y, 10)), np.searchsorted(xs, np.round(x, 10))] = z
    Z -= np.nanmin(Z)
    missing = int(np.isnan(Z).sum())
    holes = np.isnan(Z)
    Zf = Z[tuple(distance_transform_edt(holes, return_distances=False, return_indices=True))] if holes.any() else Z
    inside = int(np.sum(holes & (Zf < 5.99)))   # failed points inside the 95% region
    med = median_filter(Zf, size=3, mode="nearest")
    spikes = int(np.sum((np.abs(Zf - med) > a.spike) & (np.minimum(Zf, med) < 15)))
    tiny = 0
    for lev in (2.30, 5.99):
        lab, n = label(Zf < lev)
        sizes = np.bincount(lab.ravel())[1:]
        tiny += int(np.sum(sizes < a.tiny_cells))
    reg = Zf < 5.99
    edge = int(reg[0, :].any() or reg[-1, :].any() or reg[:, 0].any() or reg[:, -1].any())
    score = spikes + 2 * tiny + inside + 20 * edge
    rows.append((score, pair, missing, inside, spikes, tiny, edge, f"[{xs[0]:.3g},{xs[-1]:.3g}]x[{ys[0]:.3g},{ys[-1]:.3g}]"))

rows.sort(reverse=True)
print(f"{'score':>7s} {'pair':16s} {'missing':>7s} {'in95':>4s} {'spikes':>6s} {'tiny':>4s} {'edge':>4s}  box")
for r in rows[: a.top]:
    print(f"{r[0]:7.1f} {r[1]:16s} {r[2]:7d} {r[3]:4d} {r[4]:6d} {r[5]:4d} {r[6]:4d}  {r[7]}")
print(f"... {len(rows)} pairs scored; score>0: {sum(r[0] > 0 for r in rows)}")
if a.csv:
    with open(a.csv, "w") as f:
        f.write("score,pair,missing,in95,spikes,tiny,edge,box\n")
        for r in rows:
            f.write(",".join(str(v) for v in r) + "\n")
