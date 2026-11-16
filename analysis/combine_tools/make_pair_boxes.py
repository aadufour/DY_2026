#!/usr/bin/env python3
"""
Per-pair 2D scan boxes from an existing (e.g. coarse) grid scan.

For each pair, take the min/max of all grid points with 2*deltaNLL < --level
(default 5.99, the 95% CL for 2 parameters) and pad them by --pad times the
extent on each side:  [lo - pad*(hi-lo), hi + pad*(hi-lo)].
  pad 0.35 -> the 95% contour fills ~60% of the box (standalone plots)
  pad 1.0  -> box 3x the 95% extent (comparison frame, e.g. mll vs triple diff)

If the 95% region reaches the edge of the scanned box on a side, that side is
"open": the new edge is twice the old one (measured from 0) and the pair needs
another pass. Pairs in --keep (e.g. exact degeneracies) keep their current box.

Output: a copy of --metadata with a "pairs" section {pair: {op1: [lo,hi], op2: [lo,hi]}}
that runScans.py uses instead of the per-operator ranges.

Usage (in datacards/<var>, pyROOT env):
  make_pair_boxes.py --metadata metadata.json --pairs cHj3_cHu,clj3_cQe --pad 0.35 -o metadata_pairs.json
  make_pair_boxes.py ... --scan-dir ../triple_diff -o ../mll/metadata_pairs.json   # boxes for another variable
"""
import argparse, json, os
import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)

p = argparse.ArgumentParser()
p.add_argument("--metadata", default="metadata.json")
p.add_argument("--scan-dir", default=".")
p.add_argument("--pairs", default="", help="comma list; default: every pair with a merged scan file")
p.add_argument("--keep", default="", help="comma list of pairs whose box is not changed")
p.add_argument("--level", type=float, default=5.99)
p.add_argument("--pad", type=float, default=0.35)
p.add_argument("-o", "--out", default="metadata_pairs.json")
a = p.parse_args()

meta = json.load(open(a.metadata))
ops = list(meta["operators"])
old_pairs = meta.get("pairs", {})
if a.pairs:
    pairs = a.pairs.split(",")
else:
    pairs = [f"{x}_{y}" for i, x in enumerate(ops) for y in ops[i + 1:]
             if os.path.exists(os.path.join(a.scan_dir, f"higgsCombine.{x}_{y}.individual.MultiDimFit.mH125.root"))]
keep = set(filter(None, a.keep.split(",")))


def sig(v, n=3):
    return float(f"{v:.{n}g}")


new_pairs, n_open = dict(old_pairs), 0
print(f"{'pair':16s} {'status':8s} {'op1 box':>22s} {'op2 box':>22s}  95% extent (op1 | op2)")
for pair in pairs:
    o1, o2 = pair.split("_")
    fn = os.path.join(a.scan_dir, f"higgsCombine.{pair}.individual.MultiDimFit.mH125.root")
    d = ROOT.RDataFrame("limit", fn).Filter("quantileExpected > -0.5").AsNumpy([f"k_{o1}", f"k_{o2}", "deltaNLL"])
    x, y, z = d[f"k_{o1}"], d[f"k_{o2}"], 2 * (d["deltaNLL"] - d["deltaNLL"].min())
    box, ext, status = {}, {}, "closed"
    for op, v in ((o1, x), (o2, y)):
        vals = np.unique(v)
        step = np.min(np.diff(vals)) if len(vals) > 1 else 0
        vmin, vmax = vals.min(), vals.max()
        sel = v[z < a.level]
        lo, hi = sel.min(), sel.max()
        ext[op] = (lo, hi)
        w = hi - lo
        nlo, nhi = lo - a.pad * w, hi + a.pad * w
        if lo <= vmin + 1.5 * step:
            nlo, status = 2 * vmin, "open"
        if hi >= vmax - 1.5 * step:
            nhi, status = 2 * vmax, "open"
        box[op] = [sig(nlo), sig(nhi)]
    if pair in keep:
        status, box = "kept", {op: old_pairs.get(pair, meta["operators"])[op] for op in (o1, o2)}
    n_open += status == "open"
    new_pairs[pair] = box
    print(f"{pair:16s} {status:8s} {str(box[o1]):>22s} {str(box[o2]):>22s}  "
          f"[{ext[o1][0]:.3g},{ext[o1][1]:.3g}] | [{ext[o2][0]:.3g},{ext[o2][1]:.3g}]")

meta["pairs"] = new_pairs
json.dump(meta, open(a.out, "w"), indent=4)
print(f"\nwrote {a.out}: {len(pairs)} pairs, {n_open} still open (re-run them with the new box)")
