#!/usr/bin/env python3
"""Scan histograms.root for combine and report bins with zero content.

Usage:
    python3 check_zero_bins.py histograms.root
    python3 check_zero_bins.py histograms.root --tol 1e-12
"""
import argparse
import uproot

parser = argparse.ArgumentParser()
parser.add_argument("rootfile")
parser.add_argument("--tol", type=float, default=0.0,
                     help="treat |value| <= tol as zero (default: exactly 0)")
args = parser.parse_args()

f = uproot.open(args.rootfile)

# top-level keys without the ";cycle" suffix, deduplicated
names = sorted(set(k.split(";")[0] for k in f.keys()))

found_any = False
for name in names:
    if "mixed" in name.lower():
        continue
    obj = f[name]
    # only look at histogram-like objects
    if not hasattr(obj, "values"):
        continue
    vals = obj.values()
    zero_bins = [i for i, v in enumerate(vals) if abs(v) <= args.tol]
    if zero_bins:
        found_any = True
        edges = obj.axes[0].edges() if hasattr(obj, "axes") else None
        print(f"{name}: {len(zero_bins)} zero bin(s) / {len(vals)} total")
        for i in zero_bins:
            if edges is not None:
                print(f"    bin {i:3d}  [{edges[i]:.1f}, {edges[i+1]:.1f}]  value={vals[i]}")
            else:
                print(f"    bin {i:3d}  value={vals[i]}")

if not found_any:
    print("No zero bins found (excluding *mixed* distributions).")
