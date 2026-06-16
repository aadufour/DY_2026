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
parser.add_argument("--negative", action="store_true",
                     help="instead of zero bins, report bins with value < -tol")
args = parser.parse_args()

f = uproot.open(args.rootfile)

# recurse through all subdirectories, only keep actual histogram objects (TH1*)
names = sorted(set(
    k.split(";")[0] for k, cls in f.classnames(recursive=True).items()
    if cls.startswith("TH1")
))

found_any = False
for name in names:
    if "mixed" in name.lower():
        continue
    obj = f[name]
    vals = obj.values()
    if args.negative:
        bad_bins = [i for i, v in enumerate(vals) if v < -args.tol]
        label = "negative"
    else:
        bad_bins = [i for i, v in enumerate(vals) if abs(v) <= args.tol]
        label = "zero"
    if bad_bins:
        found_any = True
        edges = obj.axes[0].edges() if hasattr(obj, "axes") else None
        print(f"{name}: {len(bad_bins)} {label} bin(s) / {len(vals)} total")
        for i in bad_bins:
            if edges is not None:
                print(f"    bin {i:3d}  [{edges[i]:.1f}, {edges[i+1]:.1f}]  value={vals[i]}")
            else:
                print(f"    bin {i:3d}  value={vals[i]}")

if not found_any:
    label = "negative" if args.negative else "zero"
    print(f"No {label} bins found (excluding *mixed* distributions).")
