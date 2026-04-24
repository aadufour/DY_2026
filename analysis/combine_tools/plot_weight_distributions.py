#!/usr/bin/env python3
"""
plot_weight_distributions.py

Plot linear and quadratic EFT weight distributions (normalized to xwgt)
for all operators in the cache. Useful to spot poorly modelled operators
(peak at 0 = low statistical power).

Usage:
    plot_weight_distributions.py --cache lhe_cache.pkl --output plots/
"""

import argparse
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt

CACHE_FILE = "/grid_mnt/data__data.polcms/cms/adufour/MG5/mg5amcnlo/CACHE/lhe_cache.pkl"
OUTPUT_DIR = "weight_plots"

parser = argparse.ArgumentParser()
parser.add_argument("--cache",  default=CACHE_FILE)
parser.add_argument("--output", default=OUTPUT_DIR)
args = parser.parse_args()

os.makedirs(args.output, exist_ok=True)

print(f"Loading cache: {args.cache}")
with open(args.cache, "rb") as f:
    cache = pickle.load(f)

w_SM     = cache["w_SM"]
xwgt     = cache["xwgt"]
w_p1_all = cache["w_p1"]
w_m1_all = cache["w_m1"]
operators = sorted(w_p1_all.keys())
print(f"Found {len(operators)} operators: {operators}\n")

for op in operators:
    wp1  = np.array(w_p1_all[op])
    wm1  = np.array(w_m1_all[op])
    wsm  = np.array(w_SM)
    xwgt_arr = np.array(xwgt)

    w_lin  = (wp1 - wm1) / (2 * xwgt_arr)
    w_quad = (wp1 + wm1 - 2 * wsm) / (2 * xwgt_arr)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    fig.suptitle(op)

    for ax, weights, label in zip(axes, [w_lin, w_quad], ["linear / xwgt", "quadratic / xwgt"]):
        finite = weights[np.isfinite(weights)]
        lo, hi = np.percentile(finite, [1, 99])
        ax.hist(finite, bins=100, range=(lo, hi), histtype="step", color="steelblue")
        ax.set_xlabel(label)
        ax.set_ylabel("events")
        ax.set_yscale("log")

    plt.tight_layout()
    out = os.path.join(args.output, f"{op}.png")
    plt.savefig(out, dpi=120)
    plt.close()
    print(f"  {op} -> {out}")

print("\nDone.")
