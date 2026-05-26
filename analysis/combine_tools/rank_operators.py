#!/usr/bin/env python3
"""
rank_operators.py
=================
Extract 68% and 95% CL intervals from combine likelihood scan ROOT files
and produce a ranked summary plot + table.

Reads: higgsCombine.{op}.individual.MultiDimFit.mH125.root  (one per operator)
       higgsCombine.{op}_stat.individual.MultiDimFit.mH125.root  (stat-only, optional)

In combine's output TTree "limit":
  deltaNLL = ΔNLL = NLL - NLL_min   (NOT multiplied by 2)
  → 68% CL : 2*deltaNLL < 1    → deltaNLL < 0.5
  → 95% CL : 2*deltaNLL < 3.84 → deltaNLL < 1.92

Usage (from the datacards_morphing/inc_mm/mll/ dir, in dy_combine_morphing env):
    python3 /path/to/rank_operators.py
    python3 /path/to/rank_operators.py --outdir ranking/ --stat
"""

import argparse
import os
import json

import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep

OPERATORS = [
    "cHDD", "cHWB", "cbWRe", "cbBRe",
    "cHj1", "cHQ1", "cHj3", "cHQ3",
    "cHu",  "cHd",  "cHbq",
    "cHl1", "cHl3", "cHe",
    "cll1",
    "clj1", "clj3",
    "cQl1", "cQl3",
    "ceu",  "ced",
    "cbe",  "cje",  "cQe",
    "clu",  "cld",  "cbl",
]

CL68  = 0.5    # deltaNLL threshold for 68% CL  (Δ(-2lnL) = 1)
CL95  = 1.92   # deltaNLL threshold for 95% CL  (Δ(-2lnL) = 3.84)


def extract_interval(fname, poi):
    """Return (lo68, hi68, lo95, hi95) or None if file missing / no valid points."""
    if not os.path.isfile(fname):
        return None
    try:
        with uproot.open(fname) as f:
            t = f["limit"]
            dnll = t["deltaNLL"].array(library="np")
            vals = t[poi].array(library="np")
    except Exception as e:
        print(f"  WARNING: could not read {fname}: {e}")
        return None

    # drop the best-fit point (deltaNLL == 0) and any failed points
    mask = dnll > 0
    if mask.sum() == 0:
        return None
    dnll = dnll[mask]
    vals = vals[mask]

    def cl_interval(threshold):
        inside = vals[dnll < threshold]
        if len(inside) == 0:
            return None, None
        return float(inside.min()), float(inside.max())

    lo68, hi68 = cl_interval(CL68)
    lo95, hi95 = cl_interval(CL95)
    return lo68, hi68, lo95, hi95


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--indir",    default=".",
                        help="Directory containing higgsCombine.*.root files (default: .)")
    parser.add_argument("--outdir",   default="ranking",
                        help="Output directory for plots and table")
    parser.add_argument("--metadata", default="metadata.json",
                        help="metadata.json with operator scan ranges")
    parser.add_argument("--stat",     action="store_true",
                        help="Also overlay stat-only intervals (_stat files)")
    parser.add_argument("--operators", nargs="+", default=None,
                        help="Subset of operators (default: all 27)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    ops = args.operators or OPERATORS

    # Load metadata for scan ranges (optional, used for context)
    meta_ops = {}
    if os.path.isfile(args.metadata):
        with open(args.metadata) as f:
            meta = json.load(f)
        meta_ops = meta.get("operators", {})

    results = []
    results_stat = []

    for op in ops:
        poi  = f"k_{op}"
        fname      = os.path.join(args.indir, f"higgsCombine.{op}.individual.MultiDimFit.mH125.root")
        fname_stat = os.path.join(args.indir, f"higgsCombine.{op}_stat.individual.MultiDimFit.mH125.root")

        iv = extract_interval(fname, poi)
        if iv is None:
            print(f"  SKIP {op}: no valid scan points in {fname}")
            continue
        lo68, hi68, lo95, hi95 = iv
        width95 = (hi95 - lo95) if (lo95 is not None and hi95 is not None) else None
        results.append({"op": op, "lo68": lo68, "hi68": hi68,
                         "lo95": lo95, "hi95": hi95, "width95": width95})

        if args.stat:
            iv_s = extract_interval(fname_stat, poi)
            if iv_s is not None:
                lo68s, hi68s, lo95s, hi95s = iv_s
                w95s = (hi95s - lo95s) if (lo95s is not None and hi95s is not None) else None
                results_stat.append({"op": op, "lo68": lo68s, "hi68": hi68s,
                                      "lo95": lo95s, "hi95": hi95s, "width95": w95s})

    if not results:
        print("No results found. Check that higgsCombine.*.individual.MultiDimFit.mH125.root files exist.")
        return

    # Sort by 95% CL interval width (ascending = most sensitive first)
    results_sorted = sorted(results, key=lambda r: r["width95"] if r["width95"] else 1e9)
    stat_by_op = {r["op"]: r for r in results_stat}

    # ------------------------------------------------------------------
    # Text table
    # ------------------------------------------------------------------
    print(f"\n{'Operator':<12}  {'lo95':>9}  {'hi95':>9}  {'width95':>9}  {'lo68':>9}  {'hi68':>9}")
    print("-" * 65)
    for r in results_sorted:
        print(f"{r['op']:<12}  {r['lo95']:>9.4f}  {r['hi95']:>9.4f}  {r['width95']:>9.4f}"
              f"  {r['lo68']:>9.4f}  {r['hi68']:>9.4f}")

    # ------------------------------------------------------------------
    # Bar chart: 95% and 68% CL intervals, ranked
    # ------------------------------------------------------------------
    hep.style.use("CMS")
    fig, ax = plt.subplots(figsize=(7, max(5, 0.38 * len(results_sorted))))

    ops_ranked  = [r["op"] for r in results_sorted]
    y           = np.arange(len(ops_ranked))

    for i, r in enumerate(results_sorted):
        # 95% CL bar (light)
        ax.barh(i, r["hi95"] - r["lo95"], left=r["lo95"],
                height=0.5, color="steelblue", alpha=0.4, label="95% CL" if i == 0 else "")
        # 68% CL bar (dark)
        if r["lo68"] is not None:
            ax.barh(i, r["hi68"] - r["lo68"], left=r["lo68"],
                    height=0.5, color="steelblue", alpha=0.9, label="68% CL" if i == 0 else "")

        # stat-only overlay
        if args.stat and r["op"] in stat_by_op:
            s = stat_by_op[r["op"]]
            if s["lo95"] is not None:
                ax.barh(i, s["hi95"] - s["lo95"], left=s["lo95"],
                        height=0.15, color="crimson", alpha=0.8,
                        label="95% CL (stat)" if i == 0 else "")

    ax.axvline(0, color="black", linewidth=1.0)
    ax.set_yticks(y)
    ax.set_yticklabels(ops_ranked, fontsize=9)
    ax.set_xlabel("Wilson coefficient value")
    ax.set_title("Operator ranking by 95% CL interval width\n(narrower = more sensitive)",
                 fontsize=11)
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    ax.invert_yaxis()  # most sensitive at top
    ax.grid(axis="x", linestyle=":", alpha=0.4)
    hep.cms.label(ax=ax, data=False, label="Simulation")

    plt.tight_layout()
    outpath = os.path.join(args.outdir, "operator_ranking.pdf")
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    fig.savefig(outpath.replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nSaved: {outpath}")
    print("Done.")


if __name__ == "__main__":
    main()
