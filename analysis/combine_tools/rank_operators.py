#!/usr/bin/env python3
"""
rank_operators.py
=================
Extract 68% and 95% CL intervals from combine likelihood scan ROOT files
and produce a ranked summary plot + table.

Two supported file formats (auto-detected):

  TTree format (combine output):
    higgsCombine.{op}.individual.MultiDimFit.mH125.root  (default pattern)
    TTree "limit" with branches deltaNLL, k_{op}
    deltaNLL = ΔNLL  →  68% CL: deltaNLL < 0.5,  95% CL: deltaNLL < 1.92

  TGraph format (mkEFTScan.py output):
    scan_{op}.root  (use --pattern "scan_{op}.root")
    TGraph keys "Stat + Syst" / "Stat only"
    x = k_{op},  y = Δ(-2lnL)  →  68% CL: y < 1,  95% CL: y < 3.84

Usage:
    # combine output (default)
    python3 rank_operators.py --indir . --outdir ranking --stat

    # mkEFTScan.py output (old validation set)
    python3 rank_operators.py --indir . --outdir ranking --pattern "scan_{op}.root" \
        --tgraph-key-syst "Stat + Syst" --tgraph-key-stat "Stat only"
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


def extract_interval_ttree(fname, poi):
    """Read combine TTree output. deltaNLL = ΔNLL (not ×2)."""
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

    mask = dnll > 0   # drop best-fit point (deltaNLL==0) and failed points
    if mask.sum() == 0:
        return None
    dnll, vals = dnll[mask], vals[mask]

    def cl_interval(threshold):
        inside = vals[dnll < threshold]
        if len(inside) == 0:
            return None, None
        return float(inside.min()), float(inside.max())

    lo68, hi68 = cl_interval(CL68)
    lo95, hi95 = cl_interval(CL95)
    return lo68, hi68, lo95, hi95


def extract_interval_tgraph(fname, key):
    """Read mkEFTScan.py TGraph output. y = Δ(-2lnL) directly."""
    if not os.path.isfile(fname):
        return None
    try:
        with uproot.open(fname) as f:
            # try exact key, then strip version suffix
            if key in f:
                g = f[key]
            elif f"{key};1" in f:
                g = f[f"{key};1"]
            else:
                available = list(f.keys())
                print(f"  WARNING: key '{key}' not found in {fname}. Available: {available}")
                return None
            vals = g.member("fX")
            dnll2 = g.member("fY")   # this is Δ(-2lnL) already
    except Exception as e:
        print(f"  WARNING: could not read {fname}: {e}")
        return None

    if len(vals) == 0:
        return None

    def cl_interval(threshold):
        inside = vals[dnll2 < threshold]
        if len(inside) == 0:
            return None, None
        return float(inside.min()), float(inside.max())

    lo68, hi68 = cl_interval(1.0)    # Δ(-2lnL) < 1
    lo95, hi95 = cl_interval(3.84)   # Δ(-2lnL) < 3.84
    return lo68, hi68, lo95, hi95


def extract_interval(fname, poi, tgraph_key=None):
    """Auto-detect format and extract interval."""
    if tgraph_key is not None:
        return extract_interval_tgraph(fname, tgraph_key)
    return extract_interval_ttree(fname, poi)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--indir",    default=".",
                        help="Directory containing scan ROOT files (default: .)")
    parser.add_argument("--outdir",   default="ranking",
                        help="Output directory for plots and table")
    parser.add_argument("--metadata", default="metadata.json",
                        help="metadata.json with operator scan ranges")
    parser.add_argument("--stat",     action="store_true",
                        help="Also overlay stat-only intervals")
    parser.add_argument("--operators", nargs="+", default=None,
                        help="Subset of operators (default: all 27)")
    parser.add_argument("--pattern",  default="higgsCombine.{op}.individual.MultiDimFit.mH125.root",
                        help="Filename pattern with {op} placeholder (default: combine TTree pattern). "
                             "Use 'scan_{op}.root' for mkEFTScan.py TGraph output.")
    parser.add_argument("--pattern-stat", default=None,
                        help="Stat-only filename pattern. Default: auto-derived from --pattern "
                             "by inserting '_stat' before '.root'.")
    parser.add_argument("--tgraph-key-syst", default=None,
                        help="TGraph key for full-syst scan (e.g. 'Stat + Syst'). "
                             "If set, files are read as TGraph instead of TTree.")
    parser.add_argument("--tgraph-key-stat", default=None,
                        help="TGraph key for stat-only scan (e.g. 'Stat only').")
    parser.add_argument("--wide", action="store_true",
                        help="Wide layout: operators on x-axis, bars vertical (good for presentations).")
    parser.add_argument("--log", action="store_true",
                        help="Log scale on the constraint axis. Bars show the 95%% (and 68%%) CL "
                             "half-width — i.e. the reach on |k| — since a log scale cannot span zero.")
    parser.add_argument("--exclude", nargs="+", default=None, metavar="OP",
                        help="Operators to remove from the plot (e.g. --exclude cHDD cHWB).")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    ops = args.operators or OPERATORS

    # Derive stat pattern if not given
    pattern_stat = args.pattern_stat or args.pattern.replace(".root", "_stat.root")
    # TGraph mode if key given
    tgraph_syst = args.tgraph_key_syst
    tgraph_stat = args.tgraph_key_stat

    # Load metadata for scan ranges (optional)
    meta_ops = {}
    if os.path.isfile(args.metadata):
        with open(args.metadata) as f:
            meta = json.load(f)
        meta_ops = meta.get("operators", {})

    results = []
    results_stat = []

    for op in ops:
        poi        = f"k_{op}"
        fname      = os.path.join(args.indir, args.pattern.format(op=op))
        fname_stat = os.path.join(args.indir, pattern_stat.format(op=op))

        iv = extract_interval(fname, poi, tgraph_key=tgraph_syst)
        if iv is None:
            print(f"  SKIP {op}: no valid scan points in {fname}")
            continue
        lo68, hi68, lo95, hi95 = iv
        width95 = (hi95 - lo95) if (lo95 is not None and hi95 is not None) else None
        results.append({"op": op, "lo68": lo68, "hi68": hi68,
                         "lo95": lo95, "hi95": hi95, "width95": width95})

        if args.stat:
            iv_s = extract_interval(fname_stat, poi, tgraph_key=tgraph_stat)
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

    # Apply --exclude filter
    if args.exclude:
        excluded = set(args.exclude)
        results_sorted = [r for r in results_sorted if r["op"] not in excluded]

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

    ops_ranked = [r["op"] for r in results_sorted]
    n = len(ops_ranked)

    if args.wide:
        # ----- Wide layout: operators on x-axis, vertical bars ----------
        fig, ax = plt.subplots(figsize=(max(10, 0.55 * n), 6))
        x = np.arange(n)
        bar_w = 0.6

        if args.log:
            # Log mode: bar height = half-width of interval (reach on |k|)
            for i, r in enumerate(results_sorted):
                hw95 = r["width95"] / 2 if r["width95"] else None
                hw68 = (r["hi68"] - r["lo68"]) / 2 if (r["lo68"] is not None and r["hi68"] is not None) else None
                if hw95:
                    ax.bar(i, hw95, width=bar_w, color="steelblue", alpha=0.4,
                           label="95% CL half-width" if i == 0 else "")
                if hw68:
                    ax.bar(i, hw68, width=bar_w, color="steelblue", alpha=0.9,
                           label="68% CL half-width" if i == 0 else "")
                if args.stat and r["op"] in stat_by_op:
                    s = stat_by_op[r["op"]]
                    hw95s = s["width95"] / 2 if s["width95"] else None
                    if hw95s:
                        ax.bar(i, hw95s, width=bar_w * 0.3, color="crimson", alpha=0.8,
                               label="95% CL half-width (stat)" if i == 0 else "")
            ax.set_yscale("log")
            ax.set_ylabel(r"95% CL half-width $|k|$")
        else:
            # Linear mode: floating bar from lo to hi
            for i, r in enumerate(results_sorted):
                if r["lo95"] is not None and r["hi95"] is not None:
                    ax.bar(i, r["hi95"] - r["lo95"], bottom=r["lo95"],
                           width=bar_w, color="steelblue", alpha=0.4,
                           label="95% CL" if i == 0 else "")
                if r["lo68"] is not None and r["hi68"] is not None:
                    ax.bar(i, r["hi68"] - r["lo68"], bottom=r["lo68"],
                           width=bar_w, color="steelblue", alpha=0.9,
                           label="68% CL" if i == 0 else "")
                if args.stat and r["op"] in stat_by_op:
                    s = stat_by_op[r["op"]]
                    if s["lo95"] is not None and s["hi95"] is not None:
                        ax.bar(i, s["hi95"] - s["lo95"], bottom=s["lo95"],
                               width=bar_w * 0.3, color="crimson", alpha=0.8,
                               label="95% CL (stat)" if i == 0 else "")
            ax.axhline(0, color="black", linewidth=1.0)
            ax.set_ylabel(r"$k$")

        ax.set_xticks(x)
        ax.set_xticklabels(ops_ranked, fontsize=8, rotation=45, ha="right")
        ax.legend(frameon=False, fontsize=9, loc="upper right")
        ax.grid(axis="y", linestyle=":", alpha=0.4)
        hep.cms.label(ax=ax, data=False, label="Simulation")

    else:
        # ----- Tall layout: operators on y-axis, horizontal bars (default) -
        fig, ax = plt.subplots(figsize=(7, max(5, 0.38 * n)))
        y = np.arange(n)

        if args.log:
            # Log mode: bar width = half-width of interval (reach on |k|)
            for i, r in enumerate(results_sorted):
                hw95 = r["width95"] / 2 if r["width95"] else None
                hw68 = (r["hi68"] - r["lo68"]) / 2 if (r["lo68"] is not None and r["hi68"] is not None) else None
                if hw95:
                    ax.barh(i, hw95, height=0.5, color="steelblue", alpha=0.4,
                            label="95% CL half-width" if i == 0 else "")
                if hw68:
                    ax.barh(i, hw68, height=0.5, color="steelblue", alpha=0.9,
                            label="68% CL half-width" if i == 0 else "")
                if args.stat and r["op"] in stat_by_op:
                    s = stat_by_op[r["op"]]
                    hw95s = s["width95"] / 2 if s["width95"] else None
                    if hw95s:
                        ax.barh(i, hw95s, height=0.15, color="crimson", alpha=0.8,
                                label="95% CL half-width (stat)" if i == 0 else "")
            ax.set_xscale("log")
            ax.set_xlabel(r"95% CL half-width $|k|$")
        else:
            # Linear mode: floating bar from lo to hi
            for i, r in enumerate(results_sorted):
                ax.barh(i, r["hi95"] - r["lo95"], left=r["lo95"],
                        height=0.5, color="steelblue", alpha=0.4,
                        label="95% CL" if i == 0 else "")
                if r["lo68"] is not None:
                    ax.barh(i, r["hi68"] - r["lo68"], left=r["lo68"],
                            height=0.5, color="steelblue", alpha=0.9,
                            label="68% CL" if i == 0 else "")
                if args.stat and r["op"] in stat_by_op:
                    s = stat_by_op[r["op"]]
                    if s["lo95"] is not None:
                        ax.barh(i, s["hi95"] - s["lo95"], left=s["lo95"],
                                height=0.15, color="crimson", alpha=0.8,
                                label="95% CL (stat)" if i == 0 else "")
            ax.axvline(0, color="black", linewidth=1.0)
            ax.set_xlabel(r"$k$")

        ax.set_yticks(y)
        ax.set_yticklabels(ops_ranked, fontsize=9)
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
