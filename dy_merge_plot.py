"""
dy_merge_plot.py

Merge multiple LHE files and plot 1D DY histograms.
Variables and classes are imported from dy_plotter.py.

Two ways to specify input files:

  1. Explicit list of LHE files:
       python dy_merge_plot.py f1.lhe f2.lhe f3.lhe --out plots/

  2. Mass-binned pattern (specific drell yan case with binned mll in different folders):
       python dy_merge_plot.py --base /path/to/DY_all --bins 50 120 200 400 --out plots/
        reads /path/to/DY_all_50_120/myLHE/unweighted_events.lhe
               /path/to/DY_all_120_200/myLHE/unweighted_events.lhe
               /path/to/DY_all_200_400/myLHE/unweighted_events.lhe

Other options
-------------
  --weight SM          use a named reweight (e.g. 'SM'); default = event weight
  --max-events 5000    cap per file (useful for quick tests)
"""

import argparse
import os
import sys

# Import everything from the single-file plotter
sys.path.insert(0, os.path.dirname(__file__))
from dy_plotter import VARIABLES, fill_and_plot, parse_lhe


# ── Input file resolution ─────────────────────────────────────────────────────

def files_from_pattern(base: str, bins: list[int]) -> list[str]:
    """
    Build file paths from a base directory prefix and a list of mass edges.

    Example:
        base = "/data/DY_all"
        bins = [50, 120, 200]
        extracts ["/data/DY_all_50_120/myLHE/unweighted_events.lhe",
           "/data/DY_all_120_200/myLHE/unweighted_events.lhe"]
    """
    paths = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        path = f"{base}_{lo}_{hi}/myLHE/unweighted_events.lhe"
        paths.append(path)
    return paths


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Merge multiple LHE files and plot 1D DY observables."
    )

    # --- input: explicit files OR pattern ---
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "lhe_files", nargs="*", metavar="LHE",
        help="Explicit list of LHE files to merge.",
    )
    input_group.add_argument(
        "--base", metavar="PATH_PREFIX",
        help="Base path prefix for mass-binned samples (use with --bins).",
    )

    parser.add_argument(
        "--bins", nargs="+", type=int, metavar="EDGE",
        help="Mass bin edges in GeV, e.g. --bins 50 120 200 400 (use with --base).",
    )
    parser.add_argument(
        "--weight", default=None,
        help="Named weight key (e.g. 'SM'). Defaults to event-level weight.",
    )
    parser.add_argument(
        "--out", default="dy_plots",
        help="Output directory (default: dy_plots/).",
    )
    parser.add_argument(
        "--max-events", type=int, default=None,
        help="Max events per file (for quick tests).",
    )

    args = parser.parse_args()

    # resolve file list
    if args.base:
        if not args.bins or len(args.bins) < 2:
            parser.error("--base requires --bins with at least two edges.")
        lhe_files = files_from_pattern(args.base, args.bins)
    else:
        lhe_files = args.lhe_files
        if not lhe_files:
            parser.error("Provide at least one LHE file.")

    # read and merge
    all_events = []
    for path in lhe_files:
        if not os.path.exists(path):
            print(f"  WARNING: not found, skipping — {path}")
            continue
        print(f"Reading {path}")
        all_events.extend(parse_lhe(path, args.weight, args.max_events))

    if not all_events:
        print("No events found.")
        return

    print(f"\nTotal: {len(all_events)} events from {len(lhe_files)} file(s)")
    print(f"Plotting {len(VARIABLES)} variables to {args.out}/\n")
    fill_and_plot(all_events, VARIABLES, args.out)
    print("\nDone.")


if __name__ == "__main__":
    main()
