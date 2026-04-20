#!/usr/bin/env python3
"""
runPlots_compare.py

Same as runPlots.py but supports overlaying a second scan for comparison.

Auto-comparison (recommended) — no file paths needed:
    # Plot full-syst as main, overlay stat-only (requires prior runScans.py --stat run)
    python3 runPlots_compare.py 1 --label "Stat + Syst" --compare-stat

    # Plot stat-only as main, overlay full-syst
    python3 runPlots_compare.py 1 --stat --label "Stat only" --compare-syst

Manual overlay (escape hatch for arbitrary comparisons):
    python3 runPlots_compare.py 1 --label "Stat + Syst" \\
        --others "higgsCombine.cHDD_stat.individual.MultiDimFit.mH125.root:2:2:Stat only"

The _stat suffix in filenames is produced automatically by runScans.py --stat.
"""

from itertools import combinations
import sys
import os
import json
from optparse import OptionParser


point_mode = {
    1: 100,
    2: 10000,
    3: 100
}


def build_outname(base_name, stat, dataasim, unblind, mode, signalPOIs):
    outname = base_name
    if signalPOIs:
        outname += "__" + "_".join(signalPOIs.split(","))
    if stat:
        outname += "_stat"
    if dataasim:
        outname += "_dataasimov"
    if unblind:
        outname += "_unblind"
    if mode >= 3:
        outname += "_profiled"
    return outname


if __name__ == "__main__":

    parser = OptionParser()

    parser.add_option(
        "--unblind",
        dest="unblind",
        default=False,
        action="store_true",
        help="Run unblind fit",
    )

    parser.add_option(
        "--data-asimov",
        dest="dataasim",
        default=False,
        action="store_true",
        help="Run fit with -t -1 --toysFreq",
    )

    parser.add_option(
        "--stat",
        dest="stat",
        default=False,
        action="store_true",
        help="Main scan is stat-only (frozen nuisances)",
    )

    parser.add_option(
        "--signalPOIs",
        type="str",
        dest="signalPOIs",
        default="",
        help="Comma separated list of parameters to run the scan on --> Define POI of profiled fits",
    )

    parser.add_option(
        "--label",
        type="str",
        dest="label",
        default="",
        help="Label for the main scan curve (passed as -mll to mkEFTScan.py)",
    )

    # --- auto-comparison flags ---

    parser.add_option(
        "--compare-stat",
        dest="compare_stat",
        default=False,
        action="store_true",
        help="Auto-overlay the stat-only (_stat) scan. Main scan must be full-syst (no --stat). "
             "Requires a prior runScans.py --stat run.",
    )

    parser.add_option(
        "--compare-syst",
        dest="compare_syst",
        default=False,
        action="store_true",
        help="Auto-overlay the full-syst scan. Main scan must be stat-only (--stat). "
             "Requires a prior runScans.py run without --stat.",
    )

    parser.add_option(
        "--compare-label",
        type="str",
        dest="compare_label",
        default="",
        help="Label for the auto-comparison curve. Defaults to 'Stat only' or 'Stat + Syst'.",
    )

    parser.add_option(
        "--compare-color",
        type="int",
        dest="compare_color",
        default=2,
        help="ROOT color index for the auto-comparison curve (default: 2 = red)",
    )

    parser.add_option(
        "--compare-linestyle",
        type="int",
        dest="compare_linestyle",
        default=2,
        help="ROOT line style for the auto-comparison curve (default: 2 = dashed)",
    )

    # --- manual escape hatch ---

    parser.add_option(
        "--others",
        type="str",
        dest="others",
        default="",
        help="Manual overlay: file.root:color:linestyle:Label (passed to mkEFTScan.py --others). "
             "Use --compare-stat / --compare-syst for the common stat vs syst case.",
    )

    (options, args) = parser.parse_args()

    if options.compare_stat and options.compare_syst:
        print("ERROR: --compare-stat and --compare-syst are mutually exclusive.")
        sys.exit(1)

    if options.compare_stat and options.stat:
        print("ERROR: --compare-stat overlays the stat-only scan, but --stat says the main scan is already stat-only.")
        sys.exit(1)

    if options.compare_syst and not options.stat:
        print("ERROR: --compare-syst overlays the full-syst scan, but --stat is not set for the main scan.")
        sys.exit(1)

    mode__ = int(sys.argv[1])

    with open('metadata.json') as file:
        metadata = json.load(file)

    ops = list(metadata["operators"].keys())

    if mode__ == 3:
        mode__ = len(ops)

    suffix__ = ""
    if len(sys.argv) > 2:
        if sys.argv[2] == "lin":
            suffix__ = "_linear"

    combos = list(combinations(ops, mode__))
    print(combos)
    cmds = []

    for c in combos:
        op_pois = " ".join([f"k_{op}" for op in c])
        base_name = "_".join(c) + suffix__

        if options.signalPOIs:
            op_pois = " ".join([f"k_{op}" for op in options.signalPOIs.split(",")])

        outname = build_outname(
            base_name,
            stat=options.stat,
            dataasim=options.dataasim,
            unblind=options.unblind,
            mode=mode__,
            signalPOIs=options.signalPOIs,
        )

        fn = f"higgsCombine.{outname}.individual.MultiDimFit.mH125.root"

        cmd = f"mkEFTScan.py {fn} -p {op_pois} -maxNLL 10 -lumi 138 -cms -preliminary -o scan_{outname} -ff png pdf root"

        # Determine main label — required for legend to appear in mkEFTScan.py
        main_label = options.label
        if not main_label:
            if options.compare_stat:
                main_label = "Stat + Syst"
            elif options.compare_syst:
                main_label = "Stat only"
        if main_label:
            cmd += f' -ml "{main_label}"'

        # build auto-comparison --others string
        others_str = options.others  # start with manual value (may be empty)

        if options.compare_stat or options.compare_syst:
            other_stat = options.compare_stat  # True  → compare file has _stat
                                               # False → compare file has no _stat (full syst)
            other_outname = build_outname(
                base_name,
                stat=other_stat,
                dataasim=options.dataasim,
                unblind=options.unblind,
                mode=mode__,
                signalPOIs=options.signalPOIs,
            )
            other_fn = f"higgsCombine.{other_outname}.individual.MultiDimFit.mH125.root"

            if not os.path.isfile(other_fn):
                print(f"WARNING: comparison file not found, skipping overlay for {c}: {other_fn}")
            else:
                default_label = "Stat only" if options.compare_stat else "Stat + Syst"
                label = options.compare_label or default_label
                auto_others = f"{other_fn}:{options.compare_color}:{options.compare_linestyle}:{label}"

                # merge with manual --others if both were given
                if others_str:
                    others_str = others_str + " " + auto_others
                else:
                    others_str = auto_others

        if others_str:
            cmd += f' --others "{others_str}"'

        cmds.append(cmd)

    print(cmds)
    for command in cmds:
        os.system(command)
