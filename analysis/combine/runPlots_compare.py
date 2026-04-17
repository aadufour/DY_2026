#!/usr/bin/env python3
"""
runPlots_compare.py

Same as runPlots.py but supports overlaying a second scan via --others and --label.

Usage example (stat-only vs full systematics):
    python3 runPlots_compare.py 1 --label "Stat only" --stat \
        --others "higgsCombine.cHDD.individual.MultiDimFit.mH125.root:2:2:Syst"

The --others argument is passed directly to mkEFTScan.py:
    file.root:color:linestyle:"Label"
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
        help="Freeze all constrained nuisances",
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
        help="Label for the main scan in the plot (passed as -mll to mkEFTScan.py)",
    )

    parser.add_option(
        "--others",
        type="str",
        dest="others",
        default="",
        help="Overlay files for comparison: file.root:color:linestyle:\"Label\" (passed to mkEFTScan.py --others)",
    )

    (options, args) = parser.parse_args()

    mode__ = int(sys.argv[1])

    with open('metadata.json') as file:
        metadata = json.load(file)

    ops = metadata["operators"].keys()

    if mode__ == 3: mode__ = len(ops)

    suffix__ = ""
    if len(sys.argv) > 2:
        if sys.argv[2] == "lin":
            suffix__ = "_linear"

    combos = list(combinations(ops, mode__))
    print(combos)
    cmds = []
    for c in combos:
        ops = " ".join([f"k_{op}" for op in c])
        name = "_".join(c) + suffix__

        outname = name
        if options.signalPOIs:
            outname = outname + "__" + "_".join(options.signalPOIs.split(","))
            ops = " ".join([f"k_{op}" for op in options.signalPOIs.split(",")])
        if options.stat:
            outname += "_stat"
        if options.dataasim:
            outname += "_dataasimov"
        if options.unblind:
            outname = outname + "_unblind"
        if mode__ >= 3:
            outname += "_profiled"

        fn = f"higgsCombine.{outname}.individual.MultiDimFit.mH125.root"

        cmd = f"mkEFTScan.py {fn} -p {ops} -maxNLL 10 -lumi 138 -cms -preliminary -o scan_{outname} -ff png pdf root"
        if options.label:
            cmd += f' -mll "{options.label}"'
        if options.others:
            cmd += f' --others "{options.others}"'

        cmds.append(cmd)

    print(cmds)
    for command in cmds:
        os.system(command)
