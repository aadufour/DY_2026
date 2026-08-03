#!/usr/bin/env python3

from itertools import combinations
from concurrent.futures import ThreadPoolExecutor
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
        "--doOnly",
        type="str",
        dest="doOnly",
        default="",
        help="Comma-separated list of operators to process (e.g. 'cHDD,cHWB'). "
             "Only combinations involving exclusively these operators are run.",
    )

    parser.add_option(
        "-j", "--cores",
        type="int",
        dest="cores",
        default=0,
        help="Number of parallel worker processes (one per operator combination by default). "
             "Capped at half the machine's CPU count as a safety measure; 0 (default) means "
             "'as many as needed, up to that cap'.",
    )

    (options, args) = parser.parse_args()

    # mode indicates if we want 1d, 2d or Nd workspaces
    # with N the number of operators
    mode__ = int(sys.argv[1])

    with open('metadata.json') as file:
        metadata = json.load(file)

    ops = list(metadata["operators"].keys())

    if options.doOnly:
        only_set = set(options.doOnly.split(","))
        unknown = only_set - set(ops)
        if unknown:
            print(f"ERROR: --doOnly contains unknown operators: {', '.join(sorted(unknown))}")
            sys.exit(1)
        ops = [op for op in ops if op in only_set]

    if mode__ == 3: mode__ = len(ops)
    
    suffix__ = ""
    if len(sys.argv) > 2:
        if sys.argv[2] == "lin":
            suffix__ = "_linear"

    # create a list whose entries will contain the 
    # operators for which we want to create the workspace
    combos = list(combinations(ops, mode__))
    print(combos)
    cmds = []
    for c in combos:
        op_pois = " ".join([f"k_{op}" for op in c])
        name = "_".join(c) + suffix__

        outname = name
        if options.signalPOIs:
            outname = outname + "__" + "_".join(options.signalPOIs.split(","))
            op_pois = " ".join([f"k_{op}" for op in options.signalPOIs.split(",")])
        if options.stat:
            outname += "_stat"
        if options.dataasim:
            outname += "_dataasimov"
        if options.unblind:
            outname = outname + "_unblind"
        if mode__ >=3:
            outname += "_profiled"

        fn = f"higgsCombine.{outname}.individual.MultiDimFit.mH125.root"

        cmds.append("mkEFTScan.py " + f"{fn} -p {op_pois} -maxNLL 10 -lumi 138 -cms -preliminary -o scan_{outname} -ff png pdf root")

    print(cmds)

    total_cores = os.cpu_count() or 1
    max_safe_cores = max(1, total_cores // 2)

    if options.cores <= 0:
        n_workers = min(len(cmds), max_safe_cores) if cmds else 1
    else:
        n_workers = options.cores
        if n_workers > max_safe_cores:
            print(
                f"WARNING: requested --cores {n_workers} exceeds half of the "
                f"{total_cores} available cores; capping at {max_safe_cores} to "
                f"avoid overloading the machine."
            )
            n_workers = max_safe_cores

    print(f"Running {len(cmds)} plot job(s) with {n_workers} parallel worker(s) "
          f"(machine has {total_cores} cores, safety cap {max_safe_cores}).")

    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        list(executor.map(os.system, cmds))
