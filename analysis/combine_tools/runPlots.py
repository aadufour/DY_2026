#!/usr/bin/env python3

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

    (options, args) = parser.parse_args()
    
    # mode indicates if we want 1d, 2d or Nd workspaces
    # with N the number of operators
    mode__ = int(sys.argv[1])
    
    with open('metadata.json') as file:
        metadata = json.load(file)

    ops = metadata["operators"].keys()

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
        if mode__ >=3:
            outname += "_profiled"
            
        fn = f"higgsCombine.{outname}.individual.MultiDimFit.mH125.root"

        cmds.append("mkEFTScan.py " + f"{fn} -p {ops} -maxNLL 10 -lumi 138 -cms -preliminary -o scan_{outname} -ff png pdf root")
       
    print(cmds)
    for command in cmds:
        os.system(command)
