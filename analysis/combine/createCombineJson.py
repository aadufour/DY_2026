#!/usr/bin/env python3

"""
File: createCombineJson.py
Author: Giacomo Boldrini
Date: 2024-08-07
Description: create a json file with the bin names and operators in them
"""

from optparse import OptionParser
import json

from HiggsAnalysis.CombinedLimit.DatacardParser import *
from python.DatacardHelpers import datacardHelper

if __name__ == "__main__":
    parser = OptionParser()

    parser.add_option(
        "--datacard",
        type="string",
        dest="datacard",
        default=None,
        help="The datacard to be skimmed",
    )

    parser.add_option(
        "--output",
        type="string",
        dest="output",
        default="jsonComb.json",
        help="The output bundle folder to be created",
    )

    parser.add_option(
        "--filter",
        type="string",
        dest="filter",
        default="All",
        help="Comma separated list of operators to be kept, all others will be discarded",
    )

    (options, args) = parser.parse_args()

    DC = datacardHelper.loadDatacard(options.datacard, parser=parser)

    interesting_ops = options.filter.split(",") if options.filter != "All" else "All"

    json__ = {}

    for channel, proc, _ in DC.keyline:
        if channel not in json__:
            json__[channel] = []
        
        #if not proc.startswith("quad_"): continue 
        if "quad_" not in proc or "sm_lin_quad" in proc: continue

        op = proc.split("quad_")[1]
        if interesting_ops != "All" and op not in interesting_ops: continue
        if op not in json__[channel]: json__[channel].append(op)

    print(json__)

    with open(options.output, "w") as f:
        json.dump(json__, f, indent=4)
