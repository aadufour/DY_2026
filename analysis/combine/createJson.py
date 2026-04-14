#! /opt/anaconda3/envs/combine/bin/python3

from optparse import OptionParser
import os, re
import json


def add_parsing_opts():
    """Function with base parsing arguments used by any script"""
    parser = OptionParser(
        usage="python create_json.py [options]",
        description="Create json file from datacard inputs",
    )
    parser.add_option(
        "--datacard",
        dest="datacard",
        default="datacard.txt",
        help="Datacard",
    )
    parser.add_option(
        "--binname",
        dest="binname",
        default="sm_lin_quad",
        help="Binname used to generate the datacards.",
    )
    parser.add_option(
        "--filterOperators",
        dest="filterOperators",
        default="",
        help="Keep only these operators in the metadata file",
    )
    
    parser.add_option(
        "--ranges",
        dest="ranges",
        default="",
        help="Ranges for the operators as k_cW=-10,10:k_cHWB=-10,10:... default everything -1,1",
    )
    return parser.parse_args()


def get_operators(text, binname, ranges, filter_=[]):
    operator_matches = [
        operator.split(binname)[-1] for operator in re.findall(binname + "\S*", text)
    ]
    operator_matches = list(filter(lambda k: "mix" not in k, operator_matches))
    if filter_:
        print("Filtering")
        operator_matches = list(filter(lambda k: k in filter_, operator_matches))
    d = {}
    for op in operator_matches:
        if op in ranges:
            d[op] = [ranges[op][0], ranges[op][1]]
        else:
            d[op] = [-1.0, 1.0]

    print(d)
    return d

def get_nuisances(text):
    nuisances_matches = []
    for line in text.split("\n"):
        if "lnN" in line: 
            nuisances_matches.append( line.split("lnN")[0].strip() )
        elif "shape" in line and "shapes" not in line: 
            nuisances_matches.append( line.split("shape")[0].strip() )
        elif "rateParam" in line: 
            nuisances_matches.append( line.split("rateParam")[0].strip() )
    return nuisances_matches 



if __name__ == "__main__":
    opts, args = add_parsing_opts()

    binname = opts.binname

    opranges = [i.split("=") for i in opts.ranges.split(":")]
    opranges_dict = {}
    for entry in opranges:
        if entry != ['']:
            op = entry[0]
            r = entry[1]
            opranges_dict[entry[0]] = [float(r.split(",")[0]), float(r.split(",")[1])]

    # Add the underscore at the end just in case analyzer does not provide it.
    # This way the operators do not have an "_" in front when parsed.
    if binname[-1] != "_":
        binname += "_"

    # Open the datacard
    datacard = open(opts.datacard, "r")
    datacard_text = datacard.read()
    datacard.close()

    filterOps = [i for i in opts.filterOperators.split(",") if i != '']

    # Now open a json file
    out_json = open("metadata.json", "w")

    # Fill out the metadata
    metadata = {
        "analysis": os.getcwd().split("/")[-1],
        "card": opts.datacard,
        "operators": get_operators(datacard_text, binname, opranges_dict, filter_=filterOps),
        "nuisances": get_nuisances(datacard_text),
    }

    out_json.write(json.dumps(metadata, indent=4))
