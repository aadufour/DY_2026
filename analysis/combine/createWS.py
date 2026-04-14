#! /opt/anaconda3/envs/combine/bin/python3

from itertools import combinations
import sys 
import os 
import json


if __name__ == "__main__":

    # mode indicates if we want 1d, 2d or Nd workspaces
    # with N the number of operators
    mode__ = int(sys.argv[1])

    with open('metadata.json') as file:
        metadata = json.load(file)

    ops = metadata["operators"].keys()

    if mode__ == 3: mode__ = len(ops)

    # create a list whose entries will contain the 
    # operators for which we want to create the workspace
    combos = list(combinations(ops, mode__))

    model__ = "HiggsAnalysis.CombinedLimit.AnomalousCouplingEFTNegative_comb:analiticAnomalousCouplingEFTNegative_comb"
    suffix__ = ""
    if len(sys.argv) > 2:
        if sys.argv[2] == "quad":
            model__ = "HiggsAnalysis.CombinedLimit.AnomalousCouplingEFTNegative_comb:analiticAnomalousCouplingEFTNegative_comb"
        elif sys.argv[2] == "lin":
            model__ = "HiggsAnalysis.CombinedLimit.AnomalousCouplingLinearEFTNegative_comb:analiticAnomalousCouplingLinearEFTNegative_comb"
            suffix__ = "_linear"

    cmds = []
    for c in combos:
        operators = ",".join(c)
        name = "_".join(c) + suffix__
        t2w_quad = "text2workspace.py datacard.txt -P {} --PO reuseCompleteDatacards -o model_{}.root --PO addToCompleteOperators={} --PO fileCombination=jsonComb.json --PO eftOperators={} --X-pack-asympows --optimize-simpdf-constraints=cms --X-optimizeMHDependency=fixed --X-allow-no-signal --X-allow-no-background".format(model__, name, ",".join(ops), operators, operators)
        #t2w_lin = "text2workspace.py datacard.txt -P HiggsAnalysis.AnalyticAnomalousCoupling.AnomalousCouplingLinearEFTNegative:analiticAnomalousCouplingLinearEFTNegative --PO reuseCompleteDatacards -o model_linear_{}.root --PO eftOperators={}".format(name, operators)\

        cmds.append(t2w_quad)
        #cmds.append(t2w_lin)

    for cmd in cmds:
        print(cmd)
        os.system(cmd)

    print("---> Done")


