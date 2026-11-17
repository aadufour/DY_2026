#!/usr/bin/env python3

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
    # --ops a,b,c : only these operators are POIs (e.g. a profiled fit on a subset, mode 3);
    # addToCompleteOperators keeps ALL operators of the datacard
    sel = list(ops)
    if "--ops" in sys.argv:
        only = set(sys.argv[sys.argv.index("--ops") + 1].split(","))
        sel = [o for o in ops if o in only]

    if mode__ == 3: mode__ = len(sel)

    # create a list whose entries will contain the 
    # operators for which we want to create the workspace
    combos = list(combinations(sel, mode__))

    model__ = "HiggsAnalysis.AnalyticAnomalousCoupling.AnomalousCouplingMorphing_comb:analiticAnomalousCouplingMorphing_comb"
    suffix__ = ""
    if len(sys.argv) > 2 and sys.argv[2] != "--ops":
        if sys.argv[2] == "quad":
            model__ = "HiggsAnalysis.AnalyticAnomalousCoupling.AnomalousCouplingMorphing_comb:analiticAnomalousCouplingMorphing_comb"
        elif sys.argv[2] == "lin":
            model__ = "HiggsAnalysis.AnalyticAnomalousCoupling.AnomalousCouplingLinearEFTNegative_comb:analiticAnomalousCouplingLinearEFTNegative_comb"
            suffix__ = "_linear"

    cmds = []
    for c in combos:
        operators = ",".join(c)
        name = "_".join(c) + suffix__
        t2w_quad = "text2workspace.py datacard.txt -P {} --PO reuseCompleteDatacards -o model_{}.root --PO addToCompleteOperators={} --PO fileCombination=jsonComb.json --PO eftOperators={} --X-pack-asympows --optimize-simpdf-constraints=cms --X-optimizeMHDependency=fixed".format(model__, name, ",".join(ops), operators)
        #t2w_lin = "text2workspace.py datacard.txt -P HiggsAnalysis.AnalyticAnomalousCoupling.AnomalousCouplingLinearEFTNegative:analiticAnomalousCouplingLinearEFTNegative --PO reuseCompleteDatacards -o model_linear_{}.root --PO eftOperators={}".format(name, operators)\

        cmds.append(t2w_quad)
        #cmds.append(t2w_lin)

    for cmd in cmds:
        print(cmd)
        os.system(cmd)

    print("---> Done")


