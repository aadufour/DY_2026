import json, glob, re
ops = ",".join(json.load(open("metadata.json"))["operators"])
m = "HiggsAnalysis.AnalyticAnomalousCoupling.AnomalousCouplingMorphing_comb:analiticAnomalousCouplingMorphing_comb"
n = 0
for s in glob.glob("condor_initial/scripts/job_*.sh"):
    t = open(s).read()
    if "text2workspace" in t:
        continue
    name = re.search(r"model_(\S+)\.root --saveWorkspace", t).group(1)
    a, b = name.split("_")
    t2w = ("text2workspace.py datacard.txt -P %s --PO reuseCompleteDatacards -o model_%s.root --PO addToCompleteOperators=%s "
           "--PO fileCombination=jsonComb.json --PO eftOperators=%s,%s --X-pack-asympows --optimize-simpdf-constraints=cms "
           "--X-optimizeMHDependency=fixed\n" % (m, name, ops, a, b))
    open(s, "w").write(t.replace("combine -M MultiDimFit", t2w + "combine -M MultiDimFit", 1))
    n += 1
print(n, "scripts patched")
