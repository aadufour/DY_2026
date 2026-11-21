#!/usr/bin/env python3
"""
Re-fit single bad points of a run_2d_perpair.sh production with a more robust minimizer setting, then merge them back.

Bad points (missing / error status / stuck outlier inside 95%) come from check_2d_scans.py --list-out.
Refit setting (tested 2026-11-20 on 14 points, recovers all): strategy 2 + fallbacks Migrad0 -> Simplex -> Migrad1;
all other options as in the production command of that pair.

  refit_2d_points.py prepare <VARDIR> <badlist> [--per-job 90]   -> <VARDIR>/full/refit/{scripts,logs,joblist.txt,scan.sub}
  condor_submit <VARDIR>/full/refit/scan.sub
  refit_2d_points.py merge <VARDIR>
     per grid node keeps the lower 2dNLL of original and refit (profile = minimum over nuisances),
     writes full/refit/merged/*.root and repoints full/final/ to them (old targets in full/refit/final_before.txt)
"""
import argparse, collections, glob, os, re

p = argparse.ArgumentParser()
p.add_argument("mode", choices=["prepare", "merge"])
p.add_argument("vardir")
p.add_argument("badlist", nargs="?")
p.add_argument("--per-job", type=int, default=90)
a = p.parse_args()
vardir = os.path.abspath(a.vardir)
full = f"{vardir}/full"
R = f"{full}/refit"
FALLBACK = "--cminFallbackAlgo Minuit2,Migrad,0:0.2 --cminFallbackAlgo Minuit2,Simplex,0:0.1 --cminFallbackAlgo Minuit2,Migrad,1:0.2"
final = lambda pair: f"{full}/final/higgsCombine.{pair}.individual.MultiDimFit.mH125.root"

if a.mode == "prepare":
    pts = collections.defaultdict(list)
    for l in open(a.badlist):
        if l.strip():
            pair, n = l.split()[:2]
            pts[pair].append(int(n))
    os.makedirs(f"{R}/scripts", exist_ok=True); os.makedirs(f"{R}/logs", exist_ok=True)
    header, cmds = None, []
    for pair, ns in sorted(pts.items()):
        ps = os.path.dirname(os.path.realpath(final(pair)))          # pass2 or pass3 dir of the final scan
        line = None
        for s in glob.glob(f"{ps}/condor_scan/scripts/*.sh"):
            txt = open(s).read().splitlines()
            hit = [l for l in txt if l.startswith("combineTool.py") and f"-n .{pair}.individual.POINTS." in l]
            if hit:
                line = hit[0]
                header = header or [l for l in txt if not l.startswith(("combineTool.py", f"cd {ps}", "set -e"))]
                break
        assert line, f"no production command for {pair} in {ps}"
        assert "--cminDefaultMinimizerStrategy=0" in line and "--cminFallbackAlgo Minuit2,Migrad,0:0.2" in line, line
        line = line.replace("--cminDefaultMinimizerStrategy=0", "--cminDefaultMinimizerStrategy=2")
        line = line.replace("--cminFallbackAlgo Minuit2,Migrad,0:0.2", FALLBACK)
        line = re.sub(r"--firstPoint \d+ --lastPoint \d+ -n \S+", "", line).rstrip()
        ini = f"higgsCombine.initialFit_{pair}.MultiDimFit.mH125.root"
        if not os.path.exists(f"{R}/{ini}"):
            os.symlink(f"{vardir}/{ini}", f"{R}/{ini}")
        cmds += [f"{line} --firstPoint {n} --lastPoint {n} -n .{pair}.refit.{n} || echo REFIT_FAIL {pair} {n}" for n in sorted(ns)]
    jobs = []
    for k in range(0, len(cmds), a.per_job):
        f = f"{R}/scripts/job_{len(jobs)}.sh"
        open(f, "w").write("\n".join(header + [f"cd {R}"] + cmds[k:k + a.per_job]) + "\n"); os.chmod(f, 0o755)
        jobs.append(f)
    open(f"{R}/joblist.txt", "w").write("\n".join(jobs) + "\n")
    sub = open(f"{full}/pass2/condor_scan_bundled/scan.sub").read().replace(f"{full}/pass2/condor_scan_bundled", R)
    open(f"{R}/scan.sub", "w").write(sub)
    print(f"{len(cmds)} points, {len(pts)} pairs, {len(jobs)} jobs -> {R}/scan.sub")

else:
    import ROOT
    ROOT.gROOT.SetBatch(True)
    os.makedirs(f"{R}/merged", exist_ok=True)
    pairs = sorted({os.path.basename(f).split(".")[1] for f in glob.glob(f"{R}/higgsCombine.*.refit.*.root")})
    before = open(f"{R}/final_before.txt", "a")
    tot = collections.Counter()
    for pair in pairs:
        orig = os.path.realpath(final(pair))
        if os.path.dirname(orig) == f"{R}/merged":
            print(f"{pair}: final already merged, skipped"); continue
        o1, o2 = pair.split("_")
        ch = ROOT.TChain("limit"); ch.Add(orig)
        n_orig = ch.GetEntries()
        for f in sorted(glob.glob(f"{R}/higgsCombine.{pair}.refit.*.root")):
            ch.Add(f)
        best, keep = {}, []
        for i in range(ch.GetEntries()):
            ch.GetEntry(i)
            if ch.quantileExpected < -0.5:
                if i < n_orig:
                    keep.append(i)                                          # best fit entry of the original scan
                continue
            node = (getattr(ch, f"k_{o1}"), getattr(ch, f"k_{o2}"))
            if node not in best or ch.deltaNLL < best[node][0]:
                if node in best:
                    tot["replaced" if i >= n_orig else "?"] += 1
                elif i >= n_orig:
                    tot["recovered"] += 1
                best[node] = (ch.deltaNLL, i)
            elif i >= n_orig:
                tot["refit_not_lower"] += 1
        keep += sorted(v[1] for v in best.values())
        out = f"{R}/merged/{os.path.basename(final(pair))}"
        fo = ROOT.TFile(out, "RECREATE"); t = ch.CloneTree(0)
        for i in keep:
            ch.GetEntry(i); t.Fill()
        t.Write(); fo.Close()
        before.write(f"{pair} {orig}\n")
        os.remove(final(pair)); os.symlink(out, final(pair))
    print(f"{len(pairs)} pairs merged into final/: {dict(tot)}")
