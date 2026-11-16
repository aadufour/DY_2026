#!/usr/bin/env python3
"""
Re-run the 2D grid scans for selected operator pairs with per-pair boxes taken
from readapt_double_boundaries.py's log (the "new[op]=(lo,hi)" values), reusing
the original per-piece commands from condor_scan/scripts and only swapping the
--setParameterRanges argument. Pieces are bundled N per short-queue job.

The merged higgsCombine.<pair>.individual.MultiDimFit.mH125.root from the old box
is renamed to *.oldbox.root so the later hadd does not overwrite it.

Usage (in the datacards/<var> dir):
    python3 redo_pairs_scan.py --log readapt.log --method extrapolated --bundle 25
    condor_submit condor_redo/scan.sub
    # after completion: hadd per pair (see printed command)
"""
import argparse, glob, os, re

p = argparse.ArgumentParser()
p.add_argument("--log", default="readapt.log")
p.add_argument("--method", default="extrapolated", help="which readapt method to redo")
p.add_argument("--src", default="condor_scan/scripts")
p.add_argument("--out", default="condor_redo")
p.add_argument("--bundle", type=int, default=25)
a = p.parse_args()

pairs = {}
pat = re.compile(r"^(\w+?)_(\w+): .*method=(\S+)\s+new\[\w+\]=\(([-+\d.e]+),([-+\d.e]+)\) new\[\w+\]=\(([-+\d.e]+),([-+\d.e]+)\)")
for line in open(a.log):
    m = pat.match(line)
    if m and "68%=ok      95%=ok" not in line and m.group(3) == a.method:
        o1, o2 = m.group(1), m.group(2)
        pairs[(o1, o2)] = f"k_{o1}={float(m.group(4))},{float(m.group(5))}:k_{o2}={float(m.group(6))},{float(m.group(7))}"
print(len(pairs), "pairs:")
for k, v in pairs.items():
    print(" ", "_".join(k), v)

header, cmds = None, {k: [] for k in pairs}
for s in glob.glob(os.path.join(a.src, "*.sh")):
    lines = open(s).read().splitlines()
    for l in lines:
        m = re.search(r"--redefineSignalPOIs=k_(\w+),k_(\w+) ", l)
        if m and (m.group(1), m.group(2)) in pairs:
            if header is None:
                header = [x for x in lines if not x.startswith("combine") and x != "set -e"]
            cmds[(m.group(1), m.group(2))].append(
                re.sub(r"--setParameterRanges=\S+", "--setParameterRanges=" + pairs[(m.group(1), m.group(2))], l))

os.makedirs(os.path.join(a.out, "scripts"))  # refuse to overwrite an existing redo
os.makedirs(os.path.join(a.out, "logs"))
allc = [c for k in pairs for c in sorted(cmds[k])]
for k in pairs:
    assert len(cmds[k]) > 0, k
    print("_".join(k), len(cmds[k]), "pieces")
jl = []
for i in range(0, len(allc), a.bundle):
    sp = os.path.abspath(os.path.join(a.out, "scripts", f"bundle_{i // a.bundle}.sh"))
    with open(sp, "w") as f:
        f.write("\n".join(header + allc[i:i + a.bundle]) + "\n")
    os.chmod(sp, 0o755)
    jl.append(sp)
open(os.path.join(a.out, "joblist.txt"), "w").write("\n".join(jl) + "\n")
sub = open("condor_scan_bundled/scan.sub").read()
sub = sub.replace(os.path.abspath("condor_scan_bundled"), os.path.abspath(a.out))
open(os.path.join(a.out, "scan.sub"), "w").write(sub)

for k in pairs:
    f = f"higgsCombine.{'_'.join(k)}.individual.MultiDimFit.mH125.root"
    if os.path.exists(f):
        os.rename(f, f.replace(".root", ".oldbox.root"))
open(os.path.join(a.out, "pairs.txt"), "w").write("\n".join("_".join(k) for k in pairs) + "\n")
print(len(jl), "bundled jobs written to", a.out)
