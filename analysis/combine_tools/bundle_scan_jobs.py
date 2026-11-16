# Bundle the remaining 200-point scan pieces (condor_scan/scripts/job_N.sh) into fewer jobs:
# pieces whose output file already exists are skipped; the rest are grouped, per pair, into
# bundles of up to N_PER consecutive pieces run one after another (same output names).
import os, re, sys
N_PER = int(sys.argv[1]) if len(sys.argv) > 1 else 25
jobs = [l.strip() for l in open("condor_scan/joblist.txt")]
todo = {}
for s in jobs:
    body = open(s).read().splitlines()
    cmd = body[-1]
    name = re.search(r" -n \.(\S+)", cmd).group(1)
    a, b = re.search(r"--firstPoint (\d+) --lastPoint (\d+)", cmd).groups()
    out = f"higgsCombine.{name}.POINTS.{a}.{b}.MultiDimFit.mH125.root"
    if os.path.exists(out):
        continue
    pair = name.split(".")[0]
    todo.setdefault(pair, []).append((int(a), body))
os.makedirs("condor_scan_bundled/scripts", exist_ok=True); os.makedirs("condor_scan_bundled/logs", exist_ok=True)
paths = []
for pair, pieces in todo.items():
    pieces.sort()
    for i in range(0, len(pieces), N_PER):
        grp = pieces[i:i + N_PER]
        head = [l for l in grp[0][1][:-1] if l != "set -e"]   # env setup, without set -e: one bad piece must not stop the rest
        p = os.path.abspath(f"condor_scan_bundled/scripts/{pair}_{grp[0][0]}.sh")
        open(p, "w").write("\n".join(head + [g[1][-1] for g in grp]) + "\n")
        os.chmod(p, 0o755); paths.append(p)
open("condor_scan_bundled/joblist.txt", "w").write("\n".join(paths) + "\n")
sub = open("condor_scan/scan.sub").read()
sub = sub.replace("condor_scan/logs", "condor_scan_bundled/logs").replace("condor_scan/joblist.txt", "condor_scan_bundled/joblist.txt")
sub = re.sub(r"T3Queue = \S+", "T3Queue = short", sub)
open("condor_scan_bundled/scan.sub", "w").write(sub)
print("pieces left:", sum(len(v) for v in todo.values()), "in", len(todo), "pairs ->", len(paths), "bundled jobs")
