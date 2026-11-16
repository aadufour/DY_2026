# Bundle 200-point scan pieces written by runScans.py --condor (<src>/scripts/job_N.sh) into fewer jobs:
# pieces whose output file already exists are skipped; the rest are grouped, per pair, into
# bundles of up to N_PER consecutive pieces run one after another (same output names), short queue.
# Usage: bundle_scan_jobs.py [N_PER=25] [src=condor_scan] [dst=<src>_bundled]
import os, re, sys
N_PER = int(sys.argv[1]) if len(sys.argv) > 1 else 25
SRC = sys.argv[2].rstrip("/") if len(sys.argv) > 2 else "condor_scan"
DST = sys.argv[3].rstrip("/") if len(sys.argv) > 3 else SRC + "_bundled"
jobs = [l.strip() for l in open(f"{SRC}/joblist.txt")]
todo = {}
for s in jobs:
    body = open(s).read().splitlines()
    cmd = body[-1]
    name = re.search(r" -n \.(\S+)", cmd).group(1)   # <pair>.individual.POINTS.<a>.<b>
    a, b = re.search(r"--firstPoint (\d+) --lastPoint (\d+)", cmd).groups()
    if os.path.exists(f"higgsCombine.{name}.MultiDimFit.mH125.root"):
        continue
    pair = name.split(".")[0]
    todo.setdefault(pair, []).append((int(a), body))
os.makedirs(f"{DST}/scripts"); os.makedirs(f"{DST}/logs")   # refuse to overwrite
paths = []
for pair, pieces in todo.items():
    pieces.sort()
    for i in range(0, len(pieces), N_PER):
        grp = pieces[i:i + N_PER]
        head = [l for l in grp[0][1][:-1] if l != "set -e"]   # env setup, without set -e: one bad piece must not stop the rest
        p = os.path.abspath(f"{DST}/scripts/{pair}_{grp[0][0]}.sh")
        open(p, "w").write("\n".join(head + [g[1][-1] for g in grp]) + "\n")
        os.chmod(p, 0o755); paths.append(p)
open(f"{DST}/joblist.txt", "w").write("\n".join(paths) + "\n")
sub = open(f"{SRC}/scan.sub").read()
sub = sub.replace(os.path.abspath(f"{SRC}/logs"), os.path.abspath(f"{DST}/logs")).replace(os.path.abspath(f"{SRC}/joblist.txt"), os.path.abspath(f"{DST}/joblist.txt"))
sub = sub.replace(f"{SRC}/logs", f"{DST}/logs").replace(f"{SRC}/joblist.txt", f"{DST}/joblist.txt")
sub = re.sub(r"T3Queue = \S+", "T3Queue = short", sub)
open(f"{DST}/scan.sub", "w").write(sub)
print("pieces left:", sum(len(v) for v in todo.values()), "in", len(todo), "pairs ->", len(paths), "bundled jobs in", DST)
