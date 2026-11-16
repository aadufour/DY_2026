#!/usr/bin/env python3
"""
Copy the 2D condor job commands of selected operator pairs from one variable's
datacard dir (e.g. mll) to another's (e.g. triple_diff), swapping the working
directory. Same model/fit options and scan boxes as the source.

Writes in <dst>:
  <out>_initial/  one job per pair: text2workspace + initial fit
  <out>_scan/     scan pieces bundled N per job
Both with scan.sub (short queue, from the source's condor_scan_bundled/scan.sub).
Submit initial first, scans only after the initial fits exist.

Usage (from anywhere):
  copy_pair_jobs.py --src .../mll --dst .../triple_diff --pairs cHj3_cHu,clj3_cQe --out sel --bundle 10
"""
import argparse, glob, os, re

p = argparse.ArgumentParser()
p.add_argument("--src", required=True)
p.add_argument("--dst", required=True)
p.add_argument("--pairs", required=True)
p.add_argument("--out", default="sel")
p.add_argument("--bundle", type=int, default=10)
a = p.parse_args()
src, dst = os.path.abspath(a.src), os.path.abspath(a.dst)
pairs = [tuple(x.split("_")) for x in a.pairs.split(",")]


def collect(subdir):
    header, cmds = None, {k: [] for k in pairs}
    for s in sorted(glob.glob(os.path.join(src, subdir, "scripts", "*.sh"))):
        lines = open(s).read().splitlines()
        hit = False
        for l in lines:
            m = re.search(r"--redefineSignalPOIs=k_(\w+),k_(\w+) ", l)
            if m and (m.group(1), m.group(2)) in cmds:
                hit = True
                key = (m.group(1), m.group(2))
        if not hit:
            continue
        if header is None:
            header = [x.replace(f"cd {src}", f"cd {dst}") for x in lines
                      if not x.startswith(("combine", "text2workspace")) and x != "set -e"]
        body = [x for x in lines if x.startswith(("combine", "text2workspace"))]
        cmds[key].append(body)
    for k, v in cmds.items():
        assert v, f"no {subdir} job found for {k}"
    return header, cmds


def write(name, header, jobs, setE):
    d = os.path.join(dst, f"{a.out}_{name}")
    os.makedirs(os.path.join(d, "scripts"))  # refuse to overwrite
    os.makedirs(os.path.join(d, "logs"))
    jl = []
    for i, body in enumerate(jobs):
        sp = os.path.join(d, "scripts", f"job_{i}.sh")
        open(sp, "w").write("\n".join(header[:1] + (["set -e"] if setE else []) + header[1:] + body) + "\n")
        os.chmod(sp, 0o755)
        jl.append(sp)
    open(os.path.join(d, "joblist.txt"), "w").write("\n".join(jl) + "\n")
    sub = open(os.path.join(src, "condor_scan_bundled", "scan.sub")).read()
    sub = sub.replace(os.path.join(src, "condor_scan_bundled"), d)
    open(os.path.join(d, "scan.sub"), "w").write(sub)
    print(f"{len(jl)} jobs -> {d}")


h, ini = collect("condor_initial")
write("initial", h, [ini[k][0] for k in pairs], True)

h, scn = collect("condor_scan")
pieces = [c for k in pairs for job in scn[k] for c in job]
write("scan", h, [pieces[i:i + a.bundle] for i in range(0, len(pieces), a.bundle)], False)
for k in pairs:
    print("_".join(k), sum(len(j) for j in scn[k]), "scan pieces")
