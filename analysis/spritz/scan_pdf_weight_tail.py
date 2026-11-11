#!/usr/bin/env python3
"""
scan_pdf_weight_tail.py

Scan NanoAOD samples for events with unphysical PDF weights: for each event, take the largest
|LHEPdfWeight[i] / LHEPdfWeight[0] - 1| over the Hessian members (default 1..100).
Prints per sample the fraction of events above several thresholds, the largest effect of the
outlier events on one member's total yield, and the worst events (run:lumi:event).

Needs only uproot, awkward, numpy (+ a valid grid proxy for remote files).

Usage:
    # files from a spritz fileset.json, samples matching a pattern
    python3 scan_pdf_weight_tail.py --fileset data/fileset.json --pattern DYJetsToMuMu_M- [--max-events 300000] [--top 5]
    # or explicit files
    python3 scan_pdf_weight_tail.py --files root://.../a.root root://.../b.root
"""
import argparse
import json
import numpy as np
import awkward as ak
import uproot

p = argparse.ArgumentParser()
p.add_argument("--fileset", help="spritz fileset.json")
p.add_argument("--pattern", default="", help="only samples whose name contains this")
p.add_argument("--files", nargs="+", help="explicit list of files (instead of --fileset)")
p.add_argument("--max-events", type=int, default=300000, help="per sample (0 = all)")
p.add_argument("--members", default="1:101", help="Hessian member range, python slice (default 1:101)")
p.add_argument("--outlier", type=float, default=10, help="threshold on |r-1| used for the effect estimate")
p.add_argument("--top", type=int, default=5, help="worst events to print per sample")
p.add_argument("--prefer", default="fnal", help="preferred redirector substring when several replicas are given")
a = p.parse_args()

lo, hi = map(int, a.members.split(":"))


def open_any(replicas):
    for f in sorted(replicas, key=lambda x: 0 if a.prefer in x else 1):
        try:
            return uproot.open(f, timeout=120)["Events"]
        except Exception as e:
            print(f"   could not open {f.split('//')[1]}: {str(e)[:60]}")
    return None


samples = {}
if a.files:
    samples["files"] = [[f] for f in a.files]
else:
    d = json.load(open(a.fileset))
    for k, v in d.items():
        if a.pattern in k:
            samples[k] = [x["path"] if isinstance(x, dict) else [x] for x in v["files"]]

cols = ["run", "luminosityBlock", "event", "genWeight", "LHEPdfWeight"]
print("%-45s %8s %7s %7s %8s %6s %5s  %s" % ("sample", "events", ">1", ">10", ">100", ">1e3", ">1e4", f"effect of |r-1|>{a.outlier:g} on one member"))
for name, files in samples.items():
    parts, n = [], 0
    for reps in files:
        t = open_any(reps)
        if t is None:
            continue
        stop = None if a.max_events == 0 else a.max_events - n
        parts.append(t.arrays(cols, entry_stop=stop))
        n += len(parts[-1])
        if a.max_events and n >= a.max_events:
            break
    if not parts:
        print(name, "no file could be read")
        continue
    ev = ak.concatenate(parts)
    W = ak.to_numpy(ev.LHEPdfWeight)
    R = W[:, lo:hi] / W[:, [0]]
    m = np.abs(R - 1).max(1)
    gw = ak.to_numpy(ev.genWeight)
    eff = np.abs((gw[m > a.outlier, None] * (R[m > a.outlier] - 1)).sum(0)) / abs(gw.sum()) * 100
    frac = [100 * np.mean(m > x) for x in (1, 10, 100)]
    cnt = [int(np.sum(m > x)) for x in (1e3, 1e4)]
    print("%-45s %8d %6.2f%% %6.2f%% %7.3f%% %6d %5d  %5.1f%% (member %d)"
          % (name[:45], len(m), *frac, *cnt, eff.max(), eff.argmax() + lo))
    for j in np.argsort(-m)[: a.top]:
        k = int(np.argmax(np.abs(R[j] - 1)))
        print(f"      {int(ev.run[j])}:{int(ev.luminosityBlock[j])}:{int(ev.event[j])}  member {k + lo}  ratio {R[j, k]:+.3g}")
