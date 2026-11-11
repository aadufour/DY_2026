#!/usr/bin/env python3
"""
scan_pdf_weight_tail.py

Scan NanoAOD samples for events with unphysical PDF weights: for each event, take the largest
|LHEPdfWeight[i] / LHEPdfWeight[0] - 1| over the Hessian members (default 1..100).
reco mll = two leading muons, no ID/selection (quick check, not the analysis selection).
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
p.add_argument("--bins", default="40,60,80,100,120,140,180,220,270,350,500,700,1000,1500,3000",
               help="mll bin edges used to say which bin the worst events fall into")
p.add_argument("--prefer", default="fnal", help="preferred redirector substring when several replicas are given")
a = p.parse_args()

lo, hi = map(int, a.members.split(":"))
edges = np.array([float(x) for x in a.bins.split(",")])


def mll(pt, eta, phi, mass):
    px, py = ak.sum(pt * np.cos(phi), 1), ak.sum(pt * np.sin(phi), 1)
    pz, E = ak.sum(pt * np.sinh(eta), 1), ak.sum(np.sqrt((pt * np.cosh(eta)) ** 2 + mass**2), 1)
    return np.sqrt(np.maximum(ak.to_numpy(E**2 - px**2 - py**2 - pz**2), 0))


def which_bin(x):
    if x < 0:
        return "no 2 muons"
    i = np.searchsorted(edges, x, side="right") - 1
    return f"{edges[i]:.0f}-{edges[i + 1]:.0f}" if 0 <= i < len(edges) - 1 else ("underflow" if x < edges[0] else "overflow")


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

cols = ["run", "luminosityBlock", "event", "genWeight", "LHEPdfWeight", "Muon_pt", "Muon_eta", "Muon_phi",
        "LHEPart_pt", "LHEPart_eta", "LHEPart_phi", "LHEPart_mass", "LHEPart_pdgId"]
print("%-45s %8s %7s %7s %8s %6s %5s  %s" % ("sample", "events", ">1", ">10", ">100", ">1e3", ">1e4", f"max over members of the shift of that member's total sample yield from events with |r-1|>{a.outlier:g}"))
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
    top = np.argsort(-m)[: a.top]
    sub = ev[top]
    lep = abs(sub.LHEPart_pdgId) == 13
    m_lhe = mll(sub.LHEPart_pt[lep], sub.LHEPart_eta[lep], sub.LHEPart_phi[lep], sub.LHEPart_mass[lep])
    mu = ak.argsort(sub.Muon_pt, ascending=False)[:, :2]
    two = ak.to_numpy(ak.num(sub.Muon_pt) >= 2)
    m_reco = np.where(two, mll(sub.Muon_pt[mu], sub.Muon_eta[mu], sub.Muon_phi[mu], 0.1057 * ak.ones_like(sub.Muon_pt[mu])), -1)
    print("      %-18s %6s %10s %8s %9s  %s" % ("run:lumi:event", "member", "ratio", "LHE mll", "reco mll", "reco mll bin"))
    for i, j in enumerate(top):
        k = int(np.argmax(np.abs(R[j] - 1)))
        tag = f"{int(ev.run[j])}:{int(ev.luminosityBlock[j])}:{int(ev.event[j])}"
        print("      %-18s %6d %+10.3g %8.0f %9.0f  %s" % (tag, k + lo, R[j, k], m_lhe[i], m_reco[i], which_bin(m_reco[i])))
