#!/usr/bin/env python3
"""
Give the mixed EFT templates (w11_*) the shape-systematic variations of the SM template, bin by bin:
    w11_<syst>Up/Down = w11 * sm_<syst>Up/Down / sm      (for every shape nuisance that acts on sm)
and switch those nuisances on for the w11 processes in the datacard (column '-' -> '1').

Why: in propcorr_v3_2d the w11 templates were produced nominal only, so when the nuisances move sm/w1/wm1 but
not w11 the morphing combination is no longer a consistent EFT prediction (spurious freedom ~ c_a*c_b).

Writes a NEW folder (datacard.txt, shapes.root, jsonComb.json, metadata.json); the input is not modified.
Usage: add_w11_syst_from_sm.py <input datacard dir> <output dir>
"""
import os, shutil, sys
import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
src, dst = sys.argv[1], sys.argv[2]
if os.path.exists(dst):
    sys.exit(f"{dst} exists, refusing to overwrite")
os.makedirs(dst)
for f in ("jsonComb.json", "metadata.json"):
    if os.path.exists(f"{src}/{f}"):
        shutil.copy(f"{src}/{f}", dst)

# ---- datacard: process columns, shape nuisances acting on sm
lines = open(f"{src}/datacard.txt").read().splitlines()
proc_line = next(i for i, l in enumerate(lines) if l.split()[:1] == ["process"])
procs = lines[proc_line].split()[1:]
sm_col = procs.index("sm")
w11_cols = [i for i, p in enumerate(procs) if p.startswith("w11_")]
syst_names = []
out = []
for l in lines:
    t = l.split()
    if len(t) == len(procs) + 2 and t[1] == "shape" and t[2 + sm_col] not in ("-", "0"):
        for i in w11_cols:
            t[2 + i] = "1"
        syst_names.append(t[0])
        l = " ".join(f"{x:<20s}" for x in t)
    out.append(l)
open(f"{dst}/datacard.txt", "w").write("\n".join(out) + "\n")
print(f"{len(w11_cols)} w11 processes, shape nuisances now applied to them: {syst_names}")

# ---- shapes: copy everything, add w11 variations
fin = ROOT.TFile.Open(f"{src}/shapes.root")
fout = ROOT.TFile.Open(f"{dst}/shapes.root", "RECREATE")
keys = [k.GetName() for k in fin.GetListOfKeys()]
for k in keys:
    h = fin.Get(k); fout.cd(); h.Write(k)
sm = fin.Get("histo_sm")
n = sm.GetNbinsX()
smv = np.array([sm.GetBinContent(i) for i in range(1, n + 1)])
nadded = 0
for p in (procs[i] for i in w11_cols):
    h11 = fin.Get(f"histo_{p}")
    for s in syst_names:
        for d in ("Up", "Down"):
            hs = fin.Get(f"histo_sm_{s}{d}")
            hn = h11.Clone(f"histo_{p}_{s}{d}")
            for i in range(1, n + 1):
                r = hs.GetBinContent(i) / smv[i - 1] if smv[i - 1] != 0 else 1.0
                hn.SetBinContent(i, h11.GetBinContent(i) * r)
                hn.SetBinError(i, h11.GetBinError(i) * r)
            fout.cd(); hn.Write(hn.GetName()); nadded += 1
fout.Close()
print(f"wrote {dst}/shapes.root: {len(keys)} input histograms + {nadded} new w11 variations")
