#!/usr/bin/env python3
"""
Positivity check of the morphing templates in a shapes.root (flat histo_* keys).

Must be >= 0 in every bin, for nominal and every Up/Down variation:
  full  : histo_sm, histo_w1_<op> (c=+1), histo_wm1_<op> (c=-1), histo_w11_<a>_<b> (ca=cb=1), backgrounds
  quad  : (w1 + wm1)/2 - sm     (|M_eft|^2, exactly >= 0 event by event)
lin = (w1 - wm1)/2 is NOT checked (it can be negative).

A bin counts as negative if value < -tol_rel * max(|template|) (default 1e-6, float rounding).

Usage: check_positivity_morphing.py shapes.root [more.root ...] [--tol-rel 1e-6] [--verbose]
"""
import argparse, re
import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
p = argparse.ArgumentParser()
p.add_argument("files", nargs="+")
p.add_argument("--tol-rel", type=float, default=1e-6)
p.add_argument("--verbose", action="store_true")
a = p.parse_args()


def arr(f, n):
    h = f.Get(n)
    return np.array([h.GetBinContent(i) for i in range(1, h.GetNbinsX() + 1)])


def neg(v):
    tol = a.tol_rel * max(np.max(np.abs(v)), 1e-300)
    return np.nonzero(v < -tol)[0]


total_bad = 0
for fn in a.files:
    f = ROOT.TFile.Open(fn)
    names = [k.GetName() for k in f.GetListOfKeys() if k.GetClassName().startswith("TH1")]
    nameset = set(names)
    bad, nchecked = [], 0
    # every template (nominal + variations), except data
    for n in names:
        if n == "histo_Data":
            continue
        v = arr(f, n); nchecked += 1
        idx = neg(v)
        if len(idx):
            bad.append(("full", n, idx, v[idx]))
    # quad for every operator and variation suffix
    for n in names:
        m = re.match(r"^histo_w1_([A-Za-z0-9]+?)((?:_\w+(?:Up|Down))?)$", n)
        if not m:
            continue
        op, suf = m.group(1), m.group(2)
        wm, sm = f"histo_wm1_{op}{suf}", f"histo_sm{suf}"
        if wm not in nameset or sm not in nameset:
            continue
        q = 0.5 * (arr(f, n) + arr(f, wm)) - arr(f, sm); nchecked += 1
        idx = neg(q)
        if len(idx):
            bad.append(("quad", f"quad_{op}{suf}", idx, q[idx]))
    print(f"== {fn}: {nchecked} templates checked, {len(bad)} with negative bins")
    for kind, n, idx, vals in bad[: (None if a.verbose else 40)]:
        print(f"   {kind:4s} {n:45s} bins {list(idx + 1)} values {np.round(vals, 4).tolist()}")
    if len(bad) > 40 and not a.verbose:
        print(f"   ... {len(bad) - 40} more (--verbose)")
    total_bad += len(bad)
print("ALL POSITIVE" if total_bad == 0 else f"NEGATIVE BINS FOUND in {total_bad} templates")
