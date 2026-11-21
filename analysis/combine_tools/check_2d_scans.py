#!/usr/bin/env python3
"""
Quality check of the 2D grid scans of one run_2d_perpair.sh production (<VARDIR>/full).

Per grid point (final/ scan, i.e. pass3 if present else pass2):
  - error   : "Minimized function has error status" logged during the fit (FASTEXIT, NaN, ...), from the condor logs
  - fastexit: FASTEXIT logged (subset of error)
  - missing : point absent from the output (fit failed)
  - outlier : 2dNLL more than --excess above the median of its 8 neighbours (fit stuck above the true minimum)
"inside" = the point lies in the 95% region, judged on the neighbour median (2dNLL < 5.99), so a bad point
does not hide itself by its own wrong value.

Usage: check_2d_scans.py <VARDIR> [--excess 1] [--csv out.csv] [--list-out bad.txt]
  --list-out: one line "pair point_index flags" per bad point inside 95% (input of refit_2d_points.py)
"""
import argparse, collections, glob, os, re
import numpy as np
import ROOT
from scipy.ndimage import generic_filter

ROOT.gROOT.SetBatch(True)
p = argparse.ArgumentParser()
p.add_argument("vardir")
p.add_argument("--excess", type=float, default=1.0)
p.add_argument("--csv", default="")
p.add_argument("--list-out", default="")
a = p.parse_args()
full = os.path.join(a.vardir, "full")

# ---- log flags, keyed by (pair, round(x), round(y)); final scan of a pair comes from pass3 if it exists
final_src = {}
for f in glob.glob(f"{full}/final/higgsCombine.*.root"):
    final_src[os.path.basename(f).split(".")[1]] = os.path.basename(os.path.dirname(os.path.realpath(f)))
flags = collections.defaultdict(set)
flags_all = set()
index = {}  # (pair, x, y) -> grid point index
re_pt = re.compile(r"^Point (\d+)/\d+ .*?k_(\w+) = (\S+) k_(\w+) = (\S+)")
for ps in ("pass2", "pass3"):
    for log in glob.glob(f"{full}/{ps}/condor_scan_bundled/logs/*.out"):
        key = None
        for l in open(log, errors="ignore"):
            m = re_pt.match(l)
            if m:
                pair = f"{m.group(2)}_{m.group(4)}"
                key = (pair, round(float(m.group(3)), 6), round(float(m.group(5)), 6)) if final_src.get(pair) == ps else None
                if key:
                    flags_all.add(key); index[key] = int(m.group(1))
                continue
            if key is None:
                continue
            if "Minimized function has error status" in l:
                flags[key].add("error")
            elif "FASTEXIT" in l:
                flags[key].add("fastexit")


def nbmedian(v):
    c = np.delete(v, 4)
    c = c[~np.isnan(c)]
    return np.median(c) if len(c) else np.nan


rows = []
badlist = []
for f in sorted(glob.glob(f"{full}/final/higgsCombine.*.root")):
    pair = os.path.basename(f).split(".")[1]
    o1, o2 = pair.split("_")
    d = ROOT.RDataFrame("limit", f).Filter("quantileExpected > -0.5").AsNumpy([f"k_{o1}", f"k_{o2}", "deltaNLL"])
    x, y, z = d[f"k_{o1}"], d[f"k_{o2}"], 2 * d["deltaNLL"]
    z = z - z.min()
    # grid from the log (all points, also the failed ones); output points snapped to the nearest grid node
    lk = [k for k in flags_all if k[0] == pair]
    xs = np.unique([k[1] for k in lk]); ys = np.unique([k[2] for k in lk])
    if len(xs) < 2:
        xs, ys = np.unique(np.round(x, 5)), np.unique(np.round(y, 5))
    ix = np.abs(x[:, None] - xs[None, :]).argmin(1); iy = np.abs(y[:, None] - ys[None, :]).argmin(1)
    Z = np.full((len(xs), len(ys)), np.nan)
    Z[ix, iy] = z
    M = generic_filter(Z, nbmedian, size=3, mode="constant", cval=np.nan)
    n_expected = 10000
    inside = M < 5.99
    outl = (Z - M > a.excess) & inside
    err = np.zeros_like(inside); fe = np.zeros_like(inside)
    for (pp, kx, ky), fl in flags.items():
        if pp != pair:
            continue
        i, j = np.abs(xs - kx).argmin(), np.abs(ys - ky).argmin()
        err[i, j] = "error" in fl; fe[i, j] = "fastexit" in fl
    missing = np.isnan(Z)
    if a.list_out:
        node = {(np.abs(xs - k[1]).argmin(), np.abs(ys - k[2]).argmin()): n for k, n in index.items() if k[0] == pair}
        for i, j in zip(*np.nonzero((missing | err | outl) & inside)):
            fl = [n for n, b in (("missing", missing), ("error", err), ("outlier", outl)) if b[i, j]]
            if (i, j) in node:
                badlist.append(f"{pair} {node[(i, j)]} {','.join(fl)}")
    rows.append(dict(pair=pair, n_missing=int(missing.sum()), missing_in=int((missing & inside).sum()),
                     n_err=int(err.sum()), err_in=int((err & inside).sum()), n_fe=int(fe.sum()), fe_in=int((fe & inside).sum()),
                     outl_in=int(outl.sum()), outl_max=float(np.nanmax(np.where(outl, Z - M, np.nan))) if outl.any() else 0.0,
                     zmax_err=float(np.nanmin(np.where(err, M, np.nan))) if err.any() else np.nan))

var = os.path.basename(os.path.normpath(a.vardir))
tot = lambda k: sum(r[k] for r in rows)
npairs = lambda k: sum(1 for r in rows if r[k])
print(f"== {var}: {len(rows)} pairs")
print(f"{'':28s}{'points':>8s}{'pairs':>7s}")
for lab, k in [("missing (failed fit)", "n_missing"), ("  inside 95%", "missing_in"),
               ("error status", "n_err"), ("  inside 95%", "err_in"), ("FASTEXIT", "n_fe"), ("  inside 95%", "fe_in"),
               (f"outlier >{a.excess:g} inside 95%", "outl_in")]:
    print(f"{lab:28s}{tot(k):8d}{npairs(k):7d}")
lo = [r["zmax_err"] for r in rows if not np.isnan(r["zmax_err"])]
if lo:
    print(f"lowest 2dNLL (neighbour median) at an error point: {min(lo):.1f}")
bad = sorted([r for r in rows if r["outl_in"] or r["err_in"] or r["missing_in"]], key=lambda r: -(r["outl_in"] + r["err_in"] + r["missing_in"]))
if bad:
    print("pairs with problems inside 95%  (outliers / errors / missing, max outlier excess):")
    for r in bad:
        print(f"   {r['pair']:18s} {r['outl_in']:3d} {r['err_in']:3d} {r['missing_in']:3d}   {r['outl_max']:.1f}")
if a.list_out:
    open(a.list_out, "w").write("\n".join(badlist) + "\n")
    print(f"{len(badlist)} bad points inside 95% -> {a.list_out}")
if a.csv:
    import csv
    with open(a.csv, "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
