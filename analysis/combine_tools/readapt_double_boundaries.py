#!/usr/bin/env python3
"""
Readapt the per-operator boundaries used for 2D (double) EFT scans, based on
the actual likelihood grids that were produced with the current boundaries.

For every operator pair with a higgsCombine.<op1>_<op2>.individual.MultiDimFit.mH125.root
grid, this:
  1. Extracts the (k_op1, k_op2, deltaNLL) grid from the "limit" tree.
  2. Finds the 68% and 95% CL contours of the *joint* 2D region using the
     same chi2(2 dof) thresholds combine's own mkEFTScan.py uses:
     2*deltaNLL = 2.30 (68%), 5.99 (95%).
  3. If the 95% contour is closed and stays away from the scanned box edges,
     measures its actual reach and proposes a box of margin * reach.
     If the grid never reaches 95% CL (or the contour is clipped by the box
     edge), fits a local quadratic surface to extrapolate where 95% CL would
     be reached, and proposes margin * that estimate (flagged "extrapolated").
  4. Aggregates per operator by taking, on each side, the most demanding
     requirement across all 26 pairs that operator appears in.

Run this on LLR where the higgsCombine*.root files live (needs pyROOT).

Typical usage:

    # quick sanity check on a couple of pairs first
    python3 readapt_double_boundaries.py --metadata metadata_double.json \\
        --scan-dir . --pairs cbBRe_cbe,cHQ1_cHd

    # full run
    python3 readapt_double_boundaries.py --metadata metadata_double.json \\
        --scan-dir . --out-metadata metadata_double_new.json --report report.json
"""

import argparse
import itertools
import json
import os

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)

CL_LEVELS = (2.30, 5.99)  # 2*deltaNLL thresholds for a 2-parameter joint region, chi2(2 dof) 68%/95%


def find_scan_file(scan_dir, op1, op2):
    for a, b in [(op1, op2), (op2, op1)]:
        fp = os.path.join(scan_dir, f"higgsCombine.{a}_{b}.individual.MultiDimFit.mH125.root")
        if os.path.isfile(fp):
            return fp, a, b
    return None, op1, op2


def load_grid(filepath, op1, op2):
    f = ROOT.TFile.Open(filepath)
    if not f or f.IsZombie():
        raise IOError(f"cannot open {filepath}")
    t = f.Get("limit")
    if not t:
        f.Close()
        raise IOError(f"no 'limit' tree in {filepath}")

    xs, ys, ds = [], [], []
    for ev in t:
        d = ev.deltaNLL
        if not np.isfinite(d) or d > 1e4 or d < -1.0:
            continue
        xs.append(getattr(ev, f"k_{op1}"))
        ys.append(getattr(ev, f"k_{op2}"))
        ds.append(d)
    f.Close()

    if len(xs) < 10:
        raise ValueError(f"too few valid grid points ({len(xs)})")

    x = np.array(xs, dtype="d")
    y = np.array(ys, dtype="d")
    d = np.array(ds, dtype="d")
    z = 2.0 * (d - d.min())
    return x, y, z


def get_contours(x, y, z, levels=CL_LEVELS, npx=150, npy=150):
    g2 = ROOT.TGraph2D(len(x), x, y, z)
    g2.SetNpx(npx)
    g2.SetNpy(npy)
    hist = g2.GetHistogram()

    # Mirrors mkEFTScan.py: bins outside the Delaunay-covered region come back
    # as 0, which would look like a (false) minimum and create spurious
    # contour crossings at the edge of the covered area. Push them well above
    # both CL levels so they don't get drawn.
    for i in range(hist.GetSize()):
        if hist.GetBinContent(i + 1) == 0:
            hist.SetBinContent(i + 1, 100)

    hist.SetContour(len(levels), np.array(levels, dtype="d"))
    c = ROOT.TCanvas("c_tmp_readapt", "", 10, 10)
    hist.Draw("CONT Z LIST")
    ROOT.gPad.Update()
    conts = ROOT.gROOT.GetListOfSpecials().FindObject("contours")

    result = []
    for i in range(len(levels)):
        pieces = []
        if conts:
            level_list = conts.At(i)
            if level_list:
                for gr in level_list:
                    n = gr.GetN()
                    px = np.frombuffer(gr.GetX(), dtype="d", count=n).copy()
                    py = np.frombuffer(gr.GetY(), dtype="d", count=n).copy()
                    pieces.append((px, py))
        result.append(pieces)
    c.Close()

    xaxis, yaxis = hist.GetXaxis(), hist.GetYaxis()
    box = (xaxis.GetXmin(), xaxis.GetXmax(), yaxis.GetXmin(), yaxis.GetXmax())
    return result, box


def is_closed_interior(px, py, box, tol_frac=0.02):
    xlo, xhi, ylo, yhi = box
    tolx = tol_frac * (xhi - xlo)
    toly = tol_frac * (yhi - ylo)
    closed = abs(px[0] - px[-1]) < tolx and abs(py[0] - py[-1]) < toly
    touches_edge = (
        np.any(px <= xlo + tolx) or np.any(px >= xhi - tolx)
        or np.any(py <= ylo + toly) or np.any(py >= yhi - toly)
    )
    return closed and not touches_edge


def extrapolate_reach(x, y, z, x0, y0, level):
    """Fit a local quadratic surface around the best-fit point and use it to
    estimate how far out the `level` contour would sit, for pairs where the
    scanned grid never gets there. Approximate (symmetric about x0,y0) —
    meant to unstick a too-small box, not to be the final answer. Re-run this
    script after regenerating the scan with the enlarged range to refine."""
    dx, dy = x - x0, y - y0
    A = np.column_stack([dx**2, dx * dy, dy**2, dx, dy, np.ones_like(dx)])
    coeffs, *_ = np.linalg.lstsq(A, z, rcond=None)
    a, b, c = coeffs[0], coeffs[1], coeffs[2]
    H = np.array([[2 * a, b], [b, 2 * c]])
    try:
        Hinv = np.linalg.inv(H)
    except np.linalg.LinAlgError:
        return None
    if Hinv[0, 0] <= 0 or Hinv[1, 1] <= 0:
        return None  # not a well-behaved bowl (e.g. a flat/degenerate direction) - can't extrapolate safely
    reach_x = float(np.sqrt(level * Hinv[0, 0]))
    reach_y = float(np.sqrt(level * Hinv[1, 1]))
    return reach_x, reach_x, reach_y, reach_y


def analyze_pair(x, y, z, op1, op2, margin, fallback_factor):
    x0, y0 = float(x[np.argmin(z)]), float(y[np.argmin(z)])
    contours, box = get_contours(x, y, z)
    xlo, xhi, ylo, yhi = box

    c68, c95 = contours
    good68 = any(is_closed_interior(px, py, box) for px, py in c68)
    good95 = any(is_closed_interior(px, py, box) for px, py in c95)

    if good95:
        method = "measured"
        allx = np.concatenate([px for px, py in c95])
        ally = np.concatenate([py for px, py in c95])
        reach = (x0 - allx.min(), allx.max() - x0, y0 - ally.min(), ally.max() - y0)
    else:
        est = extrapolate_reach(x, y, z, x0, y0, CL_LEVELS[1])
        if est is None:
            method = "degenerate_fallback"
            reach = (
                fallback_factor * (x0 - xlo), fallback_factor * (xhi - x0),
                fallback_factor * (y0 - ylo), fallback_factor * (yhi - y0),
            )
        else:
            method = "extrapolated"
            reach = est

    reach_x_lo, reach_x_hi, reach_y_lo, reach_y_hi = reach
    new_x_lo = x0 - margin * reach_x_lo
    new_x_hi = x0 + margin * reach_x_hi
    new_y_lo = y0 - margin * reach_y_lo
    new_y_hi = y0 + margin * reach_y_hi

    return {
        "op1": op1, "op2": op2,
        "x0": x0, "y0": y0,
        "box": box,
        "good68": bool(good68), "good95": bool(good95),
        "method": method,
        "new_range_op1": (new_x_lo, new_x_hi),
        "new_range_op2": (new_y_lo, new_y_hi),
    }


def aggregate(results, ops):
    cand = {op: {"lo": [], "hi": [], "lo_src": [], "hi_src": []} for op in ops}
    for res in results:
        op1, op2 = res["op1"], res["op2"]
        lo1, hi1 = res["new_range_op1"]
        lo2, hi2 = res["new_range_op2"]
        cand[op1]["lo"].append(lo1); cand[op1]["lo_src"].append(op2)
        cand[op1]["hi"].append(hi1); cand[op1]["hi_src"].append(op2)
        cand[op2]["lo"].append(lo2); cand[op2]["lo_src"].append(op1)
        cand[op2]["hi"].append(hi2); cand[op2]["hi_src"].append(op1)

    new_bounds, drivers = {}, {}
    for op in ops:
        if not cand[op]["lo"]:
            continue
        lo_idx = int(np.argmin(cand[op]["lo"]))
        hi_idx = int(np.argmax(cand[op]["hi"]))
        new_bounds[op] = (cand[op]["lo"][lo_idx], cand[op]["hi"][hi_idx])
        drivers[op] = {"lo_driver": cand[op]["lo_src"][lo_idx], "hi_driver": cand[op]["hi_src"][hi_idx]}
    return new_bounds, drivers


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--metadata", default="metadata_double.json", help="Used only for the list of operator names")
    p.add_argument("--scan-dir", default=".", help="Directory containing higgsCombine.*.individual.MultiDimFit.mH125.root")
    p.add_argument("--margin", type=float, default=1.5, help="Box half-width = margin * measured/estimated 95%% CL reach")
    p.add_argument("--fallback-factor", type=float, default=10.0,
                    help="Multiplier applied to the current scanned half-width when even a quadratic extrapolation fails (degenerate direction / completely flat likelihood)")
    p.add_argument("--tolerance", type=float, default=0.10, help="Relative tolerance vs the current scanned box to call an operator 'good as is'")
    p.add_argument("--out-metadata", default="metadata_double_new.json")
    p.add_argument("--report", default="boundary_report.json")
    p.add_argument("--pairs", default="", help="Comma separated op1_op2 pairs to restrict to (for quick testing)")
    args = p.parse_args()

    with open(args.metadata) as fh:
        metadata = json.load(fh)
    ops = list(metadata["operators"].keys())

    pairs = list(itertools.combinations(ops, 2))
    if args.pairs:
        wanted = set(args.pairs.split(","))
        pairs = [(a, b) for a, b in pairs if f"{a}_{b}" in wanted or f"{b}_{a}" in wanted]

    results, missing = [], []
    for op1, op2 in pairs:
        fp, a, b = find_scan_file(args.scan_dir, op1, op2)
        if fp is None:
            missing.append(f"{op1}_{op2} (no scan file)")
            continue
        try:
            x, y, z = load_grid(fp, a, b)
            res = analyze_pair(x, y, z, a, b, args.margin, args.fallback_factor)
            results.append(res)
            print(f"{a}_{b}: 68%={'ok' if res['good68'] else 'MISSING':7s} "
                  f"95%={'ok' if res['good95'] else 'MISSING':7s} "
                  f"method={res['method']:18s} "
                  f"new[{a}]=({res['new_range_op1'][0]:+.4g},{res['new_range_op1'][1]:+.4g}) "
                  f"new[{b}]=({res['new_range_op2'][0]:+.4g},{res['new_range_op2'][1]:+.4g})")
        except Exception as e:
            print(f"[ERROR] {op1}_{op2}: {e}")
            missing.append(f"{op1}_{op2} (error: {e})")

    new_bounds, drivers = aggregate(results, ops)

    good_ops, redo_ops = [], []
    for op in ops:
        if op not in new_bounds:
            continue
        old_lo, old_hi = metadata["operators"][op]
        new_lo, new_hi = new_bounds[op]
        span = old_hi - old_lo
        if abs(new_lo - old_lo) <= args.tolerance * span and abs(new_hi - old_hi) <= args.tolerance * span:
            good_ops.append(op)
        else:
            redo_ops.append(op)

    print("\n=== GOOD AS IS ===")
    for op in good_ops:
        print(f"  {op}: {metadata['operators'][op]}")

    print("\n=== NEEDS REDO ===")
    for op in redo_ops:
        print(f"  {op}: {metadata['operators'][op]} -> "
              f"[{new_bounds[op][0]:.4g}, {new_bounds[op][1]:.4g}]  "
              f"(driven by lo:{drivers[op]['lo_driver']}, hi:{drivers[op]['hi_driver']})")

    if missing:
        print(f"\n=== MISSING / FAILED ({len(missing)}) ===")
        for m in missing:
            print(" ", m)

    new_metadata = json.loads(json.dumps(metadata))
    for op, (lo, hi) in new_bounds.items():
        new_metadata["operators"][op] = [round(lo, 6), round(hi, 6)]
    with open(args.out_metadata, "w") as fh:
        json.dump(new_metadata, fh, indent=4)
    print(f"\nWrote {args.out_metadata}")

    report = {
        "pairs": results,
        "missing": missing,
        "good_ops": good_ops,
        "redo_ops": redo_ops,
        "new_bounds": {op: list(v) for op, v in new_bounds.items()},
        "drivers": drivers,
    }
    with open(args.report, "w") as fh:
        json.dump(report, fh, indent=2, default=float)
    print(f"Wrote {args.report}")


if __name__ == "__main__":
    main()
