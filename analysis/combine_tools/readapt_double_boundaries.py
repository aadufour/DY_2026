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
     measures its actual reach (fill fraction = reach / current half-width).
     An operator is "good as is" only if this fill fraction is >= --min-fill
     (default 0.30) for every pair it appears in; otherwise a new box is
     proposed so the 95% CL fills --target-fill of it (default 0.5, i.e. box
     half-width = 2x the reach). If the grid never reaches 95% CL (or the
     contour is clipped by the box edge), a local quadratic surface is fit to
     extrapolate where 95% CL would be reached, and the same --target-fill
     sizing is applied to that estimate (flagged "extrapolated").
  4. Aggregates per operator by taking, on each side, the most demanding
     requirement across all 26 pairs that operator appears in. Operators
     that are already "good as is" keep their current boundary untouched in
     the output, so already-fine pairs don't drift run to run.

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

    # Grids produced with --doSplitPoints get hadd'd back together from N
    # sub-jobs, and each sub-job re-writes its own copy of the best-fit
    # reference point - so the true minimum ends up duplicated N times on
    # top of itself. Left in, that cluster of coincident points sits right
    # at the most important spot (the minimum) and degrades both the
    # Delaunay triangulation used for contour extraction and the local
    # quadratic fit used for extrapolation. Dedupe on (x,y), keeping the
    # lowest deltaNLL seen for each coordinate.
    best = {}
    for ev in t:
        d = ev.deltaNLL
        if not np.isfinite(d) or d > 1e4 or d < -1.0:
            continue
        x = getattr(ev, f"k_{op1}")
        y = getattr(ev, f"k_{op2}")
        key = (round(x, 8), round(y, 8))
        if key not in best or d < best[key][2]:
            best[key] = (x, y, d)
    f.Close()

    if len(best) < 10:
        raise ValueError(f"too few valid grid points ({len(best)})")

    xs = [v[0] for v in best.values()]
    ys = [v[1] for v in best.values()]
    ds = [v[2] for v in best.values()]

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


def extrapolate_reach(x, y, z, x0, y0, box, level, zcap=50.0):
    """Fit a local quadratic surface around the best-fit point and use it to
    estimate how far out the `level` contour would sit, for pairs where the
    scanned grid never gets there. Approximate (symmetric about x0,y0) —
    meant to unstick a too-small box, not to be the final answer. Re-run this
    script after regenerating the scan with the enlarged range to refine.

    x,y are rescaled by the current box half-widths before fitting: the two
    operators in a pair can have wildly different natural scales (e.g. a
    range of 0.8 vs 40), and fitting the raw coordinates leaves the
    least-squares design matrix badly conditioned, which is what was causing
    the fit to swing wildly between re-runs. Points with z above `zcap` are
    dropped too — those are usually saturated/non-converged grid points far
    from the minimum and can dominate an unweighted fit without reflecting
    the actual local curvature.
    """
    xlo, xhi, ylo, yhi = box
    sx = max(xhi - x0, x0 - xlo, 1e-12)
    sy = max(yhi - y0, y0 - ylo, 1e-12)

    mask = z <= zcap
    if mask.sum() < 6:
        mask = np.ones_like(z, dtype=bool)

    dx = (x[mask] - x0) / sx
    dy = (y[mask] - y0) / sy
    zz = z[mask]

    A = np.column_stack([dx**2, dx * dy, dy**2, dx, dy, np.ones_like(dx)])
    coeffs, *_ = np.linalg.lstsq(A, zz, rcond=None)
    a, b, c = coeffs[0], coeffs[1], coeffs[2]
    # For z = a*dx^2 + b*dx*dy + c*dy^2 = (dv)^T H (dv), H is [[a, b/2], [b/2, c]]
    # (off-diagonal is b/2, not b - matching dv^T H dv = a*dx^2 + 2*(b/2)*dx*dy + c*dy^2).
    H = np.array([[a, b / 2], [b / 2, c]])
    try:
        Hinv = np.linalg.inv(H)
    except np.linalg.LinAlgError:
        return None
    if Hinv[0, 0] <= 0 or Hinv[1, 1] <= 0:
        return None  # not a well-behaved bowl (e.g. a flat/degenerate direction) - can't extrapolate safely
    reach_x = float(np.sqrt(level * Hinv[0, 0])) * sx
    reach_y = float(np.sqrt(level * Hinv[1, 1])) * sy
    return reach_x, reach_x, reach_y, reach_y


def analyze_pair(x, y, z, op1, op2, margin, fallback_factor, zcap):
    x0, y0 = float(x[np.argmin(z)]), float(y[np.argmin(z)])
    contours, box = get_contours(x, y, z)
    xlo, xhi, ylo, yhi = box

    c68, c95 = contours
    good68 = any(is_closed_interior(px, py, box) for px, py in c68)
    good95 = any(is_closed_interior(px, py, box) for px, py in c95)

    fill_op1 = fill_op2 = (None, None)

    if good95:
        method = "measured"
        allx = np.concatenate([px for px, py in c95])
        ally = np.concatenate([py for px, py in c95])
        reach = (x0 - allx.min(), allx.max() - x0, y0 - ally.min(), ally.max() - y0)
        reach_x_lo, reach_x_hi, reach_y_lo, reach_y_hi = reach
        fill_op1 = (
            reach_x_lo / (x0 - xlo) if (x0 - xlo) > 0 else 1.0,
            reach_x_hi / (xhi - x0) if (xhi - x0) > 0 else 1.0,
        )
        fill_op2 = (
            reach_y_lo / (y0 - ylo) if (y0 - ylo) > 0 else 1.0,
            reach_y_hi / (yhi - y0) if (yhi - y0) > 0 else 1.0,
        )
    else:
        est = extrapolate_reach(x, y, z, x0, y0, box, CL_LEVELS[1], zcap=zcap)
        if est is None:
            method = "degenerate_fallback"
            reach = (
                fallback_factor * (x0 - xlo), fallback_factor * (xhi - x0),
                fallback_factor * (y0 - ylo), fallback_factor * (yhi - y0),
            )
        else:
            method = "extrapolated"
            reach = est

        # If the scan already went out to the current box edge and STILL
        # didn't find a closed 95% contour, the true reach can't be smaller
        # than that edge - clamp so a not-converged pair never proposes
        # shrinking. Without this, a quadratic fit that underestimates (e.g.
        # because the true curve steepens faster than parabolic further out)
        # can silently shrink a box that just needs to grow.
        reach = (
            max(reach[0], x0 - xlo), max(reach[1], xhi - x0),
            max(reach[2], y0 - ylo), max(reach[3], yhi - y0),
        )

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
        "fill_op1": fill_op1, "fill_op2": fill_op2,
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
    p.add_argument("--target-fill", type=float, default=0.5,
                    help="When a redo is needed, size the new box so the measured/estimated 95%% CL reach fills this fraction of the box (default 0.5, i.e. box half-width = 2x the reach)")
    p.add_argument("--min-fill", type=float, default=0.25,
                    help="An operator is 'good as is' only if, for every pair it appears in, the 95%% CL fills at least this fraction of the box (default 0.25, i.e. box at most 4x the reach). Below this the box is considered too big and gets resized toward --target-fill.")
    p.add_argument("--max-fill", type=float, default=0.90,
                    help="Upper bound on the fill fraction for 'good as is' (default 0.90, i.e. at least 10%% margin between the 95%% CL and the box edge). Above this the box is considered too tight and gets resized toward --target-fill even though the contour still technically closes.")
    p.add_argument("--fallback-factor", type=float, default=10.0,
                    help="Multiplier applied to the current scanned half-width when even a quadratic extrapolation fails (degenerate direction / completely flat likelihood)")
    p.add_argument("--fit-zcap", type=float, default=50.0,
                    help="Grid points with 2*deltaNLL above this are excluded from the local quadratic extrapolation fit (keeps saturated/non-converged points from skewing it)")
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
            margin = 1.0 / args.target_fill
            res = analyze_pair(x, y, z, a, b, margin, args.fallback_factor, args.fit_zcap)
            results.append(res)
            fill_str = ""
            if res["method"] == "measured":
                f1lo, f1hi = res["fill_op1"]
                f2lo, f2hi = res["fill_op2"]
                fill_str = f" fill[{a}]=({f1lo:.2f},{f1hi:.2f}) fill[{b}]=({f2lo:.2f},{f2hi:.2f})"
            print(f"{a}_{b}: 68%={'ok' if res['good68'] else 'MISSING':7s} "
                  f"95%={'ok' if res['good95'] else 'MISSING':7s} "
                  f"method={res['method']:18s} "
                  f"new[{a}]=({res['new_range_op1'][0]:+.4g},{res['new_range_op1'][1]:+.4g}) "
                  f"new[{b}]=({res['new_range_op2'][0]:+.4g},{res['new_range_op2'][1]:+.4g})"
                  f"{fill_str}")
        except Exception as e:
            print(f"[ERROR] {op1}_{op2}: {e}")
            missing.append(f"{op1}_{op2} (error: {e})")

    new_bounds, drivers = aggregate(results, ops)

    # An operator is only "good as is" if EVERY pair it appears in has a
    # closed 95% CL contour whose fill fraction (reach / box half-width) is
    # between --min-fill and --max-fill on both sides. Too low -> box is way
    # bigger than it needs to be; too high -> CL is too close to the edge for
    # comfort even though it technically still closes. A single unsatisfied
    # pair (not converged, or fill fraction out of band) marks the whole
    # operator for redo.
    satisfied = {op: True for op in ops}
    seen = {op: False for op in ops}
    for res in results:
        op1, op2 = res["op1"], res["op2"]
        seen[op1] = seen[op2] = True
        if not res["good95"]:
            satisfied[op1] = False
            satisfied[op2] = False
            continue
        f1lo, f1hi = res["fill_op1"]
        f2lo, f2hi = res["fill_op2"]
        if not (args.min_fill <= f1lo <= args.max_fill and args.min_fill <= f1hi <= args.max_fill):
            satisfied[op1] = False
        if not (args.min_fill <= f2lo <= args.max_fill and args.min_fill <= f2hi <= args.max_fill):
            satisfied[op2] = False

    good_ops = [op for op in ops if seen[op] and satisfied[op]]
    redo_ops = [op for op in ops if seen[op] and not satisfied[op]]

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

    # Only operators flagged for redo get their boundary rewritten - "good as
    # is" operators keep whatever is already in --metadata untouched, so
    # already-fine pairs don't drift/oscillate run to run.
    new_metadata = json.loads(json.dumps(metadata))
    for op in redo_ops:
        lo, hi = new_bounds[op]
        new_metadata["operators"][op] = [round(lo, 6), round(hi, 6)]
    with open(args.out_metadata, "w") as fh:
        json.dump(new_metadata, fh, indent=4)
    print(f"\nWrote {args.out_metadata}")

    report = {
        "pairs": results,
        "missing": missing,
        "good_ops": good_ops,
        "redo_ops": redo_ops,
        "new_bounds": {op: list(new_bounds[op]) for op in redo_ops},
        "drivers": {op: drivers[op] for op in redo_ops},
    }
    with open(args.report, "w") as fh:
        json.dump(report, fh, indent=2, default=float)
    print(f"Wrote {args.report}")


if __name__ == "__main__":
    main()
