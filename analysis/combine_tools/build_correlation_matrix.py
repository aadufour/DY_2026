#!/usr/bin/env python3
"""
Build an approximate operator-operator correlation matrix from the already-
scanned pairwise 2D likelihood grids.

This is the "pairwise, others fixed at SM" approximation: for each pair, a
local quadratic is fit to the 2*deltaNLL surface near the best-fit point,
giving a 2x2 covariance matrix V = A^-1 (A = [[a, b/2], [b/2, c]] for
z = a*dx^2 + b*dx*dy + c*dy^2), and the implied correlation coefficient is
read off as rho = V12 / sqrt(V11*V22). This is NOT the same as the
marginalized correlation from a full simultaneous 27-operator fit (which
profiles the other 25 operators instead of fixing them to 0) - see
notes/combine.md, "2D (Double) EFT Scans" chapter, for the distinction.

If --single-scan-dir is given, pairs where that quadratic fit fails
("degenerate/flat" - not enough resolvable 2D curvature) fall back to
tracing the likelihood's profiled valley (the y that minimizes z for each
scanned x, and vice versa) instead, which only needs local per-slice minima
rather than a fittable overall bowl. The resulting slope is combined with
each operator's own (independently robust) 1D-scan curvature to derive rho
= slope * sqrt(c_1d/a_1d). Every pair fit this way is clearly logged as
using the fallback, both on stdout and in the missing/fallback summary.

Usage:
    python3 build_correlation_matrix.py --metadata metadata.json --scan-dir . \\
        --single-scan-dir ../../datacards_single/inc_mm/mll \\
        --out-matrix correlation_matrix.csv --out-plot correlation_matrix.pdf
"""

import argparse
import csv
import itertools
import json
import os

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)


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

    # Same dedup as readapt_double_boundaries.py's load_grid(): --doSplitPoints
    # jobs each re-write their own copy of the best-fit reference point, so
    # hadd stacks duplicates on top of the true minimum.
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


def fit_correlation(x, y, z, zcap):
    """Fit a local quadratic near the minimum and return the implied
    correlation coefficient, or None if the fit isn't a well-behaved bowl.
    Only points with z <= zcap are used, keeping the fit inside the region
    where the quadratic (Gaussian) approximation is expected to hold."""
    x0, y0 = x[np.argmin(z)], y[np.argmin(z)]
    dx_full, dy_full = x - x0, y - y0

    sx = max(np.abs(dx_full).max(), 1e-12)
    sy = max(np.abs(dy_full).max(), 1e-12)

    mask = z <= zcap
    if mask.sum() < 6:
        mask = np.ones_like(z, dtype=bool)

    dx = dx_full[mask] / sx
    dy = dy_full[mask] / sy
    zz = z[mask]

    design = np.column_stack([dx**2, dx * dy, dy**2, dx, dy, np.ones_like(dx)])
    coeffs, *_ = np.linalg.lstsq(design, zz, rcond=None)
    a, b, c = coeffs[0], coeffs[1], coeffs[2]

    A = np.array([[a, b / 2], [b / 2, c]])
    try:
        V = np.linalg.inv(A)
    except np.linalg.LinAlgError:
        return None
    if V[0, 0] <= 0 or V[1, 1] <= 0:
        return None

    rho = V[0, 1] / np.sqrt(V[0, 0] * V[1, 1])
    return float(np.clip(rho, -1.0, 1.0))


def load_1d_grid(filepath, op):
    f = ROOT.TFile.Open(filepath)
    if not f or f.IsZombie():
        raise IOError(f"cannot open {filepath}")
    t = f.Get("limit")
    if not t:
        f.Close()
        raise IOError(f"no 'limit' tree in {filepath}")

    best = {}
    for ev in t:
        d = ev.deltaNLL
        if not np.isfinite(d) or d > 1e4 or d < -1.0:
            continue
        x = getattr(ev, f"k_{op}")
        key = round(x, 8)
        if key not in best or d < best[key][1]:
            best[key] = (x, d)
    f.Close()

    if len(best) < 6:
        raise ValueError(f"too few valid 1D grid points ({len(best)})")

    xs = np.array([v[0] for v in best.values()], dtype="d")
    ds = np.array([v[1] for v in best.values()], dtype="d")
    z = 2.0 * (ds - ds.min())
    return xs, z


def fit_1d_curvature(x, z, zcap):
    """Fit a local 1D parabola z = A*dx^2 near the minimum and return the
    curvature A in physical (non-normalized) units, or None if not a
    well-behaved bowl. This is the standalone single-operator scan, which
    doesn't suffer from the "not enough resolvable 2D curvature" problem
    that motivates the valley-slope fallback in the first place."""
    x0 = x[np.argmin(z)]
    dx_full = x - x0
    sx = max(np.abs(dx_full).max(), 1e-12)

    mask = z <= zcap
    if mask.sum() < 4:
        mask = np.ones_like(z, dtype=bool)

    dx = dx_full[mask] / sx
    zz = z[mask]

    design = np.column_stack([dx**2, dx, np.ones_like(dx)])
    coeffs, *_ = np.linalg.lstsq(design, zz, rcond=None)
    a_norm = coeffs[0]
    if a_norm <= 0:
        return None
    return float(a_norm / sx**2)  # back to physical (raw-coordinate) units


def profile_slope(u, v, zz):
    """For each distinct value of u in the grid, find the v that minimizes
    zz there (the profiled/valley trajectory), then fit a line v = m*u + q
    to those points. Needs far less resolvable curvature than fitting the
    full 2D bowl - it only needs the *location* of each slice's minimum,
    not the overall curvature magnitude."""
    ru = np.round(u, 8)
    uniq = np.unique(ru)
    pu, pv = [], []
    for uv in uniq:
        sel = ru == uv
        i = np.argmin(zz[sel])
        pu.append(uv)
        pv.append(v[sel][i])
    if len(pu) < 4:
        return None
    pu, pv = np.array(pu), np.array(pv)
    design = np.column_stack([pu, np.ones_like(pu)])
    coeffs, *_ = np.linalg.lstsq(design, pv, rcond=None)
    return float(coeffs[0])


def fit_correlation_valley_fallback(x, y, z, a_1d, c_1d, zcap):
    """Fallback for pairs where the full 2D quadratic fit fails (not enough
    resolvable curvature in some direction). Traces the likelihood's
    profiled valley instead - which only needs local per-slice minima, not
    a fittable overall bowl - to get a slope, then combines that slope with
    the operators' own (independently robust) 1D-scan curvatures to derive
    rho = slope * sqrt(c_1d/a_1d). Returns None if this also isn't usable."""
    if a_1d is None or c_1d is None or a_1d <= 0 or c_1d <= 0:
        return None

    mask = z <= zcap
    if mask.sum() < 8:
        mask = np.ones_like(z, dtype=bool)
    xs, ys, zs = x[mask], y[mask], z[mask]

    m_xy = profile_slope(xs, ys, zs)  # dy/dx from profiling y given x
    m_yx = profile_slope(ys, xs, zs)  # dx/dy from profiling x given y

    slopes = []
    if m_xy is not None:
        slopes.append(m_xy)
    if m_yx is not None and abs(m_yx) > 1e-9:
        slopes.append(1.0 / m_yx)
    if not slopes:
        return None
    m = float(np.mean(slopes))

    rho = m * np.sqrt(c_1d / a_1d)
    return float(np.clip(rho, -1.0, 1.0))


# Physical grouping from latex/table/propcorr_ops_table.tex + non_propcorr_ops_table.tex
# (the \midrule divisions in those tables): bosonic/Higgs-current, Higgs-quark
# current + dipole, Higgs-lepton current, four-lepton, four-fermion semileptonic.
PHYSICS_ORDER = [
    "cHDD", "cHWB",
    "cHj1", "cHj3", "cHQ1", "cHQ3", "cHu", "cHd", "cHbq", "cbWRe", "cbBRe",
    "cHl1", "cHl3", "cHe",
    "cll1",
    "clj1", "clj3", "cQl1", "cQl3", "ceu", "ced", "cbe", "cje", "cQe", "clu", "cld", "cbl",
]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--metadata", default="metadata.json", help="Used only for the list of operator names")
    p.add_argument("--scan-dir", default=".", help="Directory containing higgsCombine.*.individual.MultiDimFit.mH125.root")
    p.add_argument("--single-scan-dir", default=None,
                    help="Directory containing the 1D higgsCombine.<op>.individual.MultiDimFit.mH125.root scans "
                         "(e.g. .../datacards_single/inc_mm/mll). If given, pairs where the full 2D quadratic fit "
                         "fails ('degenerate/flat') fall back to a valley-slope + 1D-curvature estimate instead of "
                         "being dropped. If omitted, those pairs are just reported missing as before.")
    p.add_argument("--zcap", type=float, default=10.0, help="2*deltaNLL cap on points used for the local quadratic fit (default 10)")
    p.add_argument("--out-matrix", default="correlation_matrix.csv")
    p.add_argument("--out-plot", default="correlation_matrix.pdf")
    p.add_argument("--order", choices=["metadata", "correlation", "type"], default="type",
                    help="'type' (default) groups operators by physical type, matching "
                         "latex/table/{propcorr,non_propcorr}_ops_table.tex. 'correlation' sorts by mean "
                         "correlation with all other operators. 'metadata' keeps metadata.json's own order.")
    args = p.parse_args()

    with open(args.metadata) as fh:
        metadata = json.load(fh)
    ops = list(metadata["operators"].keys())
    n = len(ops)
    idx = {op: i for i, op in enumerate(ops)}

    corr = np.full((n, n), np.nan)
    np.fill_diagonal(corr, 1.0)

    curvature_cache = {}

    def get_1d_curvature(op):
        if op in curvature_cache:
            return curvature_cache[op]
        a_op = None
        if args.single_scan_dir:
            fp = os.path.join(args.single_scan_dir, f"higgsCombine.{op}.individual.MultiDimFit.mH125.root")
            if os.path.isfile(fp):
                try:
                    x1d, z1d = load_1d_grid(fp, op)
                    a_op = fit_1d_curvature(x1d, z1d, zcap=args.zcap)
                except Exception as e:
                    print(f"  [WARN] 1D curvature fit failed for {op}: {e}")
        curvature_cache[op] = a_op
        return a_op

    missing = []
    fallback_pairs = []
    for op1, op2 in itertools.combinations(ops, 2):
        fp, a, b = find_scan_file(args.scan_dir, op1, op2)
        if fp is None:
            missing.append(f"{op1}_{op2} (no scan file)")
            continue
        try:
            x, y, z = load_grid(fp, a, b)
            rho = fit_correlation(x, y, z, zcap=args.zcap)
            method = "quadratic"
            if rho is None and args.single_scan_dir:
                a_1d = get_1d_curvature(a)
                c_1d = get_1d_curvature(b)
                rho = fit_correlation_valley_fallback(x, y, z, a_1d, c_1d, zcap=args.zcap)
                method = "valley-slope fallback"
            if rho is None:
                missing.append(f"{op1}_{op2} (degenerate/flat fit)")
                continue
            i, j = idx[a], idx[b]
            corr[i, j] = rho
            corr[j, i] = rho
            if method == "quadratic":
                print(f"{a}_{b}: rho = {rho:+.3f}")
            else:
                print(f"{a}_{b}: rho = {rho:+.3f}  [FITTED VIA VALLEY-SLOPE FALLBACK - full 2D quadratic fit failed]")
                fallback_pairs.append(f"{a}_{b}")
        except Exception as e:
            missing.append(f"{op1}_{op2} (error: {e})")

    if args.order == "type":
        missing_from_grouping = [op for op in ops if op not in PHYSICS_ORDER]
        order_names = [op for op in PHYSICS_ORDER if op in ops] + missing_from_grouping
        if missing_from_grouping:
            print(f"\n[WARN] not in PHYSICS_ORDER, appended at the end: {missing_from_grouping}")
        order = [ops.index(op) for op in order_names]
        ops = order_names
        corr = corr[np.ix_(order, order)]
    elif args.order == "correlation":
        means = []
        for i in range(n):
            row = np.delete(corr[i], i)  # exclude self-correlation (always 1.0)
            means.append(np.nanmean(row) if np.any(np.isfinite(row)) else np.nan)
        means = np.array(means)
        # operators with no valid correlation at all sort to the end, not treated as "most anticorrelated"
        order = sorted(range(n), key=lambda i: (np.isnan(means[i]), means[i] if not np.isnan(means[i]) else 0.0))
        ops = [ops[i] for i in order]
        corr = corr[np.ix_(order, order)]
        print("\nOperator order (mean correlation with all other operators, most anticorrelated first):")
        for new_i, old_i in enumerate(order):
            label = f"{means[old_i]:+.3f}" if np.isfinite(means[old_i]) else "n/a"
            print(f"  {ops[new_i]:8s} mean rho = {label}")

    with open(args.out_matrix, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([""] + ops)
        for i, op in enumerate(ops):
            writer.writerow([op] + [f"{v:.4f}" if np.isfinite(v) else "" for v in corr[i]])
    print(f"\nWrote {args.out_matrix}")

    if fallback_pairs:
        print(f"\n{len(fallback_pairs)} pairs used the valley-slope fallback (full 2D quadratic fit failed):")
        for pp in fallback_pairs:
            print(" ", pp)

    if missing:
        total = n * (n - 1) // 2
        print(f"\n{len(missing)} / {total} pairs missing/failed:")
        for m in missing:
            print(" ", m)

    import matplotlib.pyplot as plt
    try:
        import mplhep as hep
        plt.style.use(hep.style.CMS)
    except ImportError:
        pass

    fig, ax = plt.subplots(figsize=(max(10, 0.45 * n), max(9, 0.45 * n)))
    masked = np.ma.masked_invalid(corr)
    im = ax.imshow(masked, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(n))
    ax.set_xticklabels(ops, rotation=90, fontsize=12)
    ax.set_yticks(range(n))
    ax.set_yticklabels(ops, fontsize=12)
    ax.set_xticks(np.arange(-0.5, n, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.5)
    ax.tick_params(which="minor", length=0)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("correlation coefficient")
    ax.set_title("Pairwise operator correlation (others fixed at SM)", fontsize=14)
    fig.tight_layout()
    fig.savefig(args.out_plot)
    fig.savefig(args.out_plot.rsplit(".", 1)[0] + ".png", dpi=150)
    print(f"Wrote {args.out_plot}")


if __name__ == "__main__":
    main()
