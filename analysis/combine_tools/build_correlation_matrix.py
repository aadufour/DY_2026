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

Usage:
    python3 build_correlation_matrix.py --metadata metadata.json --scan-dir . \\
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

    missing = []
    for op1, op2 in itertools.combinations(ops, 2):
        fp, a, b = find_scan_file(args.scan_dir, op1, op2)
        if fp is None:
            missing.append(f"{op1}_{op2} (no scan file)")
            continue
        try:
            x, y, z = load_grid(fp, a, b)
            rho = fit_correlation(x, y, z, zcap=args.zcap)
            if rho is None:
                missing.append(f"{op1}_{op2} (degenerate/flat fit)")
                continue
            i, j = idx[a], idx[b]
            corr[i, j] = rho
            corr[j, i] = rho
            print(f"{a}_{b}: rho = {rho:+.3f}")
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
