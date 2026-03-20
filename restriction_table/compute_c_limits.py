"""
compute_c_limits.py

For each SMEFT operator, scans the Wilson coefficient C and computes
68% / 95% CL intervals using the Poisson log-likelihood ratio.

Setup:
    data_obs  = sm  (by construction)
    expected(C, bin) = sm + C * lin + C² * quad
    lin = sm_lin_quad - sm - quad   (pure linear weight at C=1)

Test statistic:
    q(C) = -2 * sum_bins [ mu(C) - data * ln(mu(C)) ]
           + 2 * sum_bins [ mu(0) - data * ln(mu(0)) ]   <- subtract best-fit (C=0)

95% CL interval: q(C) < 3.841
68% CL interval: q(C) < 1.000

Usage:
    python3 compute_c_limits.py
    python3 compute_c_limits.py --root /path/to/histograms.root
    python3 compute_c_limits.py --c-range 5 --points 2000
"""

import argparse
import numpy as np
import uproot

OPERATORS = [
    'cHDD', 'cHWB', 'cHj1',
    'cHj3', 'cHu',  'cHd',
    'cHl1', 'cHl3', 'cHe',
    'cll1', 'clj1', 'clj3',
    'ceu',  'ced',  'cje',
    'clu',  'cld',
]
CHANNEL = "triple_DY"

parser = argparse.ArgumentParser()
parser.add_argument("--root",    default="/Users/albertodufour/code/DY2026/restriction_table/histograms.root")
parser.add_argument("--c-range", type=float, default=5.0,   help="Scan C in [-c_range, c_range]")
parser.add_argument("--points",  type=int,   default=5000,  help="Number of scan points")
parser.add_argument("--out",     default="/Users/albertodufour/code/DY2026/restriction_table/c_limits.txt")
args = parser.parse_args()

C_vals = np.linspace(-args.c_range, args.c_range, args.points)

def poisson_nll(mu, data):
    """
    -2 * sum[ mu - data*ln(mu) ]  for bins where mu > 0.
    Bins with mu <= 0 are skipped (they arise only for large |C| and are
    a sign that the EFT expansion has broken down).
    """
    safe = mu > 0
    nll = np.where(safe, mu - data * np.log(np.where(safe, mu, 1.0)), 0.0)
    return 2.0 * nll.sum()


print(f"\nOpening {args.root}")
with uproot.open(args.root) as f:
    sm_vals = f[f"{CHANNEL}/sm"].values().copy()
    data    = sm_vals   # data_obs = sm by construction

    results = {}
    print(f"\n{'Operator':<10}  {'68% CL low':>12}  {'68% CL high':>12}  {'95% CL low':>12}  {'95% CL high':>12}  {'q(C=1)':>8}")
    print("-" * 75)

    for op in OPERATORS:
        quad_vals = f[f"{CHANNEL}/quad_{op}"].values().copy()
        slq_vals  = f[f"{CHANNEL}/sm_lin_quad_{op}"].values().copy()
        lin_vals  = slq_vals - sm_vals - quad_vals   # pure linear component

        nll0 = poisson_nll(sm_vals, data)   # best-fit q (C=0, data=sm → nll0 ~ 0)

        q = np.array([
            poisson_nll(sm_vals + c * lin_vals + c**2 * quad_vals, data) - nll0
            for c in C_vals
        ])

        # 68% and 95% CL intervals (q < threshold)
        def interval(threshold):
            inside = C_vals[q < threshold]
            if len(inside) == 0:
                return (np.nan, np.nan)
            return (inside[0], inside[-1])

        lo68, hi68 = interval(1.000)
        lo95, hi95 = interval(3.841)

        # q at C=1 (how disfavoured is C=1)
        q_at_1 = float(np.interp(1.0, C_vals, q))

        results[op] = dict(lo68=lo68, hi68=hi68, lo95=lo95, hi95=hi95, q_at_1=q_at_1)
        print(f"  {op:<10}  {lo68:>12.4f}  {hi68:>12.4f}  {lo95:>12.4f}  {hi95:>12.4f}  {q_at_1:>8.2f}")

# ---- Write output table -----------------------------------------------
with open(args.out, "w") as fout:
    fout.write(f"# C limits from {args.root}\n")
    fout.write(f"# scan: C in [{-args.c_range}, {args.c_range}], {args.points} points\n")
    fout.write(f"# data_obs = sm  (Asimov at C=0)\n#\n")
    fout.write(f"# {'operator':<10}  {'68lo':>10}  {'68hi':>10}  {'95lo':>10}  {'95hi':>10}  {'q(C=1)':>8}\n")
    for op, r in results.items():
        fout.write(f"  {op:<10}  {r['lo68']:>10.4f}  {r['hi68']:>10.4f}  "
                   f"  {r['lo95']:>10.4f}  {r['hi95']:>10.4f}  {r['q_at_1']:>8.2f}\n")

print(f"\nResults written to {args.out}")
