"""
build_datacard.py

Build the combine input ROOT file and datacard.txt from the LHE cache.
Supports one or multiple EFT operators simultaneously.

Usage:
    python3 build_datacard.py --op cHDD
    python3 build_datacard.py --op cHDD --C 0.5 --lumi 59740
    python3 build_datacard.py --op cHDD --unrolled
    python3 build_datacard.py --op cHDD cHW --C 1.0 0.5 --lumi 59740
    python3 build_datacard.py --all_op --lumi 59740
    python3 build_datacard.py --op cHDD --output histograms.root --datacard datacard.txt --cache lhe_cache.pkl
"""

import argparse
import os
import pickle
from itertools import combinations

import boost_histogram as bh
import numpy as np
import uproot

# ---- Config --------------------------------------------------

CACHE_FILE    = "/Users/albertodufour/code/DY2026/analysis/lhe_cache.pkl"
OUTPUT_FILE   = "/Users/albertodufour/code/DY2026/analysis/histograms.root"
DATACARD_FILE = "/Users/albertodufour/code/DY2026/analysis/datacard.txt"
CHANNEL       = "triple_DY"

MLL_EDGES   = np.array([50, 70, 90, 110, 200, 800, 1400, 2000, 2400, 3000], dtype=float)
RAP_EDGES   = np.array([0.0, 0.5, 1.0, 2.5], dtype=float)
CSTAR_EDGES = np.array([-1.0, -0.5, 0.0, 0.5, 1.0], dtype=float)

# ---- Args ----------------------------------------------------

parser = argparse.ArgumentParser(description="Build combine ROOT file and datacard from LHE cache.")
op_group = parser.add_mutually_exclusive_group(required=True)
op_group.add_argument("--op",     nargs="+",
                      help="One or more operator names (must be in cache, e.g. --op cHDD cHW)")
op_group.add_argument("--all_op", action="store_true",
                      help="Use all operators found in cache")
parser.add_argument("--C",        type=float, nargs="+", default=[1.0],
                    help="Reference Wilson coefficient(s) for sm_lin_quad template. "
                         "One value is broadcast to all operators; "
                         "or pass one per operator (default: 1.0)")
parser.add_argument("--lumi",     type=float, default=1.0,
                    help="Luminosity in pb^-1 (default: 1.0 = raw weights in pb)")
parser.add_argument("--unrolled", action="store_true",
                    help="Use unrolled 3D (mll x rap x cstar) histogram instead of 1D mll")
parser.add_argument("--cache",    default=CACHE_FILE, help="Path to lhe_cache.pkl")
parser.add_argument("--output",   default=OUTPUT_FILE, help="Output ROOT file path")
parser.add_argument("--datacard", default=DATACARD_FILE, help="Output datacard path")
args = parser.parse_args()

LUMI = args.lumi

# ---- Load cache ----------------------------------------------

if not os.path.exists(args.cache):
    raise FileNotFoundError(
        f"Cache not found: {args.cache}\n"
        "Run build_cache.py first to generate it."
    )

print("Loading cache ...")
with open(args.cache, "rb") as f:
    cache = pickle.load(f)

mll_arr   = cache["mll"]
rap_arr   = cache["rap"]
cstar_arr = cache["cstar"]
w_SM      = cache["w_SM"]
w_p1_all  = cache["w_p1"]
w_m1_all  = cache["w_m1"]
w_pp_all  = cache.get("w_pp", {})
print(f"  {len(mll_arr):,} events loaded\n")

# ---- Resolve operator list -----------------------------------

if args.all_op:
    OPERATORS = sorted(w_p1_all.keys())
    print(f"--all_op: found {len(OPERATORS)} operators in cache: {OPERATORS}")
else:
    OPERATORS = args.op

# Validate all requested operators exist in cache
missing = [op for op in OPERATORS if op not in w_p1_all]
if missing:
    available = sorted(w_p1_all.keys())
    raise KeyError(
        f"Operator(s) not found in cache: {missing}\n"
        f"Available: {available}\n"
        f"Add them to build_cache.py's OPERATORS list and rebuild."
    )

# ---- Resolve C values ----------------------------------------

if len(args.C) == 1:
    C_values = {op: args.C[0] for op in OPERATORS}
elif len(args.C) == len(OPERATORS):
    C_values = {op: c for op, c in zip(OPERATORS, args.C)}
else:
    raise ValueError(
        f"--C expects either 1 value (broadcast) or {len(OPERATORS)} values "
        f"(one per operator), got {len(args.C)}."
    )

# ---- Print summary -------------------------------------------

print(f"Operators  : {OPERATORS}")
print(f"C_ref      : { {op: C_values[op] for op in OPERATORS} }")
print(f"Luminosity : {LUMI} pb^-1")
print(f"Observable : {'unrolled 3D' if args.unrolled else 'mll 1D'}")
print(f"Cache      : {args.cache}")
print(f"Output     : {args.output}")
print(f"Datacard   : {args.datacard}")
print()

# ---- Histogram builders --------------------------------------

def _x_axis(label):
    return bh.axis.Variable(MLL_EDGES, metadata=label)


def _unrolled_axis():
    n = (len(MLL_EDGES)-1) * (len(RAP_EDGES)-1) * (len(CSTAR_EDGES)-1)
    return bh.axis.Regular(n, 0, n, metadata="Unrolled bin (m_{ll} #times y_{ll} #times cos#theta*)")


def make_mll_1d(weights, proc_label):
    h = bh.Histogram(_x_axis("m_{ll} [GeV]"), storage=bh.storage.Weight())
    h.fill(mll_arr, weight=weights * LUMI)
    h.metadata = proc_label
    return h


def _make_3d_filled(weights):
    h = bh.Histogram(
        bh.axis.Variable(MLL_EDGES),
        bh.axis.Variable(RAP_EDGES),
        bh.axis.Variable(CSTAR_EDGES),
        storage=bh.storage.Weight(),
    )
    h.fill(mll_arr, rap_arr, cstar_arr, weight=weights * LUMI)
    return h


def unroll_3d(h3d, proc_label):
    v   = h3d.view()
    n   = v["value"].size
    h1d = bh.Histogram(_unrolled_axis(), storage=bh.storage.Weight())
    h1d.view()["value"]    = v["value"].T.flatten()
    h1d.view()["variance"] = v["variance"].T.flatten()
    h1d.metadata = proc_label
    return h1d


def make_hist(weights, proc_label):
    if args.unrolled:
        return unroll_3d(_make_3d_filled(weights), proc_label)
    return make_mll_1d(weights, proc_label)

# ---- Weight decomposition and histograms ---------------------

# sm and data_obs are shared across all operators
histograms = {}
histograms["sm"]       = make_hist(w_SM, f"SM -- {CHANNEL}")
histograms["data_obs"] = make_hist(w_SM, f"data_obs (= SM) -- {CHANNEL}")

for op in OPERATORS:
    C   = C_values[op]
    wp1 = w_p1_all[op]
    wm1 = w_m1_all[op]

    w_lin         = 0.5 * (wp1 - wm1)
    w_quad        = 0.5 * (wp1 + wm1) - w_SM
    w_slq         = w_SM + C * w_lin + C**2 * w_quad
    w_quad_scaled = C**2 * w_quad

    if w_slq.sum() < 0:
        print(f"WARNING: [{op}] full cross section is negative at C={C} -- "
              "C may be outside EFT validity.\n")

    histograms[f"quad_{op}"]        = make_hist(w_quad_scaled, f"quad {op}  C={C} -- {CHANNEL}")
    histograms[f"sm_lin_quad_{op}"] = make_hist(w_slq,         f"SM+lin+quad {op}  C={C} -- {CHANNEL}")

# ---- Interference (mixed) terms for operator pairs -----------

OP_PAIRS = list(combinations(OPERATORS, 2))
if len(OPERATORS) >= 2 and not w_pp_all:
    print("WARNING: cache has no 'w_pp' cross-term weights -- skipping interference shapes.\n"
          "         Rebuild the cache with the updated build_cache.py to include them.")

for op1, op2 in OP_PAIRS:
    pair = (op1, op2)
    if pair not in w_pp_all:
        pair = (op2, op1)
    if pair not in w_pp_all:
        print(f"WARNING: no cross-weight found for ({op1}, {op2}), skipping.")
        continue

    C1 = C_values[op1]
    C2 = C_values[op2]

    # w_inter = coefficient of C1*C2 per event
    # Derived from: w_pp = w_SM + A1 + A2 + B11 + B22 + B12
    #               w_p1[op] = w_SM + A + B  =>  w_inter = w_pp - w_p1[op1] - w_p1[op2] + w_SM
    w_inter        = w_pp_all[pair] - w_p1_all[op1] - w_p1_all[op2] + w_SM
    w_inter_scaled = C1 * C2 * w_inter

    name = f"sm_lin_quad_mixed_{op1}_{op2}"
    histograms[name] = make_hist(w_inter_scaled, f"mixed {op1}x{op2}  C1={C1} C2={C2} -- {CHANNEL}")

# ---- Print histogram summary ---------------------------------

print(f"{'Process':<32}  {'Integral':>14}  {'sqrt(SumW2)':>14}")
print("-" * 64)
for name, h in histograms.items():
    intg = h.values().sum()
    unc  = np.sqrt(h.variances().sum())
    print(f"  {name:<30}  {intg:>14.4e}  {unc:>14.4e}")

n_bins = histograms["sm"].values().size
print(f"\nHistogram bins : {n_bins}")
if args.unrolled:
    n_m = len(MLL_EDGES)   - 1
    n_r = len(RAP_EDGES)   - 1
    n_c = len(CSTAR_EDGES) - 1
    print(f"  = {n_m} mll x {n_r} rap x {n_c} cstar")
print()

# ---- Write ROOT file -----------------------------------------

print(f"Writing {args.output} ...")
with uproot.recreate(args.output) as rf:
    for name, h in histograms.items():
        key = f"{CHANNEL}/{name}"
        rf[key] = h
        print(f"  wrote {key}")

print(f"\nOutput : {os.path.abspath(args.output)}")

# ---- Build process list and rates ----------------------------

# Order: sm, then for each operator: quad_{op}, sm_lin_quad_{op},
#        then for each pair: sm_lin_quad_mixed_{op1}_{op2}
processes = ["sm"]
for op in OPERATORS:
    processes.append(f"quad_{op}")
    processes.append(f"sm_lin_quad_{op}")
for op1, op2 in OP_PAIRS:
    name = f"sm_lin_quad_mixed_{op1}_{op2}"
    if name in histograms:
        processes.append(name)

rates       = {p: histograms[p].values().sum() for p in processes}
observation = histograms["data_obs"].values().sum()

# Process indices (all EFT templates get negative indices):
#   sm                           -> -1
#   quad_{OPi}                   -> -(2i+2)
#   sm_lin_quad_{OPi}            -> -(2i+3)
#   sm_lin_quad_mixed_{OPi_OPj}  -> -(2*N+2 + k)  (k = 0-based pair index)
n_ops = len(OPERATORS)
proc_indices = {"sm": -1}
for i, op in enumerate(OPERATORS):
    proc_indices[f"quad_{op}"]        = -(2*i + 2)
    proc_indices[f"sm_lin_quad_{op}"] = -(2*i + 3)
for k, (op1, op2) in enumerate(OP_PAIRS):
    name = f"sm_lin_quad_mixed_{op1}_{op2}"
    if name in histograms:
        proc_indices[name] = -(2*n_ops + 2 + k)

# ---- Write datacard ------------------------------------------

root_abs  = os.path.abspath(args.output)
n_proc    = len(processes)

sep_long  = "-" * 130
sep_short = "-" * 45

lines = [
    sep_long,
    "imax 1",
    f"jmax {n_proc - 1}",
    "kmax 1",
    sep_long,
    f"shapes *         {CHANNEL}    {root_abs} $CHANNEL/$PROCESS $CHANNEL/$PROCESS_$SYSTEMATIC",
    f"shapes data_obs  {CHANNEL}    {root_abs} $CHANNEL/data_obs",
    sep_long,
    f"bin          {CHANNEL}",
    f"observation  {observation:.4g}",
    sep_long,
]

# bin / process / process index / rate rows
bin_row  = "bin         " + "".join(f"  {CHANNEL:<25}" for _ in processes)
proc_row = "process     " + "".join(f"  {p:<25}" for p in processes)
idx_row  = "process     " + "".join(f"  {proc_indices[p]:<25}" for p in processes)
rate_row = "rate        " + "".join(f"  {rates[p]:<25.4g}" for p in processes)

lines += [bin_row, proc_row, idx_row, rate_row, sep_short]

# systematics
lumi_row = "lumi  lnN   " + "".join(f"  {1.02:<25}" for _ in processes)
lines.append(lumi_row)
lines.append("")

datacard_text = "\n".join(lines)

print(f"\nWriting {args.datacard} ...")
with open(args.datacard, "w") as f:
    f.write(datacard_text)

print(f"Datacard : {os.path.abspath(args.datacard)}")
print()
print(datacard_text)
