#!/usr/bin/env python3
"""
build_datacard_new.py
=====================
Simplified datacard builder for the LHE pipeline.

Differences vs build_datacard_reco_bins.py:
  - Reads from lhe_cache_new.pkl (produced by build_cache_new.py)
  - No --pdf-flavour option: hardcoded to 5-flavour NNPDF31 (pdf_325300)
  - MUF-only scale envelope (MUR=1 fixed — no effect at LO)
  - RECO mll binning (34 bins, 50-3000 GeV)
  - N_gen normalization fix applied

Usage:
    python3 build_datacard_new.py --all_op --lumi 59740
    python3 build_datacard_new.py --op cHDD cHWB --lumi 59740
"""

import argparse
import os
import pickle
from itertools import combinations

import boost_histogram as bh
import numpy as np
import uproot

# ---- Config ------------------------------------------------------------------

CACHE_FILE    = "/grid_mnt/data__data.polcms/cms/adufour/LHE/SYST_slc7/CACHE/lhe_cache_parallel.pkl"
OUTPUT_FILE   = "histograms.root"
DATACARD_FILE = "datacard.txt"
CHANNEL       = "triple_DY"
PDF_LABEL     = "NNPDF31_nnlo_as_0118_mc_hessian_pdfas (325300, 5-flavour, 103 members)"
PDF_KEY       = "pdf_325300"

# Broader binning -- fine RECO binning is too statistics-limited for this sample.
# Based on the triple_diff MLL_EDGES convention in analysis/spritz/plot_eft.py,
# starting at 50 GeV (cache's actual lower bound) with extra high-mass bins
# extending to the full 3000 GeV range.
MLL_EDGES = np.array(
    [50, 60, 80, 100, 120, 140, 180, 220, 270, 350, 500, 750, 1000, 3000], dtype=float
)
RAP_EDGES   = np.array([0.0, 0.5, 1.0, 2.5], dtype=float)
CSTAR_EDGES = np.array([-1.0, -0.5, 0.0, 0.5, 1.0], dtype=float)

# ---- Args --------------------------------------------------------------------

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
op_group = parser.add_mutually_exclusive_group(required=True)
op_group.add_argument("--op",     nargs="+")
op_group.add_argument("--all_op", action="store_true")
parser.add_argument("--C",        type=float, nargs="+", default=[1.0])
parser.add_argument("--lumi",     type=float, default=59740.,
                    help="Luminosity in pb^-1 (default: 59740 = 2018 full)")
parser.add_argument("--unrolled", action="store_true")
parser.add_argument("--cache",    default=CACHE_FILE)
parser.add_argument("--output",   default=OUTPUT_FILE)
parser.add_argument("--datacard", default=DATACARD_FILE)
args = parser.parse_args()

LUMI = args.lumi

# ---- Load cache --------------------------------------------------------------

if not os.path.exists(args.cache):
    raise FileNotFoundError(f"Cache not found: {args.cache}\nRun build_cache_new.py first.")

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

w_scale_all   = cache.get("w_scale",       {})
w_pdf_central = cache.get("w_pdf_central", None)
pdf_arr       = cache.get(PDF_KEY,         None)

# ---- Normalization fix -------------------------------------------------------
# The cache merges 7 mll-binned samples each with N_GEN_PER_SAMPLE events.
# Each sample's weights sum to σ_sample × N_GEN_PER_SAMPLE (not ÷ N_gen).
# Dividing by N_GEN_PER_SAMPLE gives sum(w)/N_GEN_PER_SAMPLE × LUMI = σ × LUMI.
# DO NOT divide by len(w_SM) (= 7 × N_GEN_PER_SAMPLE) — that underestimates by 7×.

N_GEN_PER_SAMPLE = 100_000   # events generated per mll-binned gridpack
N_gen = len(w_SM)
print(f"  N_gen total                    : {N_gen:,}")
print(f"  N_gen per mll-bin sample       : {N_GEN_PER_SAMPLE:,}")
print(f"  sum(w_SM) before fix           : {w_SM.sum():.4e} pb")
w_SM          = w_SM / N_GEN_PER_SAMPLE
w_p1_all      = {op: w / N_GEN_PER_SAMPLE for op, w in w_p1_all.items()}
w_m1_all      = {op: w / N_GEN_PER_SAMPLE for op, w in w_m1_all.items()}
w_pp_all      = {k:  w / N_GEN_PER_SAMPLE for k,  w in w_pp_all.items()}
w_scale_all   = {k:  w / N_GEN_PER_SAMPLE for k,  w in w_scale_all.items()}
if w_pdf_central is not None:
    w_pdf_central = w_pdf_central / N_GEN_PER_SAMPLE
if pdf_arr is not None:
    pdf_arr = pdf_arr / N_GEN_PER_SAMPLE
print(f"  sum(w_SM) after fix            : {w_SM.sum():.4e} pb")
print(f"  N_expected at {LUMI:.0f} pb^-1       : {w_SM.sum() * LUMI:.4e}")
print(f"  (cross-check: sum(xwgt)/N_GEN_PER_SAMPLE*L = {cache.get('xwgt', w_SM*N_GEN_PER_SAMPLE).sum()/N_GEN_PER_SAMPLE*LUMI:.4e})")
print()

has_scale = bool(w_scale_all)
has_pdf   = pdf_arr is not None and w_pdf_central is not None

print(f"  {len(mll_arr):,} events loaded")
print(f"  MUF scale variation keys : {len(w_scale_all)}"
      + ("" if has_scale else "  (none — skipped)"))
if has_pdf:
    print(f"  PDF set : {PDF_LABEL} ({pdf_arr.shape[1]} members)")
else:
    print(f"  PDF ({PDF_KEY}) : not found in cache — skipped")
print()

# ---- Operators ---------------------------------------------------------------

if args.all_op:
    OPERATORS = sorted(w_p1_all.keys())
    print(f"--all_op: {len(OPERATORS)} operators found in cache")
else:
    OPERATORS = args.op

missing = [op for op in OPERATORS if op not in w_p1_all]
if missing:
    raise KeyError(f"Operators not in cache: {missing}\nAvailable: {sorted(w_p1_all.keys())}")

if len(args.C) == 1:
    C_values = {op: args.C[0] for op in OPERATORS}
elif len(args.C) == len(OPERATORS):
    C_values = {op: c for op, c in zip(OPERATORS, args.C)}
else:
    raise ValueError(f"--C: expected 1 or {len(OPERATORS)} values, got {len(args.C)}")

print(f"Operators  : {OPERATORS}")
print(f"Luminosity : {LUMI} pb^-1")
print(f"Observable : {'unrolled 3D' if args.unrolled else '1D mll'}")
syst_list = ["lumi"] + (["muf_scale"] if has_scale else []) + (["pdf"] if has_pdf else [])
print(f"Systematics: {', '.join(syst_list)}")
print()

# ---- Histogram helpers -------------------------------------------------------

def make_hist(weights, label=""):
    if args.unrolled:
        h3 = bh.Histogram(
            bh.axis.Variable(MLL_EDGES),
            bh.axis.Variable(RAP_EDGES),
            bh.axis.Variable(CSTAR_EDGES),
            storage=bh.storage.Weight())
        h3.fill(mll_arr, rap_arr, cstar_arr, weight=weights * LUMI)
        v = h3.view()
        n = (len(MLL_EDGES)-1) * (len(RAP_EDGES)-1) * (len(CSTAR_EDGES)-1)
        ax = bh.axis.Regular(n, 0, n, metadata="Unrolled (mll x yZ x cstar)")
        h1 = bh.Histogram(ax, storage=bh.storage.Weight())
        h1.view()["value"]    = v["value"].T.flatten()
        h1.view()["variance"] = v["variance"].T.flatten()
        h1.metadata = label
        return h1
    else:
        h = bh.Histogram(
            bh.axis.Variable(MLL_EDGES, metadata="m_{ll} [GeV]"),
            storage=bh.storage.Weight())
        h.fill(mll_arr, weight=weights * LUMI)
        h.metadata = label
        return h


def _pdf_updown(w_nom, h_nom):
    """Asymmetric PDF uncertainty: quadrature sum over eigenvector deviations."""
    pdf_central_set = pdf_arr[:, 0]
    nominal  = h_nom.values().flatten()
    rep_vals = np.array([
        make_hist(w_nom * (pdf_arr[:, i] / pdf_central_set)).values().flatten()
        for i in range(1, pdf_arr.shape[1])
    ])
    sigma_up   = np.zeros_like(nominal)
    sigma_down = np.zeros_like(nominal)
    for b in range(len(nominal)):
        dev    = rep_vals[:, b] - nominal[b]
        up_dev = dev[dev > 0]
        dn_dev = dev[dev < 0]
        if len(up_dev) > 0: sigma_up[b]   = np.sqrt(np.sum(up_dev**2))
        if len(dn_dev) > 0: sigma_down[b] = np.sqrt(np.sum(dn_dev**2))
    shape = h_nom.values().shape
    h_up = h_nom.copy(); h_up.view()["value"] = (nominal + sigma_up).reshape(shape)
    h_dn = h_nom.copy(); h_dn.view()["value"] = np.maximum(nominal - sigma_down, 0).reshape(shape)
    return h_up, h_dn

# ---- Build nominal histograms -------------------------------------------------

histograms   = {}
proc_weights = {}
histograms["sm"]       = make_hist(w_SM, "SM")
histograms["data_obs"] = make_hist(w_SM, "data_obs (Asimov = SM)")
proc_weights["sm"]     = w_SM

OP_PAIRS = list(combinations(OPERATORS, 2))

for op in OPERATORS:
    C      = C_values[op]
    wp1    = w_p1_all[op]
    wm1    = w_m1_all[op]
    w_lin  = 0.5 * (wp1 - wm1)
    w_quad = 0.5 * (wp1 + wm1) - w_SM
    w_slq  = w_SM + C * w_lin + C**2 * w_quad
    if w_slq.sum() < 0:
        print(f"WARNING: [{op}] total cross section negative at C={C}")
    histograms[f"quad_{op}"]        = make_hist(C**2 * w_quad, f"quad {op}")
    histograms[f"sm_lin_quad_{op}"] = make_hist(w_slq,         f"SM+lin+quad {op}")
    proc_weights[f"quad_{op}"]        = C**2 * w_quad
    proc_weights[f"sm_lin_quad_{op}"] = w_slq

for op1, op2 in OP_PAIRS:
    pair = (op1, op2) if (op1, op2) in w_pp_all else (op2, op1)
    if pair not in w_pp_all:
        print(f"WARNING: no cross-weight for ({op1},{op2}), skipping.")
        continue
    C1, C2  = C_values[op1], C_values[op2]
    w_inter = w_pp_all[pair] - w_p1_all[op1] - w_p1_all[op2] + w_SM
    w_mixed = C1 * C2 * w_inter
    histograms[f"sm_lin_quad_mixed_{op1}_{op2}"]  = make_hist(w_mixed)
    proc_weights[f"sm_lin_quad_mixed_{op1}_{op2}"] = w_mixed

# ---- Build systematic shape variants -----------------------------------------

scale_keys = list(w_scale_all.keys())   # MUF-only variations

nominal_procs = [k for k in histograms if k != "data_obs"]
for proc in nominal_procs:
    h     = histograms[proc]
    w_nom = proc_weights[proc]
    if has_scale:
        scale_hists = [make_hist(w_nom * (w_scale_all[k] / w_pdf_central)) for k in scale_keys]
        all_vals    = np.array([s.values() for s in scale_hists])
        h_up = scale_hists[0].copy(); h_up.view()["value"] = all_vals.max(axis=0)
        h_dn = scale_hists[0].copy(); h_dn.view()["value"] = all_vals.min(axis=0)
        histograms[f"{proc}_muf_scaleUp"]   = h_up
        histograms[f"{proc}_muf_scaleDown"] = h_dn
    if has_pdf:
        histograms[f"{proc}_pdfUp"], histograms[f"{proc}_pdfDown"] = _pdf_updown(w_nom, h)

# ---- Print summary -----------------------------------------------------------

print(f"{'Process':<42}  {'Integral':>14}  {'sqrt(SumW2)':>14}")
print("-" * 74)
for name in nominal_procs + ["data_obs"]:
    h    = histograms[name]
    intg = h.values().sum()
    unc  = np.sqrt(h.variances().sum())
    print(f"  {name:<40}  {intg:>14.4e}  {unc:>14.4e}")
print()

# ---- Write ROOT file ---------------------------------------------------------

print(f"Writing {args.output} ...")
with uproot.recreate(args.output) as rf:
    for name, h in histograms.items():
        rf[f"{CHANNEL}/{name}"] = h
print(f"  {len(histograms)} histograms written\n")

# ---- Datacard ----------------------------------------------------------------

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

n_ops        = len(OPERATORS)
proc_indices = {"sm": -1}
for i, op in enumerate(OPERATORS):
    proc_indices[f"quad_{op}"]        = -(2*i + 2)
    proc_indices[f"sm_lin_quad_{op}"] = -(2*i + 3)
for k, (op1, op2) in enumerate(OP_PAIRS):
    name = f"sm_lin_quad_mixed_{op1}_{op2}"
    if name in histograms:
        proc_indices[name] = -(2*n_ops + 2 + k)

n_syst   = 1 + (1 if has_scale else 0) + (1 if has_pdf else 0)
col      = lambda s: f"  {str(s):<28}"
sep_long = "-" * 130
sep_short= "-" * 45

lines = [
    sep_long,
    "imax 1",
    f"jmax {len(processes) - 1}",
    f"kmax {n_syst}",
    sep_long,
    (f"shapes *         {CHANNEL}  {args.output}"
     f"  $CHANNEL/$PROCESS  $CHANNEL/$PROCESS_$SYSTEMATIC"),
    (f"shapes data_obs  {CHANNEL}  {args.output}"
     f"  $CHANNEL/data_obs"),
    sep_long,
    f"bin          {CHANNEL}",
    f"observation  {observation:.6g}",
    sep_long,
    "bin         " + "".join(col(CHANNEL)           for _ in processes),
    "process     " + "".join(col(p)                 for p in processes),
    "process     " + "".join(col(proc_indices[p])   for p in processes),
    "rate        " + "".join(col(f"{rates[p]:.6g}") for p in processes),
    sep_short,
    "lumi        lnN  " + "".join(col("1.02") for _ in processes),
]
if has_scale:
    lines.append("muf_scale   shape" + "".join(col("1") for _ in processes))
if has_pdf:
    lines.append("pdf         shape" + "".join(col("1") for _ in processes))
lines.append("")

datacard_text = "\n".join(lines)
os.makedirs(os.path.dirname(os.path.abspath(args.datacard)), exist_ok=True)
print(f"Writing {args.datacard} ...")
with open(args.datacard, "w") as f:
    f.write(datacard_text)

print(f"\n{'-'*72}")
print(datacard_text)
