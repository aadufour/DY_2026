"""
build_datacard.py

Build the Combine input ROOT file and datacard.txt from the LHE cache.
Supports one or multiple EFT operators simultaneously.

Systematics implemented
-----------------------------------
  lumi       lnN    2% flat on all processes
  qcd_scale  shape  Envelope of MUR/MUF variations (factorised from SM)
  pdf        shape  RMS over PDF replicas (NNPDF MC method, factorised from SM)

Factorisation assumption for EFT templates:
  The relative shape variation is taken from the SM template and applied
  multiplicatively to all EFT components (quad, sm_lin_quad, mixed).

Usage:
    python3 build_datacard.py --op cHDD
    python3 build_datacard.py --op cHDD --C 0.5 --lumi 59740
    python3 build_datacard.py --op cHDD --unrolled
    python3 build_datacard.py --op cHDD cHWB --C 1.0 --lumi 59740
    python3 build_datacard.py --all_op --lumi 59740
"""

import argparse
import os
import pickle
from itertools import combinations

import boost_histogram as bh
import numpy as np
import uproot

# --------- Config ---------------------------------------

CACHE_FILE    = "/grid_mnt/data__data.polcms/cms/adufour/MG5/mg5amcnlo/CACHE/lhe_cache.pkl"
OUTPUT_FILE   = "/grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/combine/histograms.root"
DATACARD_FILE = "/grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/combine/datacard.txt"
CHANNEL       = "triple_DY"

MLL_EDGES   = np.array([50, 70, 90, 110, 200, 800, 1400, 2000, 2400, 3000], dtype=float)
RAP_EDGES   = np.array([0.0, 0.5, 1.0, 2.5], dtype=float)
CSTAR_EDGES = np.array([-1.0, -0.5, 0.0, 0.5, 1.0], dtype=float)

# ---------- Args ----------------------------------------------

parser = argparse.ArgumentParser()
op_group = parser.add_mutually_exclusive_group(required=True)
op_group.add_argument("--op",     nargs="+")
op_group.add_argument("--all_op", action="store_true")
parser.add_argument("--C",        type=float, nargs="+", default=[1.0])
parser.add_argument("--lumi",     type=float, default=1.0,
                    help="Luminosity in pb^-1 (default: 1.0 = raw weights in pb)")
parser.add_argument("--unrolled", action="store_true")
parser.add_argument("--cache",    default=CACHE_FILE)
parser.add_argument("--output",   default=OUTPUT_FILE)
parser.add_argument("--datacard", default=DATACARD_FILE)
args = parser.parse_args()

LUMI = args.lumi

# --------------- Load cache ----------------------------------------

if not os.path.exists(args.cache):
    raise FileNotFoundError(f"Cache not found: {args.cache}\nRun build_cache.py first.")

print("Loading cache ...")
with open(args.cache, "rb") as f:
    cache = pickle.load(f)

mll_arr     = cache["mll"]
rap_arr     = cache["rap"]
cstar_arr   = cache["cstar"]
w_SM        = cache["w_SM"]
w_p1_all    = cache["w_p1"]
w_m1_all    = cache["w_m1"]
w_pp_all    = cache.get("w_pp",    {})
w_scale_all     = cache.get("w_scale",       {})
w_pdf_all       = cache.get("w_pdf",         {})
w_pdf_central   = cache.get("w_pdf_central", None)

has_scale = bool(w_scale_all)
has_pdf   = bool(w_pdf_all) and w_pdf_central is not None

print(f"  {len(mll_arr):,} events loaded")
print(f"  Scale variation keys : {len(w_scale_all)}"
      + ("" if has_scale else "  (none — will be skipped)"))
print(f"  PDF replica keys     : {len(w_pdf_all)}"
      + ("" if has_pdf else "  (none — will be skipped)"))
print()

# -------------- Operator list ---------------------------------------

if args.all_op:
    OPERATORS = sorted(w_p1_all.keys())
    print(f"--all_op: {len(OPERATORS)} operators found in cache")
else:
    OPERATORS = args.op

missing = [op for op in OPERATORS if op not in w_p1_all]
if missing:
    raise KeyError(f"Operators not in cache: {missing}\nAvailable: {sorted(w_p1_all.keys())}")

# -------------- C values -----------------------------------------

if len(args.C) == 1:
    C_values = {op: args.C[0] for op in OPERATORS}
elif len(args.C) == len(OPERATORS):
    C_values = {op: c for op, c in zip(OPERATORS, args.C)}
else:
    raise ValueError(f"--C: expected 1 or {len(OPERATORS)} values, got {len(args.C)}")

# -------------------- Summary ----------------------------------------

print(f"Operators  : {OPERATORS}")
print(f"C_ref      : { {op: C_values[op] for op in OPERATORS} }")
print(f"Luminosity : {LUMI} pb^-1")
print(f"Observable : {'unrolled 3D (mll x rap x cstar)' if args.unrolled else '1D mll'}")
syst_list = ["lumi"] + (["qcd_scale"] if has_scale else []) + (["pdf"] if has_pdf else [])
print(f"Systematics: {', '.join(syst_list)}")
print()

# ------------------ Histograms ----------------------------------

def make_hist(weights, label=""):
    if args.unrolled:
        h3 = bh.Histogram(
            bh.axis.Variable(MLL_EDGES),
            bh.axis.Variable(RAP_EDGES),
            bh.axis.Variable(CSTAR_EDGES),
            storage=bh.storage.Weight())
        h3.fill(mll_arr, rap_arr, cstar_arr, weight=weights * LUMI)
        v  = h3.view()
        n  = (len(MLL_EDGES)-1) * (len(RAP_EDGES)-1) * (len(CSTAR_EDGES)-1)
        ax = bh.axis.Regular(n, 0, n,
            metadata="Unrolled bin (m_{ll} x y_{ll} x cos(theta*))")
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

def _apply_ratio(h_nom, ratio):
    """Return copy of h_nom with bins multiplied by ratio array."""
    h = h_nom.copy()
    h.view()["value"]    = h_nom.view()["value"]    * ratio
    h.view()["variance"] = h_nom.view()["variance"] * ratio**2
    return h

def _safe_ratio(var_vals, nom_vals):
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(nom_vals > 0, var_vals / nom_vals, 1.0)

def _pdf_updown(w_nom, h_nom):
    """Asymmetric RMS envelope over PDF replicas for a given process weight array."""
    central_vals = make_hist(w_nom * w_pdf_central).values().flatten()
    rep_vals     = np.array([make_hist(w_nom * w_pdf_all[k]).values().flatten()
                             for k in w_pdf_all])
    sigma_up   = np.zeros_like(central_vals)
    sigma_down = np.zeros_like(central_vals)
    for b in range(central_vals.shape[0]):
        dev = rep_vals[:, b] - central_vals[b]
        up_dev = dev[dev > 0]
        dn_dev = dev[dev < 0]
        if len(up_dev) > 0: sigma_up[b]   = np.sqrt(np.mean(up_dev**2))
        if len(dn_dev) > 0: sigma_down[b] = np.sqrt(np.mean(dn_dev**2))
    shape = h_nom.values().shape
    up_ratio = _safe_ratio((central_vals + sigma_up).reshape(shape),          h_nom.values())
    dn_ratio = _safe_ratio(np.maximum(central_vals - sigma_down, 0).reshape(shape), h_nom.values())
    return _apply_ratio(h_nom, up_ratio), _apply_ratio(h_nom, dn_ratio)

# ----------- Pre-compute scale variation keys (exclude central) ------------------------------

scale_keys = []
if has_scale:
    all_scale_keys = list(w_scale_all.keys())
    central_idx = next(
        (i for i, k in enumerate(all_scale_keys)
         if "1.0" in k and k.lower().count("1.") >= 2),
        None)
    scale_keys = [k for i, k in enumerate(all_scale_keys) if i != central_idx]
    print(f"QCD scale: {len(scale_keys)} non-central variations\n")

# ---- Build nominal process histograms -----------------------------------------------------------

histograms   = {}
proc_weights = {}
histograms["sm"]       = make_hist(w_SM, "SM")
histograms["data_obs"] = make_hist(w_SM, "data_obs (Asimov = SM)")
proc_weights["sm"]     = w_SM

for op in OPERATORS:
    C      = C_values[op]
    wp1    = w_p1_all[op]
    wm1    = w_m1_all[op]
    w_lin  = 0.5 * (wp1 - wm1)
    w_quad = 0.5 * (wp1 + wm1) - w_SM
    w_slq  = w_SM + C * w_lin + C**2 * w_quad
    if w_slq.sum() < 0:
        print(f"WARNING: [{op}] total cross section negative at C={C}")
    histograms[f"quad_{op}"]        = make_hist(C**2 * w_quad, f"quad {op} C={C}")
    histograms[f"sm_lin_quad_{op}"] = make_hist(w_slq,         f"SM+lin+quad {op} C={C}")
    proc_weights[f"quad_{op}"]        = C**2 * w_quad
    proc_weights[f"sm_lin_quad_{op}"] = w_slq

OP_PAIRS = list(combinations(OPERATORS, 2))
for op1, op2 in OP_PAIRS:
    pair = (op1, op2) if (op1, op2) in w_pp_all else (op2, op1)
    if pair not in w_pp_all:
        print(f"WARNING: no cross-weight for ({op1},{op2}), skipping.")
        continue
    C1, C2  = C_values[op1], C_values[op2]
    w_inter = w_pp_all[pair] - w_p1_all[op1] - w_p1_all[op2] + w_SM
    w_mixed = C1 * C2 * w_inter
    histograms[f"sm_lin_quad_mixed_{op1}_{op2}"]  = make_hist(w_mixed, f"mixed {op1}x{op2}")
    proc_weights[f"sm_lin_quad_mixed_{op1}_{op2}"] = w_mixed

# ---- Build Up/Down shape variants for all processes (except data_obs) ----------------------

nominal_procs = [k for k in histograms if k != "data_obs"]

for proc in nominal_procs:
    h     = histograms[proc]
    w_nom = proc_weights[proc]
    if has_scale:
        scale_var_vals = np.array([
            make_hist(w_nom * w_scale_all[k]).values() for k in scale_keys
        ])
        scale_up_ratio   = _safe_ratio(scale_var_vals.max(axis=0), h.values())
        scale_down_ratio = _safe_ratio(scale_var_vals.min(axis=0), h.values())
        histograms[f"{proc}_qcd_scaleUp"]   = _apply_ratio(h, scale_up_ratio)
        histograms[f"{proc}_qcd_scaleDown"] = _apply_ratio(h, scale_down_ratio)
    if has_pdf:
        histograms[f"{proc}_pdfUp"], histograms[f"{proc}_pdfDown"] = _pdf_updown(w_nom, h)

    #implementing asymmetric pdf 
    # if has_pdf:
    #     pdf_var_vals = np.array([make_hist(w_pdf_all[k]).values() for k in w_pdf_all])
        
    #     up_sum   = np.zeros_like(sm_nom_vals)
    #     dn_sum   = np.zeros_like(sm_nom_vals)
    #     up_count = np.zeros_like(sm_nom_vals)
    #     dn_count = np.zeros_like(sm_nom_vals)
        
    #     for rep in pdf_var_vals:
    #         above = rep > sm_nom_vals
    #         up_sum   += np.where(above, (rep - sm_nom_vals)**2, 0)
    #         dn_sum   += np.where(~above, (rep - sm_nom_vals)**2, 0)
    #         up_count += above
    #         dn_count += ~above
        
    #     sigma_up = np.sqrt(np.where(up_count > 0, up_sum / up_count, 0))
    #     sigma_dn = np.sqrt(np.where(dn_count > 0, dn_sum / dn_count, 0))
        
    #     pdf_up_ratio   = _safe_ratio(sm_nom_vals + sigma_up, sm_nom_vals)
    #     pdf_down_ratio = _safe_ratio(np.maximum(sm_nom_vals - sigma_dn, 0), sm_nom_vals)

# ---- Print nominal summary --------------------------------------------------------

print(f"{'Process':<42}  {'Integral':>14}  {'sqrt(SumW2)':>14}")
print("-" * 74)
for name in nominal_procs + ["data_obs"]:
    h    = histograms[name]
    intg = h.values().sum()
    unc  = np.sqrt(h.variances().sum())
    print(f"  {name:<40}  {intg:>14.4e}  {unc:>14.4e}")
print()

# ---- Write ROOT file --------------------------------------

os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
print(f"Writing {args.output} ...")
with uproot.recreate(args.output) as rf:
    for name, h in histograms.items():
        rf[f"{CHANNEL}/{name}"] = h
print(f"  {len(histograms)} histograms written\n")

# ---- Process list, indices, rates ----------------------

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

# ---- Write datacard -------------------------

n_syst    = 1 + (1 if has_scale else 0) + (1 if has_pdf else 0)
col       = lambda s: f"  {str(s):<28}"
sep_long  = "-" * 130
sep_short = "-" * 45

lines = [
    sep_long,
    "imax 1",
    f"jmax {len(processes) - 1}",
    f"kmax {n_syst}",
    sep_long,
    (f"shapes *         {CHANNEL}  {os.path.abspath(args.output)}"
     f"  $CHANNEL/$PROCESS  $CHANNEL/$PROCESS_$SYSTEMATIC"),
    (f"shapes data_obs  {CHANNEL}  {os.path.abspath(args.output)}"
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
    lines.append("qcd_scale   shape" + "".join(col("1") for _ in processes))
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
