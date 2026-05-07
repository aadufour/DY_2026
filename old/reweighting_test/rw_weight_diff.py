"""
rw_weight_diff.py

For each event in the reweighted LHE file, compute:

    delta_w(op) = w_op - w_SM

where w_SM  = event.weights['SM']   (all C = 0)
      w_op  = event.weights[op_name] (one operator C = 1, rest 0)

Saves one plot per operator into OUT_DIR.
"""

import os
import warnings
import pylhe
import numpy as np
import boost_histogram as bh
import mplhep as hep
import matplotlib.pyplot as plt

# ── Config ────────────────────────────────────────────────────────────────────
LHE_FILE = (
    "/Users/albertodufour/MG5_2_9_18/mg5amcnlo/DY_all/myLHE/unweighted_events.lhe"
)
OUT_DIR = "/Users/albertodufour/code/DY2026/reweighting_test/weight_diff_plots"
os.makedirs(OUT_DIR, exist_ok=True)

OPERATORS = [
    'cHDD', 'cHWB', 'cHj1',
    'cHj3', 'cHu', 'cHd',
    'cHl1', 'cHl3', 'cHe',
    'cll1', 'clj1', 'clj3',
    'ceu',  'ced',  'cje',
    'clu',  'cld',
    'minuscHDD', 'minuscHWB', 'minuscHj1',
    'minuscHj3', 'minuscHu', 'minuscHd',
    'minuscHl1', 'minuscHl3', 'minuscHe',
    'minuscll1', 'minusclj1', 'minusclj3',
    'minusceu',  'minusced',  'minuscje',
    'minusclu',  'minuscld',
]

N_BINS   = 100
BIN_RANGE = None   # set to (lo, hi) to fix range; None = auto per operator

# ── Read LHE and collect weight differences ───────────────────────────────────
print(f"Reading {LHE_FILE}")

# dict: op_name -> list of delta_w values
deltas = {op: [] for op in OPERATORS}

with warnings.catch_warnings():
    warnings.simplefilter("ignore", DeprecationWarning)
    events = list(pylhe.read_lhe_with_attributes(LHE_FILE))

for i, event in enumerate(events):
    if not (i + 1) % 5000:
        print(f"Processed {i + 1} events")

    w_sm = event.weights['SM']

    for op in OPERATORS:
        deltas[op].append(np.abs(event.weights[op] - w_sm))
 

   
# uniq = np.unique(deltas["ced"])

# for i in uniq:
#     v = 0
#     for j in deltas["ced"]:
#         if j == i: v+=1

#     print(f"val {i} {v}")

# print(np.unique(deltas["ced"]))   

   
n_events = i + 1
print(f"Done: {n_events} events")

# ── Fill boost histograms ──────────────────────────────────
hists = {}
for op in OPERATORS:
    arr = np.array(deltas[op])
    lo, hi = arr.min(), arr.max()
    # small padding so extreme values land inside the axis
    # pad = 0.05 * (hi - lo) if hi != lo else 1.0
    # axis = bh.axis.Regular(N_BINS, lo - pad, hi + pad)
    axis = bh.axis.Regular(N_BINS, 0. ,0.1)
    h = bh.Histogram(axis, storage=bh.storage.Double())
    h.fill(arr)
    hists[op] = h


# ── Plotting ──────────────────────────────────────────────
hep.style.use("CMS")

for op in OPERATORS:
    fig, ax = plt.subplots(figsize=(6, 5))


    hep.histplot(
        hists[op],
        ax=ax,
        color="steelblue",
        alpha=0.7,
        histtype="fill",
        edgecolor="steelblue",
        linewidth=0.8,
    )

    ax.axvline(0, color="red", linewidth=1.0, linestyle="--", label="SM")
    ax.set_xlabel(r"$w_\mathrm{op} - w_\mathrm{SM}$")
    ax.set_ylabel("Events")
    ax.legend()
    # ax.set_xscale("log")
    ax.set_yscale("log")
    hep.cms.label(f"{op}", ax=ax)

    out_path = os.path.join(OUT_DIR, f"weight_diff_{op}.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out_path}")
