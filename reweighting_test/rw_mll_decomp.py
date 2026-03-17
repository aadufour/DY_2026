"""
rw_mll_decomp.py

Plot the m_ll distribution decomposed into:
  - SM          weight = w_SM
  - Linear      weight = C   * 0.5*(w[op] - w[minus_op])
  - Quadratic   weight = C^2 * (0.5*(w[op] + w[minus_op]) - w_SM)
  - Full        weight = w_SM + C*w_lin + C^2*w_quad

Usage:
  python3 rw_mll_decomp.py                      # all 17 operators, C=1
  python3 rw_mll_decomp.py --cHDD               # only cHDD
  python3 rw_mll_decomp.py --cHDD --changeC 2.5 # cHDD with C=2.5
"""

import os
import warnings
import argparse
import numpy as np
import boost_histogram as bh
import mplhep as hep
import matplotlib.pyplot as plt
import pylhe

# ── Config ────────────────────────────────────────────────────────────────────
LHE_FILE = (
    "/Users/albertodufour/MG5_2_9_18/mg5amcnlo/DY_all/myLHE/unweighted_events.lhe"
)
OUT_DIR = "/Users/albertodufour/code/DY2026/reweighting_test/mll_decomp_plots"

OPERATORS = [
    'cHDD', 'cHWB', 'cHj1',
    'cHj3', 'cHu',  'cHd',
    'cHl1', 'cHl3', 'cHe',
    'cll1', 'clj1', 'clj3',
    'ceu',  'ced',  'cje',
    'clu',  'cld',
]

# m_ll binning
N_BINS    = 100
MLL_LO    = 50.0   # GeV
MLL_HI    = 1000.0  # GeV

# ── Argument parsing ──────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
for op in OPERATORS:
    parser.add_argument(f'--{op}', action='store_true', default=False)
parser.add_argument('--changeC', type=float, default=1.0,
                    help='Wilson coefficient value C (default: 1.0)')
args = parser.parse_args()

C = args.changeC
requested = [op for op in OPERATORS if getattr(args, op)]
ops_to_plot = requested if requested else OPERATORS
print(f"Operators to plot: {ops_to_plot}")
print(f"C = {C}")

# ── Read LHE — single pass, collect arrays ────────────────────────────────────
print(f"\nReading {LHE_FILE}")

mll_vals = []
w_SM     = []
w_p1     = {op: [] for op in ops_to_plot}   # C = +1
w_m1     = {op: [] for op in ops_to_plot}   # C = -1

with warnings.catch_warnings():
    warnings.simplefilter("ignore", DeprecationWarning)
    events = pylhe.read_lhe_with_attributes(LHE_FILE)

    for i, event in enumerate(events):
        if (i + 1) % 5000 == 0:
            print(f"  processed {i + 1} events")

        # ── compute m_ll ──────────────────────────────────────────────────────
        # MG5 LHE: we know the final state is exactly two leptons
        l1, l2 = [p for p in event.particles if int(p.status) == 1]
        px = l1.px + l2.px
        py = l1.py + l2.py
        pz = l1.pz + l2.pz
        e  = l1.e  + l2.e
        mll = np.sqrt(max(e**2 - px**2 - py**2 - pz**2, 0.0))

        if not (MLL_LO <= mll <= MLL_HI):
            continue

        mll_vals.append(mll)
        w_SM.append(event.weights['SM'])
        for op in ops_to_plot:
            w_p1[op].append(event.weights[op])
            w_m1[op].append(event.weights[f'minus{op}'])

n_ev = len(mll_vals)
print(f"Done — {n_ev} events in [{MLL_LO}, {MLL_HI}] GeV\n")

# ── Convert to numpy ──────────────────────────────────────────────────────────
mll_arr = np.array(mll_vals)
wSM_arr = np.array(w_SM)

axis = bh.axis.Regular(N_BINS, MLL_LO, MLL_HI)

# SM histogram (shared across all operator plots)
h_SM = bh.Histogram(axis, storage=bh.storage.Weight())
h_SM.fill(mll_arr, weight=wSM_arr)

# ── Fill one pair of histograms per operator ───────────────────────────────────
os.makedirs(OUT_DIR, exist_ok=True)
hep.style.use("CMS")

for op in ops_to_plot:
    wp1 = np.array(w_p1[op])
    wm1 = np.array(w_m1[op])

    w_lin  = 0.5 * (wp1 - wm1)
    w_quad = 0.5 * (wp1 + wm1) - wSM_arr
    w_full = wSM_arr + C * w_lin + C**2 * w_quad

    h_lin  = bh.Histogram(axis, storage=bh.storage.Weight())
    h_quad = bh.Histogram(axis, storage=bh.storage.Weight())
    h_full = bh.Histogram(axis, storage=bh.storage.Weight())
    h_lin.fill(mll_arr,  weight= C * w_lin)
    h_quad.fill(mll_arr, weight= C**2 * w_quad)
    h_full.fill(mll_arr, weight=w_full)

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig, (ax, ax_ratio) = plt.subplots(
        2, 1, figsize=(7, 7),
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
        sharex=True,
    )

    hep.histplot(h_SM,   ax=ax, label="SM",   color="black",       linewidth=1.5, histtype="step")
    hep.histplot(h_lin,  ax=ax, label="Lin",  color="steelblue",   linewidth=1.5, histtype="step", linestyle="--")
    hep.histplot(h_quad, ax=ax, label="Quad", color="firebrick",   linewidth=1.5, histtype="step", linestyle=":")
    hep.histplot(h_full, ax=ax, label="Full", color="forestgreen", linewidth=1.5, histtype="step", linestyle="-.")

    # ── Ratio panel: Full / SM ─────────────────────────────────────────────────
    sm_vals   = h_SM.values()
    full_vals = h_full.values()
    sm_var    = h_SM.variances()
    full_var  = h_full.variances()

    safe      = sm_vals != 0
    ratio     = np.where(safe, full_vals / np.where(safe, sm_vals, 1), np.nan)
    ratio_err = np.where(
        safe & (full_vals != 0),
        ratio * np.sqrt(
            full_var / np.where(full_vals != 0, full_vals**2, 1) +
            sm_var   / np.where(safe,            sm_vals**2,   1)
        ),
        np.nan,
    )

    bin_centers = h_SM.axes[0].centers
    ax_ratio.errorbar(
        bin_centers, ratio, yerr=ratio_err,
        fmt='o', color="forestgreen", markersize=3, linewidth=0.8,
    )
    ax_ratio.axhline(1.0, color="black", linewidth=0.8, linestyle="--")
    ax_ratio.set_ylabel("Full / SM")
    ax_ratio.set_ylim(0.5, 1.5)
    ax_ratio.set_xlabel(r"$m_{\ell\ell}$ [GeV]")

    ax.set_ylabel("Weighted events")
    ax.legend(frameon=False, loc="upper right", fontsize=9)
    ax.loglog()
    hep.cms.label(f"{op}, $C={C}$", ax=ax, loc=0, data=True, lumi=None)

    out_path = os.path.join(OUT_DIR, f"mll_decomp_{op}_C{C}.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out_path}")

print("\nAll done.")
