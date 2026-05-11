"""
plot_smeft_zpeak.py

Plots mll in the Z peak region (50-100 GeV) for:
  - SM (LHEReweightingWeight[0])
  - SM + EFT at varying coupling values

Usage (inside apptainer):
    python plot_smeft_zpeak.py --pkl condor/results_merged_new.pkl \
                                --region inc_mm \
                                --operator 1 \
                                --cvals 0.5 1.0 2.0

Operator index mapping (1..27 = op_i at +1):
    Use build_cache.py to map indices to operator names.
    Default: operator 1 (first EFT operator).
"""

import argparse
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
from spritz.framework.framework import read_chunks, add_dict_iterable

plt.style.use(hep.style.CMS)

# ── Argument parsing ──────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--pkl",      default="condor/results_merged_new.pkl",
                    help="Path to merged pkl from spritz-merge")
parser.add_argument("--region",   default="inc_mm",
                    choices=["inc_ee", "inc_mm", "inc_em"],
                    help="Region to plot")
parser.add_argument("--operator", type=int, default=1,
                    help="Operator index (1..27) for EFT reweighting")
parser.add_argument("--cvals",    type=float, nargs="+", default=[0.5, 1.0, 2.0],
                    help="Wilson coefficient values to plot")
parser.add_argument("--mll-lo",   type=float, default=50.0)
parser.add_argument("--mll-hi",   type=float, default=100.0)
parser.add_argument("--nbins",    type=int,   default=50)
parser.add_argument("--outdir",   default="plots_zpeak")
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

MLL_BINS = [
    "50_120", "120_200", "200_400", "400_600",
    "600_800", "800_1000", "1000_3000",
]
DATASETS = [f"DYSMEFTsim_LO_mll_{b}" for b in MLL_BINS]

# ── Load merged pkl ───────────────────────────────────────────────────────────
print(f"Loading {args.pkl} ...")
chunks = read_chunks(args.pkl)

# Accumulate per-event arrays across all chunks and datasets
mll_all    = []
weight_all = []
w_sm_all   = []
w_kp_all   = []   # operator k at +1
w_km_all   = []   # operator k at -1  (index = k + 27)

op_idx_p = args.operator          # 1..27
op_idx_m = args.operator + 27     # 28..54

n_events_total = 0

for chunk in chunks:
    result = chunk.get("result", {})
    if not result:
        continue
    real = result.get("real_results", {})
    for dataset in DATASETS:
        if dataset not in real:
            continue
        ev = real[dataset]["events"].get(args.region, {})
        if not ev or len(ev.get("mll", [])) == 0:
            continue

        mll    = np.array(ev["mll"],              dtype=np.float64)
        weight = np.array(ev["weight"],           dtype=np.float64)
        w_sm   = np.array(ev[f"w_0"],             dtype=np.float64)
        w_kp   = np.array(ev[f"w_{op_idx_p}"],   dtype=np.float64)
        w_km   = np.array(ev[f"w_{op_idx_m}"],   dtype=np.float64)

        # mask to Z peak region
        mask = (mll >= args.mll_lo) & (mll <= args.mll_hi)
        mll_all.append(mll[mask])
        weight_all.append(weight[mask])
        w_sm_all.append(w_sm[mask])
        w_kp_all.append(w_kp[mask])
        w_km_all.append(w_km[mask])
        n_events_total += mask.sum()

if n_events_total == 0:
    print(f"No events found in region {args.region} with mll in [{args.mll_lo}, {args.mll_hi}]")
    exit(1)

mll    = np.concatenate(mll_all)
weight = np.concatenate(weight_all)
w_sm   = np.concatenate(w_sm_all)
w_kp   = np.concatenate(w_kp_all)
w_km   = np.concatenate(w_km_all)

print(f"Total events after selection: {len(mll)}")

# ── EFT decomposition ─────────────────────────────────────────────────────────
w_lin  = 0.5 * (w_kp - w_km)
w_quad = 0.5 * (w_kp + w_km) - w_sm

# ── Histogram bins ────────────────────────────────────────────────────────────
edges = np.linspace(args.mll_lo, args.mll_hi, args.nbins + 1)

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, (ax, ax_ratio) = plt.subplots(
    2, 1, sharex=True,
    gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
    figsize=(8, 7), dpi=150,
)

# SM
nom_w = weight * w_sm
h_sm, _ = np.histogram(mll, bins=edges, weights=nom_w)
h_sm2, _ = np.histogram(mll, bins=edges, weights=nom_w**2)

hep.histplot(h_sm, edges, ax=ax, label="SM", color="black", linewidth=2, histtype="step")

colors = plt.cm.plasma(np.linspace(0.2, 0.8, len(args.cvals)))

for c, color in zip(args.cvals, colors):
    w_eft = w_sm + c * w_lin + c**2 * w_quad
    eft_w = weight * w_eft
    h_eft, _ = np.histogram(mll, bins=edges, weights=eft_w)
    label = f"SM+EFT op{args.operator}, c={c}"
    hep.histplot(h_eft, edges, ax=ax, label=label, color=color, linewidth=1.5,
                 histtype="step", linestyle="--")

    # ratio EFT/SM
    ratio = np.where(h_sm > 0, h_eft / h_sm, np.nan)
    x = 0.5 * (edges[:-1] + edges[1:])
    ax_ratio.plot(x, ratio, color=color, linewidth=1.5, linestyle="--")

ax_ratio.axhline(1.0, color="black", linewidth=1, linestyle=":")
ax_ratio.set_ylabel("EFT / SM")
ax_ratio.set_ylim(0.5, 1.5)
ax_ratio.set_xlabel(r"$m_{\ell\ell}$ [GeV]")

ax.set_ylabel("Weighted events")
ax.set_yscale("log")
ax.legend(fontsize=9, frameon=False)
hep.cms.label("Private", data=False, ax=ax, year="2018")

outpath = os.path.join(args.outdir, f"zpeak_{args.region}_op{args.operator}.png")
fig.savefig(outpath, bbox_inches="tight")
plt.close()
print(f"Saved {outpath}")
