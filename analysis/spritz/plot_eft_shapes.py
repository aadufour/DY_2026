#!/usr/bin/env python3
"""
plot_eft_shapes.py — EFT shape comparison vs SM
Adapted from Giacomo's check_shapes.py.

Reads histos.root from spritz-postproc.
For each operator, plots SM (black), EFT at c=+1 (orange), EFT at c=-1 (blue)
reconstructed from the lin/quad components stored in histos.root.

Usage (from dy_smeftsim_v4/, in analysis_venv):
    python3 /path/to/plot_eft_shapes.py
    python3 /path/to/plot_eft_shapes.py --region inc_ee
    python3 /path/to/plot_eft_shapes.py --operators cHDD cHWB
"""

import argparse
import os

import uproot
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep

# --------------------------------------------------
# Operator index → real name (LHE reweighting order)
# --------------------------------------------------

OPERATORS = [
    "cHDD", "cHWB", "cbWRe", "cbBRe",
    "cHj1", "cHQ1", "cHj3",  "cHQ3",
    "cHu",  "cHd",  "cHbq",
    "cHl1", "cHl3", "cHe",
    "cll1",
    "clj1", "clj3",
    "cQl1", "cQl3",
    "ceu",  "ced",
    "cbe",  "cje",  "cQe",
    "clu",  "cld",  "cbl",
]

# --------------------------------------------------
# Configuration
# --------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("--root",      default="histos.root")
parser.add_argument("--region",    default="inc_mm", choices=["inc_ee", "inc_mm", "inc_em"])
parser.add_argument("--variable",  default="mll")
parser.add_argument("--outdir",    default="check")
parser.add_argument("--operators", nargs="+", default=None,
                    help="Subset of operators by real name (e.g. cHDD cHWB). Default: all.")
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

root_file = uproot.open(args.root)

# --------------------------------------------------
# Helper
# --------------------------------------------------

def read_hist(sample):
    path = f"{args.region}/{args.variable}/nominal/histo_{sample}"
    h = root_file[path]
    return h.values(), h.axes[0].edges()

# --------------------------------------------------
# Load SM once
# --------------------------------------------------

vals_sm, edges = read_hist("DYSMEFTsim_SM")
centers = 0.5 * (edges[:-1] + edges[1:])

# --------------------------------------------------
# Loop over operators
# --------------------------------------------------

ops_to_plot = [
    (i + 1, op) for i, op in enumerate(OPERATORS)
    if args.operators is None or op in args.operators
]

for idx, op_name in ops_to_plot:
    tag = f"op{idx:02d}"

    try:
        vals_lin,  _ = read_hist(f"DYSMEFTsim_{tag}_lin")
        vals_quad, _ = read_hist(f"DYSMEFTsim_{tag}_quad")
    except KeyError:
        print(f"  WARNING: {tag} ({op_name}) not found in ROOT file, skipping.")
        continue

    # Reconstruct EFT at c=+1 and c=-1 from lin/quad components
    # (equivalent to Giacomo's raw op and op_m1 histograms)
    vals_cp1 = vals_sm + vals_lin + vals_quad   # c = +1
    vals_cm1 = vals_sm - vals_lin + vals_quad   # c = -1

    hist_names = {
        "SM":              vals_sm,
        f"{op_name} c=+1": vals_cp1,
        f"{op_name} c=-1": vals_cm1,
    }

    colors = {
        "SM":              "black",
        f"{op_name} c=+1": "orange",
        f"{op_name} c=-1": "blue",
    }

    # --------------------------------------------------
    # Style
    # --------------------------------------------------

    hep.style.use("CMS")

    # --------------------------------------------------
    # Figure with ratio panel
    # --------------------------------------------------

    fig, (ax, rax) = plt.subplots(
        2, 1,
        figsize=(8, 8),
        gridspec_kw={"height_ratios": [3, 1]},
        sharex=True,
    )

    # --------------------------------------------------
    # Upper pad: shapes
    # --------------------------------------------------

    for label, values in hist_names.items():
        hep.histplot(
            values,
            edges,
            label=label,
            ax=ax,
            histtype="step",
            linewidth=2,
            color=colors[label],
        )

    ax.legend()
    ax.set_ylabel("Events")
    ax.set_yscale("log")
    hep.cms.label(ax=ax, data=False)

    # --------------------------------------------------
    # Ratio pad
    # --------------------------------------------------

    for label, values in hist_names.items():
        if label == "SM":
            continue

        ratio = np.divide(
            values,
            vals_sm,
            out=np.zeros_like(values, dtype=float),
            where=vals_sm != 0,
        )

        rax.step(
            edges[:-1],
            ratio,
            where="post",
            linewidth=2,
            label=f"{label}/SM",
            color=colors[label],
        )

    rax.axhline(1.0, color="black", linestyle="--")
    rax.set_ylabel("Ratio")
    rax.set_xlabel(r"$m_{\ell\ell}$ [GeV]")
    rax.set_ylim(0.5, 1.5)

    # --------------------------------------------------
    # Save
    # --------------------------------------------------

    plt.tight_layout()
    plt.savefig(f"{args.outdir}/mll_shapes_{op_name}.png", dpi=300)
    plt.savefig(f"{args.outdir}/mll_shapes_{op_name}.pdf")
    plt.close()

    print(f"Saved: {args.outdir}/mll_shapes_{op_name}.png")
