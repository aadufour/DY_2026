"""
plot_eft_shapes.py  —  EFT shape comparison vs SM
===================================================
Adapted from Giacomo's check_shapes.py.
Reads histos.root produced by spritz-postproc and plots:
  - One PNG per operator: SM vs EFT at coupling c, with ratio panel
  - One summary PNG: EFT/SM ratio for all operators overlaid

The ROOT file structure is:
  {region}/{variable}/nominal/histo_DYSMEFTsim_SM
  {region}/{variable}/nominal/histo_DYSMEFTsim_op01_lin
  {region}/{variable}/nominal/histo_DYSMEFTsim_op01_quad
  ...

EFT reconstruction: H(c) = H_SM + c*H_lin + c^2*H_quad

Usage (analysis_venv, from dy_smeftsim_v4/):
    python3 plot_eft_shapes.py
    python3 plot_eft_shapes.py --root histos.root --region inc_ee --cval 2.0
    python3 plot_eft_shapes.py --summary-only
"""

import argparse
import os

import uproot
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep

plt.style.use(hep.style.CMS)

# ── Operator index → real name mapping (LHE reweighting order) ────────────────
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
]  # op01=OPERATORS[0], op02=OPERATORS[1], ...

# ── Argument parsing ──────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--root",     default="histos.root", help="Path to histos.root")
parser.add_argument("--outdir",   default="plots_eft",   help="Output directory")
parser.add_argument("--region",   default="inc_mm",      choices=["inc_ee", "inc_mm", "inc_em"])
parser.add_argument("--variable", default="mll",         help="Variable to plot")
parser.add_argument("--cval",     type=float, default=1.0, help="Wilson coefficient value")
parser.add_argument("--operators", nargs="+", default=None,
                    help="Subset of operators by real name (e.g. cHDD cHWB). Default: all.")
parser.add_argument("--summary-only", action="store_true",
                    help="Skip per-operator plots, only make summary ratio plot")
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

# ── Open ROOT file ────────────────────────────────────────────────────────────
print(f"Opening {args.root} ...")
root_file = uproot.open(args.root)

def get_hist(sample):
    """Read histogram values and edges for a given sample name."""
    path = f"{args.region}/{args.variable}/nominal/histo_{sample}"
    try:
        h = root_file[path]
        return h.values(), h.axes[0].edges()
    except KeyError:
        return None, None

# ── Load SM ───────────────────────────────────────────────────────────────────
vals_sm, edges = get_hist("DYSMEFTsim_SM")
if vals_sm is None:
    raise RuntimeError(f"SM histogram not found. Check region/variable/root file.")

x = 0.5 * (edges[:-1] + edges[1:])
print(f"SM total yield: {vals_sm.sum():.1f}")

# ── Select operators to plot ──────────────────────────────────────────────────
if args.operators:
    ops_to_plot = [(i+1, op) for i, op in enumerate(OPERATORS) if op in args.operators]
else:
    ops_to_plot = list(enumerate(OPERATORS, start=1))

# ── Per-operator plots ────────────────────────────────────────────────────────
ratio_all = {}   # op_name → ratio array (for summary)

for idx, op_name in ops_to_plot:
    tag = f"op{idx:02d}"

    vals_lin,  _ = get_hist(f"DYSMEFTsim_{tag}_lin")
    vals_quad, _ = get_hist(f"DYSMEFTsim_{tag}_quad")

    if vals_lin is None or vals_quad is None:
        print(f"  WARNING: {tag} ({op_name}) not found, skipping.")
        continue

    c = args.cval
    vals_eft = vals_sm + c * vals_lin + c**2 * vals_quad
    ratio = np.where(vals_sm > 0, vals_eft / vals_sm, np.nan)
    ratio_all[op_name] = ratio

    if args.summary_only:
        continue

    # ── Individual plot ───────────────────────────────────────────────────────
    fig, (ax, ax_r) = plt.subplots(
        2, 1, sharex=True, figsize=(10, 8),
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
    )

    hep.histplot(vals_sm,  edges, ax=ax, label="SM",
                 color="black", linewidth=2, histtype="step")
    hep.histplot(vals_eft, edges, ax=ax, label=rf"SM + {op_name}, $c={c}$",
                 color="#e41a1c", linewidth=1.8, linestyle="--", histtype="step")
    hep.histplot(vals_sm + c * vals_lin, edges, ax=ax,
                 label=rf"SM + lin only ($c={c}$)",
                 color="#377eb8", linewidth=1.2, linestyle=":", histtype="step")

    ax.set_ylabel("Events")
    ax.set_yscale("log")
    ax.legend(fontsize=10, frameon=False)
    hep.cms.label("Private", data=False, ax=ax, year="2018", lumi=59.8)

    ax_r.axhline(1.0, color="black", linewidth=1, linestyle=":")
    ax_r.step(edges[:-1], ratio, where="post", color="#e41a1c", linewidth=1.8)
    ax_r.set_ylim(0.5, 1.5)
    ax_r.set_ylabel("EFT / SM", fontsize=11)
    ax_r.set_xlabel(r"$m_{\ell\ell}$ [GeV]")

    outpath = os.path.join(args.outdir, f"{args.variable}_{args.region}_{op_name}_c{c}.png")
    fig.savefig(outpath, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outpath}")

# ── Summary plot ──────────────────────────────────────────────────────────────
if ratio_all:
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.axhline(1.0, color="black", linewidth=1.5, linestyle=":")

    colors = plt.cm.tab20(np.linspace(0, 1, len(ratio_all)))
    for (op_name, ratio), color in zip(ratio_all.items(), colors):
        ax.step(edges[:-1], ratio, where="post",
                label=op_name, color=color, linewidth=1.2, alpha=0.85)

    ax.set_ylabel("EFT / SM", fontsize=12)
    ax.set_ylim(0.5, 1.5)
    ax.set_xlabel(r"$m_{\ell\ell}$ [GeV]")
    ax.set_xscale("log")
    ax.legend(fontsize=8, frameon=True, ncol=3, loc="upper right")
    hep.cms.label("Private", data=False, ax=ax, year="2018", lumi=59.8,
                  label=rf"EFT/SM, $c={args.cval}$, {args.region}")

    outpath = os.path.join(args.outdir,
                           f"summary_{args.variable}_{args.region}_c{args.cval}.png")
    fig.savefig(outpath, bbox_inches="tight")
    plt.close()
    print(f"\nSummary saved: {outpath}")

print("Done.")
