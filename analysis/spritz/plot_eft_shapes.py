"""
plot_eft_shapes.py  —  EFT shape comparison vs SM
===================================================
Reads the merged pkl from spritz-merge and produces:
  - One plot per operator: SM vs EFT at c=1 (lin+quad), with ratio panel
  - A summary plot: EFT/SM ratio for all operators overlaid

Usage (inside apptainer, from dy_smeftsim_v4/):
    python /grid_mnt/.../plot_eft_shapes.py --pkl condor/results_merged.pkl
    python /grid_mnt/.../plot_eft_shapes.py --pkl condor/results_merged.pkl \\
        --region inc_ee --variable mll --cval 2.0 --outdir plots_eft
"""

import argparse
import os
import pickle

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np

plt.style.use(hep.style.CMS)

# ── Operator names (order matches LHE reweighting indices 1..27) ──────────────
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

MLL_BINS = ["50_120", "120_200", "200_400", "400_600", "600_800", "800_1000", "1000_3000"]

# ── Argument parsing ──────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--pkl",      required=True,        help="Path to results_merged.pkl")
parser.add_argument("--outdir",   default="plots_eft",  help="Output directory")
parser.add_argument("--region",   default="inc_mm",     choices=["inc_ee", "inc_mm", "inc_em"])
parser.add_argument("--variable", default="mll",        help="Variable to plot")
parser.add_argument("--cval",     type=float, default=1.0, help="Wilson coefficient value")
parser.add_argument("--operators", nargs="+", default=None,
                    help="Subset of operators to plot (default: all 27)")
parser.add_argument("--summary-only", action="store_true",
                    help="Only produce the summary ratio plot, skip per-operator plots")
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

ops_to_plot = args.operators if args.operators else OPERATORS

# ── Load merged pkl ───────────────────────────────────────────────────────────
print(f"Loading {args.pkl} ...")
with open(args.pkl, "rb") as f:
    results = pickle.load(f)
print(f"  {len(results)} keys found.")

# ── Helper: sum histogram across all mll bins ─────────────────────────────────
def get_h(subsample, variable=args.variable, region=args.region):
    """Sum hist.Hist across all 7 mll bins for a given subsample."""
    h_total = None
    for b in MLL_BINS:
        key = f"DYSMEFTsim_LO_mll_{b}_{subsample}"
        if key not in results:
            continue
        h = results[key]["histos"][variable]
        try:
            h_reg = h[region, "nominal"]
        except Exception:
            try:
                h_reg = h[region]
            except Exception as e:
                print(f"  WARNING: cannot slice {key}: {e}")
                continue
        h_total = h_reg if h_total is None else h_total + h_reg
    return h_total


def eft_values(op, c, vals_sm, vals_lin, vals_quad):
    return vals_sm + c * vals_lin + c**2 * vals_quad


# ── Load SM once ──────────────────────────────────────────────────────────────
h_sm = get_h("SM")
if h_sm is None:
    raise RuntimeError("SM histogram not found — check pkl keys and region name.")

vals_sm = h_sm.values()
edges   = h_sm.axes[0].edges
x       = 0.5 * (edges[:-1] + edges[1:])

print(f"SM total yield: {vals_sm.sum():.1f}")

# ── Per-operator plots ────────────────────────────────────────────────────────
ratio_all   = {}   # op_name → ratio array (for summary plot)
missing_ops = []

for op in ops_to_plot:
    h_lin  = get_h(f"{op}_lin")
    h_quad = get_h(f"{op}_quad")

    if h_lin is None or h_quad is None:
        print(f"  WARNING: {op} not found in pkl, skipping.")
        missing_ops.append(op)
        continue

    vals_lin  = h_lin.values()
    vals_quad = h_quad.values()
    vals_eft  = eft_values(op, args.cval, vals_sm, vals_lin, vals_quad)

    ratio = np.where(vals_sm > 0, vals_eft / vals_sm, np.nan)
    ratio_all[op] = ratio

    if args.summary_only:
        continue

    # ── Individual operator plot ──────────────────────────────────────────────
    fig, (ax, ax_r) = plt.subplots(
        2, 1, sharex=True,
        figsize=(10, 8),
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
    )

    # SM
    hep.histplot(vals_sm, edges, ax=ax,
                 label="SM", color="black", linewidth=2, histtype="step")

    # EFT
    hep.histplot(vals_eft, edges, ax=ax,
                 label=f"SM + {op}, c={args.cval}",
                 color="#e41a1c", linewidth=1.8, linestyle="--", histtype="step")

    # Linear term only (shows sign and magnitude of interference)
    vals_lin_only = vals_sm + args.cval * vals_lin
    hep.histplot(vals_lin_only, edges, ax=ax,
                 label=f"SM + lin only (c={args.cval})",
                 color="#377eb8", linewidth=1.2, linestyle=":", histtype="step")

    ax.set_ylabel("Events")
    ax.set_yscale("log")
    ax.legend(fontsize=10, frameon=False)
    hep.cms.label("Private", data=False, ax=ax, year="2018", lumi=59.8)

    # Ratio
    ax_r.axhline(1.0, color="black", linewidth=1, linestyle=":")
    ax_r.step(edges[:-1], ratio, where="post", color="#e41a1c", linewidth=1.8)
    ax_r.set_ylabel("EFT / SM", fontsize=11)
    ax_r.set_ylim(0.5, 1.5)
    ax_r.set_xlabel(rf"$m_{{\ell\ell}}$ [GeV]" if variable == "mll" else args.variable)

    outpath = os.path.join(args.outdir, f"{args.variable}_{args.region}_{op}_c{args.cval}.png")
    fig.savefig(outpath, bbox_inches="tight")
    plt.close()
    print(f"  Saved {outpath}")

# ── Summary plot: all operators' ratios overlaid ──────────────────────────────
if ratio_all:
    fig, ax = plt.subplots(figsize=(12, 6))

    cmap = plt.cm.tab20(np.linspace(0, 1, len(ratio_all)))
    ax.axhline(1.0, color="black", linewidth=1.5, linestyle=":")

    for (op, ratio), color in zip(ratio_all.items(), cmap):
        # Only plot where SM > 0
        ax.step(edges[:-1], ratio, where="post",
                label=op, color=color, linewidth=1.2, alpha=0.85)

    ax.set_ylabel("EFT / SM", fontsize=12)
    ax.set_ylim(0.5, 1.5)
    ax.set_xlabel(rf"$m_{{\ell\ell}}$ [GeV]" if args.variable == "mll" else args.variable)
    ax.set_xscale("log")
    ax.legend(fontsize=7, frameon=True, ncol=3, loc="upper right")
    hep.cms.label("Private", data=False, ax=ax, year="2018", lumi=59.8,
                  label=f"EFT/SM, c={args.cval}, {args.region}")

    outpath = os.path.join(args.outdir, f"summary_{args.variable}_{args.region}_c{args.cval}.png")
    fig.savefig(outpath, bbox_inches="tight")
    plt.close()
    print(f"\nSummary plot saved: {outpath}")

if missing_ops:
    print(f"\nMissing operators (not in pkl): {missing_ops}")

print("Done.")
