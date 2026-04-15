#!/usr/bin/env python3
"""
plot_pdf_syst.py

Two-panel plot per process:
  TOP    — nominal mll distribution
  BOTTOM — PDF uncertainty per bin (+σ_up / -σ_down) in %

Reads nominal/up/down directly from histograms.root (built by build_datacard.py).
"""

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import uproot

# ── Args ────────────────────────────────────────────────────────────────────

parser = argparse.ArgumentParser()
parser.add_argument("--input",   default="histograms.root")
parser.add_argument("--channel", default="triple_DY")
parser.add_argument("--syst",    default="pdf", help="Systematic name (pdf or qcd_scale)")
parser.add_argument("--procs",   nargs="+", default=None,
                    help="Processes to plot (default: all nominals)")
parser.add_argument("--outdir",  default="plots")
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

# ── Load ROOT file ───────────────────────────────────────────────────────────

f = uproot.open(args.input)
keys = [k for k, v in f.classnames().items() if v.startswith("TH")]

# find nominal processes (no Up/Down suffix, not data_obs)
channel = args.channel
syst    = args.syst

if args.procs:
    nominal_procs = args.procs
else:
    nominal_procs = sorted([
        k.split("/")[1] for k in keys
        if k.startswith(f"{channel}/")
        and not k.endswith("Up")
        and not k.endswith("Down")
        and "data_obs" not in k
    ])

print(f"Plotting systematics '{syst}' for: {nominal_procs}")

# ── Plot one figure per process ──────────────────────────────────────────────

for proc in nominal_procs:
    nom_key = f"{channel}/{proc}"
    up_key  = f"{channel}/{proc}_{syst}Up"
    dn_key  = f"{channel}/{proc}_{syst}Down"

    if up_key not in [k for k in f.classnames()]:
        print(f"  Skipping {proc} — no {syst}Up/Down found")
        continue

    h_nom = f[nom_key]
    h_up  = f[up_key]
    h_dn  = f[dn_key]

    nom_vals = h_nom.values().flatten()
    up_vals  = h_up.values().flatten()
    dn_vals  = h_dn.values().flatten()

    # bin edges: use axis from histogram
    try:
        edges = h_nom.axis().edges()
    except Exception:
        edges = np.arange(len(nom_vals) + 1, dtype=float)

    bin_centres = 0.5 * (edges[:-1] + edges[1:])

    # fractional uncertainty in %
    safe_nom  = np.where(nom_vals > 0, nom_vals, np.nan)
    frac_up   =  (up_vals - nom_vals) / safe_nom * 100
    frac_down = -(nom_vals - dn_vals) / safe_nom * 100

    # ── Figure ───────────────────────────────────────────────────────────────

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(9, 7), sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05}
    )

    def draw_step(ax, vals, **kwargs):
        ax.step(edges, np.append(vals, vals[-1]), where="post", **kwargs)

    # top: nominal + up/down
    draw_step(ax_top, nom_vals, color="black",     lw=1.5, label="nominal")
    draw_step(ax_top, up_vals,  color="tomato",    lw=1.0, ls="--", label=f"{syst} up")
    draw_step(ax_top, dn_vals,  color="steelblue", lw=1.0, ls="--", label=f"{syst} down")
    ax_top.fill_between(
        np.repeat(edges, 2)[1:-1],
        np.repeat(dn_vals, 2),
        np.repeat(up_vals, 2),
        alpha=0.15, color="gray"
    )
    ax_top.set_yscale("log")
    ax_top.set_ylabel("Events", fontsize=12)
    ax_top.set_title(f"{proc}  —  {syst} systematic", fontsize=12)
    ax_top.legend(fontsize=10)

    # bottom: fractional uncertainty
    draw_step(ax_bot, frac_up,   color="tomato",    lw=1.5, label=r"$+\sigma$")
    draw_step(ax_bot, frac_down, color="steelblue", lw=1.5, label=r"$-\sigma$")
    ax_bot.axhline(0, color="black", lw=0.8, ls="--")
    ax_bot.set_ylabel(f"{syst} / nominal  [%]", fontsize=11)
    ax_bot.set_xlabel("Unrolled bin" if len(edges) > 15 else r"$m_{\ell\ell}$  [GeV]", fontsize=12)
    ax_bot.legend(fontsize=10)

    # x-axis: log scale only if it looks like mll (few bins with large range)
    if len(edges) <= 15 and edges[-1] / edges[0] > 5:
        ax_bot.set_xscale("log")
    ax_bot.set_xlim(edges[0], edges[-1])

    out = os.path.join(args.outdir, f"{proc}_{syst}.pdf")
    fig.savefig(out, bbox_inches="tight")
    print(f"  Saved → {out}")
    plt.close(fig)

print("Done.")
