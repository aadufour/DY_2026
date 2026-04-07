#!/usr/bin/env python3
"""
plot_pdf_updown.py

Giacomo's prescription — produces 3 distributions + ratio panel:

Step 1: Fill 100 histograms, one per PDF replica.
        Each event contributes with weight = ev.weights[replica_id]
        (which is already XWGTUP * PDFWeight[i] as stored by MG5)

Step 2: For each bin b, split the 100 replica values into:
          Uppp  = { h_i(b) : h_i(b) > h_nominal(b) }
          Downnn = { h_i(b) : h_i(b) < h_nominal(b) }

Step 3: For each bin b:
          UP(b)   = RMS( Uppp  )  = sqrt( mean( h_i(b)^2 for i in Uppp   ) )
          DOWN(b) = RMS( Downnn ) = sqrt( mean( h_i(b)^2 for i in Downnn ) )

Step 4: Build the UP and DOWN distributions by setting bin content = UP(b), DOWN(b)

Plot:
  Top panel   — three step histograms: nominal, UP, DOWN
  Bottom panel — ratio: UP/nominal and DOWN/nominal vs mll
"""

import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pylhe

# ── Config ───────────────────────────────────────────────────────────────────

LHE_FILE = os.path.join(os.path.dirname(__file__), "sm_pdf_test_200k.lhe")
PLOT_DIR  = os.path.join(os.path.dirname(__file__), "plots_pdf_study")
os.makedirs(PLOT_DIR, exist_ok=True)

CENTRAL_ID = '45'
PDF_IDS    = [str(i) for i in range(46, 146)]   # replicas 1-100

MLL_EDGES   = np.logspace(np.log10(60), np.log10(800), 40)
BIN_CENTRES = 0.5 * (MLL_EDGES[:-1] + MLL_EDGES[1:])
MLL_MIN, MLL_MAX = MLL_EDGES[0], MLL_EDGES[-1]
N_BINS = len(MLL_EDGES) - 1

# ── Kinematic helper ──────────────────────────────────────────────────────────

def mll_of(p1, p2):
    p = np.array(p1) + np.array(p2)
    return np.sqrt(max(p[3]**2 - p[0]**2 - p[1]**2 - p[2]**2, 0.0))

# ── Read LHE ──────────────────────────────────────────────────────────────────

buf_mll     = []
buf_central = []
buf_pdf     = [[] for _ in PDF_IDS]

print(f"Reading {LHE_FILE} ...")
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    for i, ev in enumerate(pylhe.read_lhe_with_attributes(LHE_FILE)):
        if (i + 1) % 10000 == 0:
            print(f"  {i + 1} events processed")

        leptons = [p for p in ev.particles
                   if int(p.status) == 1 and abs(int(p.id)) in {11, 13}]
        if len(leptons) < 2:
            continue

        v = [[p.px, p.py, p.pz, p.e] for p in leptons[:2]]
        m = mll_of(v[0], v[1])
        if not (MLL_MIN <= m <= MLL_MAX):
            continue

        buf_mll.append(m)
        buf_central.append(ev.weights[CENTRAL_ID])
        for j, k in enumerate(PDF_IDS):
            buf_pdf[j].append(ev.weights[k])

mll_arr   = np.array(buf_mll)
w_central = np.array(buf_central)
w_pdf     = np.array(buf_pdf, dtype=np.float64)   # (100, N_events)
print(f"Events kept: {len(mll_arr):,}\n")

# ── Step 1: fill 100 replica histograms + nominal ─────────────────────────────

def fill(weights):
    vals, _ = np.histogram(mll_arr, bins=MLL_EDGES, weights=weights)
    return vals.astype(float)

h_nominal = fill(w_central)                                  # (N_bins,)
h_replicas = np.array([fill(w_pdf[i]) for i in range(100)]) # (100, N_bins)

# ── Steps 2 & 3: asymmetric RMS per bin ──────────────────────────────────────

h_up   = np.zeros(N_BINS)
h_down = np.zeros(N_BINS)

for b in range(N_BINS):
    vals = h_replicas[:, b]      # 100 replica values for this bin
    nom  = h_nominal[b]

    up_vals   = vals[vals > nom]
    down_vals = vals[vals < nom]

    # RMS of the replica values themselves (not deviations)
    h_up[b]   = np.sqrt(np.mean(up_vals**2))   if len(up_vals)   > 0 else nom
    h_down[b] = np.sqrt(np.mean(down_vals**2)) if len(down_vals) > 0 else nom

# ── Plot ──────────────────────────────────────────────────────────────────────

def step(ax, edges, vals, **kw):
    ax.step(np.append(edges[:-1], edges[-1]),
            np.append(vals, vals[-1]),
            where='post', **kw)

safe_nominal = np.where(h_nominal > 0, h_nominal, np.nan)

fig, (ax_top, ax_bot) = plt.subplots(
    2, 1, figsize=(8, 7),
    sharex=True,
    gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.05}
)

# ── top: three distributions ──────────────────────────────────────────────────

step(ax_top, MLL_EDGES, h_nominal, color='k',         lw=1.5, label='Nominal')
step(ax_top, MLL_EDGES, h_up,      color='steelblue', lw=1.2, ls='--', label='PDF up')
step(ax_top, MLL_EDGES, h_down,    color='tomato',    lw=1.2, ls='--', label='PDF down')

ax_top.set_yscale('log')
ax_top.set_ylabel('Events  [a.u.]', fontsize=12)
ax_top.legend(fontsize=10)
ax_top.set_xlim(MLL_MIN, MLL_MAX)

# ── bottom: ratio UP/nominal and DOWN/nominal ─────────────────────────────────

step(ax_bot, MLL_EDGES, h_up   / safe_nominal, color='steelblue', lw=1.2, ls='--', label='up / nominal')
step(ax_bot, MLL_EDGES, h_down / safe_nominal, color='tomato',    lw=1.2, ls='--', label='down / nominal')
ax_bot.axhline(1.0, color='k', lw=1.0, ls='-')

ax_bot.set_xscale('log')
ax_bot.set_xlabel(r'$m_{\ell\ell}$  [GeV]', fontsize=12)
ax_bot.set_ylabel('variation / nominal', fontsize=11)
ax_bot.legend(fontsize=9)
ax_bot.set_ylim(0.9, 1.1)
ax_bot.grid(axis='y', alpha=0.3)

fig.tight_layout()
fig.savefig(f"{PLOT_DIR}/plot_pdf_updown.pdf")
print("Saved plot_pdf_updown.pdf")
