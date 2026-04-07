#!/usr/bin/env python3
"""
plot_scale_updown.py

QCD scale uncertainty using the 7-point envelope prescription.

Weight IDs (from initrwgt header, MUR=1/MUF=1/PDF=central is the nominal ID 45):
  ID  1 : MUR=0.5  MUF=0.5
  ID  6 : MUR=0.5  MUF=1.0
  ID 11 : MUR=0.5  MUF=2.0
  ID 16 : MUR=1.0  MUF=0.5
  ID 25 : MUR=1.0  MUF=2.0
  ID 35 : MUR=2.0  MUF=1.0
  ID 40 : MUR=2.0  MUF=2.0

Prescription (standard CMS envelope):
  For each bin b:
    UP(b)   = max over 7 variations
    DOWN(b) = min over 7 variations

  Unlike PDFs, scale points are NOT a statistical ensemble — there is no
  "average" to take. The envelope captures the largest excursion in either
  direction and is the standard CMS/LHC treatment for QCD scale uncertainty.

Plot:
  Top panel    — nominal + UP + DOWN as step histograms
  Bottom panel — (UP - nominal)/nominal  [positive]
                 (DOWN - nominal)/nominal [negative]
                 directly shows the fractional scale uncertainty per bin
"""

import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pylhe

# ── Config ────────────────────────────────────────────────────────────────────

LHE_FILE = os.path.join(os.path.dirname(__file__), "sm_pdf_test_200k.lhe")
PLOT_DIR  = os.path.join(os.path.dirname(__file__), "plots_pdf_study")
os.makedirs(PLOT_DIR, exist_ok=True)

CENTRAL_ID = '45'

# All 8 non-nominal static scale variation IDs (no dynamic scale choice)
# MG5 computes all 9 combinations of (MUR, MUF) in {0.5, 1, 2}^2.
# The nominal (MUR=1, MUF=1) is ID 45 = CENTRAL_ID.
# We take the envelope over all remaining 8 points — no convention applied.
SCALE_IDS = {
    '1':  (0.5, 0.5),
    '6':  (0.5, 1.0),
    '11': (0.5, 2.0),
    '16': (1.0, 0.5),
    '25': (1.0, 2.0),
    '30': (2.0, 0.5),
    '35': (2.0, 1.0),
    '40': (2.0, 2.0),
}

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
buf_scale   = {k: [] for k in SCALE_IDS}

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
        for k in SCALE_IDS:
            buf_scale[k].append(ev.weights[k])

mll_arr   = np.array(buf_mll)
w_central = np.array(buf_central)
w_scale   = {k: np.array(buf_scale[k]) for k in SCALE_IDS}
print(f"Events kept: {len(mll_arr):,}\n")

# ── Fill histograms ───────────────────────────────────────────────────────────

def fill(weights):
    vals, _ = np.histogram(mll_arr, bins=MLL_EDGES, weights=weights)
    return vals.astype(float)

h_nominal = fill(w_central)                                         # (N_bins,)
h_vars    = np.array([fill(w_scale[k]) for k in SCALE_IDS])        # (7, N_bins)

# ── Envelope prescription ─────────────────────────────────────────────────────

h_up   = h_vars.max(axis=0)   # per-bin maximum across 7 variations
h_down = h_vars.min(axis=0)   # per-bin minimum across 7 variations

safe_nominal = np.where(h_nominal > 0, h_nominal, np.nan)

delta_up   =  (h_up   - h_nominal) / safe_nominal   # positive
delta_down =  (h_down - h_nominal) / safe_nominal   # negative

# ── Plot ──────────────────────────────────────────────────────────────────────

def step(ax, edges, vals, **kw):
    ax.step(np.append(edges[:-1], edges[-1]),
            np.append(vals, vals[-1]),
            where='post', **kw)

fig, (ax_top, ax_bot) = plt.subplots(
    2, 1, figsize=(8, 7),
    sharex=True,
    gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.05}
)

# ── Top panel: nominal + envelope ─────────────────────────────────────────────

step(ax_top, MLL_EDGES, h_nominal, color='k',       lw=1.5,       label='Nominal  (MUR=MUF=1)')
step(ax_top, MLL_EDGES, h_up,      color='tomato',  lw=1.2, ls='--', label='Scale up   (envelope max)')
step(ax_top, MLL_EDGES, h_down,    color='tomato',  lw=1.2, ls=':',  label='Scale down (envelope min)')

ax_top.set_yscale('log')
ax_top.set_ylabel('Events  [a.u.]', fontsize=12)
ax_top.legend(fontsize=10)
ax_top.set_xlim(MLL_MIN, MLL_MAX)
ax_top.set_title('QCD scale uncertainty', fontsize=12)

# ── Bottom panel: fractional deviation ────────────────────────────────────────

step(ax_bot, MLL_EDGES, delta_up,   color='tomato', lw=1.5, ls='--', label='(up $-$ nom) / nom')
step(ax_bot, MLL_EDGES, delta_down, color='tomato', lw=1.5, ls=':',  label='(down $-$ nom) / nom')
ax_bot.axhline(0.0, color='k', lw=1.0)
ax_bot.fill_between(BIN_CENTRES, delta_down, delta_up,
                    step='mid', alpha=0.15, color='tomato')

ax_bot.set_xscale('log')
ax_bot.set_xlabel(r'$m_{\ell\ell}$  [GeV]', fontsize=12)
ax_bot.set_ylabel(r'$\delta$ / nominal', fontsize=11)
ax_bot.legend(fontsize=9)
ax_bot.grid(axis='y', alpha=0.3)

fig.tight_layout()
fig.savefig(f"{PLOT_DIR}/plot_scale_updown.pdf")
print("Saved plot_scale_updown.pdf")
