#!/usr/bin/env python3
"""
pdf_weights_study.py

Read sm_pdf_test LHE (with PDF + scale systematics) and produce 4 diagnostic plots:

  1. mll spectrum: central + PDF +-1sigma band + scale envelope
  2. Ratio (variation / central) vs mll  [PDF replicas + scale variations]
  3. Per-event PDF replica ratio distribution (histogram)
  4. PDF uncertainty fraction per mll bin

Weight ID map (from initrwgt header):
  ID 45             : central  (MUR=1, MUF=1, PDF member 0)  ← denominator
  IDs 46 - 145      : PDF replicas 1 - 100
  IDs 1,6,11,16,25,35,40 : 7-point scale envelope (fixed scale choice only)
      1  : (muR=0.5, muF=0.5)
      6  : (muR=0.5, muF=1.0)
      11 : (muR=0.5, muF=2.0)
      16 : (muR=1.0, muF=0.5)
      25 : (muR=1.0, muF=2.0)
      35 : (muR=2.0, muF=1.0)
      40 : (muR=2.0, muF=2.0)

Usage:
    conda activate dy_analysis
    python3 pdf_weights_study.py
"""

import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pylhe

# ── Config ─────────────────────────────────────────────────────────────────

LHE_FILE = os.path.join(os.path.dirname(__file__), "sm_pdf_test_200k.lhe")
PLOT_DIR  = os.path.join(os.path.dirname(__file__), "plots_pdf_study")
os.makedirs(PLOT_DIR, exist_ok=True)

CENTRAL_ID = '45'
PDF_IDS    = [str(i) for i in range(46, 146)]   # replicas 1-100

SCALE_IDS  = {             # id : (muR, muF) label
    '1':  (0.5, 0.5),
    '6':  (0.5, 1.0),
    '11': (0.5, 2.0),
    '16': (1.0, 0.5),
    '25': (1.0, 2.0),
    '35': (2.0, 1.0),
    '40': (2.0, 2.0),
}

# 200k events: extend to 800 GeV, still log-spaced for good resolution
# near the Z peak without drowning in empty high-mll bins.
MLL_EDGES   = np.logspace(np.log10(60), np.log10(800), 40)
BIN_CENTRES = 0.5 * (MLL_EDGES[:-1] + MLL_EDGES[1:])
N_BINS      = len(MLL_EDGES) - 1

MLL_MIN, MLL_MAX = MLL_EDGES[0], MLL_EDGES[-1]

# ── Kinematic helper ────────────────────────────────────────────────────────

def mll_of(p1, p2):
    p = np.array(p1) + np.array(p2)
    return np.sqrt(max(p[3]**2 - p[0]**2 - p[1]**2 - p[2]**2, 0.0))

# ── Read LHE ────────────────────────────────────────────────────────────────

buf_mll     = []
buf_central = []
buf_pdf     = [[] for _ in PDF_IDS]   # list per replica
buf_scale   = {k: [] for k in SCALE_IDS}

print(f"Reading {LHE_FILE} ...")
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    for i, ev in enumerate(pylhe.read_lhe_with_attributes(LHE_FILE)):
        if (i + 1) % 2000 == 0:
            print(f"  {i + 1} events processed")

        leptons = [p for p in ev.particles
                   if int(p.status) == 1 and abs(int(p.id)) in {11, 13}]
        if len(leptons) < 2:
            continue

        v = [[p.px, p.py, p.pz, p.e] for p in leptons[:2]]
        m = mll_of(v[0], v[1])
        if not (MLL_MIN <= m <= MLL_MAX):
            continue

        w_c = ev.weights[CENTRAL_ID]
        buf_mll.append(m)
        buf_central.append(w_c)
        for j, k in enumerate(PDF_IDS):
            buf_pdf[j].append(ev.weights[k])
        for k in SCALE_IDS:
            buf_scale[k].append(ev.weights[k])

mll_arr   = np.array(buf_mll)
w_central = np.array(buf_central)
w_pdf     = np.array(buf_pdf, dtype=np.float64)   # (100, N_events)
w_scale   = {k: np.array(v) for k, v in buf_scale.items()}

print(f"Events kept: {len(mll_arr):,}\n")


# ── Helper ──────────────────────────────────────────────────────────────────

def hist(weights):
    vals, _ = np.histogram(mll_arr, bins=MLL_EDGES, weights=weights)
    return vals.astype(float)


def step(ax, edges, vals, **kw):
    """Draw a histogram as a step line."""
    ax.step(np.append(edges[:-1], edges[-1]),
            np.append(vals, vals[-1]),
            where='post', **kw)


# ── Build histograms ─────────────────────────────────────────────────────────

h_central = hist(w_central)
h_pdf     = np.array([hist(w_pdf[i]) for i in range(100)])  # (100, N_bins)
h_scale   = {k: hist(w_scale[k]) for k in SCALE_IDS}

h_pdf_mean = h_pdf.mean(axis=0)
h_pdf_std  = h_pdf.std(axis=0)

scale_stack = np.stack(list(h_scale.values()))   # (7, N_bins)
h_scale_lo  = scale_stack.min(axis=0)
h_scale_hi  = scale_stack.max(axis=0)


# ── Plot 1: mll spectrum with bands ─────────────────────────────────────────

fig, ax = plt.subplots(figsize=(8, 5))

ax.fill_between(BIN_CENTRES,
                h_central - h_pdf_std, h_central + h_pdf_std,
                alpha=0.7, color='steelblue', label=r'PDF $\pm 1 \sigma$ (RMS of 100 replicas)')
# ax.fill_between(BIN_CENTRES,
#                 h_scale_lo, h_scale_hi,
#                 alpha=0.35, color='tomato', label='Scale envelope  (7-point)')
step(ax, MLL_EDGES, h_central, color='k', lw=1.0
    #  ,label='Central  (MUR=MUF=1, member 0)'
     )

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(MLL_MIN, MLL_MAX)
ax.set_xlabel(r'$m_{\ell\ell}$  [GeV]', fontsize=12)
ax.set_ylabel('Events  [a.u.]', fontsize=12)
# ax.set_title('SM DY — mll with PDF and scale uncertainties')
ax.legend()
fig.tight_layout()
fig.savefig(f"{PLOT_DIR}/plot1_mll_bands.pdf")
print("Saved plot1_mll_bands.pdf")


# # ── Plot 2: ratio (variation / central) vs mll ──────────────────────────────

# fig, (ax_pdf, ax_scale) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

# # top: all 100 PDF replicas
# safe_central = np.where(h_central > 0, h_central, np.nan)
# for i in range(100):
#     ax_pdf.plot(BIN_CENTRES, h_pdf[i] / safe_central,
#                 color='steelblue', alpha=0.18, lw=0.8)
# # overlay ±1σ band
# ax_pdf.fill_between(BIN_CENTRES,
#                     (h_central - h_pdf_std) / safe_central,
#                     (h_central + h_pdf_std) / safe_central,
#                     alpha=0.35, color='steelblue', label='±1σ RMS')
# ax_pdf.axhline(1, color='k', lw=1.2, ls='--')
# ax_pdf.set_ylabel('replica / central', fontsize=11)
# ax_pdf.set_title('PDF replica ratios vs mll')
# ax_pdf.set_ylim(0.7, 1.3)
# ax_pdf.legend(fontsize=9)

# # bottom: 7-point scale variations
# colors = plt.cm.tab10(np.linspace(0, 0.7, len(SCALE_IDS)))
# for (k, (mur, muf)), col in zip(SCALE_IDS.items(), colors):
#     ratio = h_scale[k] / safe_central
#     ax_scale.plot(BIN_CENTRES, ratio, color=col, lw=1.8,
#                   label=rf'$\mu_R$={mur}, $\mu_F$={muf}')
# ax_scale.axhline(1, color='k', lw=1.2, ls='--')
# ax_scale.set_ylabel('variation / central', fontsize=11)
# ax_scale.set_xscale('log')
# ax_scale.set_xlabel(r'$m_{\ell\ell}$  [GeV]', fontsize=12)
# ax_scale.set_title('Scale variation ratios vs mll')
# ax_scale.legend(fontsize=8, ncol=2)

# fig.tight_layout()
# fig.savefig(f"{PLOT_DIR}/plot2_ratios.pdf")
# print("Saved plot2_ratios.pdf")


# # ── Plot 3: distribution of per-event PDF ratios ─────────────────────────────

# ratios_flat = (w_pdf / w_central[np.newaxis, :]).flatten()

# fig, ax = plt.subplots(figsize=(7, 4))
# ax.hist(ratios_flat, bins=120, range=(0.5, 1.5), density=True,
#         color='steelblue', alpha=0.75, edgecolor='none')
# ax.axvline(1.0, color='k', lw=1.5, ls='--', label='central')
# ax.axvline(ratios_flat.mean(), color='tomato', lw=1.5, ls='-',
#            label=f'mean = {ratios_flat.mean():.4f}')
# ax.set_xlabel(r'$w_{\rm replica}\,/\,w_{\rm central}$  per event', fontsize=12)
# ax.set_ylabel('Density', fontsize=12)
# ax.set_title('Per-event PDF replica ratios  (all 100 replicas × all events)')
# ax.legend()
# fig.tight_layout()
# fig.savefig(f"{PLOT_DIR}/plot3_ratio_distribution.pdf")
# print("Saved plot3_ratio_distribution.pdf")


# # ── Plot 4: PDF uncertainty fraction per mll bin ─────────────────────────────

# pdf_frac = h_pdf_std / safe_central

# fig, ax = plt.subplots(figsize=(7, 4))
# step(ax, MLL_EDGES, pdf_frac * 100, color='steelblue', lw=2.0)
# ax.set_xscale('log')
# ax.set_xlabel(r'$m_{\ell\ell}$  [GeV]', fontsize=12)
# ax.set_ylabel(r'$\sigma_{\rm PDF}\,/\,\sigma_{\rm central}$  [%]', fontsize=12)
# ax.set_title('PDF uncertainty fraction per mll bin')
# ax.set_ylim(0)
# ax.grid(axis='y', alpha=0.3)
# fig.tight_layout()
# fig.savefig(f"{PLOT_DIR}/plot4_pdf_fraction.pdf")
# print("Saved plot4_pdf_fraction.pdf")
