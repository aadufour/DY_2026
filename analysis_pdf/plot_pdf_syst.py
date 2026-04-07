#!/usr/bin/env python3
"""
plot_pdf_syst.py

Two-panel plot:
  TOP    — nominal mll distribution
  BOTTOM — PDF uncertainty per bin:
               +σ_up   (RMS of replica deviations above nominal)
               -σ_down (RMS of replica deviations below nominal, plotted negative)

Weight IDs (from initrwgt header of sm_pdf_test):
  '45'      → central / nominal  (MUR=1, MUF=1, PDF member 0)
  '46'-'145' → PDF replicas 1-100
"""

import warnings
import numpy as np
import matplotlib.pyplot as plt
import pylhe

# ── Config ──────────────────────────────────────────────────────────────────

LHE_FILE   = "/Users/albertodufour/code/DY2026/analysis_pdf/sm_pdf_test_200k.lhe"
OUT_FILE   = "/Users/albertodufour/code/DY2026/analysis_pdf/plots/pdf_syst.pdf"

CENTRAL_ID = '45'
PDF_IDS    = [str(i) for i in range(46, 146)]   # replicas 1-100

MLL_EDGES  = np.logspace(np.log10(70), np.log10(1000), 25)
BIN_CENTRES = 0.5 * (MLL_EDGES[:-1] + MLL_EDGES[1:])

# ── Read LHE ────────────────────────────────────────────────────────────────

buf_mll     = []
buf_central = []
buf_pdf     = [[] for _ in range(100)]   # one list per replica

print("Reading LHE...")
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    for n, ev in enumerate(pylhe.read_lhe_with_attributes(LHE_FILE)):
        if (n + 1) % 50000 == 0:
            print(f"  {n+1} events")

        # find the two final-state leptons
        leptons = [p for p in ev.particles
                   if int(p.status) == 1 and abs(int(p.id)) in {11, 13}]
        if len(leptons) < 2:
            continue

        # compute mll
        px = sum(p.px for p in leptons[:2])
        py = sum(p.py for p in leptons[:2])
        pz = sum(p.pz for p in leptons[:2])
        e  = sum(p.e  for p in leptons[:2])
        mll = np.sqrt(max(e**2 - px**2 - py**2 - pz**2, 0.0))

        if not (MLL_EDGES[0] <= mll <= MLL_EDGES[-1]):
            continue

        buf_mll.append(mll)
        buf_central.append(ev.weights[CENTRAL_ID])
        for i, pid in enumerate(PDF_IDS):
            buf_pdf[i].append(ev.weights[pid])

mll_arr  = np.array(buf_mll)
w_central = np.array(buf_central)
w_pdf     = np.array([np.array(buf_pdf[i]) for i in range(100)])  # (100, N_events)

print(f"Events in range: {len(mll_arr)}")

# ── Fill histograms ──────────────────────────────────────────────────────────

def fill(weights):
    h, _ = np.histogram(mll_arr, bins=MLL_EDGES, weights=weights)
    return h.astype(float)

h_nominal = fill(w_central)                              # shape (N_bins,)
h_pdf     = np.array([fill(w_pdf[i]) for i in range(100)])  # shape (100, N_bins)

# ── Compute σ_up and σ_down per bin ─────────────────────────────────────────
#
#   For each bin b:
#     deviations[i, b] = h_pdf[i, b] - h_nominal[b]
#
#     σ_up[b]   = RMS of deviations[i, b] for replicas i where deviation > 0
#     σ_down[b] = RMS of deviations[i, b] for replicas i where deviation < 0

deviations = h_pdf - h_nominal[np.newaxis, :]   # (100, N_bins)

N_bins   = len(MLL_EDGES) - 1
sigma_up   = np.zeros(N_bins)
sigma_down = np.zeros(N_bins)

for b in range(N_bins):
    dev = deviations[:, b]                  # 100 numbers for this bin

    up_devs   = dev[dev > 0]
    down_devs = dev[dev < 0]

    if len(up_devs) > 0:
        sigma_up[b]   = np.sqrt(np.mean(up_devs**2))
    if len(down_devs) > 0:
        sigma_down[b] = np.sqrt(np.mean(down_devs**2))

# ── Plot ─────────────────────────────────────────────────────────────────────

fig, (ax_top, ax_bot) = plt.subplots(
    2, 1,
    figsize=(8, 7),
    sharex=True,
    gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.05}
)

# helper: draw a step histogram given edges and bin values
def draw_step(ax, edges, values, **kwargs):
    ax.step(edges, np.append(values, values[-1]), where='post', **kwargs)

# ── Top panel: nominal ───────────────────────────────────────────────────────

draw_step(ax_top, MLL_EDGES, h_nominal, color='black', lw=1.5
        #   , label='Nominal (PDF member 0)'
          )
ax_top.set_yscale('log')
ax_top.set_ylabel('Events', fontsize=12)
ax_top.legend(fontsize=11)
# ax_top.set_title('SM DY   |   PDF uncertainty from NNPDF31 (100 replicas)', fontsize=12)

# ── Bottom panel: ±σ / nominal  (in %) ───────────────────────────────────────

# avoid division by zero in empty bins
safe_nominal = np.where(h_nominal > 0, h_nominal, np.nan)

frac_up   =  sigma_up   / safe_nominal * 100   # positive
frac_down = -sigma_down / safe_nominal * 100   # negative

draw_step(ax_bot, MLL_EDGES, frac_up,   color='tomato',    lw=1.5, label=r'$+\sigma_{\rm up}$')
draw_step(ax_bot, MLL_EDGES, frac_down, color='steelblue', lw=1.5, label=r'$-\sigma_{\rm down}$')
ax_bot.axhline(0, color='black', lw=0.8, ls='--')

ax_bot.set_xlabel(r'$m_{\ell\ell}$  [GeV]', fontsize=12)
ax_bot.set_ylabel(r'$\sigma_{\rm PDF}$ / nominal  [%]', fontsize=12)
ax_bot.legend(fontsize=10)

# shared x-axis
ax_bot.set_xscale('log')
ax_bot.set_xlim(MLL_EDGES[0], MLL_EDGES[-1])

import os
os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)
fig.savefig(OUT_FILE, bbox_inches='tight')
print(f"Saved → {OUT_FILE}")
plt.show()
