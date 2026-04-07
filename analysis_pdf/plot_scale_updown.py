#!/usr/bin/env python3
"""
plot_scale_updown.py


Weight IDs (MUR=1/MUF=1/PDF=central is the nominal ID 45):
  ID  1 : MUR=0.5  MUF=0.5
  ID  6 : MUR=0.5  MUF=1.0
  ID 11 : MUR=0.5  MUF=2.0
  ID 16 : MUR=1.0  MUF=0.5
  ID 25 : MUR=1.0  MUF=2.0
  ID 30 : MUR=2.0  MUF=0.5
  ID 35 : MUR=2.0  MUF=1.0
  ID 40 : MUR=2.0  MUF=2.0

Prescription (standard CMS envelope):
  For each bin b:
    UP(b)   = max over 7 variations
    DOWN(b) = min over 7 variations


Plot:
  Top panel    - nominal + UP + DOWN as step histograms
  Bottom panel - (UP - nominal)/nominal  [positive]
                 (DOWN - nominal)/nominal [negative]
                 directly shows the fractional scale uncertainty per bin
"""

import os
import pickle

import matplotlib.pyplot as plt
import numpy as np

# ---- Config ----------------------

CACHE_FILE = os.path.join(os.path.dirname(__file__), "cache_pdf.pkl")
PLOT_DIR   = os.path.join(os.path.dirname(__file__), "plots_pdf_study")
os.makedirs(PLOT_DIR, exist_ok=True)

MLL_EDGES        = np.logspace(np.log10(60), np.log10(3000), 40)
BIN_CENTRES      = 0.5 * (MLL_EDGES[:-1] + MLL_EDGES[1:])
N_BINS           = len(MLL_EDGES) - 1
MLL_MIN, MLL_MAX = MLL_EDGES[0], MLL_EDGES[-1]

# --- Load cache --------------------------------

print(f"Loading {CACHE_FILE} ...")
with open(CACHE_FILE, 'rb') as f:
    cache = pickle.load(f)

mll_arr   = cache['mll']        # (N,)
w_central = cache['w_central']  # (N,)
w_scale   = cache['w_scale']    # (8, N) - rows ordered as cache['scale_ids']
scale_ids = cache['scale_ids']  # ['1','6','11','16','25','30','35','40']
print(f"Events loaded: {len(mll_arr):,}\n")

# ---- Fill histograms ----------------------------

def fill(weights):
    vals, _ = np.histogram(mll_arr, bins=MLL_EDGES, weights=weights)
    return vals.astype(float)

h_nominal = fill(w_central)                                              # (N_bins,)
h_vars    = np.array([fill(w_scale[i]) for i in range(len(scale_ids))])  # (8, N_bins)

# ----- Envelope prescription --------------

h_up   = h_vars.max(axis=0)   
h_down = h_vars.min(axis=0)   

# check
safe_nominal = np.where(h_nominal > 0, h_nominal, np.nan)

delta_up   =  (h_up   - h_nominal) / safe_nominal   
delta_down =  (h_down - h_nominal) / safe_nominal  

# ------ Plot -----------------------------------

def step(ax, edges, vals, **kw):
    ax.step(np.append(edges[:-1], edges[-1]),
            np.append(vals, vals[-1]),
            where='post', **kw)

fig, (ax_top, ax_bot) = plt.subplots(
    2, 1, figsize=(8, 7),
    sharex=True,
    gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.05}
)

# --------------- Top panel---------------------------

step(ax_top, MLL_EDGES, h_nominal, color='k',       lw=1.5,       label='Nominal  (MUR=MUF=1)')
step(ax_top, MLL_EDGES, h_up,      color='tomato',  lw=1.2, ls='--', label='Scale up   (envelope max)')
step(ax_top, MLL_EDGES, h_down,    color='tomato',  lw=1.2, ls=':',  label='Scale down (envelope min)')

ax_top.set_yscale('log')
ax_top.set_ylabel('Events  [a.u.]', fontsize=12)
ax_top.legend(fontsize=10)
ax_top.set_xlim(MLL_MIN, MLL_MAX)
ax_top.set_title('QCD scale uncertainty', fontsize=12)

# ----- Bottom panel --------------------------

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




# ------ Individual scale variations --------------------------


# Map weight ID to MUR, MUF
SCALE_LABELS = {
    '1':  (0.5, 0.5),
    '6':  (0.5, 1.0),
    '11': (0.5, 2.0),
    '16': (1.0, 0.5),
    '25': (1.0, 2.0),
    '30': (2.0, 0.5),
    '35': (2.0, 1.0),
    '40': (2.0, 2.0),
}

for i, sid in enumerate(scale_ids):
    mur, muf = SCALE_LABELS.get(sid, ('?', '?'))
    h_var = h_vars[i]
    delta  = (h_var - h_nominal) / safe_nominal

    fig_i, (at, ab) = plt.subplots(
        2, 1, figsize=(8, 7),
        sharex=True,
        gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.05}
    )

    # Top panel
    step(at, MLL_EDGES, h_nominal, color='k',      lw=1.5,       label='Nominal  (MUR=MUF=1)')
    step(at, MLL_EDGES, h_var,     color='tomato', lw=1.2, ls='--',
         label=f'ID {sid}  MUR={mur}  MUF={muf}')

    at.set_yscale('log')
    at.set_ylabel('Events  [a.u.]', fontsize=12)
    at.legend(fontsize=10)
    at.set_xlim(MLL_MIN, MLL_MAX)
    at.set_title(f'QCD scale variation - ID {sid}  |  MUR={mur}  MUF={muf}', fontsize=12)

    # Bottom panel
    step(ab, MLL_EDGES, delta, color='tomato', lw=1.5, ls='--',
         label=r'(var $-$ nom) / nom')
    ab.axhline(0.0, color='k', lw=1.0)
    ab.fill_between(BIN_CENTRES, 0, delta, step='mid', alpha=0.15, color='tomato')

    ab.set_xscale('log')
    ab.set_xlabel(r'$m_{\ell\ell}$  [GeV]', fontsize=12)
    ab.set_ylabel(r'$\delta$ / nominal', fontsize=11)
    ab.legend(fontsize=9)
    ab.grid(axis='y', alpha=0.3)

    fig_i.tight_layout()
    out = f"{PLOT_DIR}/plot_scale_var_ID{sid}_MUR{mur}_MUF{muf}.pdf"
    fig_i.savefig(out)
    plt.close(fig_i)
    print(f"Saved {out}")

# ---------- Merged canvas ---------------------


from matplotlib.gridspec import GridSpec

n_vars = len(scale_ids)   # 8
n_cols = 4
n_macro_rows = (n_vars + n_cols - 1) // n_cols   # 2

fig_all = plt.figure(figsize=(n_cols * 5, n_macro_rows * 6))
gs = GridSpec(
    n_macro_rows * 2, n_cols,
    figure=fig_all,
    height_ratios=[3, 1] * n_macro_rows,
    hspace=0.05,
    wspace=0.35,
)

for i, sid in enumerate(scale_ids):
    mur, muf = SCALE_LABELS.get(sid, ('?', '?'))
    h_var = h_vars[i]
    delta  = (h_var - h_nominal) / safe_nominal

    col       = i % n_cols
    macro_row = i // n_cols

    at = fig_all.add_subplot(gs[macro_row * 2,     col])
    ab = fig_all.add_subplot(gs[macro_row * 2 + 1, col], sharex=at)

    # Top panel
    step(at, MLL_EDGES, h_nominal, color='k',      lw=1.5,       label='Nominal')
    step(at, MLL_EDGES, h_var,     color='tomato', lw=1.2, ls='--',
         label=f'ID {sid}')

    at.set_yscale('log')
    at.set_ylabel('Events [a.u.]', fontsize=9)
    at.legend(fontsize=8)
    at.set_xlim(MLL_MIN, MLL_MAX)
    at.set_title(f'ID {sid}  MUR={mur}  MUF={muf}', fontsize=10)
    plt.setp(at.get_xticklabels(), visible=False)

    # Bottom panel
    step(ab, MLL_EDGES, delta, color='tomato', lw=1.5, ls='--')
    ab.axhline(0.0, color='k', lw=1.0)
    ab.fill_between(BIN_CENTRES, 0, delta, step='mid', alpha=0.15, color='tomato')

    ab.set_xscale('log')
    ab.set_xlabel(r'$m_{\ell\ell}$ [GeV]', fontsize=9)
    ab.set_ylabel(r'$\delta$/nom', fontsize=9)
    ab.grid(axis='y', alpha=0.3)

fig_all.suptitle('QCD scale variations', fontsize=13, y=1.01)
fig_all.tight_layout()
out_all = f"{PLOT_DIR}/plot_scale_all_variations.pdf"
fig_all.savefig(out_all, bbox_inches='tight')
plt.close(fig_all)
print(f"Saved {out_all}")

