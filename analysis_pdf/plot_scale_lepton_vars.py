#!/usr/bin/env python3
"""
plot_scale_lepton_vars.py

QCD scale uncertainty using the 7-point envelope prescription,
for single-lepton kinematic variables:
  ptl1  -- leading-lepton transverse momentum [GeV]
  etal1 -- leading-lepton pseudorapidity
  phil1 -- leading-lepton azimuthal angle [rad]

"Leading lepton" = the one with highest pT in the event.

These quantities must be present in cache_pdf.pkl.
If they are missing, rebuild the cache with build_cache_pdf.py
(which now stores them).

For each variable the script produces:
  1. One envelope plot  (nominal + UP envelope + DOWN envelope, top+ratio panels)
  2. Eight individual scale-variation plots  (one per MUR/MUF point)
  3. One merged canvas with all eight individual plots

Weight IDs used (same as plot_scale_updown.py):
  ID  1 : MUR=0.5  MUF=0.5
  ID  6 : MUR=0.5  MUF=1.0
  ID 11 : MUR=0.5  MUF=2.0
  ID 16 : MUR=1.0  MUF=0.5
  ID 25 : MUR=1.0  MUF=2.0
  ID 30 : MUR=2.0  MUF=0.5
  ID 35 : MUR=2.0  MUF=1.0
  ID 40 : MUR=2.0  MUF=2.0
"""

import os
import pickle

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

# ── Config ────────────────────────────────────────────────────────────────────

CACHE_FILE = os.path.join(os.path.dirname(__file__), "cache_pdf.pkl")
PLOT_DIR   = os.path.join(os.path.dirname(__file__), "plots_pdf_study")
os.makedirs(PLOT_DIR, exist_ok=True)

# ── Variable definitions ──────────────────────────────────────────────────────
#
# Each entry:
#   key    -- key in the cache dict
#   edges  -- bin edges (numpy array)
#   xlabel -- x-axis label (LaTeX, all ASCII source)
#   title  -- base title string for plot titles
#   xscale -- 'log' or 'linear'
#   xlim   -- (xmin, xmax) passed to set_xlim (None = use edges extremes)

VARIABLES = [
    {
        'key':    'ptl1',
        'edges':  np.logspace(np.log10(20.0), np.log10(1500.0), 40),
        'xlabel': r'$p_{\rm T}^{\ell_1}$  [GeV]',
        'title':  'Leading-lepton $p_{\\rm T}$',
        'xscale': 'log',
        'xlim':   None,
    },
    {
        'key':    'etal1',
        'edges':  np.linspace(-5.0, 5.0, 40),
        'xlabel': r'$\eta^{\ell_1}$',
        'title':  'Leading-lepton pseudorapidity',
        'xscale': 'linear',
        'xlim':   (-5.0, 5.0),
    },
    {
        'key':    'phil1',
        'edges':  np.linspace(-np.pi, np.pi, 40),
        'xlabel': r'$\phi^{\ell_1}$  [rad]',
        'title':  'Leading-lepton azimuthal angle',
        'xscale': 'linear',
        'xlim':   (-np.pi, np.pi),
    },
]

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

# ── Load cache ────────────────────────────────────────────────────────────────

print(f"Loading {CACHE_FILE} ...")
with open(CACHE_FILE, 'rb') as f:
    cache = pickle.load(f)

# Verify required keys are present
required = ['w_central', 'w_scale', 'scale_ids'] + [v['key'] for v in VARIABLES]
missing  = [k for k in required if k not in cache]
if missing:
    raise KeyError(
        f"Cache is missing keys: {missing}\n"
        "Rebuild cache_pdf.pkl with the updated build_cache_pdf.py."
    )

w_central = cache['w_central']   # (N,)
w_scale   = cache['w_scale']     # (8, N)
scale_ids = cache['scale_ids']   # ['1','6','11','16','25','30','35','40']
print(f"Events loaded: {len(w_central):,}\n")

# ── Helpers ───────────────────────────────────────────────────────────────────

def fill(obs_arr, edges, weights):
    """Weighted histogram of obs_arr into edges."""
    vals, _ = np.histogram(obs_arr, bins=edges, weights=weights)
    return vals.astype(float)

def step(ax, edges, vals, **kw):
    """Draw a step histogram from pre-computed bin values."""
    ax.step(
        np.append(edges[:-1], edges[-1]),
        np.append(vals, vals[-1]),
        where='post', **kw
    )

# ── Loop over variables ───────────────────────────────────────────────────────

for var in VARIABLES:
    key    = var['key']
    edges  = var['edges']
    xlabel = var['xlabel']
    title  = var['title']
    xscale = var['xscale']
    xlim   = var['xlim'] if var['xlim'] is not None else (edges[0], edges[-1])

    obs_arr     = cache[key]                              # (N,)
    bin_centres = 0.5 * (edges[:-1] + edges[1:])

    h_nominal = fill(obs_arr, edges, w_central)
    h_vars    = np.array([fill(obs_arr, edges, w_scale[i])
                          for i in range(len(scale_ids))])  # (8, N_bins)

    # Envelope
    h_up   = h_vars.max(axis=0)
    h_down = h_vars.min(axis=0)

    safe_nom   = np.where(h_nominal > 0, h_nominal, np.nan)
    delta_up   = (h_up   - h_nominal) / safe_nom
    delta_down = (h_down - h_nominal) / safe_nom

    # ── 1. Envelope plot ──────────────────────────────────────────────────────

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(8, 7),
        sharex=True,
        gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.05}
    )

    step(ax_top, edges, h_nominal, color='k',      lw=1.5,
         label='Nominal  (MUR=MUF=1)')
    step(ax_top, edges, h_up,      color='tomato', lw=1.2, ls='--',
         label='Scale up   (envelope max)')
    step(ax_top, edges, h_down,    color='tomato', lw=1.2, ls=':',
         label='Scale down (envelope min)')

    ax_top.set_yscale('log')
    ax_top.set_ylabel('Events  [a.u.]', fontsize=12)
    ax_top.legend(fontsize=10)
    ax_top.set_xlim(*xlim)
    # ax_top.set_title(f'QCD scale uncertainty  --  {title}', fontsize=12)

    step(ax_bot, edges, delta_up,   color='tomato', lw=1.5, ls='--',
         label='(up $-$ nom) / nom')
    step(ax_bot, edges, delta_down, color='tomato', lw=1.5, ls=':',
         label='(down $-$ nom) / nom')
    ax_bot.axhline(0.0, color='k', lw=1.0)
    ax_bot.fill_between(bin_centres, delta_down, delta_up,
                        step='mid', alpha=0.15, color='tomato')

    ax_bot.set_xscale(xscale)
    ax_bot.set_xlabel(xlabel, fontsize=12)
    ax_bot.set_ylabel(r'$\delta$ / nominal', fontsize=11)
    ax_bot.legend(fontsize=9)
    ax_bot.grid(axis='y', alpha=0.3)

    fig.tight_layout()
    out = f"{PLOT_DIR}/plot_scale_updown_{key}.pdf"
    fig.savefig(out)
    plt.close(fig)
    print(f"Saved {out}")

    # ── 2. Individual scale-variation plots ───────────────────────────────────

    for i, sid in enumerate(scale_ids):
        mur, muf = SCALE_LABELS.get(sid, ('?', '?'))
        h_var  = h_vars[i]
        delta  = (h_var - h_nominal) / safe_nom

        fig_i, (at, ab) = plt.subplots(
            2, 1, figsize=(8, 7),
            sharex=True,
            gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.05}
        )

        step(at, edges, h_nominal, color='k',      lw=1.5,
             label='Nominal  (MUR=MUF=1)')
        step(at, edges, h_var,     color='tomato', lw=1.2, ls='--',
             label=f'ID {sid}  MUR={mur}  MUF={muf}')

        at.set_yscale('log')
        at.set_ylabel('Events  [a.u.]', fontsize=12)
        at.legend(fontsize=10)
        at.set_xlim(*xlim)
        at.set_title(
            f'QCD scale var. - {title}  |  ID {sid}  MUR={mur}  MUF={muf}',
            fontsize=11
        )

        step(ab, edges, delta, color='tomato', lw=1.5, ls='--',
             label=r'(var $-$ nom) / nom')
        ab.axhline(0.0, color='k', lw=1.0)
        ab.fill_between(bin_centres, 0, delta, step='mid', alpha=0.15,
                        color='tomato')

        ab.set_xscale(xscale)
        ab.set_xlabel(xlabel, fontsize=12)
        ab.set_ylabel(r'$\delta$ / nominal', fontsize=11)
        ab.legend(fontsize=9)
        ab.grid(axis='y', alpha=0.3)

        fig_i.tight_layout()
        out_i = f"{PLOT_DIR}/plot_scale_var_{key}_ID{sid}_MUR{mur}_MUF{muf}.pdf"
        fig_i.savefig(out_i)
        plt.close(fig_i)
        print(f"Saved {out_i}")

    # ── 3. Merged canvas: all 8 individual variations ─────────────────────────

    n_scale  = len(scale_ids)   # 8
    n_cols   = 4
    n_mrows  = (n_scale + n_cols - 1) // n_cols   # 2

    fig_all = plt.figure(figsize=(n_cols * 5, n_mrows * 6))
    gs = GridSpec(
        n_mrows * 2, n_cols,
        figure=fig_all,
        height_ratios=[3, 1] * n_mrows,
        hspace=0.05,
        wspace=0.35,
    )

    for i, sid in enumerate(scale_ids):
        mur, muf = SCALE_LABELS.get(sid, ('?', '?'))
        h_var  = h_vars[i]
        delta  = (h_var - h_nominal) / safe_nom

        col      = i % n_cols
        mrow     = i // n_cols

        at = fig_all.add_subplot(gs[mrow * 2,     col])
        ab = fig_all.add_subplot(gs[mrow * 2 + 1, col], sharex=at)

        step(at, edges, h_nominal, color='k',      lw=1.5, label='Nominal')
        step(at, edges, h_var,     color='tomato', lw=1.2, ls='--',
             label=f'ID {sid}')

        at.set_yscale('log')
        at.set_ylabel('Events [a.u.]', fontsize=9)
        at.legend(fontsize=8)
        at.set_xlim(*xlim)
        at.set_title(f'ID {sid}  MUR={mur}  MUF={muf}', fontsize=10)
        plt.setp(at.get_xticklabels(), visible=False)

        step(ab, edges, delta, color='tomato', lw=1.5, ls='--')
        ab.axhline(0.0, color='k', lw=1.0)
        ab.fill_between(bin_centres, 0, delta, step='mid', alpha=0.15,
                        color='tomato')

        ab.set_xscale(xscale)
        ab.set_xlabel(xlabel, fontsize=9)
        ab.set_ylabel(r'$\delta$/nom', fontsize=9)
        ab.grid(axis='y', alpha=0.3)

    # fig_all.suptitle(
    #     f'QCD scale variations -- all 8 points  ({title})',
    #     fontsize=13, y=1.01
    # )
    fig_all.tight_layout()
    out_all = f"{PLOT_DIR}/plot_scale_all_variations_{key}.pdf"
    fig_all.savefig(out_all, bbox_inches='tight')
    plt.close(fig_all)
    print(f"Saved {out_all}")
