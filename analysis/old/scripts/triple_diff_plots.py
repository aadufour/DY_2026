"""
rw_triple_diff.py

Triple-differential SMEFT decomposition in (m_ll, y_ll, cos_theta*).

Fills 3D boost_histogram, then unrolls to 1D.

Observables
-----------
  mll    : dilepton invariant mass [GeV]
  rap    : dilepton rapidity
  cstar  : cos(theta*) (angle of l- vs z in rest frame)

Weight decomposition (per operator)
------------------------------------
  SM     : w_SM
  Lin    : C    * 0.5*(w[op] - w[minus_op])
  Quad   : C^2  * (0.5*(w[op] + w[minus_op]) - w_SM)
  Full   : w_SM + C*w_lin + C^2*w_quad

Cache
-----
  Parsed LHE data is stored in lhe_cache.pkl (same directory) to avoid
  re-reading the LHE files on every run.  The cache is created either by
  running build_cache.py once, or automatically by this script on the first
  run.  Pass --rebuild-cache to force regeneration from the LHE files.

Usage
-----
  python3 rw_triple_diff.py                           # all 17 operators, C=1
  python3 rw_triple_diff.py --cHDD                    # only cHDD
  python3 rw_triple_diff.py --cHDD --changeC 2.0      # cHDD with C=2
  python3 rw_triple_diff.py --cHDD --plot1d           # also save 1D projections
"""

import os
import pickle
import warnings
import argparse
import numpy as np
import boost_histogram as bh
import mplhep as hep
import matplotlib.pyplot as plt
import pylhe

# ── Config ────────────────────────────────────────────────────────────────────
MLL_BIN_EDGES = [50, 120, 200, 400, 600, 800, 1000, 3000]
LHE_FILES = [
    f"/Users/albertodufour/MG5_2_9_18/mg5amcnlo/DY_all_{lo}_{hi}/myLHE/unweighted_events.lhe"
    for lo, hi in zip(MLL_BIN_EDGES[:-1], MLL_BIN_EDGES[1:])
]
OUT_DIR = "/Users/albertodufour/code/DY2026/analysis/triple_diff_plots"

OPERATORS = [
    'cHDD', 'cHWB', 'cHj1',
    'cHj3', 'cHu',  'cHd',
    'cHl1', 'cHl3', 'cHe',
    'cll1', 'clj1', 'clj3',
    'ceu',  'ced',  'cje',
    'clu',  'cld',
]

# ── Binning ───────────────────────────────────────────────────────────────────
# mll: defining custom bin edges variable-width edges capturing the Z peak and high-mass tail (from GB)
# MLL_EDGES   = np.array([50, 64, 76, 82, 86, 90, 98, 103, 121, 127, 130, 133, 148, 151, 154, 157, 163, 166, 172, 178, 184, 205, 210, 220, 235, 240, 260, 265, 325, 345, 500, 530, 570, 618, 654, 708, 3000], dtype=float)
MLL_EDGES   = np.array([50, 70, 90, 110, 200, 800, 1400, 2000, 2400, 3000], dtype=float)
# MLL_EDGES   = np.array([50, 90,  200, 3000], dtype=float)
RAP_EDGES   = np.array([0.0, 0.5, 1.0, 2.5], dtype=float)
CSTAR_EDGES = np.array([-1, -0.5, 0.0, 0.5, 1], dtype=float) 

N_MLL   = len(MLL_EDGES)   - 1
N_RAP   = len(RAP_EDGES)   - 1 
N_CSTAR = len(CSTAR_EDGES) - 1  
N_TOT   = N_MLL * N_RAP * N_CSTAR 


# ── Kinematic observables ─────────────────────────────────────────────────────
# 4-vector convention: [px, py, pz, E]

def mll(p1, p2):
    p = np.array(p1) + np.array(p2)
    return np.sqrt(max(p[3]**2 - sum(p[i]**2 for i in range(3)), 0.0))


def rap(p1, p2):
    p  = np.array(p1) + np.array(p2)
    E  = p[3]
    pz = p[2]
    return np.abs(0.5 * np.log((E + pz) / (E - pz)))


def cstar(p1, p2):
    """cos(theta*): angle of e- vs z-axis in dilepton rest frame."""
    p1 = np.array(p1)
    p2 = np.array(p2)
    p  = p1 + p2
    E  = p[3]
    pz = p[2]
    mass  = mll(p1, p2)
    beta  = pz / E
    gamma = E / mass
    # boost e- along z into rest frame
    pz1_boosted = gamma * (p1[2] - beta * p1[3])
    p1_mag      = np.sqrt(p1[0]**2 + p1[1]**2 + pz1_boosted**2)
    return pz1_boosted / p1_mag


# ── Unrolling ─────────────────────────────────────────────────────────────────
def unroll_3d(h3d):
    """
    Unroll a 3D bh.Histogram (with Weight storage) to 1D.

    For shape (n_mll, n_rap, n_cstar), the 1D bin index is:
        i_cstar * (n_rap * n_mll) + i_rap * n_mll + i_mll
    """
    v   = h3d.view()                    # structured array, shape (nx, ny, nz)
    n   = v['value'].size
    ax  = bh.axis.Regular(n, 0, n)
    h1d = bh.Histogram(ax, storage=bh.storage.Weight())
    h1d.view()['value']    = v['value'].T.flatten()
    h1d.view()['variance'] = v['variance'].T.flatten()
    return h1d


def make_3d():
    return bh.Histogram(
        bh.axis.Variable(MLL_EDGES),
        bh.axis.Variable(RAP_EDGES),
        bh.axis.Variable(CSTAR_EDGES),
        storage=bh.storage.Weight(),
    )


# ── Argument parsing ──────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
for op in OPERATORS:
    parser.add_argument(f'--{op}', action='store_true', default=False)
parser.add_argument('--changeC', type=float, default=1.0,
                    help='Wilson coefficient value C (default: 1.0)')
parser.add_argument('--plot1d', action='store_true', default=False,
                    help='Also save 1D projection plots (mll, rap, cstar)')
parser.add_argument('--plot-components', action='store_true', default=False,
                    help='Also save individual PDF canvases for Lin, Quad, Full (no log scale)')
parser.add_argument('--out-dir', type=str, default=None,
                    help='Output directory for plots (default: analysis/triple_diff_plots)')
parser.add_argument('--rebuild-cache', action='store_true', default=False,
                    help='Force re-reading LHE files even if cache exists')
args   = parser.parse_args()
C      = args.changeC
requested   = [op for op in OPERATORS if getattr(args, op)]
ops_to_plot = requested if requested else OPERATORS
print(f"Operators : {ops_to_plot}")
print(f"C         : {C}")
print(f"Unrolled bins: {N_TOT} = ({N_MLL} mll * {N_RAP} rap * {N_CSTAR} cstar)")


# ── Load or build cache ───────────────────────────────────────────────────────
# The cache (lhe_cache.pkl) stores pre-parsed LHE event arrays so that
# subsequent runs skip the slow pylhe parsing step.
#
# The cache can be produced in two ways:
#   1. Standalone:  run build_cache.py once (recommended for first-time setup
#                   or when the LHE files change).
#   2. On-the-fly:  if lhe_cache.pkl is absent, this script reads the LHE
#                   files directly and writes the cache itself before
#                   proceeding to the histogram filling.
#
# To force a rebuild from LHE files, pass --rebuild-cache on the command line.
CACHE_FILE = os.path.join(os.path.dirname(__file__), "/Users/albertodufour/code/DY2026/analysis/lhe_cache.pkl")

if os.path.exists(CACHE_FILE) and not args.rebuild_cache:
    print(f"Loading cache from {CACHE_FILE}")
    with open(CACHE_FILE, 'rb') as f:
        cache = pickle.load(f)
    mll_arr   = cache['mll']
    rap_arr   = cache['rap']
    cstar_arr = cache['cstar']
    wSM_arr   = cache['w_SM']
    w_p1_all  = cache['w_p1']
    w_m1_all  = cache['w_m1']
    print(f"  {len(mll_arr)} events loaded\n")
else:
    if args.rebuild_cache:
        print("--rebuild-cache set, re-reading LHE files.")
    else:
        print(f"No cache found at {CACHE_FILE}. Reading LHE files.")
        print("  Tip: run build_cache.py once to avoid this on future runs.\n")

    mll_vals   = []
    rap_vals   = []
    cstar_vals = []
    w_SM       = []
    w_p1_all   = {op: [] for op in OPERATORS}
    w_m1_all   = {op: [] for op in OPERATORS}

    for lhe_file in LHE_FILES:
        print(f"\nReading {lhe_file}")
        if not os.path.exists(lhe_file):
            print(f"  WARNING: file not found, skipping.")
            continue

        n_before = len(mll_vals)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            events = pylhe.read_lhe_with_attributes(lhe_file)

            for i, event in enumerate(events):
                if (i + 1) % 5000 == 0:
                    print(f"  processed {i + 1} events")

                leptons = [
                    p for p in event.particles
                    if int(p.status) == 1 and abs(int(p.id)) in {11, 13}
                ]
                if len(leptons) < 2:
                    continue

                lm = next((p for p in leptons if int(p.id) > 0), leptons[0])
                lp = next((p for p in leptons if int(p.id) < 0), leptons[1])
                v_lm = [lm.px, lm.py, lm.pz, lm.e]
                v_lp = [lp.px, lp.py, lp.pz, lp.e]

                m  = mll(v_lm, v_lp)
                y  = rap(v_lm, v_lp)
                cs = cstar(v_lm, v_lp)

                if not (MLL_EDGES[0] <= m <= MLL_EDGES[-1]):
                    continue
                if not (RAP_EDGES[0] <= y <= RAP_EDGES[-1]):
                    continue
                if not (-1.0 <= cs <= 1.0):
                    continue

                mll_vals.append(m)
                rap_vals.append(y)
                cstar_vals.append(cs)
                w_SM.append(event.weights['SM'])
                for op in OPERATORS:
                    w_p1_all[op].append(event.weights[op])
                    w_m1_all[op].append(event.weights[f'minus{op}'])

        print(f"  {len(mll_vals) - n_before} events passed from this file")

    mll_arr   = np.array(mll_vals)
    rap_arr   = np.array(rap_vals)
    cstar_arr = np.array(cstar_vals)
    wSM_arr   = np.array(w_SM)
    w_p1_all  = {op: np.array(v) for op, v in w_p1_all.items()}
    w_m1_all  = {op: np.array(v) for op, v in w_m1_all.items()}

    print(f"\nTotal — {len(mll_arr)} events passed kinematic cuts")
    print(f"Saving cache to {CACHE_FILE}")
    with open(CACHE_FILE, 'wb') as f:
        pickle.dump({
            'mll': mll_arr, 'rap': rap_arr, 'cstar': cstar_arr,
            'w_SM': wSM_arr, 'w_p1': w_p1_all, 'w_m1': w_m1_all,
        }, f)
    print(f"  Cache saved ({os.path.getsize(CACHE_FILE)/1e6:.1f} MB)\n")

# SM 3D histogram (shared)
h3d_SM = make_3d()
h3d_SM.fill(mll_arr, rap_arr, cstar_arr, weight=wSM_arr)
h1d_SM = unroll_3d(h3d_SM)

# ── 1D projection ──────────────────────────────────────────────────────
def plot_1d_projections(h3d_SM, h3d_lin, h3d_quad, h3d_full, op, C, out_dir):
    """
    Project the 3D histograms onto each axis independently (summing over the
    other two) and save a 3-panel figure: mll | rap | cstar.
    bh slicing: sum = integrate out that axis.
    """
    # axis 0 = mll, axis 1 = rap, axis 2 = cstar
    proj = {
        'mll':   (0, r"$m_{\ell\ell}$ [GeV]"),
        'rap':   (1, r"$y_{\ell\ell}$"),
        'cstar': (2, r"$\cos\theta^*$"),
    }

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.subplots_adjust(wspace=0.35)

    for ax, (name, (axis_idx, xlabel)) in zip(axes, proj.items()):
        # sum over the other two axes
        other = tuple(i for i in range(3) if i != axis_idx)
        p_SM   = h3d_SM.project(axis_idx)
        p_lin  = h3d_lin.project(axis_idx)
        p_quad = h3d_quad.project(axis_idx)
        p_full = h3d_full.project(axis_idx)

        hep.histplot(p_SM,   ax=ax, label="SM",   color="black",       linewidth=1.5, histtype="step")
        hep.histplot(p_lin,  ax=ax, label="Lin",  color="steelblue",   linewidth=1.5, histtype="step", linestyle="--")
        hep.histplot(p_quad, ax=ax, label="Quad", color="firebrick",   linewidth=1.5, histtype="step", linestyle=":")
        hep.histplot(p_full, ax=ax, label="Full", color="forestgreen", linewidth=1.5, histtype="step", linestyle="-.")

        ax.set_xlabel(xlabel)
        # ax.set_ylabel("Weighted events")
        ax.legend(frameon=False, loc="upper right", fontsize=8)
        ax.semilogy()

    hep.style.use("CMS")
    hep.cms.label(llabel="", rlabel=f"{op}={C}", ax=axes[0])
    # hep.cms.label(f"", ax=axes[0], loc=0, data=True, lumi=None)
    # hep.cms.label(f"{op}={C}", ax=axes[1], loc=0, data=True, lumi=None)
    
    # hep.cms.label(data=True, lumi=None)

    out_path = os.path.join(out_dir, f"1d_proj_{op}_C{C}.pdf")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out_path}")


# ── Individual component plots (no log scale) ─────────────────────────────────
def plot_component_canvases(h1d_SM, h1d_lin, h1d_quad, h1d_full, op, C, out_dir):
    """
    Save separate PDF canvases for Lin, Quad, and Full (no log scale).
    Each figure shows the component alongside SM for reference.
    """
    components = [
        ("Lin",  h1d_lin,  "steelblue",   "--"),
        ("Quad", h1d_quad, "firebrick",   ":"),
        ("Full", h1d_full, "forestgreen", "-."),
    ]
    for name, h_comp, color, ls in components:
        fig, ax = plt.subplots(figsize=(12, 5))

        hep.histplot(h_comp, ax=ax, label=name, color=color, linewidth=1.5, histtype="step", linestyle=ls)
        ax.axhline(0.0, color="black", linewidth=0.8, linestyle="-")

        # slice separators
        for i_slice in range(1, N_RAP * N_CSTAR):
            ax.axvline(i_slice * N_MLL, color="gray", linewidth=0.4, linestyle="--", alpha=0.4)

        # slice labels
        for i_cs, cs_lo, cs_hi in zip(range(N_CSTAR), CSTAR_EDGES[:-1], CSTAR_EDGES[1:]):
            for i_rap, rap_lo, rap_hi in zip(range(N_RAP), RAP_EDGES[:-1], RAP_EDGES[1:]):
                i_slice  = i_cs * N_RAP + i_rap
                x_center = (i_slice + 0.5) * N_MLL
                label    = f"y[{rap_lo:.1f},{rap_hi:.1f}]\ncs[{cs_lo:.1f},{cs_hi:.1f}]"
                ax.text(x_center, 1.01, label, transform=ax.get_xaxis_transform(),
                        ha='center', va='bottom', fontsize=4, color='gray')

        ax.set_xlabel("Unrolled bins")
        ax.set_ylabel("Weighted events")
        ax.legend(frameon=False, loc="upper right", fontsize=9)
        hep.cms.label(f"{op}={C}  [{name}]", ax=ax, loc=0, data=True, lumi=None)

        out_path = os.path.join(out_dir, f"component_{name.lower()}_{op}_C{C}.pdf")
        fig.savefig(out_path, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved {out_path}")


# ── Per-operator: fill, unroll, plot ──────────────────────────────────────────
if args.out_dir:
    OUT_DIR = args.out_dir
os.makedirs(OUT_DIR, exist_ok=True)
hep.style.use("CMS")
hep.cms.label(data=True, lumi=None)

for op in ops_to_plot:
    wp1 = w_p1_all[op]
    wm1 = w_m1_all[op]

    w_lin  = 0.5 * (wp1 - wm1)
    w_quad = 0.5 * (wp1 + wm1) - wSM_arr
    w_full = wSM_arr + C * w_lin + C**2 * w_quad

    # sanity check
    print(f"  {op}: "
          f"lin={w_lin.sum():.3e}  "
          f"quad={w_quad.sum():.3e}  "
          f"full={w_full.sum():.3e}")
    if w_full.sum() < 0:
        print(f"    !!! full cross section is negative : C={C} may be outside EFT validity")

    # fill 3D
    h3d_lin  = make_3d()
    h3d_quad = make_3d()
    h3d_full = make_3d()
    h3d_lin.fill( mll_arr, rap_arr, cstar_arr, weight=C      * w_lin)
    h3d_quad.fill(mll_arr, rap_arr, cstar_arr, weight=C**2   * w_quad)
    h3d_full.fill(mll_arr, rap_arr, cstar_arr, weight=w_full)

    # unroll to 1D
    h1d_lin  = unroll_3d(h3d_lin)
    h1d_quad = unroll_3d(h3d_quad)
    h1d_full = unroll_3d(h3d_full)

    # ── Plot unrolled ──────────────────────────────────────────────────────────
    fig, (ax, ax_ratio) = plt.subplots(
        2, 1, figsize=(12, 7),
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
        sharex=True,
    )

    hep.histplot(h1d_SM,   ax=ax, label="SM",   color="black",       linewidth=1.5, histtype="step")
    hep.histplot(h1d_lin,  ax=ax, label="Lin",  color="steelblue",   linewidth=1.5, histtype="step", linestyle="--")
    hep.histplot(h1d_quad, ax=ax, label="Quad", color="firebrick",   linewidth=1.5, histtype="step", linestyle=":")
    hep.histplot(h1d_full, ax=ax, label="Full", color="forestgreen", linewidth=1.5, histtype="step", linestyle="-.")

    # ── Ratio panel: Full / SM ─────────────────────────────────────────────────
    sm_vals   = h1d_SM.values()
    full_vals = h1d_full.values()
    sm_var    = h1d_SM.variances()
    full_var  = h1d_full.variances()

    # avoid division by zero
    safe      = sm_vals != 0
    ratio     = np.where(safe, full_vals / np.where(safe, sm_vals, 1), np.nan)
    # error propagation: sigma_ratio = ratio * sqrt(sigma_full**2/full**2 + sigma_SM**2/SM**2)
    ratio_err = np.where(
        safe & (full_vals != 0),
        ratio * np.sqrt(
            full_var / np.where(full_vals != 0, full_vals**2, 1) +
            sm_var   / np.where(safe,            sm_vals**2,   1)
        ),
        np.nan,
    )

    bin_centers = h1d_SM.axes[0].centers
    ax_ratio.errorbar(
        bin_centers, ratio, yerr=ratio_err,
        fmt='o', color="forestgreen", markersize=2, linewidth=0.8, label="Full / SM",
    )
    ax_ratio.axhline(1.0, color="black", linewidth=0.8, linestyle="--")
    ax_ratio.set_ylabel("Full / SM")
    ax_ratio.set_ylim(0, 2)

    # slice separators on both panels
    for i_slice in range(1, N_RAP * N_CSTAR):
        for a in (ax, ax_ratio):
            a.axvline(i_slice * N_MLL, color="gray", linewidth=0.4, linestyle="--", alpha=0.4)

    # label each (rap, cstar) slice at the top of main panel
    for i_cs, cs_lo, cs_hi in zip(range(N_CSTAR), CSTAR_EDGES[:-1], CSTAR_EDGES[1:]):
        for i_rap, rap_lo, rap_hi in zip(range(N_RAP), RAP_EDGES[:-1], RAP_EDGES[1:]):
            i_slice  = i_cs * N_RAP + i_rap
            x_center = (i_slice + 0.5) * N_MLL
            label    = f"y[{rap_lo:.1f},{rap_hi:.1f}]\ncs[{cs_lo:.1f},{cs_hi:.1f}]"
            ax.text(x_center, 1.01, label, transform=ax.get_xaxis_transform(),
                    ha='center', va='bottom', fontsize=4, color='gray')

    ax.set_ylabel("Weighted events")
    ax.legend(frameon=False, loc="upper right", fontsize=9)
    ax.semilogy()
    # ax_ratio.set_xlabel(f"Unrolled bin  — mll fastest ({N_MLL}), rap×cstar ({N_RAP}×{N_CSTAR}) — {N_TOT} total")
    ax_ratio.set_xlabel(f"Unrolled bins")
    hep.cms.label(f"{op}={C}", ax=ax, loc=0, data=True, lumi=None)

    out_path = os.path.join(OUT_DIR, f"triple_diff_{op}_C{C}.pdf")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out_path}")

    if args.plot1d:
        plot_1d_projections(h3d_SM, h3d_lin, h3d_quad, h3d_full, op, C, OUT_DIR)

    if args.plot_components:
        plot_component_canvases(h1d_SM, h1d_lin, h1d_quad, h1d_full, op, C, OUT_DIR)

print("\nAll done.")
