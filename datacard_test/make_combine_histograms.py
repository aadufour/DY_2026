"""
make_combine_histograms.py

Produce histograms.root for CMS combine shape analysis from the LHE weight
cache built by rw_build_cache.py.

------------------------------------------------------------

Output ROOT structure:
    histograms.root
    triple_DY/
        data_obs       (= sm)
        sm             (pure SM component)
        quad_<op>      (C^2 * w_quad -- pure quadratic EFT)
        sm_lin_quad_<op>  (full prediction at C_ref)

Weight decomposition (matches rw_triple_diff.py):
    w_lin       = 0.5 * (w[+op] - w[-op])
    w_quad      = 0.5 * (w[+op] + w[-op]) - w_SM
    sm_lin_quad = w_SM + C*w_lin + C^2*w_quad

Usage:
    python3 make_combine_histograms.py
    python3 make_combine_histograms.py --op cHDD
    python3 make_combine_histograms.py --op cHDD --C 0.5 --lumi 59740
    python3 make_combine_histograms.py --op cHDD --unrolled   # 3D unrolled
    python3 make_combine_histograms.py --op cHDD --no-plot    # skip PNG
"""

import argparse
import os
import pickle

import boost_histogram as bh
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import uproot

# ---- Config --------------------------------------------------

CACHE_FILE  = "/Users/albertodufour/code/DY2026/triple_diff_test/lhe_cache.pkl"
OUTPUT_FILE = "histograms.root"
CHANNEL     = "triple_DY"

# mll binning  (matches rw_triple_diff.py)
MLL_EDGES   = np.array([50, 70, 90, 110, 200, 800, 1400, 2000, 2400, 3000], dtype=float)

# 3D binning   (only used with --unrolled)
RAP_EDGES   = np.array([0.0, 0.5, 1.0, 2.5], dtype=float)
CSTAR_EDGES = np.array([-1.0, -0.5, 0.0, 0.5, 1.0], dtype=float)


# ---- Args --------------------------------------------------------

parser = argparse.ArgumentParser(description="Build combine input ROOT file from LHE cache.")
parser.add_argument("--op",      default="cHDD",
                    help="Operator name (must be in cache, e.g. cHDD)")
parser.add_argument("--C",       type=float, default=1.0,
                    help="Reference Wilson coefficient value for sm_lin_quad template (default: 1.0)")
parser.add_argument("--lumi",    type=float, default=1.0,
                    help="Luminosity in pb^-1 -- scales all yields (default: 1.0 = raw weights in pb)")
parser.add_argument("--unrolled", action="store_true",
                    help="Use unrolled 3D (mll x rap x cstar) histogram instead of 1D mll")
parser.add_argument("--cache",   default=CACHE_FILE, help="Path to lhe_cache.pkl")
parser.add_argument("--output",  default=OUTPUT_FILE, help="Output ROOT file path")
parser.add_argument("--no-plot", action="store_true",  help="Skip saving bh PNG")
args = parser.parse_args()

OP   = args.op
C    = args.C
LUMI = args.lumi

# Get plot output path from ROOT file path
plot_path = os.path.splitext(args.output)[0] + f"_{OP}.png"

print(f"Operator   : {OP}")
print(f"C_ref      : {C}")
print(f"Luminosity : {LUMI} pb^-1")
print(f"Observable : {'unrolled 3D' if args.unrolled else 'mll 1D'}")
print(f"Cache      : {args.cache}")
print(f"Output     : {args.output}")
print()

# -- Load cache -----------------------------------------------

if not os.path.exists(args.cache):
    raise FileNotFoundError(
        f"Cache not found: {args.cache}\n"
        "Run rw_build_cache.py first to generate it."
    )

print("Loading cache ...")
with open(args.cache, "rb") as f:
    cache = pickle.load(f)

mll_arr   = cache["mll"]
rap_arr   = cache["rap"]
cstar_arr = cache["cstar"]
w_SM      = cache["w_SM"]
w_p1_all  = cache["w_p1"]   # weights with operator set to +1
w_m1_all  = cache["w_m1"]   # weights with operator set to -1
print(f"  {len(mll_arr):,} events loaded\n")

if OP not in w_p1_all:
    available = sorted(w_p1_all.keys())
    raise KeyError(
        f"Operator '{OP}' not found in cache.\n"
        f"Available: {available}\n"
        f"Add '{OP}' to rw_build_cache.py's OPERATORS list and rebuild."
    )

# ---- Weight decomposition ------------------------------------------

wp1 = w_p1_all[OP]
wm1 = w_m1_all[OP]

w_lin         = 0.5 * (wp1 - wm1)
w_quad        = 0.5 * (wp1 + wm1) - w_SM
w_slq         = w_SM + C * w_lin + C**2 * w_quad   # sm_lin_quad at C_ref
w_quad_scaled = C**2 * w_quad

if w_slq.sum() < 0:
    print(f"WARNING: full cross section is negative at C={C} -- C may be outside EFT validity.\n")

# ── Histogram builders ─────────────────────────────────────────────────────────
# axis metadata : ROOT axis title  (uproot reads bh axis.metadata as x-title)
# hist metadata  : ROOT hist title (shown in TBrowser title bar)

def _x_axis(label):
    """Return a Variable axis with the given ROOT-style title."""
    return bh.axis.Variable(MLL_EDGES, metadata=label)


def _unrolled_axis():
    n = (len(MLL_EDGES)-1) * (len(RAP_EDGES)-1) * (len(CSTAR_EDGES)-1)
    return bh.axis.Regular(n, 0, n, metadata="Unrolled bin (m_{ll} #times y_{ll} #times cos#theta*)")


def make_mll_1d(weights, proc_label):
    """1D m_ll histogram with axis / histogram titles set for TBrowser."""
    h = bh.Histogram(
        _x_axis("m_{ll} [GeV]"),
        storage=bh.storage.Weight(),
    )
    h.fill(mll_arr, weight=weights * LUMI)
    h.metadata = f"{proc_label} -- {CHANNEL} -- {OP}  C={C}"
    return h


def _make_3d_filled(weights):
    h = bh.Histogram(
        bh.axis.Variable(MLL_EDGES),
        bh.axis.Variable(RAP_EDGES),
        bh.axis.Variable(CSTAR_EDGES),
        storage=bh.storage.Weight(),
    )
    h.fill(mll_arr, rap_arr, cstar_arr, weight=weights * LUMI)
    return h


def unroll_3d(h3d, proc_label):
    v   = h3d.view()
    n   = v["value"].size
    h1d = bh.Histogram(_unrolled_axis(), storage=bh.storage.Weight())
    h1d.view()["value"]    = v["value"].T.flatten()
    h1d.view()["variance"] = v["variance"].T.flatten()
    h1d.metadata = f"{proc_label} -- {CHANNEL} -- {OP}  C={C}"
    return h1d


def make_hist(weights, proc_label):
    if args.unrolled:
        return unroll_3d(_make_3d_filled(weights), proc_label)
    return make_mll_1d(weights, proc_label)

# ---- Fill ---------------------------------------------------------------

process_weights = {
    "sm":                (w_SM,          "SM"),
    f"quad_{OP}":        (w_quad_scaled, f"quad  {OP}"),
    f"sm_lin_quad_{OP}": (w_slq,         f"SM+lin+quad  {OP}"),
}

histograms = {name: make_hist(w, lbl) for name, (w, lbl) in process_weights.items()}

# data_obs: clone of sm
h_data = make_hist(w_SM, "data_obs (= SM)")
histograms["data_obs"] = h_data

# ---- Print summary ------------------------------------------------

print(f"{'Process':<32}  {'Integral':>14}  {'sqrt(SumW2)':>14}")
print("-" * 64)
for name, h in histograms.items():
    intg = h.values().sum()
    unc  = np.sqrt(h.variances().sum())
    print(f"  {name:<30}  {intg:>14.4e}  {unc:>14.4e}")

n_bins = histograms["sm"].values().size
print(f"\nHistogram bins : {n_bins}")
if args.unrolled:
    n_m = len(MLL_EDGES)   - 1
    n_r = len(RAP_EDGES)   - 1
    n_c = len(CSTAR_EDGES) - 1
    print(f"  = {n_m} mll x {n_r} rap x {n_c} cstar")
print()

# ---- Write ROOT file ------------------------------------------------

print(f"Writing {args.output} ...")
with uproot.recreate(args.output) as rf:
    for name, h in histograms.items():
        key = f"{CHANNEL}/{name}"
        rf[key] = h
        print(f"  wrote {key}")

print(f"\nOutput : {os.path.abspath(args.output)}")

# ---- Readback check: print what actually ended up in the file ----------------
print("\nReadback check:")
with uproot.open(args.output) as rf:
    for name in rf[CHANNEL].keys():
        h    = rf[f"{CHANNEL}/{name}"]
        print(f"  {name:<30}  title='{h.title}'  "
              f"xaxis='{h.member('fXaxis').member('fTitle')}'  "
              f"nbins={h.member('fXaxis').member('fNbins')}")

# ----- PNG plot ------------------------------------------------

def save_plot(histograms, out_path, logy=True):
    """Overlay plot (CMS style) with a Full/SM ratio panel."""
    hep.style.use("CMS")

    # ax: main histogram panel (3x taller); ax_r: ratio panel below
    fig, (ax, ax_r) = plt.subplots(
        2, 1, figsize=(9, 7),
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
        sharex=True,
    )

    h_sm   = histograms["sm"]
    h_full = histograms[f"sm_lin_quad_{OP}"]
    edges  = h_sm.axes[0].edges
    cx     = 0.5 * (edges[:-1] + edges[1:])   # bin centres for errorbar

    # SM
    hep.histplot(h_sm.values(), edges, ax=ax,
                 label=f"SM  (integral={h_sm.values().sum():.3g})",
                 color="#1f77b4", linestyle="-", linewidth=1.8,
                 histtype="step", yerr=np.sqrt(h_sm.variances()))

    # quad
    h_quad = histograms[f"quad_{OP}"]
    hep.histplot(h_quad.values(), edges, ax=ax,
                 label=f"quad  (integral={h_quad.values().sum():.3g})",
                 color="#d62728", linestyle="--", linewidth=1.8,
                 histtype="step", yerr=np.sqrt(h_quad.variances()))

    # sm_lin_quad (full)
    hep.histplot(h_full.values(), edges, ax=ax,
                 label=f"sm_lin_quad  (integral={h_full.values().sum():.3g})",
                 color="#2ca02c", linestyle="-.", linewidth=1.8,
                 histtype="step", yerr=np.sqrt(h_full.variances()))

    # data_obs: drawn as errorbar points, not a stepped line
    h_data = histograms["data_obs"]
    ax.errorbar(cx, h_data.values(), yerr=np.sqrt(h_data.variances()),
                fmt="o", color="black",
                label=f"data_obs  (integral={h_data.values().sum():.3g})",
                markersize=5, linewidth=1.2, zorder=5)

    # ---- Ratio panel: full / SM ------------------------------------------------
    sm_v    = h_sm.values()
    full_v  = h_full.values()
    sm_e2   = h_sm.variances()
    full_e2 = h_full.variances()
    safe    = sm_v > 0
    ratio   = np.where(safe, full_v / np.where(safe, sm_v, 1), np.nan)
    rerr    = np.where(
        safe & (full_v > 0),
        ratio * np.sqrt(
            full_e2 / np.where(full_v > 0, full_v**2, 1)
            + sm_e2 / np.where(safe,        sm_v**2,  1)
        ),
        np.nan,
    )
    ax_r.errorbar(cx, ratio, yerr=rerr,
                  fmt="o", color="#2ca02c", markersize=4, linewidth=0.9)

    ax_r.axhline(1.0, color="black", linewidth=0.8, linestyle="--")
    ax_r.set_ylim(0.5, 1.7)
    ax_r.set_ylabel("Full / SM", fontsize=11)
    ax_r.set_xlabel(
        "Unrolled bin" if args.unrolled else r"$m_{\ell\ell}$ [GeV]",
        fontsize=13,
    )
    if not args.unrolled:
        # ax_r.set_xscale("log")
        pass

    # ----- Main panel decoration ------------------------------------------------
    ax.set_ylabel("Events / bin", fontsize=13)
    if logy:
        ax.set_yscale("log")
        ax.set_ylim(bottom=1e-6 * h_sm.values()[h_sm.values() > 0].max())
    ax.legend(
        frameon=True, framealpha=0.8, loc="upper right",
        fontsize=9, ncol=2,
    )
    if not args.unrolled:
        # ax.set_xscale("log")
        pass

    hep.cms.label(
        ax=ax,
        llabel="",
        rlabel=rf"$\mathcal{{L}}$ = {LUMI:.0f} pb$^{{-1}}$   $C_{{{OP}}}$ = {C}",
        loc=0,
    )

    suffix = "_log" if logy else "_lin"
    fname  = out_path.replace(".png", f"{suffix}.png")
    fig.savefig(fname, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {fname}")


if not args.no_plot:
    print("\nSaving bh plots ...")
    save_plot(histograms, plot_path, logy=True)
    save_plot(histograms, plot_path, logy=False)

# ---- ROOT macro ------------------------

macro_path = os.path.splitext(args.output)[0] + "_view.C"
macro = f"""\
// Auto-generated by make_combine_histograms.py
// Load in ROOT with:  root -l {macro_path}
//
void {os.path.splitext(os.path.basename(args.output))[0]}_view() {{
    TFile *f = TFile::Open("{os.path.abspath(args.output)}");
    TDirectory *dir = (TDirectory*)f->Get("{CHANNEL}");

    const char* procs[] = {{"sm","quad_{OP}","sm_lin_quad_{OP}","data_obs"}};
    Color_t cols[]      = {{4, 2, 8, 1}};
    Style_t lines[]     = {{1, 2, 9, 1}};
    Width_t widths[]    = {{2, 2, 2, 1}};

    TCanvas *c = new TCanvas("c_{CHANNEL}","{CHANNEL}  {OP}  C={C}",900,600);
    c->SetLogy();
    gStyle->SetOptStat(0);
    TLegend *leg = new TLegend(0.62,0.65,0.92,0.88);

    TH1 *href = nullptr;
    for (int i = 0; i < 4; i++) {{
        TH1 *h = (TH1*)dir->Get(procs[i]);
        if (!h) continue;
        h->SetLineColor(cols[i]);
        h->SetLineStyle(lines[i]);
        h->SetLineWidth(widths[i]);
        h->SetFillStyle(0);
        leg->AddEntry(h, procs[i], "l");
        if (!href) {{ h->Draw("HIST"); href = h; }}
        else        h->Draw("HIST SAME");
    }}
    leg->Draw();
    c->SaveAs("{os.path.splitext(args.output)[0]}_{OP}.pdf");
}}
"""
with open(macro_path, "w") as fmac:
    fmac.write(macro)
print(f"\nROOT macro : {macro_path}")
print(f"  run with : root -l {macro_path}")
