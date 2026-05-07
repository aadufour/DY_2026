import pylhe
import numpy as np
import boost_histogram as bh
import mplhep as hep
import matplotlib.pyplot as plt

# ── LHE files ────────────────────────────────────────────────────────────────
# weight_id=None  → use nominal event weight (eventinfo.weight)  [C=+1 sample]
# weight_id='rwgt_1' → use MG5 reweight block weight             [C=0, C=-1]
LHE_FILES = {
    "C=+1":    ("/Users/albertodufour/MG5_2_9_18/mg5amcnlo/DY_all/myLHE/unweighted_events.lhe",   None),
    "SM (C=0)":("/Users/albertodufour/MG5_2_9_18/mg5amcnlo/DY_all/myLHE/weighted_events_C0.lhe",  "rwgt_1"),
    "C=-1":    ("/Users/albertodufour/MG5_2_9_18/mg5amcnlo/DY_all/myLHE/weighted_events_C-1.lhe", None),
}

# ── Kinematic helpers (same convention as test_pylhe3.py) ────────────────────
def mll(p1, p2):
    p = np.array(p1) + np.array(p2)
    return np.sqrt(p[3] ** 2 - sum(p[l] ** 2 for l in range(3)))


def rap(p1, p2):
    p = np.array(p1) + np.array(p2)
    E = p[3]
    pz = p[2]
    return 0.5 * np.log((E + pz) / (E - pz))


def cstar(p1, p2):
    """cos θ* of the lepton in the dilepton rest frame."""
    p1 = np.array(p1)
    p2 = np.array(p2)

    p = p1 + p2
    E = p[3]
    pz = p[2]
    mass = mll(p1, p2)

    beta = pz / E
    gamma = E / mass

    pz1_boosted = gamma * (p1[2] - beta * p1[3])
    p1_mag = np.sqrt(p1[0] ** 2 + p1[1] ** 2 + pz1_boosted ** 2)

    return pz1_boosted / p1_mag


# ── Histogram axes ────────────────────────────────────────────────────────────
ax_mll = bh.axis.Regular(200, 60, 120, metadata="m_{ll} [GeV]")
ax_rap = bh.axis.Regular(100, -5, 5, metadata="y_{ll}")
ax_cs = bh.axis.Regular(50, -1, 1, metadata="cos θ*")


def make_hists():
    return (
        bh.Histogram(ax_mll, storage=bh.storage.Double()),
        bh.Histogram(ax_rap, storage=bh.storage.Double()),
        bh.Histogram(ax_cs, storage=bh.storage.Double()),
    )


# ── Fill histograms ---------------------------------
def fill_from_lhe(path, label, weight_id=None):
    """
    weight_id=None    : nominal weight from event header (C=+1 file)
    weight_id='rwgt_1': MG5 reweight block (C=0 / C=-1 files, same
                        phase-space points reweighted to new parameters)
    """
    vals_mll, vals_rap, vals_cs, vals_w = [], [], [], []

    for i, event in enumerate(pylhe.read_lhe_with_attributes(path)):
        if not (i + 1) % 5000:
            print(f"  [{label}] processed {i + 1} events")

        if weight_id is None:
            w = event.eventinfo.weight
        else:
            w = event.weights[weight_id]

        final = [p for p in event.particles if p.status == 1.0]
        pmu = [[final[k].px, final[k].py, final[k].pz, final[k].e] for k in range(2)]
        vals_mll.append(mll(pmu[0], pmu[1]))
        vals_rap.append(rap(pmu[0], pmu[1]))
        vals_cs.append(cstar(pmu[0], pmu[1]))
        vals_w.append(w)

    print(f"  [{label}] done — {i + 1} events total, filling histograms...")

    w_arr = np.array(vals_w)
    h_mll, h_rap, h_cs = make_hists()
    h_mll.fill(np.array(vals_mll), weight=w_arr)
    h_rap.fill(np.array(vals_rap), weight=w_arr)
    h_cs.fill(np.array(vals_cs),   weight=w_arr)

    return h_mll, h_rap, h_cs


# ── PLOTTING ──────────────────────────────────────────────────────────────
COLORS = ["tab:red", "tab:blue", "tab:green"]
STYLES = ["-", "--", ":"]

OUT_DIR = "/Users/albertodufour/code/DY2026/reweighting_test"


def plot_observable(hists_dict, obs_label, out_file, logy=False):
    hep.style.use("CMS")
    fig, ax = plt.subplots()
    for (label, h), color, ls in zip(hists_dict.items(), COLORS, STYLES):
        hep.histplot(
            h,
            ax=ax,
            label=label,
            color=color,
            linestyle=ls,
            density=True,
            histtype="step",
            linewidth=1.5,
        )
    ax.set_xlabel(obs_label)
    ax.set_ylabel("Normalised events / bin")
    if logy:
        ax.set_yscale("log")
    ax.legend()
    hep.cms.label("MG5", ax=ax)
    fig.savefig(f"{OUT_DIR}/{out_file}", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out_file}")


##---------------------------------------------------------------------------------------
if __name__ == "__main__":
    results = {}
    for label, (path, weight_id) in LHE_FILES.items():
        print(f"\nReading {label} : {path}")
        results[label] = fill_from_lhe(path, label, weight_id=weight_id)

    # unpack per-observable dicts
    mll_hists = {lbl: results[lbl][0] for lbl in results}
    rap_hists = {lbl: results[lbl][1] for lbl in results}
    cs_hists = {lbl: results[lbl][2] for lbl in results}


    plot_observable(mll_hists, r"$m_{\ell\ell}$ [GeV]", "mll_SMEFT.png")
    plot_observable(rap_hists, r"$y_{\ell\ell}$", "rap_SMEFT.png")
    plot_observable(cs_hists, r"$\cos\theta^*$", "cstar_SMEFT.png")
