# ruff: noqa: E501
"""
config_dy_smeft_eft.py  —  DY SMEFTsim LO, EFT histogram approach
=================================================================
Uses runner_eft.py + subsamples to store EFT weight components as
separate histograms (no per-event arrays). After spritz-merge the pkl
contains, for each dataset × subsample:
  results["DYSMEFTsim_LO_mll_50_120_SM"]["histos"]["mll"]
  results["DYSMEFTsim_LO_mll_50_120_op01_lin"]["histos"]["mll"]
  results["DYSMEFTsim_LO_mll_50_120_op01_quad"]["histos"]["mll"]
  ...

To reconstruct the distribution at coupling c for operator k:
  H(c) = H_SM + c * H_opk_lin + c**2 * H_opk_quad

LHEReweightingWeight indexing (406 total):
  0        = SM
  1..27    = op_k at c=+1
  28..54   = op_k at c=-1  (k+27)
  55..405  = pairs (op_i, op_j)
"""

import json

import hist
import numpy as np
from spritz.framework.framework import cmap_petroff, get_fw_path

fw_path = get_fw_path()
with open(f"{fw_path}/data/common/lumi.json") as file:
    lumis = json.load(file)

year = "Full2018v9"
lumi = lumis[year]["tot"] / 1000
plot_label = "DY SMEFT LO EFT"
year_label = "2018"
njobs = 1000

runner = f"{fw_path}/src/spritz/runners/runner_3DY_trees_singleTriggers.py"

special_analysis_cfg = {
    "do_theory_variations": False,
}

MLL_BINS = [
    "50_120",
    "120_200",
    "200_400",
    "400_600",
    "600_800",
    "800_1000",
    "1000_3000",
]

# ── EFT subsamples ────────────────────────────────────────────────────────────
# One subsample per EFT component: SM + lin_k + quad_k for k = 1..27
# These string expressions are eval()'d inside the runner.
rwgt = "events.LHEReweightingWeight"

subsamples_eft = {
    "SM": f"{rwgt}[:, 0]",
}
for _k in range(1, 28):
    _km = _k + 27
    subsamples_eft[f"op{_k:02d}_lin"]  = f"0.5 * ({rwgt}[:, {_k}] - {rwgt}[:, {_km}])"
    subsamples_eft[f"op{_k:02d}_quad"] = f"0.5 * ({rwgt}[:, {_k}] + {rwgt}[:, {_km}] - 2 * {rwgt}[:, 0])"

# ── Datasets ──────────────────────────────────────────────────────────────────
datasets = {}
for b in MLL_BINS:
    name = f"DYSMEFTsim_LO_mll_{b}"
    datasets[name] = {
        "files": name,
        "task_weight": 8,
        "read_form": "mc",
        "subsamples": subsamples_eft,
    }

# ── Samples (for spritz-plot / spritz-postproc) ───────────────────────────────
# Each dataset × subsample becomes a separate entry in the merged pkl.
# We define only SM here; add operators as needed for plotting.
samples = {}
colors = {}
for b in MLL_BINS:
    key = f"DYSMEFTsim_LO_mll_{b}_SM"
    samples[key] = {"samples": [key]}
    colors[key] = cmap_petroff[0]

# ── Regions ───────────────────────────────────────────────────────────────────
preselections = lambda events: (events.mll > 50)  # noqa E731

regions = {
    "inc_ee": {
        "func": lambda events: preselections(events) & events["ee"],
        "mask": 0,
    },
    "inc_mm": {
        "func": lambda events: preselections(events) & events["mm"],
        "mask": 0,
    },
    "inc_em": {
        "func": lambda events: preselections(events) & events["em"],
        "mask": 0,
    },
}

# ── Binning ───────────────────────────────────────────────────────────────────
mll_bins = [
    *range(50, 76, 5),    # 50–75:  5 GeV steps
    *range(76, 106, 2),   # 76–105: 2 GeV steps (Z peak)
    *range(106, 120, 5),  # 106–119: 5 GeV steps
    120, 150, 200, 250, 300, 400, 600, 800, 1000, 1500, 3000,
]
costheta_bins = [-1, -0.6, -0.2, 0.2, 0.6, 1]
yZ_bins       = [-3.0, -1.5, 0.0, 1.5, 3.0]


def cos_theta_star(l1, l2):
    get_sign = lambda nr: nr / abs(nr)  # noqa E731
    return (
        2 * get_sign((l1 + l2).pz) / (l1 + l2).mass
        * get_sign(l1.pdgId)
        * (l2.pz * l1.energy - l1.pz * l2.energy)
        / np.sqrt(((l1 + l2).mass) ** 2 + ((l1 + l2).pt) ** 2)
    )


# ── Variables (histograms, no per-event arrays) ───────────────────────────────
variables = {
    "mll": {
        "func": lambda events: (events.Lepton[:, 0] + events.Lepton[:, 1]).mass,
        "axis": hist.axis.Variable(mll_bins, name="mll"),
    },
    "costhetastar_bins": {
        "func": lambda events: cos_theta_star(events.Lepton[:, 0], events.Lepton[:, 1]),
        "axis": hist.axis.Variable(costheta_bins, name="costhetastar_bins"),
    },
    "yZ_bins": {
        "func": lambda events: (events.Lepton[:, 0] + events.Lepton[:, 1]).eta,
        "axis": hist.axis.Variable(yZ_bins, name="yZ_bins"),
    },
}

nuisances = {}
nuisances["stat"] = {
    "type": "auto",
    "maxPoiss": "10",
    "includeSignal": "0",
    "samples": {},
}

check_weights = {}
