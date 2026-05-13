# ruff: noqa: E501
"""
config_dy_smeft_v5.py  —  DY SMEFTsim LO, EFT raw weights (Giacomo's approach)
=================================================================================
Stores raw LHE reweighting weights at c=+1 and c=-1 for each operator,
using real operator names. Lin/quad decomposition is done at plotting time.

This fixes the postproc normalisation issue of the lin/quad approach (v4):
  sumw_p1 ≈ sumw_SM  →  ratio h_op/h_SM is correct after postproc.

After spritz-merge the pkl contains, for each dataset x subsample:
  results["DYSMEFTsim_LO_mll_50_120_SM"]["histos"]["mll"]
  results["DYSMEFTsim_LO_mll_50_120_cHDD"]["histos"]["mll"]    # c=+1
  results["DYSMEFTsim_LO_mll_50_120_cHDD_m1"]["histos"]["mll"] # c=-1
  ...

At plotting time:
  vals_lin  = 0.5 * (vals_p1 - vals_m1)
  vals_quad = 0.5 * (vals_p1 + vals_m1) - vals_sm
  H(c) = H_SM + c * H_lin + c**2 * H_quad

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

runner = "/grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/spritz/runner_dy_smeft.py"

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

# -- EFT subsamples ------------------------------------------------------------
# One subsample per EFT component: SM + lin_k + quad_k for each operator.
# Operator names and indices extracted from LHE reweighting block.
# Format: (mask_expression, weight_expression) — both eval()'d in the runner.
# All events pass the mask; weight is the EFT component for that subsample.
rwgt = "events.LHEReweightingWeight"
all_mask = "ak.ones_like(events.run) == 1"

# LHE weight index → operator name (indices 1..27 = c=+1, 28..54 = c=-1)
OPERATORS = [
    "cHDD", "cHWB", "cbWRe", "cbBRe",
    "cHj1", "cHQ1", "cHj3",  "cHQ3",
    "cHu",  "cHd",  "cHbq",
    "cHl1", "cHl3", "cHe",
    "cll1",
    "clj1", "clj3",
    "cQl1", "cQl3",
    "ceu",  "ced",
    "cbe",  "cje",  "cQe",
    "clu",  "cld",  "cbl",
]  # 27 operators, matching LHE indices 1..27

# Raw weights at c=+1 and c=-1 — like Giacomo's approach.
# sumw_p1 ≈ sumw_SM so postproc normalisation is consistent.
# Lin/quad decomposition is done at plotting time.
subsamples_eft = {
    "SM": (all_mask, f"{rwgt}[:, 0]"),
}
for _k, _name in enumerate(OPERATORS, start=1):
    _km = _k + 27
    subsamples_eft[f"{_name}"]    = (all_mask, f"{rwgt}[:, {_k}]")   # c=+1
    subsamples_eft[f"{_name}_m1"] = (all_mask, f"{rwgt}[:, {_km}]")  # c=-1

# -- Datasets ------------------------------------------------------------------
datasets = {}
for b in MLL_BINS:
    name = f"DYSMEFTsim_LO_mll_{b}"
    datasets[name] = {
        "files": name,
        "task_weight": 8,
        "read_form": "mc",
        "subsamples": subsamples_eft,
    }

# -- Samples (for spritz-plot / spritz-postproc) -------------------------------
# One entry per subsample (SM + lin + quad for each operator), summing all mll bins.
# This writes every histogram into histos.root for clean uproot-based plotting.
samples = {}
colors = {}

# SM
samples["DYSMEFTsim_SM"] = {
    "samples": [f"DYSMEFTsim_LO_mll_{b}_SM" for b in MLL_BINS]
}
colors["DYSMEFTsim_SM"] = cmap_petroff[0]

# One entry per operator: raw weight at c=+1 and c=-1, real operator names.
for _i, _op in enumerate(OPERATORS, start=1):
    samples[f"DYSMEFTsim_{_op}"] = {
        "samples": [f"DYSMEFTsim_LO_mll_{b}_{_op}" for b in MLL_BINS]
    }
    samples[f"DYSMEFTsim_{_op}_m1"] = {
        "samples": [f"DYSMEFTsim_LO_mll_{b}_{_op}_m1" for b in MLL_BINS]
    }
    colors[f"DYSMEFTsim_{_op}"]    = cmap_petroff[(_i - 1) % len(cmap_petroff)]
    colors[f"DYSMEFTsim_{_op}_m1"] = cmap_petroff[(_i - 1) % len(cmap_petroff)]

# -- Regions -------------------------------------------------------------------
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

# -- Binning -------------------------------------------------------------------
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


# -- Variables (histograms, no per-event arrays) -------------------------------
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
