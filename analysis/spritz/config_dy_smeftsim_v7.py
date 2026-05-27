# ruff: noqa: E501
"""
config_dy_smeftsim_v7.py  —  DY SMEFTsim LO, morphing-convention names + theory systematics
==============================================================================================
Changes from v6 (config_dy_smeft_v5.py):

1. Sample/subsample names use morphing convention directly:
     SM       → "sm"
     c=+1 op  → "w1_{op}"
     c=-1 op  → "wm1_{op}"
   This means histos.root from spritz-postproc already has the correct histogram
   names expected by AnomalousCouplingMorphing.py. build_shapes_morphing.py no
   longer needs to translate names — it just copies and adds histo_Data.

2. do_theory_variations = True  →  runner fills PDF and QCD scale variation slices
   in the "syst" axis of each histogram.

3. nuisances["QCDscale"] and nuisances["PDF"] declared with lheScaleWeight /
   lhePdfWeight types so spritz-datacard produces the corresponding shape
   nuisance lines in the datacard.

Runner: runner_dy_smeft_v7.py  (theory variations unlocked)
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

runner = "/grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/spritz/runner_dy_smeft_v7.py"

special_analysis_cfg = {
    "do_theory_variations": True,
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
# Names use morphing convention so histos.root is directly compatible with
# AnomalousCouplingMorphing.py:
#   sm      = SM template        (LHE weight index 0)
#   w1_{op} = c=+1 raw template  (LHE weight index k,   k = 1..27)
#   wm1_{op}= c=-1 raw template  (LHE weight index k+27)
rwgt = "events.LHEReweightingWeight"
all_mask = "ak.ones_like(events.run) == 1"

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

subsamples_eft = {
    "sm": (all_mask, f"{rwgt}[:, 0]"),
}
for _k, _name in enumerate(OPERATORS, start=1):
    _km = _k + 27
    subsamples_eft[f"w1_{_name}"]  = (all_mask, f"{rwgt}[:, {_k}]")    # c=+1
    subsamples_eft[f"wm1_{_name}"] = (all_mask, f"{rwgt}[:, {_km}]")   # c=-1

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
# Morphing-convention names → histos.root has histo_sm, histo_w1_cHDD, etc.
samples = {}
colors = {}

# SM
samples["sm"] = {
    "samples": [f"DYSMEFTsim_LO_mll_{b}_sm" for b in MLL_BINS]
}
colors["sm"] = cmap_petroff[0]

# One entry per operator: raw weight at c=+1 and c=-1
for _i, _op in enumerate(OPERATORS, start=1):
    samples[f"w1_{_op}"] = {
        "samples": [f"DYSMEFTsim_LO_mll_{b}_w1_{_op}" for b in MLL_BINS]
    }
    samples[f"wm1_{_op}"] = {
        "samples": [f"DYSMEFTsim_LO_mll_{b}_wm1_{_op}" for b in MLL_BINS]
    }
    colors[f"w1_{_op}"]  = cmap_petroff[(_i - 1) % len(cmap_petroff)]
    colors[f"wm1_{_op}"] = cmap_petroff[(_i - 1) % len(cmap_petroff)]

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


# -- Variables (histograms) ----------------------------------------------------
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

# -- Nuisances -----------------------------------------------------------------
nuisances = {}
nuisances["stat"] = {
    "type": "auto",
    "maxPoiss": "10",
    "includeSignal": "0",
    "samples": {},
}
# QCDScale: 8 variations (first 4 + last 4 of LHEScaleWeight, see theory_unc.py)
# → envelope (max/min deviation from nominal)
_qcd_vars = [f"QCDScale_{i}" for i in range(8)]

# PDFWeight: 103 entries (NNPDF31_nnlo_as_0118: 0=central, 1..100=replicas, 101..102=alphas vars)
# → square (RMS across all replicas)
_pdf_vars = [f"PDFWeight_{i}" for i in range(103)]

nuisances["QCDscale"] = {
    "name": "QCDscale",
    "type": "shape",
    "kind": "lheScaleWeight_envelope",
    "samples": {s: _qcd_vars for s in samples},
}
nuisances["PDF"] = {
    "name": "PDF",
    "type": "shape",
    "kind": "lhePdfWeight_square",
    "samples": {s: _pdf_vars for s in samples},
}

check_weights = {}
