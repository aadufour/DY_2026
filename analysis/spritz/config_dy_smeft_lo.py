# ruff: noqa: E501

import json

import hist
import numpy as np
from spritz.framework.framework import cmap_petroff, get_fw_path


fw_path = get_fw_path()
with open(f"{fw_path}/data/common/lumi.json") as file:
    lumis = json.load(file)

year = "Full2018v9"
lumi = lumis[year]["tot"] / 1000  # All of 2018
plot_label = "DY SMEFT LO"
year_label = "2018"
njobs = 1000

runner = f"{fw_path}/src/spritz/runners/runner_3DY_trees_singleTriggers.py"

special_analysis_cfg = {
    "do_theory_variations": False
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

datasets = {}
for b in MLL_BINS:
    name = f"DYSMEFTsim_LO_mll_{b}"
    datasets[name] = {
        "files": name,
        "task_weight": 8,
        "read_form": "mc",
    }

samples = {}
for b in MLL_BINS:
    name = f"DYSMEFTsim_LO_mll_{b}"
    samples[name] = {"samples": [name]}

colors = {}

# regions
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

gen_mll_bins   = [50, 120, 200, 400, 600, 800, 1000, 3000]
costheta_bins  = [-1, -0.6, -0.2, 0.2, 0.6, 1]
yZ_bins        = [-3.0, -1.5, 0.0, 1.5, 3.0]


def cos_theta_star(l1, l2):
    get_sign = lambda nr: nr / abs(nr)  # noqa E731
    return (
        2 * get_sign((l1 + l2).pz) / (l1 + l2).mass
        * get_sign(l1.pdgId)
        * (l2.pz * l1.energy - l1.pz * l2.energy)
        / np.sqrt(((l1 + l2).mass) ** 2 + ((l1 + l2).pt) ** 2)
    )


variables = {
    "mll": {
        "func": lambda events: (events.Lepton[:, 0] + events.Lepton[:, 1]).mass,
        "axis": hist.axis.Regular(60, 50, 3000, name="mll"),
        "save_events": True,
    },
    "costhetastar_bins": {
        "func": lambda events: cos_theta_star(events.Lepton[:, 0], events.Lepton[:, 1]),
        "axis": hist.axis.Variable(costheta_bins, name="costhetastar_bins"),
        "save_events": True,
    },
    "yZ_bins": {
        "func": lambda events: (events.Lepton[:, 0] + events.Lepton[:, 1]).eta,
        "axis": hist.axis.Variable(yZ_bins, name="yZ_bins"),
        "save_events": True,
    },
    "Gen_mll": {
        "func": lambda events: events.Gen_mll,
        "save_events": True,
    },
}

# Save all 406 LHEReweightingWeights by index.
# Ordering: 0=SM, 1..27=op_i(+1), 28..54=op_i(-1), 55..405=pairs(op_i,op_j)
# Use build_cache.py to map indices to operator names.
for _i in range(406):
    variables[f"w_{_i}"] = {
        "func": lambda events, i=_i: events.LHEReweightingWeight[:, i],
        "save_events": True,
    }

nuisances = {}

nuisances["stat"] = {
    "type": "auto",
    "maxPoiss": "10",
    "includeSignal": "0",
    "samples": {},
}

check_weights = {}
