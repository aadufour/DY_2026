# ruff: noqa: E501
"""
config_syst.py  —  Fabian's backgrounds + DY SMEFTsim LO EFT signal + theory systematics
==========================================================================================
Same as config.py but with systematics enabled:
  - do_theory_variations = True  → QCDscale, PDFweight, alphaS, PSWeight from NanoAOD
  - do_variations = True         → all variation histograms kept (not stripped to nom-only)
  - do_rochester_variations = True → Rochester stat (100) + syst (3 sets)
  - noStat = True on all EFT samples → zeroes bin variances so combine autoMCStats
                                        does not apply to EFT templates
  - stat autoMCStats nuisance added for backgrounds only (includeSignal=0)

Experimental systematics (mu_*, PU, prefire, tt_ptrw) commented out — pending Giacomo.

Runner: runner.py (unchanged — flags drive behaviour via special_analysis_cfg)
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
plot_label = "DY EFT"
year_label = "2018"
njobs = 3000

runner = "/grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/spritz/runner.py"

special_analysis_cfg = {
    "do_theory_variations": True,       # QCDscale + PDF + PSWeight from LHEScaleWeight/LHEPdfWeight/PSWeight
    "do_rochester_variations": True,    # Rochester stat (100 toys) + syst (3 sets)
    "do_variations": True,              # keep all variation histograms (not just nom)
    "invert_one_isolation_loose": False,
    "invert_one_isolation_control": False,
    "skip_genmatching": False,
    "reweight_fakes": False,
}

# ---- EFT setup ---------------------------------------------------------------
MLL_BINS = [
    "50_120", "120_200", "200_400", "400_600",
    "600_800", "800_1000", "1000_3000",
]

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
]  # 27 operators, matching LHE weight indices 1..27

rwgt = "events.LHEReweightingWeight"
all_mask = "ak.ones_like(events.run) == 1"

subsamples_eft = {"sm": (all_mask, f"{rwgt}[:, 0]")}
for _k, _name in enumerate(OPERATORS, start=1):
    _km = _k + 27
    subsamples_eft[f"w1_{_name}"]  = (all_mask, f"{rwgt}[:, {_k}]")
    subsamples_eft[f"wm1_{_name}"] = (all_mask, f"{rwgt}[:, {_km}]")

# ---- Datasets ----------------------------------------------------------------
datasets = {}

# -- Fabian's backgrounds (DY MiNNLO, TT, VV, single top, GGToLL) --
datasets["DYmm_M-10to50"]    = {"files": "DYJetsToMuMu_M-10to50",     "task_weight": 8, "max_weight": 1e9, "read_form": "mc"}
datasets["DYmm_M-50to100"]   = {"files": "DYJetsToMuMu",              "task_weight": 8, "max_weight": 1e9, "read_form": "mc"}
datasets["DYmm_M-100to200"]  = {"files": "DYJetsToMuMu_M-100to200",   "task_weight": 8, "max_weight": 1e9, "read_form": "mc"}
datasets["DYmm_M-200to400"]  = {"files": "DYJetsToMuMu_M-200to400",   "task_weight": 8, "max_weight": 1e9, "read_form": "mc"}
datasets["DYmm_M-400to500"]  = {"files": "DYJetsToMuMu_M-400to500",   "task_weight": 8, "max_weight": 1e9, "read_form": "mc"}
datasets["DYmm_M-500to700"]  = {"files": "DYJetsToMuMu_M-500to700",   "task_weight": 8, "max_weight": 1e9, "read_form": "mc"}
datasets["DYmm_M-700to800"]  = {"files": "DYJetsToMuMu_M-700to800",   "task_weight": 8, "max_weight": 1e9, "read_form": "mc"}
datasets["DYmm_M-800to1000"] = {"files": "DYJetsToMuMu_M-800to1000",  "task_weight": 8, "max_weight": 1e9, "read_form": "mc"}
datasets["DYmm_M-1000to1500"]= {"files": "DYJetsToMuMu_M-1000to1500", "task_weight": 8, "max_weight": 1e9, "read_form": "mc"}
datasets["DYmm_M-1500to2000"]= {"files": "DYJetsToMuMu_M-1500to2000", "task_weight": 8, "max_weight": 1e9, "read_form": "mc"}
datasets["DYmm_M-2000toInf"] = {"files": "DYJetsToMuMu_M-2000toInf",  "task_weight": 8, "max_weight": 1e9, "read_form": "mc"}
datasets["DYtt"]             = {"files": "DYJetsToTauTau",             "task_weight": 8, "max_weight": 1e9, "read_form": "mc"}

datasets["ST_s-channel"]           = {"files": "ST_s-channel",           "task_weight": 8, "read_form": "mc"}
datasets["ST_t-channel_top_5f"]    = {"files": "ST_t-channel_top_5f",    "task_weight": 8, "read_form": "mc"}
datasets["ST_t-channel_antitop_5f"]= {"files": "ST_t-channel_antitop_5f","task_weight": 8, "read_form": "mc"}
datasets["ST_tW_top_noHad"]        = {"files": "ST_tW_top_noHad",        "task_weight": 8, "read_form": "mc"}
datasets["ST_tW_antitop_noHad"]    = {"files": "ST_tW_antitop_noHad",    "task_weight": 8, "read_form": "mc"}

datasets["TTTo2L2Nu"] = {"files": "TTTo2L2Nu", "task_weight": 8, "top_pt_rwgt": True, "read_form": "mc"}
datasets["WWTo2L2Nu"] = {"files": "WWTo2L2Nu", "task_weight": 8, "read_form": "mc"}
datasets["WZ"]        = {"files": "WZ_TuneCP5_13TeV-pythia8", "task_weight": 8, "read_form": "mc"}
datasets["ZZ"]        = {"files": "ZZ_TuneCP5_13TeV-pythia8", "task_weight": 8, "read_form": "mc"}

for _tag, _m in [
    ("M-10to30",    ["El-El", "Inel-El_El-Inel", "Inel-Inel"]),
    ("M-30to50",    ["El-El", "Inel-El_El-Inel", "Inel-Inel"]),
    ("M-50to200",   ["El-El", "Inel-El_El-Inel", "Inel-Inel"]),
    ("M-200to1500", ["El-El", "Inel-El_El-Inel", "Inel-Inel"]),
    ("M-1500toInf", ["El-El", "Inel-El_El-Inel", "Inel-Inel"]),
]:
    for _mode in _m:
        _name = f"GGToMuMu_{_tag}_{_mode}"
        datasets[_name] = {"files": _name, "task_weight": 8, "read_form": "mc"}

# -- EFT signal (DY SMEFTsim LO, 7 mll bins, 55 subsamples each) --
# EFT: 406 → runner filters events to those with exactly 406 LHEReweightingWeights
for b in MLL_BINS:
    name = f"DYSMEFTsim_LO_mll_{b}"
    datasets[name] = {
        "files": name,
        "task_weight": 8,
        "read_form": "mc",
        "EFT": 406,
        "subsamples": subsamples_eft,
    }

# -- Data (SingleMuon 2018 A-D) --
DataRun = [
    ["A", "Run2018A-UL2018-v1"],
    ["B", "Run2018B-UL2018-v1"],
    ["C", "Run2018C-UL2018-v1"],
    ["D", "Run2018D-UL2018-v1"],
]

samples_data = []
for era, sd in DataRun:
    tag = f"SingleMuon_{sd}"
    if "Muon" in sd or "Run2018" in sd:
        tag = tag.replace("v1", "GT36")
    datasets[f"SingleMuon_{era}"] = {
        "files": tag,
        "trigger_sel": "events.SingleMu",
        "read_form": "data",
        "is_data": True,
        "era": f"UL2018{era}",
    }
    samples_data.append(f"SingleMuon_{era}")

# ---- Samples -----------------------------------------------------------------
samples = {
    "Data": {"samples": samples_data, "is_data": True},
    "GGToLL": {
        "samples": [
            f"GGToMuMu_{m}_{mode}"
            for m in ["M-10to30","M-30to50","M-50to200","M-200to1500","M-1500toInf"]
            for mode in ["El-El","Inel-El_El-Inel","Inel-Inel"]
        ]
    },
    "Single Top": {
        "samples": [
            "ST_s-channel", "ST_t-channel_top_5f", "ST_t-channel_antitop_5f",
            "ST_tW_top_noHad", "ST_tW_antitop_noHad",
        ]
    },
    "TT":   {"samples": ["TTTo2L2Nu"]},
    "WW":   {"samples": ["WWTo2L2Nu"]},
    "WZ":   {"samples": ["WZ"]},
    "ZZ":   {"samples": ["ZZ"]},
    "DYtt": {"samples": ["DYtt"]},
    "DYll": {
        "samples": [
            "DYmm_M-10to50", "DYmm_M-50to100", "DYmm_M-100to200",
            "DYmm_M-200to400", "DYmm_M-400to500", "DYmm_M-500to700",
            "DYmm_M-700to800", "DYmm_M-800to1000", "DYmm_M-1000to1500",
            "DYmm_M-1500to2000", "DYmm_M-2000toInf",
        ],
        "is_signal": True,
    },
    # EFT signal subsamples (morphing convention)
    # noStat=True: zeroes bin variances in shapes.root → autoMCStats does not
    # apply to EFT templates (bin errors are correlated across sm/w1/wm1, not independent Poisson)
    "sm": {"samples": [f"DYSMEFTsim_LO_mll_{b}_sm" for b in MLL_BINS], "noStat": True},
}
for _op in OPERATORS:
    samples[f"w1_{_op}"]  = {"samples": [f"DYSMEFTsim_LO_mll_{b}_w1_{_op}"  for b in MLL_BINS], "noStat": True}
    samples[f"wm1_{_op}"] = {"samples": [f"DYSMEFTsim_LO_mll_{b}_wm1_{_op}" for b in MLL_BINS], "noStat": True}

# ---- Colors ------------------------------------------------------------------
colors = {}
colors["Fakes"]      = cmap_petroff[0]
colors["GGToLL"]     = cmap_petroff[1]
colors["Single Top"] = cmap_petroff[2]
colors["TT"]         = cmap_petroff[3]
colors["WW"]         = cmap_petroff[4]
colors["WZ"]         = cmap_petroff[5]
colors["ZZ"]         = "#e42536"
colors["DYtt"]       = "#964a8b"
colors["DYll"]       = "#f89c20"
colors["sm"]         = "#000000"
for _i, _op in enumerate(OPERATORS):
    colors[f"w1_{_op}"]  = cmap_petroff[_i % len(cmap_petroff)]
    colors[f"wm1_{_op}"] = cmap_petroff[_i % len(cmap_petroff)]

# ---- Regions -----------------------------------------------------------------
preselections = lambda events: (events.mll > 50)  # noqa E731

regions = {
    "inc_mm": {
        "func": lambda events: preselections(events) & events["mm"],
        "mask": 0,
    },
    "inc_mm_ss": {
        "func": lambda events: preselections(events) & events["mm_ss"],
        "mask": 0,
    },
}

# ---- Variables ---------------------------------------------------------------
# Fine mll binning for EFT combine fits (34 bins, 50–3000 GeV)
mll_bins = [
    *range(50, 76, 5),    # 50–75:  5 GeV steps
    *range(76, 106, 2),   # 76–105: 2 GeV steps (Z peak)
    *range(106, 120, 5),  # 106–119: 5 GeV steps
    120, 150, 200, 250, 300, 400, 600, 800, 1000, 1500, 3000,
]


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
        "axis": hist.axis.Variable(mll_bins, name="mll"),
        "label": "$m_{\\ell\\ell}$",
        "unit": "GeV",
    },
    "costhetastar": {
        "func": lambda events: cos_theta_star(events.Lepton[:, 0], events.Lepton[:, 1]),
        "axis": hist.axis.Regular(50, -1, 1, name="costhetastar"),
        "label": "$\\cos\\theta^*$",
    },
    "rapll_abs": {
        "func": lambda events: abs((events.Lepton[:, 0] + events.Lepton[:, 1]).rapidity),
        "axis": hist.axis.Regular(50, 0, 2.5, name="rapll_abs"),
        "label": "$|y_{\\ell\\ell}|$",
    },
}

# ---- Nuisances ---------------------------------------------------------------
mc_samples = [s for s in samples if not samples[s].get("is_data", False)]
bkg_samples = [s for s in mc_samples if not samples[s].get("noStat", False) and not samples[s].get("is_signal", False)]
theory_samples = ["DYll", "DYtt", "Single Top", "TT", "WW"]

nuisances = {
    "lumi": {
        "name": "lumi",
        "type": "lnN",
        "samples": {s: "1.0084" for s in mc_samples},
    },
    # autoMCStats: MC statistical uncertainty via Barlow-Beeston lite.
    # includeSignal=0: exclude signal (DYll, index=0). EFT templates additionally
    # have noStat=True so their bin variances are zeroed in shapes.root.
    "stat": {
        "type": "auto",
        "maxPoiss": "10",
        "includeSignal": "0",
        "samples": {},
    },
    # ---- Theory systematics (from NanoAOD branches) --------------------------
    "QCDscale": {
        "name": "QCDScale",
        "type": "shape",
        "kind": "envelope",
        "samples": (
            {k: [f"QCDScale_{i}" for i in [0,1,3,4,5,7,8]] for k in ["Single Top", "TT", "WW"]}
            | {k: [(f"QCDScale_{2*i}", f"QCDScale_{i}") for i in [0,1,3,4,5,7,8]] for k in ["DYll", "DYtt"]}
        ),
        "is_theory_unc": True,
    },
    "PDFweight": {
        "name": "PDFweight",
        "type": "shape",
        "kind": "square",
        "samples": {k: [f"PDFWeight_{i}" for i in range(101)] for k in theory_samples},
        "is_theory_unc": True,
    },
    "alphaS": {
        "name": "alphaS",
        "type": "shape",
        "kind": "envelope",
        "samples": {k: [f"PDFWeight_{i}" for i in [101, 102]] for k in ["DYll", "DYtt"]},
        "is_theory_unc": True,
    },
    "PSWeight": {
        "name": "PSWeight",
        "type": "shape",
        "kind": "envelope",
        "samples": {k: [f"PSWeight_{i}" for i in range(4)] for k in ["DYll", "DYtt", "Single Top", "TT", "WW", "WZ", "ZZ"]},
        "is_theory_unc": True,
    },
    # ---- Experimental systematics (pending Giacomo — enable after audit) ------
    # "mu_reco": {
    #     "name": "mu_reco",
    #     "type": "shape",
    #     "samples": {k: None for k in bkg_samples},
    #     "kind": "weight",
    # },
    # "mu_idiso": {
    #     "name": "mu_idiso",
    #     "type": "shape",
    #     "samples": {k: None for k in bkg_samples},
    #     "kind": "weight",
    # },
    # "mu_trig": {
    #     "name": "mu_trig",
    #     "type": "shape",
    #     "samples": {k: None for k in bkg_samples},
    #     "kind": "weight",
    # },
    # "PU": {
    #     "name": "PU",
    #     "type": "shape",
    #     "samples": {k: None for k in bkg_samples},
    #     "kind": "weight",
    # },
    # "prefireWeight": {
    #     "name": "prefireWeight",
    #     "type": "shape",
    #     "samples": {k: None for k in bkg_samples},
    #     "kind": "weight",
    # },
    # "tt_ptrw": {
    #     "name": "tt_ptrw",
    #     "type": "shape",
    #     "samples": {"TT": None},
    #     "kind": "weight",
    # },
    # "rochester_stat": {
    #     "name": "rochester_stat",
    #     "type": "shape",
    #     "kind": "stdev",
    #     "samples": {k: [f"rochester_stat{i}" for i in range(100)] for k in bkg_samples},
    # },
    # "rochester_syst": {
    #     "name": "rochester_syst",
    #     "type": "shape",
    #     "kind": "square",
    #     "samples": {k: [f"rochester_{s}" for s in ["set2", "set3", "set4"]] for k in bkg_samples},
    # },
}

check_weights = {}
