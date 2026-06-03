# ruff: noqa: E501

import json

import awkward as ak
import hist
import numpy as np
from spritz.framework.framework import cmap_pastel, cmap_petroff, get_fw_path

fw_path = get_fw_path()

year = "Full2018v9"
runner = f"{fw_path}/src/spritz/runners/runner_3DY.py"

with open(f"{fw_path}/data/common/lumi.json") as file:
    lumis = json.load(file)

lumi = lumis[year]["tot"] / 1000  # All of 2018
plot_label = "DY"
year_label = "2018"
njobs = 500

special_analysis_cfg = {
    "do_theory_variations": True,
    "do_rochester_variations": True,
    "do_variations": True,
    "invert_one_isolation_loose": False,
    "invert_one_isolation_control": False,
    "skip_genmatching": False,
    "reweight_fakes": True
}

bins = {
    "mll": np.linspace(40, 200, 64),
}

datasets = {
    # DYmm and DYtt commented out — xrootd auth failures, rerun with fresh proxy
    # "DYmm_M-10to50": { "files": "DYJetsToMuMu_M-10to50", "task_weight": 8, "max_weight": 1e9 },
    # "DYmm_M-50to100": { "files": "DYJetsToMuMu", "task_weight": 8, "max_weight": 1e9 },
    # "DYmm_M-100to200": { "files": "DYJetsToMuMu_M-100to200", "task_weight": 8, "max_weight": 1e9 },
    # "DYmm_M-200to400": { "files": "DYJetsToMuMu_M-200to400", "task_weight": 8, "max_weight": 1e9 },
    # "DYmm_M-400to500": { "files": "DYJetsToMuMu_M-400to500", "task_weight": 8, "max_weight": 1e9 },
    # "DYmm_M-500to700": { "files": "DYJetsToMuMu_M-500to700", "task_weight": 8, "max_weight": 1e9 },
    # "DYmm_M-700to800": { "files": "DYJetsToMuMu_M-700to800", "task_weight": 8, "max_weight": 1e9 },
    # "DYmm_M-800to1000": { "files": "DYJetsToMuMu_M-800to1000", "task_weight": 8, "max_weight": 1e9 },
    # "DYmm_M-1000to1500": { "files": "DYJetsToMuMu_M-1000to1500", "task_weight": 8, "max_weight": 1e9 },
    # "DYmm_M-1500to2000": { "files": "DYJetsToMuMu_M-1500to2000", "task_weight": 8, "max_weight": 1e9 },
    # "DYmm_M-2000toInf": { "files": "DYJetsToMuMu_M-2000toInf", "task_weight": 8, "max_weight": 1e9 },
    # "DYtt": { "files": "DYJetsToTauTau", "task_weight": 8, "max_weight": 1e9 },
    "ST_s-channel": {
        "files": "ST_s-channel",
        "task_weight": 8,
    },
    "ST_t-channel_top_5f": {
        "files": "ST_t-channel_top_5f",
        "task_weight": 8,
    },
    "ST_t-channel_antitop_5f": {
        "files": "ST_t-channel_antitop_5f",
        "task_weight": 8,
    },
    "ST_tW_top_noHad": {
        "files": "ST_tW_top_noHad",
        "task_weight": 8,
    },
    "ST_tW_antitop_noHad": {
        "files": "ST_tW_antitop_noHad",
        "task_weight": 8,
    },
    "TTTo2L2Nu": {
        "files": "TTTo2L2Nu",
        "task_weight": 8,
        "top_pt_rwgt": True,
    },
    # "TTToSemiLeptonic": {
    #     "files": "TTToSemiLeptonic",
    #     "task_weight": 8,
    #     "top_pt_rwgt": True,
    # },
    "WWTo2L2Nu": {
        "files": "WWTo2L2Nu",
        "task_weight": 8,
    },
    "WZ": {
        "files": "WZ_TuneCP5_13TeV-pythia8",
        "task_weight": 8,
    },
    "ZZ": {
        "files": "ZZ_TuneCP5_13TeV-pythia8",
        "task_weight": 8,
    },
    # GGToMuMu commented out — xrootd auth failures, rerun with fresh proxy
    # "GGToMuMu_M-10to30_El-El": { "files": "GGToMuMu_M-10to30_El-El", "task_weight": 8 },
    # "GGToMuMu_M-10to30_Inel-El_El-Inel": { "files": "GGToMuMu_M-10to30_Inel-El_El-Inel", "task_weight": 8 },
    # "GGToMuMu_M-10to30_Inel-Inel": { "files": "GGToMuMu_M-10to30_Inel-Inel", "task_weight": 8 },
    # "GGToMuMu_M-30to50_El-El": { "files": "GGToMuMu_M-30to50_El-El", "task_weight": 8 },
    # "GGToMuMu_M-30to50_Inel-El_El-Inel": { "files": "GGToMuMu_M-30to50_Inel-El_El-Inel", "task_weight": 8 },
    # "GGToMuMu_M-30to50_Inel-Inel": { "files": "GGToMuMu_M-30to50_Inel-Inel", "task_weight": 8 },
    # "GGToMuMu_M-50to200_El-El": { "files": "GGToMuMu_M-50to200_El-El", "task_weight": 8 },
    # "GGToMuMu_M-50to200_Inel-El_El-Inel": { "files": "GGToMuMu_M-50to200_Inel-El_El-Inel", "task_weight": 8 },
    # "GGToMuMu_M-50to200_Inel-Inel": { "files": "GGToMuMu_M-50to200_Inel-Inel", "task_weight": 8 },
    # "GGToMuMu_M-200to1500_El-El": { "files": "GGToMuMu_M-200to1500_El-El", "task_weight": 8 },
    # "GGToMuMu_M-200to1500_Inel-El_El-Inel": { "files": "GGToMuMu_M-200to1500_Inel-El_El-Inel", "task_weight": 8 },
    # "GGToMuMu_M-200to1500_Inel-Inel": { "files": "GGToMuMu_M-200to1500_Inel-Inel", "task_weight": 8 },
    # "GGToMuMu_M-1500toInf_El-El": { "files": "GGToMuMu_M-1500toInf_El-El", "task_weight": 8 },
    # "GGToMuMu_M-1500toInf_Inel-El_El-Inel": { "files": "GGToMuMu_M-1500toInf_Inel-El_El-Inel", "task_weight": 8 },
    # "GGToMuMu_M-1500toInf_Inel-Inel": { "files": "GGToMuMu_M-1500toInf_Inel-Inel", "task_weight": 8 },
    # "WJetsToLNu_0J": {
    #     "files": "WJetsToLNu_0J",
    #     "task_weight": 8,
    # },
    # "WJetsToLNu_1J": {
    #     "files": "WJetsToLNu_1J",
    #     "task_weight": 8,
    # },
    # "WJetsToLNu_2J": {
    #     "files": "WJetsToLNu_2J",
    #     "task_weight": 8,
    # },
}


for dataset in datasets:
    datasets[dataset]["read_form"] = "mc"


DataRun = [
    ["A", "Run2018A-UL2018-v1"],
    ["B", "Run2018B-UL2018-v1"],
    ["C", "Run2018C-UL2018-v1"],
    ["D", "Run2018D-UL2018-v1"],
]

DataSets = ["SingleMuon"]

DataTrig = {
    "SingleMuon": "events.SingleMu",
}


samples_data = []
for era, sd in DataRun:
    for pd in DataSets:
        tag = pd + "_" + sd

        if "Run2018" in sd and "Muon" in pd:
            tag = tag.replace("v1","GT36")

        datasets[f"{pd}_{era}"] = {
            "files": tag,
            "trigger_sel": DataTrig[pd],
            "read_form": "data",
            "is_data": True,
            "era": f"UL2018{era}",
        }
        samples_data.append(f"{pd}_{era}")


samples = {
    "Data": {
        "samples": samples_data,
        "is_data": True,
    },
    # "W+Jets": {
    #     "samples": [
    #         "WJetsToLNu_0J",
    #         "WJetsToLNu_1J",
    #         "WJetsToLNu_2J",
    #    ]
    # },
    # "GGToLL": { ... },  # commented out — xrootd auth failures
    "Single Top": {
        "samples": [
            "ST_s-channel",
            "ST_t-channel_top_5f",
            "ST_t-channel_antitop_5f",
            "ST_tW_top_noHad",
            "ST_tW_antitop_noHad",
        ]
    },
    "TT": {
        "samples": [
            "TTTo2L2Nu",
        ]
    },
    # "TTToSemiLeptonic": {
    #     "samples": [
    #         "TTToSemiLeptonic"
    #     ]
    # },
    "WW": {
        "samples": [
            "WWTo2L2Nu",
        ]
    },
    "WZ": {
        "samples": [
            "WZ",
        ]
    },
    "ZZ": {
        "samples": [
            "ZZ"
        ]
    },
    # "DYtt": { ... },   # commented out — xrootd auth failures
    # "DYll": { ... },   # commented out — xrootd auth failures
}

colors = {}
#colors["W+Jets"] = cmap_petroff[0]
colors["Fakes"] = cmap_petroff[0]
# colors["GGToLL"] = cmap_petroff[1]  # commented out
colors["Single Top"] = cmap_petroff[2]
colors["TT"] = cmap_petroff[3]
colors["WW"] = cmap_petroff[4]
colors["WZ"] = cmap_petroff[5]
# colors["ZZ"] = cmap_petroff[6]
# colors["DYtt"] = cmap_petroff[8]
# colors["DYll"] = cmap_petroff[9]

#ugly fix for this cmap
colors["ZZ"]   = "#e42536"
# colors["DYtt"] = "#964a8b"  # commented out
# colors["DYll"] = "#f89c20"  # commented out

# regions

preselections = lambda events: (events.mll > 40)

regions = {
    "inc_mm": {
        "func": lambda events: preselections(events) & events["mm"] & (events.mll < 500),
        "mask": 0
    },
    # "veto_mm": {
    #     "func": lambda events: preselections(events) & ((events.mll < 60) | (events.mll > 120)) & (events.mll < 500) & events["mm"],
    #     "mask": 0
    # },
    "inc_mm_ss": {
        "func": lambda events: preselections(events) & (events.mll < 500) & events["mm_ss"],
        "mask": 0
    },
    # "veto_mm_ss": {
    #     "func": lambda events: preselections(events) & ((events.mll < 60) | (events.mll > 120)) & (events.mll < 500) & events["mm_ss"],
    #     "mask": 0
    # },
}

def cos_theta_star(l1, l2):
    get_sign = lambda nr: nr/abs(nr)
    return 2*get_sign((l1+l2).pz)/(l1+l2).mass * get_sign(l1.pdgId)*(l2.pz*l1.energy-l1.pz*l2.energy)/np.sqrt(((l1+l2).mass)**2+((l1+l2).pt)**2)

def transverse_mass(l, nu):
    return np.sqrt(2*l.pt*nu.pt*(1-np.cos(l.phi-nu.phi)))

def iso_transverse_mass(l1, l2, nu):
    return ak.where(
        l1.pfRelIso04_all < l2.pfRelIso04_all,
        transverse_mass(l1, nu),
        transverse_mass(l2, nu)
    )

variables = {
    #############
    # Dilepton
    #############
    "mll": {
        "func": lambda events: (events.Lepton[:, 0] + events.Lepton[:, 1]).mass,
        "axis": hist.axis.Regular(64, 40, 200, name="mll"),
        "label": "$m_{\\ell\\ell}$",
        "unit": "GeV"
    },
    "mll_medium": {
        "func": lambda events: (events.Lepton[:, 0] + events.Lepton[:, 1]).mass,
        "axis": hist.axis.Variable([40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,
            130,140,150,160,170,180,190,200,210,220,230,240,255,270,285,300,325,350,
            375,400,450,500], name="mll_medium"),
        "label": "$m_{\\ell\\ell}$",
        "unit": "GeV",
        "xlog": True
    },
    "ptll": {
        "func": lambda events: (events.Lepton[:, 0] + events.Lepton[:, 1]).pt,
        "axis": hist.axis.Regular(40, 0, 400, name="ptll"),
        "label": "$p_{T}^{\\ell\\ell}$",
        "unit": "GeV"
    },
    "costhetastar": {
        "func": lambda events: cos_theta_star(events.Lepton[:, 0], events.Lepton[:, 1]),
        "axis": hist.axis.Regular(50, -1, 1, name="costhetastar"),
        "label": "$cos\\,\\theta^{\\ast}$"
    },
    "rapll_abs": {
        "func": lambda events: abs((events.Lepton[:, 0] + events.Lepton[:, 1]).rapidity),
        "axis": hist.axis.Regular(50, 0, 2.5, name="rapll_abs"),
        "label": "$|y_{\\ell\\ell}|$"
    },
    #############
    # Single lepton
    #############
    "ptl1": {
        "func": lambda events: events.Lepton[:, 0].pt,
        "axis": hist.axis.Regular(50, 30, 280, name="ptl1"),
        "label": "$p_{T}^{\\ell_{1}}$",
        "unit": "GeV"
    },
    "etal1": {
        "func": lambda events: events.Lepton[:, 0].eta,
        "axis": hist.axis.Regular(50, -2.5, 2.5, name="etal1"),
        "label": "$\\eta_{\\ell_{1}}$"
    },
    "ptl2": {
        "func": lambda events: events.Lepton[:, 1].pt,
        "axis": hist.axis.Regular(50, 15, 165, name="ptl2"),
        "label": "$p_{T}^{\\ell_{2}}$",
        "unit": "GeV"
    },
    "etal2": {
        "func": lambda events: events.Lepton[:, 1].eta,
        "axis": hist.axis.Regular(50, -2.5, 2.5, name="etal2"),
        "label": "$\\eta_{\\ell_{2}}$"
    },
    #############
    # Multi-differential
    #############
    "triple_diff": {
        "axis": [
            hist.axis.Variable([40,60,80,100,120,140,180,220,270,350,500], name="mll"),
            hist.axis.Variable([-1.0,-0.6,-0.2,0.2,0.6,1.0], name="costhetastar"),
            hist.axis.Variable([0.0,0.48,0.96,1.44,2.4], name="rapll_abs"),
        ],
        "label": ["$m_{\\ell\\ell}$", "$cos\\,\\theta^{\\ast}$", "$|y_{\\ell\\ell}|$"],
        "unit": ["GeV","",""],
        "xlog": True
    },
}

mc_samples = [skey for skey in samples if not samples[skey].get('is_data',False)]

nuisances = {
    "lumi": {
        "name": "lumi",
        "type": "lnN",
        "samples": dict((skey, "1.0084") for skey in mc_samples)
    },
    ## Use the following if you want to apply the automatic combine MC stat nuisances
    # "stat": {
    #     "type": "auto",
    #     "maxPoiss": "10",
    #     "includeSignal": "0",
    #     "samples": {}
    # },
    # "mu_reco": {
    #     "name": "mu_reco",
    #     "type": "shape",
    #     "samples": mc_samples,
    #     "kind": "weight"
    # },
    # "mu_idiso": {
    #     "name": "mu_idiso",
    #     "type": "shape",
    #     "samples": mc_samples,
    #     "kind": "weight"
    # },
    # "mu_trig": {
    #     "name": "mu_trig",
    #     "type": "shape",
    #     "samples": mc_samples,
    #     "kind": "weight"
    # },
    # "PU": {
    #     "name": "PU",
    #     "type": "shape",
    #     "samples": mc_samples,
    #     "kind": "weight"
    # },
    # "prefireWeight": {
    #     "name": "prefireWeight",
    #     "type": "shape",
    #     "samples": mc_samples,
    #     "kind": "weight"
    # },
    # "tt_ptrw": {
    #     "name": "tt_ptrw",
    #     "type": "shape",
    #     "samples": ['TT'],
    #     "kind": "weight"
    # },
    # "rochester_stat": {
    #     "name": "rochester_stat",
    #     "type": "shape",
    #     "kind": "stdev",
    #     "samples": { k: [f"rochester_stat{i}" for i in range(100)] for k in samples},
    # },
    # "rochester_syst": {
    #     "name": "rochester_syst",
    #     "type": "shape",
    #     "kind": "square",
    #     "samples": { k: [f"rochester_{set_i}" for set_i in ["set2","set3","set4"]] for k in samples},
    # },
    # "QCDscale": {
    #     "name": "QCDScale",
    #     "type": "shape",
    #     "kind": "envelope",
    #     "samples": ({ k: [f"QCDScale_{i}" for i in [0,1,3,4,5,7,8]] for k in ['Single Top', 'TT', 'WW'] }
    #         | { k: [(f"QCDScale_{2*i}", f"QCDScale_{i}") for i in [0,1,3,4,5,7,8]] for k in ['DYll', 'DYtt'] }),
    #     "is_theory_unc": True
    # },
    # "PDFweight": {
    #     "name": "PDFweight",
    #     "type": "shape",
    #     "kind": "square",
    #     "samples": { k: [f"PDFWeight_{i}" for i in range(101)] for k in ['DYll', 'DYtt', 'Single Top', 'TT', 'WW'] },
    #     "is_theory_unc": True
    # },
    # "alphaS": {
    #     "name": "alphaS",
    #     "type": "shape",
    #     "kind": "envelope",
    #     "samples": { k: [f"PDFWeight_{i}" for i in [101,102]] for k in ['DYll', 'DYtt'] },
    #     "is_theory_unc": True
    # },
    # "PSWeight": {
    #     "name": "PSWeight",
    #     "type": "shape",
    #     "kind": "envelope",
    #     "samples": { k: [f"PSWeight_{i}" for i in range(4)] for k in ['DYll', 'DYtt', 'Single Top', 'TT', 'WW', 'WZ', 'ZZ'] },
    #     "is_theory_unc": True
    # },
}
