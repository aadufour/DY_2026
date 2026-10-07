import json
import hist
import numpy as np
from itertools import combinations_with_replacement
from spritz.framework.framework import cmap_pastel, cmap_petroff, get_fw_path, get_rw_idx, get_eft_points

fw_path = get_fw_path()
with open(f"{fw_path}/data/common/lumi.json") as file:
    lumis = json.load(file)

year = "Full2018v9"
lumi = lumis[year]["tot"] / 1000
lumi_unc = lumis[year]["rel_unc"]
plot_label = "2018 EFT"
year_label = "2018"
njobs = 1000

runner = f"{fw_path}/src/spritz/runners/runner_3DY_eft.py"

special_analysis_cfg = {
    "do_variations": True,
    "do_theory_variations": True, # 116 variations
    "do_rochester_stat_variations": False, # 100 variations
    "do_jet_variations": False, # 24 variations
    "invert_one_isolation_loose": False,
    "invert_one_isolation_control": False,
    "reweight_fakes": True,
}

# Higher order corrections
n3lo_qcd = {
    "file": f"{fw_path}/data/common/kfactor_ewscheme3_3D_N3LO_N3LL_NNLO_NNLL.root",
    "object": "ratio_N3LO+N3LL_over_NNLO+NNLL",
    "name": "N3LO_QCD",
}
nlo_ew = {
    "file": f"{fw_path}/data/common/powheg_ew_ratio.root",
    "object": "h_ratio",
    "name": "NLO_EW",
}

# SMEFT samples
file_key = "DYMuMu_LO_EFT_SMEFTsim_propcorr_mll%s_Photos_startingOne"
eft_bins = ["50_120","120_200","200_400","400_600"]#,"600_800","800_1000","1000_3000"]
operators = ["clj1", "clj3", "ceu", "ced", "cje", "clu", "cld"]

eft_points = get_eft_points(operators)
cov_pairs = list(combinations_with_replacement(eft_points, 2))

eft_datasets = {
    f"DYmm_mll{b}": {
        "files": file_key % b,
        "task_weight": 8,
        "eft_reweighting": {
            "weight_branch": "LHEReweightingWeight",
            "points": {p: get_rw_idx(file_key % b, p) for p in eft_points},
            "covariance_pairs": cov_pairs
        }
    } for b in eft_bins
}

# SM samples
dy_bins = ["10to50","50to100","100to200","200to400","400to500","500to700"]#,"700to800",
    #"800to1000","1000to1500","1500to2000","2000toInf"]
ggll_bins = ["10to30","30to50","50to200","200to1500"]#,"1500toInf"]

datasets = eft_datasets | { 
    **{
        f"DYmm_MiNNLO_M-{b}": {
            "files": f"DYJetsToMuMu_M-{b}" if b!="50to100" else "DYJetsToMuMu",
            "task_weight": 8,
            "max_weight": 1e9,
            "ho_corrections": [n3lo_qcd, nlo_ew] 
        } for b in dy_bins 
    },
    "DYtt": {
        "files": "DYJetsToTauTau",
        "task_weight": 8,
        "max_weight": 1e9,
        "ho_corrections": [nlo_ew]
    },
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
    "TTToSemiLeptonic": {
        "files": "TTToSemiLeptonic",
        "task_weight": 8,
        "top_pt_rwgt": True,
        "genmatching_nlep": 1,
    },
    "WWTo2L2Nu": {
        "files": "WWTo2L2Nu",
        "task_weight": 8,
    },
    "WZTo3LNu": {
        "files": "WZTo3LNu",
        "task_weight": 8,
    },
    "WZTo2Q2L": {
        "files": "WZTo2Q2L",
        "task_weight": 8,
    },
    "ZZTo4L": {
        "files": "ZZTo4L",
        "task_weight": 8,
    },
    "ZZTo2L2Nu": {
        "files": "ZZTo2L2Nu",
        "task_weight": 8,
    },
    "ZZTo2Q2L": {
        "files": "ZZTo2Q2L",
        "task_weight": 8,
    },
    **{
        f"GGToMuMu_M-{b}_{case}": {
            "files": f"GGToMuMu_M-{b}_{case}",
            "task_weight": 8 
        } for b in ggll_bins 
        for case in ["El-El", "Inel-El_El-Inel", "Inel-Inel"] 
    }
}

for dataset in datasets:
    datasets[dataset]["read_form"] = "mc"

# Data
samples_data = []
for era in ["A", "B", "C", "D"]:
    datasets[f"SingleMuon_{era}"] = {
        "files": f"SingleMuon_Run{year_label}{era}-UL{year_label}-GT36",
        "trigger_sel": "events.SingleMu",
        "read_form": "data",
        "is_data": True,
        "era": f"UL{year_label}{era}"
    }
    samples_data.append(f"SingleMuon_{era}")

# Merge samples
samples = {
    "Data": {
        "samples": samples_data,
        "is_data": True,
    },
    "GGToLL": { 
        "samples": [
            f"GGToMuMu_M-{b}_{case}" for b in ggll_bins 
            for case in ["El-El", "Inel-El_El-Inel", "Inel-Inel"]
        ] 
    },
    "Single_Top": {
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
            "TTToSemiLeptonic"
        ]
    },
    "VV": {
        "samples": [
            "WWTo2L2Nu",
            "WZTo3LNu",
            "WZTo2Q2L",
            "ZZTo4L",
            "ZZTo2L2Nu",
            "ZZTo2Q2L"
        ]
    },
    "DYtt": {
        "samples": [
            "DYtt"
        ]
    },
    "DYmm_MiNNLO": {
        "samples": [
            f"DYmm_MiNNLO_M-{b}" for b in dy_bins
        ],
        "exclude_from_datacard": True
    },
    **{
        f"DYmm_{point}": {
            "samples": [f"{dataset}_{point}" for dataset in eft_datasets],
            "is_smeft": True,
            "is_signal": point != "sm" 
        } for point in eft_points
    },
    **{
        f"DYmm_cov_{a}_{b}": {
            "samples": [f"{dataset}_cov_{a}_{b}" for dataset in eft_datasets],
            "covariance_of": (f"DYmm_{a}", f"DYmm_{b}"),
            "is_variance": True,
            "exclude_from_datacard": True,
        } for a, b in cov_pairs
    }
}

colors = {
    "Fakes": cmap_petroff[0],
    "GGToLL": cmap_petroff[1],
    "Single_Top": cmap_petroff[2],
    "TT": cmap_petroff[3],
    "VV": cmap_petroff[4],
    "DYtt": cmap_petroff[8],
    "DYmm_MiNNLO": cmap_petroff[9],
    "DYmm": cmap_petroff[9]
}

# Renormalize samples from "reference" to "target"
renorm_samples = {
    "target": "DYmm_MiNNLO",
    "reference": "DYmm_sm",
    "samples": [f"DYmm_{point}" for point in eft_points]+[f"DYmm_cov_{a}_{b}" for a, b in cov_pairs]
}

# Define regions
preselections = lambda events: (50 < events.mll) & (events.mll < 500)

regions = {
    "bveto_mm": {
        "func": lambda events: preselections(events) & events.mm & events.bveto,
        "mask": 0,
    },
    "bveto_mm_ss": {
        "func": lambda events: preselections(events) & events.mm_ss & events.bveto,
        "mask": 0,
    },
}

# Define variables and histograms
def cos_theta_star(l1, l2):
    get_sign = lambda nr: nr / abs(nr)  # noqa: E731
    return (
        2 * get_sign((l1 + l2).pz) / (l1 + l2).mass
        * get_sign(l1.pdgId)
        * (l2.pz * l1.energy - l1.pz * l2.energy)
        / np.sqrt(((l1 + l2).mass) ** 2 + ((l1 + l2).pt) ** 2)
    )

variables = {
    "mll": {
        "func": lambda events: (events.Lepton[:, 0] + events.Lepton[:, 1]).mass,
        "axis": hist.axis.Regular(60, 50, 200, name="mll"),
        "label": "$m_{\\ell\\ell}$",
        "unit": "GeV"
    },
    "mll_medium": {
        "func": lambda events: (events.Lepton[:, 0] + events.Lepton[:, 1]).mass,
        "axis": hist.axis.Variable([50,55,60,65,70,75,80,85,90,95,100,105,110,
            115,120,130,140,150,160,170,180,190,200,220,240,260,280,300,325,350,375,
            400,450,500], name="mll_medium"),
        "label": "$m_{\\ell\\ell}$",
        "unit": "GeV",
        "xlog": True
    },
    "costhetastar": {
        "func": lambda events: cos_theta_star(events.Lepton[:, 0], events.Lepton[:, 1]),
        "axis": hist.axis.Regular(50, -1, 1, name="costhetastar"),
        "label": "$cos\\,\\theta^{\\ast}$"
    },
    "rapll_abs": {
        "func": lambda events: abs((events.Lepton[:, 0] + events.Lepton[:, 1]).rapidity),
        "axis": hist.axis.Regular(48, 0, 2.4, name="rapll_abs"),
        "label": "$|y_{\\ell\\ell}|$"
    },
    "triple_diff": {
        "axis": [
            hist.axis.Variable([50,70,80,90,100,110,120,150,200,300,500], name="mll"),
            hist.axis.Variable([-1.0,-0.5,0.0,0.5,1.0], name="costhetastar"),
            hist.axis.Variable([0.0,0.48,0.96,1.44,2.4], name="rapll_abs"),
        ],
        "label": ["$m_{\\ell\\ell}$", "$cos\\,\\theta^{\\ast}$", "$|y_{\\ell\\ell}|$"],
        "unit": ["GeV", "", ""],
        "xlog": True,
    },
}

cards_regions = ["bveto_mm"]
cards_variables = ["triple_diff"]
covariance_file = "covariance.root"

all_samples = [s for s in samples if not samples[s].get("is_variance")]
mc_samples = [s for s in all_samples if not samples[s].get("is_data")]

# Fakes
fakes_dict = {
    "regions": [
        {"target": "bveto_mm", "source": "bveto_mm_ss"}
    ],
    "subtract_mc": [
        s for s in mc_samples if not (samples[s].get("is_smeft"))
    ],
    "nuisances": [
        "Fakes transfer factor: Fit", "Fakes transfer factor: Model"
    ]
}

nuisances = {
    "lumi": {
        "name": "lumi",
        "type": "lnN",
        "samples": {s: str(lumi_unc) for s in mc_samples}
    },
    ## Use the following if you want to apply the automatic combine MC stat nuisances
    "stat": {
        "type": "auto",
        "maxPoiss": "10",
        "includeSignal": "0",
        "samples": {}
    },
    "Pile-up corr.": {
        "name": "PU",
        "type": "shape",
        "samples": mc_samples,
        "kind": "weight"
    },
    "L1 pre-firing corr.": {
        "name": "prefireWeight",
        "type": "shape",
        "samples": mc_samples,
        "kind": "weight"
    },
    #############
    # Leptons
    #############
    "Trigger SF": {
        "name": "mu_trig",
        "type": "shape",
        "samples": mc_samples,
        "kind": "weight"
    },
    "Muon Reconstruction SF": {
        "name": "mu_reco",
        "type": "shape",
        "samples": mc_samples,
        "kind": "weight"
    },
    "Muon ID SF": {
        "name": "mu_id",
        "type": "shape",
        "samples": mc_samples,
        "kind": "weight"
    },
    "Muon Isolation SF": {
        "name": "mu_iso",
        "type": "shape",
        "samples": mc_samples,
        "kind": "weight"
    },
    "Rochester corr. (syst)": {
        "name": "rochester_syst",
        "type": "shape",
        "kind": "square",
        "samples": all_samples,
        "variations": [
            {"label": "Rochester corr. set2", "tag": "rochester_set2"},
            {"label": "Rochester corr. set3", "tag": "rochester_set3"},
            {"label": "Rochester corr. set4", "tag": "rochester_set4"}
        ]
    },
    #############
    # Theory
    #############
    "NLO EW correction": {
        "name": "NLO_EW",
        "type": "shape",
        "samples": ["DYmm_MiNNLO", "DYtt", *[f"DYmm_{point}" for point in eft_points]],
        "kind": "weight"
    },
    "N3LO QCD correction": {
        "name": "N3LO_QCD",
        "type": "shape",
        "samples": ["DYmm_MiNNLO", *[f"DYmm_{point}" for point in eft_points]],
        "kind": "weight"
    },
    "Top $p_{T}$ corr.": {
        "name": "tt_ptrw",
        "type": "shape",
        "samples": ["TT"],
        "kind": "weight"
    },
    #############
    # b-tagging
    #############
    "puidSF": {
        "name": "puidSF",
        "type": "shape",
        "samples": mc_samples,
        "kind": "weight"
    },
    "btagSF_sf": {
        "name": "btagSF_sf",
        "type": "shape",
        "samples": mc_samples,
        "kind": "weight"
    },
    "btagSF_eff": {
        "name": "btagSF_eff",
        "type": "shape",
        "samples": mc_samples,
        "kind": "weight"
    },
    #############
    # Fakes
    #############
    "Fakes transfer factor: Fit": {
        "name": "fakes_param",
        "type": "shape",
        "kind": "weight",
        "samples": all_samples,
    },
    "Fakes transfer factor: Model": {
        "name": "fakes_model",
        "type": "shape",
        "kind": "envelope",
        "samples": all_samples,
        "variations": [
            {"label": "fakes_model", "tag": "fakes_model"}
        ],
    },
}

corrections = {
    "Pile-up corr.": { 
        "name": "PU",
        "samples": mc_samples 
    },
    "L1 pre-firing corr.": {
        "name": "prefireWeight",
        "samples": mc_samples 
    },
    "Trigger SF": { 
        "name": "mu_trig",
        "samples": mc_samples 
    },
    "Muon Reconstruction SF": { 
        "name": "mu_reco",
        "samples": mc_samples 
    },
    "Muon ID SF": { 
        "name": "mu_id",
        "samples": mc_samples 
    },
    "Muon Isolation SF": { 
        "name": "mu_iso",
        "samples": mc_samples 
    },
    "Rochester corr.": { 
        "name": "rochester",
        "samples": all_samples, 
        "related_nuisances": ["Rochester corr. (syst)"]
    },
    "NLO EW correction": { 
        "name": "NLO_EW",
        "samples": ["DYmm_MiNNLO", "DYtt", *[f"DYmm_{point}" for point in eft_points]] 
    },
    "N3LO QCD correction": { 
        "name": "N3LO_QCD",
        "samples": ["DYmm_MiNNLO", *[f"DYmm_{point}" for point in eft_points]] 
    },
    "Top $p_{T}$ corr.": { 
        "name": "tt_ptrw",
        "samples": ["TT"] 
    },
    "puidSF": { 
        "name": "puidSF",
        "samples": mc_samples 
    },
    "btagSF": { 
        "name": "btagSF",
        "samples": mc_samples,
        "related_nuisances": ["btagSF_sf", "btagSF_eff"] 
    },
    "Fakes transfer factor": { 
        "name": "fakes",
        "samples": all_samples,
        "related_nuisances": ["Fakes transfer factor: Fit", "Fakes transfer factor: Model"] 
    },
}
