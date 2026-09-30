{
  "flags": [
    "goodVertices",
    "globalSuperTightHalo2016Filter",
    "HBHENoiseFilter",
    "HBHENoiseIsoFilter",
    "EcalDeadCellTriggerPrimitiveFilter",
    "BadPFMuonFilter",
    "BadPFMuonDzFilter",
    "hfNoisyHitsFilter",
    "eeBadScFilter",
    "ecalBadCalibFilter"
  ],
  "leptonSF": "RPLME_PATH_FW/data/Full2018v9/clib/lepton_sf.json.gz",
  "puWeightsKey": "Collisions18_UltraLegacy_goldenJSON",
  "run_to_era": "RPLME_PATH_FW/data/Full2018v9/clib/run_to_era.json.gz",
  "era": "Full2018v9",
  "year": "2018",
  "leptonsWP": {
    "muWP": "tightId"
  },
  "muRecoSfKey": "NUM_TrackerMuons_DEN_genTracks",
  "muIdSfKey": "NUM_TightID_DEN_TrackerMuons",
  "muIsoSfKey": "NUM_TightRelIso_DEN_TightIDandIPCut",
  "muTrigSfKey": "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight",
  "bTag": {
    "btagLoose": 0.049,
    "btagMedium": 0.2783,
    "btagTight": 0.71
  },
  "bVeto": {
    "wp": "Loose",
    "pt": 30
  },
  "jme": {
    "jet_veto_tag": "Summer19UL18_V1",
    "lvl_compound": "L1L2L3Res",
    "jet_algo": "AK4PFchs",
    "jer_tag": "Summer19UL18_JRV2_MC",
    "jec_tag": {
      "mc": "Summer19UL18_V5_MC",
      "data": {
        "UL2018A": "Summer19UL18_RunA_V5_DATA",
        "UL2018B": "Summer19UL18_RunB_V5_DATA",
        "UL2018C": "Summer19UL18_RunC_V5_DATA",
        "UL2018D": "Summer19UL18_RunD_V5_DATA"
      }
    },
    "jes": [
      "Absolute",
      "Absolute_2018",
      "BBEC1",
      "BBEC1_2018",
      "EC2",
      "EC2_2018",
      "FlavorQCD",
      "HF",
      "HF_2018",
      "RelativeBal",
      "RelativeSample_2018"
    ]
  },
  "jer_smear": "RPLME_PATH_FW/data/Full2018v9/clib/jer_smear.json.gz",
  "jet_jerc": "RPLME_PATH_FW/data/Full2018v9/clib/jet_jerc.json.gz",
  "puidSF": "RPLME_PATH_FW/data/Full2018v9/clib/jmar.json.gz",
  "jetvetomaps": "RPLME_PATH_FW/data/Full2018v9/clib/jetvetomaps.json.gz",
  "met": "RPLME_PATH_FW/data/Full2018v9/clib/met.json.gz",
  "btagSF": "RPLME_PATH_FW/data/Full2018v9/clib/btagging.json.gz",
  "btagEfficiency": "RPLME_PATH_FW/data/Full2018v9/clib/btagging_eff.json.gz",
  "puWeights": "RPLME_PATH_FW/data/Full2018v9/clib/puWeights.json.gz",
  "lumiMask": "RPLME_PATH_FW/data/Full2018v9/lumimask/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt",
  "rochester_file": "RPLME_PATH_FW/data/Full2018v9/rochester/RoccoR2018UL.txt",
  "fakesRW": "RPLME_PATH_FW/data/Full2018v9/clib/fakes_rw.json",
  "nloewRW": "RPLME_PATH_FW/data/common/nlo_ew.json"
}