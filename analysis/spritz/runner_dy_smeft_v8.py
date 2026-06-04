"""
runner_dy_smeft_v8.py
=====================
v8 of the DY SMEFTsim runner. Changes from v7:

1. check_weights axis removed from histograms.
   v7 produced 4-axis histograms (variable, category, syst, check_weights).
   spritz postproc (post_process_eft.py) expects exactly 3 axes (variable,
   category, syst). v7 required a manual pkl preprocessing step to drop the
   4th axis before running postproc. v8 never adds it, so postproc works
   directly on the merged pkl without any preprocessing.

   The check_weights loop is kept (runs once with cwgt="nominal") but the
   axis is not declared and check_weights= is not passed to fill().

Inherited from v7:
- Theory variations controlled by special_analysis_cfg["do_theory_variations"]
- Nom-only filter conditional on do_theory_variations
- histogram filling uses events[f"weight_{dataset_name}"][mask] for correct
  subsample EFT weights
- "1000_3000" mll LHE filter
- removed onnxruntime / DNN imports
"""
import gc
import json
import sys
import traceback as tb
from copy import deepcopy

import awkward as ak
import correctionlib
import hist
import numpy as np
import spritz.framework.variation as variation_module
import uproot
import vector
from spritz.framework.framework import (
    big_process,
    get_analysis_dict,
    get_fw_path,
    read_chunks,
    write_chunks,
)
from spritz.modules.basic_selections import (
    LumiMask,
    lumi_mask,
    pass_flags,
    pass_trigger,
    pass_weightfilter
)
from spritz.modules.btag_sf import btag_sf
from spritz.modules.gen_analysis import gen_analysis
from spritz.modules.jet_sel import jetSel
from spritz.modules.jme import (
    correct_jets_data,
    correct_jets_mc,
    jet_veto,
    remove_jets_HEM_issue,
)
from spritz.modules.lepton_sel import createLepton, leptonSel
from spritz.modules.lepton_sf import lepton_sf
from spritz.modules.prompt_gen import prompt_gen_match_leptons
from spritz.modules.puid_sf import puid_sf
from spritz.modules.puweight import puweight_sf
from spritz.modules.rochester import correctRochester, getRochester
try:
    from spritz.modules.run_assign import assign_run
except ImportError:
    from spritz.modules.run_assign import assign_run_period as assign_run
from spritz.modules.theory_unc import theory_unc
from spritz.modules.trigger_sf import (
    trigger_sf,
    match_trigger_object
)

# trigger_sf_latinos is not available in all apptainer versions and is not used
try:
    from spritz.modules.trigger_sf_latinos import trigger_sf_latinos
except ImportError:
    trigger_sf_latinos = None

vector.register_awkward()

print("uproot version", uproot.__version__)
print("awkward version", ak.__version__)

path_fw = get_fw_path()
with open("cfg.json") as file:
    txt = file.read()
    txt = txt.replace("RPLME_PATH_FW", path_fw)
    cfg = json.loads(txt)

ceval_puWeight = correctionlib.CorrectionSet.from_file(cfg["puWeights"])
ceval_lepton_sf = correctionlib.CorrectionSet.from_file(cfg["leptonSF"])
ceval_assign_run = correctionlib.CorrectionSet.from_file(cfg["run_to_era"])

cset_trigger = correctionlib.CorrectionSet.from_file(cfg["triggerSF"])
rochester = getRochester(cfg)

analysis_path = sys.argv[1]
analysis_cfg = get_analysis_dict(analysis_path)
special_analysis_cfg = analysis_cfg["special_analysis_cfg"]


def ensure_not_none(arr):
    if ak.any(ak.is_none(arr)):
        raise Exception("There are some None in branch", arr[ak.is_none(arr)])
    return ak.fill_none(arr, -9999.9)


def process(events, **kwargs):
    dataset = kwargs["dataset"]
    trigger_sel = kwargs.get("trigger_sel", "")
    isData = kwargs.get("is_data", False)
    era = kwargs.get("era", None)
    isData = kwargs.get("is_data", False)
    subsamples = kwargs.get("subsamples", {})
    special_weight = eval(kwargs.get("weight", "1.0"))

    print("SumW and NEvents before cuts")
    n_rwgt = events["nLHEReweightingWeight"] if "nLHEReweightingWeight" in ak.fields(events) else "N/A"
    print(len(events), n_rwgt)

    variations = variation_module.Variation()
    variations.register_variation([], "nom")

    if isData:
        events["weight"] = ak.ones_like(events.run)
    else:
        events["weight"] = events.genWeight
    if "EFT" in kwargs.keys() and kwargs["EFT"]:
        neft_rwgts = kwargs["EFT"]
        print("neft_rwgts--> ", neft_rwgts)
        events = events[ak.num(events.LHEReweightingWeight) == neft_rwgts]
        events["rwgt"] = ak.pad_none(
            events.LHEReweightingWeight, neft_rwgts, clip=True, axis=1
        )
        events["rwgt"] = ak.fill_none(events.rwgt, 0.0)
    if isData:
        lumimask = LumiMask(cfg["lumiMask"])
        events = lumi_mask(events, lumimask)
    else:
        events = pass_weightfilter(events, kwargs.get("max_weight", None))
        events = events[events.pass_weightfilter]

    sumw = ak.sum(events.weight)
    nevents_raw = ak.num(events.weight, axis=0)
    nevents = int(np.array(nevents_raw).sum())
    print("SumW and NEvents before cuts")
    print(sumw, nevents)
    print('-----')

    if special_weight != 1.0:
        print(f"Using special weight for {dataset}: {special_weight}")

    events["weight"] = events.weight * special_weight

    # pass trigger and flags
    events = assign_run(events, isData, cfg, ceval_assign_run)
    events = pass_trigger(events, cfg["era"])
    events = pass_flags(events, cfg["flags"])

    events = events[events.pass_flags & events.pass_trigger]

    print("Number of events after trigger and flags: ", len(events))

    if isData:
        events = events[eval(trigger_sel)]

    events = createLepton(events)

    events = leptonSel(events, cfg)
    events["Lepton"] = events.Lepton[events.Lepton.isLoose]
    events = events[ak.num(events.Lepton) >= 2]
    events = events[events.Lepton[:, 0].pt >= 20]
    events = events[events.Lepton[:, 1].pt >= 10]

    print("Number of events after basic Lepton sel: ", len(events))

    if not isData:
        events = prompt_gen_match_leptons(events)

    events = events[events.PV.npvsGood > 0]

    print("Number of events after prompt match and NPV>0: ", len(events))

    if kwargs.get("top_pt_rwgt", False):
        top_particle_mask = (events.GenPart.pdgId == 6) & ak.values_astype(
            (events.GenPart.statusFlags >> 13) & 1, bool
        )
        toppt = ak.fill_none(
            ak.mask(events, ak.num(events.GenPart[top_particle_mask]) >= 1)
            .GenPart[top_particle_mask][:, -1]
            .pt,
            0.0,
        )

        atop_particle_mask = (events.GenPart.pdgId == -6) & ak.values_astype(
            (events.GenPart.statusFlags >> 13) & 1, bool
        )
        atoppt = ak.fill_none(
            ak.mask(events, ak.num(events.GenPart[atop_particle_mask]) >= 1)
            .GenPart[atop_particle_mask][:, -1]
            .pt,
            0.0,
        )

        events['topPtWeight'] = (toppt * atoppt > 0.0) * np.sqrt(
            (0.103*np.exp(-0.0118*toppt) - 0.000134*toppt + 0.973)
            * (0.103*np.exp(-0.0118*atoppt) - 0.000134*atoppt + 0.973)
        ) + (toppt * atoppt <= 0.0)
    else:
        events['topPtWeight'] = ak.ones_like(events.weight)

    # Correct Muons with rochester
    try:
        events = correctRochester(events, isData, rochester, s=5)           # Giacomo's signature
    except TypeError:
        events, variations = correctRochester(events, variations, isData, rochester, s=5)  # Fabian's signature
    events = match_trigger_object(events, cfg)

    doTheoryVariations = special_analysis_cfg.get("do_theory_variations", False)

    if not isData:
        events, variations = puweight_sf(events, variations, ceval_puWeight, cfg)
        events, variations = trigger_sf(events, variations, ceval_lepton_sf, cfg)
        events, variations = lepton_sf(events, variations, ceval_lepton_sf, cfg)

        if "L1PreFiringWeight" in ak.fields(events):
            events["prefireWeight"] = events.L1PreFiringWeight.Nom
            events["prefireWeight_up"] = events.L1PreFiringWeight.Up
            events["prefireWeight_down"] = events.L1PreFiringWeight.Dn

            variations.register_variation(
                columns=["prefireWeight"],
                variation_name="prefireWeight_up",
                format_rule=lambda _, var_name: var_name,
            )
            variations.register_variation(
                columns=["prefireWeight"],
                variation_name="prefireWeight_down",
                format_rule=lambda _, var_name: var_name,
            )
        else:
            events["prefireWeight"] = ak.ones_like(events.weight)

        if doTheoryVariations:
            events, variations = theory_unc(events, variations)
    else:
        pass

    # Regions definitions
    regions = deepcopy(analysis_cfg["regions"])
    variables = deepcopy(analysis_cfg["variables"])
    check_weights = deepcopy(analysis_cfg["check_weights"])

    def fun(events, sample):
        return events[f"weight_{sample}"]

    check_weights["nominal"] = {
        "func": lambda events, sample: fun(events, sample)
    }

    # Only keep nominal variation (unless theory variations are enabled)
    if not doTheoryVariations:
        variations.variations_dict = {
            k: v for k, v in variations.variations_dict.items() if k == "nom"
        }

    # 3-axis histograms: (variable, category, syst)
    # check_weights axis deliberately omitted — postproc expects exactly 3 axes.
    default_axis = [
        hist.axis.StrCategory(
            [region for region in regions],
            name="category",
        ),
        hist.axis.StrCategory(
            sorted(list(variations.get_variations_all())),
            name="syst"
        ),
    ]

    results = {dataset: {"sumw": sumw, "nevents": nevents, "events": 0, "histos": 0}}
    if subsamples != {}:
        results = {}
        for subsample in subsamples:
            results[f"{dataset}_{subsample}"] = {
                "sumw": sumw,
                "nevents": nevents,
                "events": 0,
                "histos": 0,
            }

    for dataset_name in results:
        _events = {}
        histos = {}
        for variable in variables:
            _events[variable] = ak.Array([])

            if "axis" in variables[variable]:
                if isinstance(variables[variable]["axis"], list):
                    histos[variable] = hist.Hist(
                        *variables[variable]["axis"],
                        *default_axis,
                        hist.storage.Weight(),
                    )
                else:
                    histos[variable] = hist.Hist(
                        variables[variable]["axis"],
                        *default_axis,
                        hist.storage.Weight(),
                    )

        results[dataset_name]["histos"] = histos
        results[dataset_name]["events"] = _events

    originalEvents = ak.copy(events)
    jet_pt_backup = ak.copy(events.Jet.pt)

    print("Doing variations")
    for variation in sorted(variations.get_variations_all()):
        events = ak.copy(originalEvents)

        print(variation)
        for switch in variations.get_variation_subs(variation):
            if len(switch) == 2:
                variation_dest, variation_source = switch
                events[variation_dest] = events[variation_source]

        # resort Leptons
        lepton_sort = ak.argsort(events[("Lepton", "pt")], ascending=False, axis=1)
        events["Lepton"] = events.Lepton[lepton_sort]

        events = events[(ak.num(events.Lepton, axis=1) >= 2)]

        eleWP = cfg["leptonsWP"]["eleWP"]
        muWP = cfg["leptonsWP"]["muWP"]

        comb = ak.ones_like(events.run) == 1.0
        for ilep in range(2):
            comb = comb & (
                events.Lepton[:, ilep]["isTightElectron_" + eleWP]
                | events.Lepton[:, ilep]["isTightMuon_" + muWP]
            )
        events["l2Tight"] = ak.copy(comb)
        events = events[events.l2Tight]

        print("l2Tight: ", len(events))

        if len(events) == 0:
            continue

        # Define lepton flavour categories
        events["ll"] = (ak.ones_like(events.weight, dtype=bool))
        events["ee"] = (events.Lepton[:, 0].pdgId * events.Lepton[:, 1].pdgId) == -11 * 11
        events["mm"] = (events.Lepton[:, 0].pdgId * events.Lepton[:, 1].pdgId) == -13 * 13
        events["tt"] = (events.Lepton[:, 0].pdgId * events.Lepton[:, 1].pdgId) == -15 * 15
        events["em"] = (events.Lepton[:, 0].pdgId * events.Lepton[:, 1].pdgId) == -11 * 13
        events["em_ss"] = (events.Lepton[:, 0].pdgId * events.Lepton[:, 1].pdgId) == 11 * 13
        events["ee_ss"] = (events.Lepton[:, 0].pdgId * events.Lepton[:, 1].pdgId) == 11 * 11
        events["mm_ss"] = (events.Lepton[:, 0].pdgId * events.Lepton[:, 1].pdgId) == 13 * 13

        if not isData:
            events["prompt_gen_match_2l"] = (
                events.Lepton[:, 0].promptgenmatched
                & events.Lepton[:, 1].promptgenmatched
            )
            events = events[events.prompt_gen_match_2l]

        # Analysis level cuts
        leptoncut = (events.ee | events.mm | events.em | events.ee_ss | events.mm_ss | events.em_ss)

        # third lepton veto
        leptoncut = leptoncut & (
            ak.fill_none(
                ak.mask(
                    ak.all(events.Lepton[:, 2:].pt < 10, axis=1),
                    ak.num(events.Lepton) >= 3,
                ),
                True,
                axis=0,
            )
        )

        leptoncut = (
            leptoncut & (events.Lepton[:, 1].pt > 15) & (
                ((events.mm | events.mm_ss) & (events.Lepton[:, 0].pt > 30)) |
                ((events.ee | events.em | events.ee_ss | events.em_ss) & (events.Lepton[:, 0].pt > 38))
            )
        )

        events = events[leptoncut]

        if len(events) == 0:
            continue

        if not isData:
            events["RecoSF"] = events.Lepton[:, 0].RecoSF * events.Lepton[:, 1].RecoSF
            events["TightSF"] = events.Lepton[:, 0].TightSF * events.Lepton[:, 1].TightSF

            events["weight"] = (
                events.weight
                * events.puWeight
                * events.topPtWeight
                * events.RecoSF
                * events.TightSF
                * events.prefireWeight
                * events.TriggerSFweight_2l
            )

        # Gen level variables
        if not isData:
            events = gen_analysis(events, dataset)
        else:
            events["Gen_mll"] = ak.full_like(events.weight, -10.0)
            events["Gen_ptll"] = ak.full_like(events.weight, -10.0)
            events["Gen_dphill"] = ak.full_like(events.weight, -10.0)

        # Variable definitions
        for variable in variables:
            if "func" in variables[variable]:
                events[variable] = variables[variable]["func"](events)

        events[dataset] = ak.ones_like(events.run) == 1.0

        # LHE-level mll filters to avoid double counting between mll-binned samples
        if not isData and not dataset.startswith("GG") and dataset not in ["WZ", "ZZ"]:
            lhe_leptons_mask = (events.LHEPart.status == 1) & (
                (abs(events.LHEPart.pdgId) == 11)
                | (abs(events.LHEPart.pdgId) == 13)
                | (abs(events.LHEPart.pdgId) == 15)
            )
            lhe_leptons = events.LHEPart[lhe_leptons_mask]
            if ak.all(ak.num(lhe_leptons) == 2):
                lhe_mll = (lhe_leptons[:, 0] + lhe_leptons[:, 1]).mass

                if "50_100" in dataset or dataset in ["DYmm_M-50", "DYee_M-50"]:
                    events = events[(lhe_mll >= 50) & (lhe_mll <= 100)]
                    print("mass < 100 > 50 GeV: ", len(events))
                if "50_120" in dataset or "50to120" in dataset:
                    events = events[(lhe_mll >= 50) & (lhe_mll <= 120)]
                    print("mass < 120 > 50 GeV: ", len(events))
                if "100_200" in dataset or "100to200" in dataset:
                    events = events[(lhe_mll > 100) & (lhe_mll <= 200)]
                    print("mass < 200 > 100 GeV: ", len(events))
                if "120_200" in dataset or "120to200" in dataset:
                    events = events[(lhe_mll > 120) & (lhe_mll <= 200)]
                    print("mass < 200 > 120 GeV: ", len(events))
                if "200_400" in dataset or "200to400" in dataset:
                    events = events[(lhe_mll > 200) & (lhe_mll <= 400)]
                    print("mass < 400 > 200 GeV: ", len(events))
                if "400_600" in dataset or "400to600" in dataset:
                    events = events[(lhe_mll > 400) & (lhe_mll <= 600)]
                    print("mass < 600 > 400 GeV: ", len(events))
                if "400_500" in dataset or "400to500" in dataset:
                    events = events[(lhe_mll > 400) & (lhe_mll <= 500)]
                    print("mass < 500 > 400 GeV: ", len(events))
                if "500_700" in dataset or "500to700" in dataset:
                    events = events[(lhe_mll > 500) & (lhe_mll <= 700)]
                    print("mass < 700 > 500 GeV: ", len(events))
                if "700_800" in dataset or "700to800" in dataset:
                    events = events[(lhe_mll > 700) & (lhe_mll <= 800)]
                    print("mass < 800 > 700 GeV: ", len(events))
                if "600_800" in dataset or "600to800" in dataset:
                    events = events[(lhe_mll > 600) & (lhe_mll <= 800)]
                    print("mass < 800 > 600 GeV: ", len(events))
                if "800_1000" in dataset or "800to1000" in dataset:
                    events = events[(lhe_mll > 800) & (lhe_mll <= 1000)]
                    print("mass < 1000 > 800 GeV: ", len(events))
                if "1000_1500" in dataset or "1000to1500" in dataset:
                    events = events[(lhe_mll > 1000) & (lhe_mll <= 1500)]
                    print("mass < 1500 > 1000 GeV: ", len(events))
                if "1000_3000" in dataset or "1000to3000" in dataset:
                    events = events[(lhe_mll > 1000) & (lhe_mll <= 3000)]
                    print("mass < 3000 > 1000 GeV: ", len(events))
                if "1500_2000" in dataset or "1500to2000" in dataset:
                    events = events[(lhe_mll > 1500) & (lhe_mll <= 2000)]
                    print("mass < 2000 > 1500 GeV: ", len(events))
                if "1500_inf" in dataset or "1500toinf" in dataset:
                    events = events[(lhe_mll > 1500)]
                    print("mass > 1500 GeV: ", len(events))
                if "2000_Inf" in dataset or "2000toInf" in dataset:
                    events = events[(lhe_mll > 2000)]
                    print("mass > 2000 GeV: ", len(events))

        events[f"mask_{dataset}"] = ak.ones_like(events.run) == 1.0
        events[f"weight_{dataset}"] = events.weight

        if subsamples != {}:
            for subsample in subsamples:
                subsample_val = subsamples[subsample]

                if isinstance(subsample_val, str):
                    subsample_mask = eval(subsample_val)
                    subsample_weight = 1.0
                elif (
                    isinstance(subsample_val, tuple) or isinstance(subsample_val, list)
                ) and len(subsample_val) == 2:
                    subsample_mask = eval(subsample_val[0])
                    subsample_weight = eval(subsample_val[1])
                else:
                    raise Exception(
                        "subsample value can either be a str (mask) or tuple/list of "
                        "len 2 (mask, weight)"
                    )

                events[f"mask_{dataset}_{subsample}"] = subsample_mask
                events[f"weight_{dataset}_{subsample}"] = (
                    events.weight * subsample_weight
                )

        for region in regions:
            regions[region]["mask"] = regions[region]["func"](events)

        # Fill histograms
        for dataset_name in results:
            results[dataset_name]["events"] = {}
            for region in regions:
                results[dataset_name]["events"][region] = {}
                for cwgt in check_weights:
                    events[cwgt] = check_weights[cwgt]["func"](events, dataset_name) if not isData else events.weight
                    mask = regions[region]["mask"] & events[f"mask_{dataset_name}"]
                    if len(events[mask]) == 0:
                        continue
                    for variable in variables:
                        if ("save_events" in variables[variable].keys()) and (variables[variable]["save_events"]):
                            results[dataset_name]["events"][region][variable] = events[variable][mask]
                            if "weight" not in results[dataset_name]["events"][region].keys():
                                results[dataset_name]["events"][region]["weight"] = events[f"weight_{dataset_name}"][mask]

                        if "axis" in variables[variable].keys():
                            if isinstance(variables[variable]["axis"], list):
                                var_names = [k.name for k in variables[variable]["axis"]]
                                vals = {
                                    var_name: events[var_name][mask] for var_name in var_names
                                }
                                results[dataset_name]["histos"][variable].fill(
                                    **vals,
                                    category=region,
                                    syst=variation,
                                    weight=events[f"weight_{dataset_name}"][mask],  # FIX: use subsample weight
                                )
                            else:
                                var_name = variables[variable]["axis"].name
                                results[dataset_name]["histos"][variable].fill(
                                    events[var_name][mask],
                                    category=region,
                                    syst=variation,
                                    weight=events[f"weight_{dataset_name}"][mask],  # FIX: use subsample weight
                                )

    gc.collect()
    return results


if __name__ == "__main__":
    chunks_readable = False
    new_chunks = read_chunks("chunks_job.pkl", readable=chunks_readable)
    print("N chunks to process", len(new_chunks))

    results = {}
    errors = []
    processed = []

    for i in range(len(new_chunks)):
        new_chunk = new_chunks[i]

        if new_chunk["result"] != {}:
            print(
                "Skip chunk",
                {k: v for k, v in new_chunk["data"].items() if k != "read_form"},
                "was already processed",
            )
            continue

        print(new_chunk["data"]["dataset"])

        try:
            new_chunks[i]["result"] = big_process(process=process, **new_chunk["data"])
            new_chunks[i]["error"] = ""
        except Exception as e:
            print("\n\nError for chunk", new_chunk, file=sys.stderr)
            nice_exception = "".join(tb.format_exception(None, e, e.__traceback__))
            print(nice_exception, file=sys.stderr)
            new_chunks[i]["result"] = {}
            new_chunks[i]["error"] = nice_exception

        print(f"Done {i+1}/{len(new_chunks)}")

    write_chunks(new_chunks, "results.pkl", readable=chunks_readable)
