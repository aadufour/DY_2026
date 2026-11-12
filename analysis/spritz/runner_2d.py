# runner_2d.py: runner_pdfcut.py (PDF CUT) + batched filling of EFT subsamples (see EFT BATCH).
# EFT subsamples given as (mask_expr, weight_expr) are filled into ONE histogram per variable with an extra
# IntCategory "subsample" axis (megahisto idea, spritz_giacomo runner_3DY_eft_full_morphing_megahisto.py),
# in batches of EFT_FILL_BATCH subsamples, and unpacked at the end into the usual per-subsample results
# entries ({dataset}_{subsample}): the output format is identical to runner.py, so merge/postproc/cards
# are unchanged. Needed for the 406 points (incl. 351 cross terms) of the 2D production.
# ruff: noqa: E501
"""
runner_combined.py — Fabian's runner_3DY.py + EFT subsample support
=====================================================================
Base: spritz_fabian/src/spritz/runners/runner_3DY.py (μμ channel, 2018 UL)
EFT additions (4 surgical changes, marked ### EFT ###):
  1. EFT filter: keep only events with exactly N LHEReweightingWeights
  2. LHE mll filter: avoid double-counting between mll-binned EFT samples
  3. Subsample handling: parse (mask_expr, weight_expr) tuples for EFT weights
  4. Histogram fill: use per-subsample EFT weight instead of global events["weight"]

All of Fabian's logic (lepton sel, SFs, Rochester, fakes, MET) is unchanged.
"""
import gc

EFT_FILL_BATCH = 50  # subsamples per fill call (memory: n_sel * batch * ~5 arrays)
import json
import sys
import traceback as tb

import awkward as ak
import correctionlib
import hist
import numpy as np
import spritz.framework.variation as variation_module
import uproot
import vector
from copy import deepcopy
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
    pass_weightfilter,
)
from spritz.modules.fake_leptons import reweightFakeLep
from spritz.modules.jme import correct_met
from spritz.modules.lepton_sel import createLepton, leptonSel
from spritz.modules.lepton_sf import lepton_sf
from spritz.modules.prompt_gen import prompt_gen_match_leptons
from spritz.modules.puweight import puweight_sf
from spritz.modules.rochester import (
    correctRochester,
    getRochester,
    varyRochester,
)
from spritz.modules.run_assign import assign_run_period
from spritz.modules.theory_unc import theory_unc
from spritz.modules.trigger_sf import (
    match_trigger_object,
    trigger_sf,
)
from spritz.modules.tt_reweight import tt_reweight

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
ceval_met = correctionlib.CorrectionSet.from_file(cfg["met"])

rochester = getRochester(cfg)

analysis_path = sys.argv[1]
analysis_cfg = get_analysis_dict(analysis_path)
special_analysis_cfg = analysis_cfg["special_analysis_cfg"]


def process(events, **kwargs):
    dataset = kwargs["dataset"]
    trigger_sel = kwargs.get("trigger_sel", "")
    isData = kwargs.get("is_data", False)
    era = kwargs.get("era", None)
    subsamples = kwargs.get("subsamples", {})
    special_weight = eval(kwargs.get("weight", "1.0"))

    variations = variation_module.Variation()
    variations.register_variation([], "nom")

    if isData:
        events["weight"] = ak.ones_like(events.run)
    else:
        events["weight"] = events.genWeight

    ### EFT 1: filter events with exactly N LHEReweightingWeights ###
    if not isData and "EFT" in kwargs and kwargs["EFT"]:
        neft_rwgts = kwargs["EFT"]
        if "LHEReweightingWeight" in ak.fields(events):
            events = events[ak.num(events.LHEReweightingWeight) == neft_rwgts]

    if isData:
        lumimask = LumiMask(cfg["lumiMask"])
        events = lumi_mask(events, lumimask)
    else:
        events = pass_weightfilter(events, kwargs.get("max_weight", None))
        events = events[events.pass_weightfilter]

    sumw = ak.sum(events.weight)
    nevents = ak.num(events.weight, axis=0)

    if special_weight != 1.0:
        print(f"Using special weight for {dataset}: {special_weight}")

    events["weight"] = events.weight * special_weight

    # pass trigger and flags
    events = assign_run_period(events, isData, cfg, ceval_assign_run)
    events = pass_trigger(events, cfg["era"])
    events = pass_flags(events, cfg["flags"])

    events = events[events.pass_flags & events.pass_trigger]

    if isData:
        events = events[eval(trigger_sel)]

    events = createLepton(events)
    events = leptonSel(events, cfg)
    events["Lepton"] = events.Lepton[events.Lepton.isLoose]

    events = events[ak.num(events.Lepton) >= 2]
    events = events[events.Lepton[:, 0].pt >= 24]
    events = events[events.Lepton[:, 1].pt >= 10]

    # LHE mll filter for DY MiNNLO (Fabian's original)
    if dataset in ["DYmm_M-50to100", "DYmm_M-50"]:
        outgoing_mask = (events.LHEPart.status == 1)
        lepton_mask = (abs(events.LHEPart.pdgId) == 13)
        lhe_leptons = events.LHEPart[outgoing_mask & lepton_mask]
        assert ak.all(ak.num(lhe_leptons) == 2)
        lhe_mll = (lhe_leptons[:, 0] + lhe_leptons[:, 1]).mass
        events = events[(50 < lhe_mll) & (lhe_mll < 100)]

    ### EFT 2: LHE mll filter for mll-binned EFT samples (avoid double-counting) ###
    if not isData and "DYSMEFTsim" in dataset:
        lhe_lep_mask = (events.LHEPart.status == 1) & (
            (abs(events.LHEPart.pdgId) == 11)
            | (abs(events.LHEPart.pdgId) == 13)
            | (abs(events.LHEPart.pdgId) == 15)
        )
        lhe_leps = events.LHEPart[lhe_lep_mask]
        if ak.all(ak.num(lhe_leps) == 2):
            lhe_mll = (lhe_leps[:, 0] + lhe_leps[:, 1]).mass
            if "50_120" in dataset:
                events = events[(lhe_mll >= 50) & (lhe_mll <= 120)]
            elif "120_200" in dataset:
                events = events[(lhe_mll > 120) & (lhe_mll <= 200)]
            elif "200_400" in dataset:
                events = events[(lhe_mll > 200) & (lhe_mll <= 400)]
            elif "400_600" in dataset:
                events = events[(lhe_mll > 400) & (lhe_mll <= 600)]
            elif "600_800" in dataset:
                events = events[(lhe_mll > 600) & (lhe_mll <= 800)]
            elif "800_1000" in dataset:
                events = events[(lhe_mll > 800) & (lhe_mll <= 1000)]
            elif "1000_3000" in dataset:
                events = events[(lhe_mll > 1000) & (lhe_mll <= 3000)]

    if not isData:
        events = prompt_gen_match_leptons(events)

    events = events[events.PV.npvsGood > 0]

    if kwargs.get("top_pt_rwgt", False):
        events, variations = tt_reweight(events, variations)
    else:
        events['topPtWeight'] = ak.ones_like(events.weight)

    if special_analysis_cfg.get("do_rochester_variations", False):
        events, variations = varyRochester(events, variations, isData, rochester)
    events, variations = correctRochester(events, variations, isData, rochester)

    events = match_trigger_object(events, cfg)
    events = correct_met(events, ceval_met, isData)

    if special_analysis_cfg.get("reweight_fakes", False):
        events, variations = reweightFakeLep(events, variations)
    else:
        events["fakeLepWeight"] = ak.ones_like(events.weight)

    if not isData:
        events, variations = puweight_sf(events, variations, ceval_puWeight, cfg)
        events, variations = trigger_sf(events, variations, ceval_lepton_sf, cfg)
        events, variations = lepton_sf(events, variations, ceval_lepton_sf, cfg)

        if "L1PreFiringWeight" in ak.fields(events):
            events["prefireWeight"] = events.L1PreFiringWeight.Nom
            events["prefireWeight_up"] = events.L1PreFiringWeight.Up
            events["prefireWeight_down"] = events.L1PreFiringWeight.Dn
            events["prefireWeight_before"] = ak.ones_like(events.L1PreFiringWeight.Nom)

            for tag in ["up", "down", "before"]:
                variations.register_variation(
                    columns=["prefireWeight"],
                    variation_name=f"prefireWeight_{tag}",
                    format_rule=lambda _, var_name: var_name,
                )
        else:
            events["prefireWeight"] = ak.ones_like(events.weight)

        if special_analysis_cfg.get("do_theory_variations", False):
            # PDF CUT: unphysical Hessian member weights (|LHEPdfWeight[i]/LHEPdfWeight[0]| > 100)
            # are reset to the nominal member; the event and its nominal weight are kept.
            if "LHEPdfWeight" in ak.fields(events):
                W = events.LHEPdfWeight
                bad = abs(W / W[:, 0]) > 100
                print("PDF cut: resetting", int(ak.sum(bad)), "member weights in", int(ak.sum(ak.any(bad, axis=1))), "events")
                events["LHEPdfWeight"] = ak.where(bad, W[:, 0], W)
            events, variations = theory_unc(events, variations)

    regions = deepcopy(analysis_cfg["regions"])
    variables = deepcopy(analysis_cfg["variables"])

    if not special_analysis_cfg.get("do_variations", False):
        variations.variations_dict = {
            k: v for k, v in variations.variations_dict.items() if k == "nom"
        }

    default_axis = [
        hist.axis.StrCategory(
            [region for region in regions],
            name="category",
        ),
        hist.axis.StrCategory(
            sorted(list(variations.get_variations_all())),
            name="syst"
        )
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

    ### EFT BATCH 1: one histogram per variable for all (mask, weight) subsamples ###
    eft_names = [s for s in subsamples if isinstance(subsamples[s], (tuple, list)) and len(subsamples[s]) == 2]
    eft_hists = {}
    for variable in variables:
        if eft_names and "axis" in variables[variable]:
            ax = variables[variable]["axis"]
            eft_hists[variable] = hist.Hist(
                *(ax if isinstance(ax, list) else [ax]),
                hist.axis.IntCategory(list(range(len(eft_names))), name="subsample"),
                *default_axis,
                hist.storage.Weight(),
            )
    eft_keys = {f"{dataset}_{s}" for s in eft_names}

    for dataset_name in results:
        if dataset_name in eft_keys:
            continue
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
        muWP = cfg["leptonsWP"]["muWP"]

        comb = (
            events.Lepton[:, 0]["isTightMuon_" + muWP] & events.Lepton[:, 1]["isTightMuon_" + muWP]
        )
        if special_analysis_cfg.get("invert_one_isolation", False):
            comb = comb & (
                (events.Lepton[:, 0]["isTightMuon_RelIso"] & ~events.Lepton[:, 1]["isTightMuon_RelIso"])
                | (~events.Lepton[:, 0]["isTightMuon_RelIso"] & events.Lepton[:, 1]["isTightMuon_RelIso"])
            )
        elif special_analysis_cfg.get("invert_one_isolation_loose", False):
            comb = comb & (
                (events.Lepton[:, 0]["isTightMuon_RelIso"] & ~events.Lepton[:, 1]["isTightMuon_RelIso_loose"])
                | (~events.Lepton[:, 0]["isTightMuon_RelIso_loose"] & events.Lepton[:, 1]["isTightMuon_RelIso"])
            )
        elif special_analysis_cfg.get("invert_one_isolation_control", False):
            comb = comb & (
                (events.Lepton[:, 0]["isTightMuon_RelIso"] & ~events.Lepton[:, 1]["isTightMuon_RelIso"] & events.Lepton[:, 1]["isTightMuon_RelIso_loose"])
                | (~events.Lepton[:, 0]["isTightMuon_RelIso"] & events.Lepton[:, 0]["isTightMuon_RelIso_loose"] & events.Lepton[:, 1]["isTightMuon_RelIso"])
            )
        elif special_analysis_cfg.get("invert_both_isolation", False):
            comb = comb & (
                ~events.Lepton[:, 0]["isTightMuon_RelIso"] & ~events.Lepton[:, 1]["isTightMuon_RelIso"]
            )
        else:
            comb = comb & (
                events.Lepton[:, 0]["isTightMuon_RelIso"] & events.Lepton[:, 1]["isTightMuon_RelIso"]
            )
        events["l2Tight"] = ak.copy(comb)
        events = events[events.l2Tight]

        if len(events) == 0:
            continue

        events["mm"] = (
            events.Lepton[:, 0].pdgId * events.Lepton[:, 1].pdgId
        ) == -13 * 13
        events["mm_ss"] = (
            events.Lepton[:, 0].pdgId * events.Lepton[:, 1].pdgId
        ) == 13 * 13

        if not isData and not special_analysis_cfg.get("skip_genmatching", False):
            events["prompt_gen_match_2l"] = (
                events.Lepton[:, 0].promptgenmatched & events.Lepton[:, 1].promptgenmatched
            )
            events = events[events.prompt_gen_match_2l]

        leptoncut = (events.mm | events.mm_ss)

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

        leptoncut = leptoncut & (events.Lepton[:, 0].pt > 29) & (events.Lepton[:, 1].pt > 15)
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

        events["fakeLepWeight"] = ak.where(
            events.mm_ss, events.fakeLepWeight, ak.ones_like(events.weight)
        )
        events["weight"] = events.weight * events.fakeLepWeight

        for variable in variables:
            if "func" in variables[variable]:
                events[variable] = variables[variable]["func"](events)

        events[dataset] = ak.ones_like(events.run) == 1.0

        ### EFT 3: subsample handling — support (mask_expr, weight_expr) tuples ###
        if subsamples != {}:
            for subsample in subsamples:
                val = subsamples[subsample]
                if isinstance(val, (tuple, list)) and len(val) == 2:
                    continue  # EFT BATCH: weights built below as one matrix
                else:
                    # Standard subsample: mask only
                    events[f"{dataset}_{subsample}"] = eval(val)
                    events[f"weight_{dataset}_{subsample}"] = events.weight

        for region in regions:
            regions[region]["mask"] = regions[region]["func"](events)
            if isData and regions[region].get("data_mll_max") is not None:
                regions[region]["mask"] = regions[region]["mask"] & (events.mll < regions[region]["data_mll_max"])

        ### EFT BATCH 2: (n_events, n_subsamples) weight matrix, mask folded into the weight ###
        if eft_names:
            w_nom = ak.to_numpy(events.weight)
            wmat = np.empty((len(w_nom), len(eft_names)))
            for k, s_ in enumerate(eft_names):
                m_expr, w_expr = subsamples[s_]
                wmat[:, k] = w_nom * ak.to_numpy(eval(w_expr)) * ak.to_numpy(eval(m_expr))
            for region in regions:
                mask_np = ak.to_numpy(regions[region]["mask"])
                n_sel = int(mask_np.sum())
                if n_sel == 0:
                    continue
                w_sel = wmat[mask_np]
                for variable, h in eft_hists.items():
                    ax = variables[variable]["axis"]
                    vnames = [a.name for a in ax] if isinstance(ax, list) else [ax.name]
                    vals = {vn: ak.to_numpy(events[vn])[mask_np] for vn in vnames}
                    for start in range(0, len(eft_names), EFT_FILL_BATCH):
                        b = min(EFT_FILL_BATCH, len(eft_names) - start)
                        h.fill(
                            **{vn: np.tile(v, b) for vn, v in vals.items()},
                            subsample=np.repeat(np.arange(start, start + b), n_sel),
                            category=region,
                            syst=variation,
                            weight=w_sel[:, start:start + b].T.reshape(-1),
                        )
            del wmat

        # Fill histograms
        for dataset_name in results:
            if dataset_name in eft_keys:
                continue
            for region in regions:
                mask = regions[region]["mask"] & events[dataset_name]

                if len(events[mask]) == 0:
                    continue

                ### EFT 4: use per-subsample weight (EFT weight for EFT samples, ###
                ###         regular weight for backgrounds)                        ###
                wkey = f"weight_{dataset_name}"
                fill_weight = (
                    events[wkey][mask]
                    if wkey in ak.fields(events)
                    else events["weight"][mask]
                )

                for variable in results[dataset_name]["histos"]:
                    if isinstance(variables[variable]["axis"], list):
                        var_names = [k.name for k in variables[variable]["axis"]]
                        vals = {
                            var_name: events[var_name][mask] for var_name in var_names
                        }
                        results[dataset_name]["histos"][variable].fill(
                            **vals,
                            category=region,
                            syst=variation,
                            weight=fill_weight,
                        )
                    else:
                        var_name = variables[variable]["axis"].name
                        results[dataset_name]["histos"][variable].fill(
                            events[var_name][mask],
                            category=region,
                            syst=variation,
                            weight=fill_weight,
                        )

    ### EFT BATCH 3: unpack into the usual per-subsample entries ###
    for k, s_ in enumerate(eft_names):
        results[f"{dataset}_{s_}"]["histos"] = {v: h[{"subsample": k}] for v, h in eft_hists.items()}
        results[f"{dataset}_{s_}"]["events"] = {v: ak.Array([]) for v in variables}
    del eft_hists

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
