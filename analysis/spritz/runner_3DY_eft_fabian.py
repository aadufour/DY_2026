import gc
import json
import sys
import traceback as tb
import awkward as ak
import numpy as np
import correctionlib
import hist
import vector
from copy import deepcopy
from spritz.framework.framework import big_process, get_analysis_dict, get_fw_path, read_chunks, write_chunks
import spritz.framework.variation as variation_module
from spritz.modules.basic_selections import LumiMask, lumi_mask, pass_flags, pass_trigger, pass_weightfilter
from spritz.modules.btag_sf import btag_sf
from spritz.modules.fake_leptons import get_fake_weights, fakes_reweight
from spritz.modules.h2erratum import h2erratum_reweight
from spritz.modules.ho_reweight import HO_reweight
from spritz.modules.jet_sel import clean_jet, jet_sel
from spritz.modules.jme import correct_jets_data, correct_jets_mc, jet_veto, remove_jets_HEM_issue
from spritz.modules.lepton_sel import create_lepton, lepton_sel
from spritz.modules.lepton_sf import lepton_sf
from spritz.modules.prefireweight import prefireweight
from spritz.modules.prompt_gen import prompt_gen_match_leptons
from spritz.modules.puid_sf import puid_sf
from spritz.modules.puweight import puweight_sf
from spritz.modules.rochester import correct_rochester, get_rochester, vary_rochester
from spritz.modules.run_assign import assign_run_period
from spritz.modules.theory_unc import theory_unc
from spritz.modules.trigger_sf import match_trigger_object, trigger_sf
from spritz.modules.tt_reweight import tt_reweight

vector.register_awkward()
print("awkward version", ak.__version__)

EFT_BATCH_SIZE = 2000

def cov_name(name_i, name_j):
    return f"cov_{name_i}_{name_j}"

def eft_points_key(dataset):
    return f"{dataset}__eft_points"

def eft_cov_key(dataset):
    return f"{dataset}__eft_covariances"

def fill_histos(histos, weights, events, regions, variables, variation, mask=None, names=[], batch_size=1):
    if mask is None:
        mask = ak.ones_like(events.run)==1.0

    for region in regions:
        mask_ = regions[region]["mask"] & mask
        nevents_sel = ak.sum(mask_)
        if nevents_sel == 0: 
            continue

        weights_sel = weights[mask_]

        for variable in histos:
            if isinstance(variables[variable]["axis"], list):
                var_names = [k.name for k in variables[variable]["axis"]]
            else:
                var_names = [variables[variable]["axis"].name]
            vals = {vn: events[vn][mask_] for vn in var_names}

            if names:
                for start in range(0, len(names), batch_size):
                    batch_idx = start // batch_size
                    b = len(names[start: start+batch_size])
                    vals_tiled = {vn: np.tile(val,b) for vn,val in vals.items()}
                    weights_tiled = np.array(weights_sel[:, start:start + b]).T.reshape(-1)
                    idx_tiled = np.repeat(np.arange(b, dtype=np.intp), nevents_sel)

                    histos[variable][batch_idx].fill(
                        **vals_tiled, subsample=idx_tiled, category=region, 
                        syst=variation, weight=weights_tiled,
                    )
            else:
                histos[variable].fill(
                    **vals, category=region, syst=variation, weight=weights_sel,
                )

##################################################

path_fw = get_fw_path()
with open("cfg.json") as file:
    txt = file.read()
    txt = txt.replace("RPLME_PATH_FW", path_fw)
    cfg = json.loads(txt)

ceval_puid = correctionlib.CorrectionSet.from_file(cfg["puidSF"])
ceval_btag = correctionlib.CorrectionSet.from_file(cfg["btagSF"])
ceval_btageff = correctionlib.CorrectionSet.from_file(cfg["btagEfficiency"])
ceval_puWeight = correctionlib.CorrectionSet.from_file(cfg["puWeights"])
ceval_lepton_sf = correctionlib.CorrectionSet.from_file(cfg["leptonSF"])
ceval_assign_run = correctionlib.CorrectionSet.from_file(cfg["run_to_era"])
rochester = get_rochester(cfg)

analysis_cfg = get_analysis_dict(sys.argv[1])
regions = deepcopy(analysis_cfg["regions"])
variables = deepcopy(analysis_cfg["variables"])

special_analysis_cfg = analysis_cfg.get("special_analysis_cfg", {})
reweight_fakes = special_analysis_cfg.get("reweight_fakes", False)
do_variations = special_analysis_cfg.get("do_variations", True)
do_rochester_stat_variations = special_analysis_cfg.get("do_rochester_stat_variations", False)
do_jet_variations = special_analysis_cfg.get("do_jet_variations", False)
do_theory_variations = special_analysis_cfg.get("do_theory_variations", False)
invert_one_isolation = special_analysis_cfg.get("invert_one_isolation", False)
invert_one_isolation_loose = special_analysis_cfg.get("invert_one_isolation_loose", False)
invert_one_isolation_control = special_analysis_cfg.get("invert_one_isolation_control", False)
invert_both_isolation = special_analysis_cfg.get("invert_both_isolation", False)

##################################################

def process(events, **kwargs):
    dataset = kwargs["dataset"]
    trigger_sel = kwargs.get("trigger_sel", "")
    isData = kwargs.get("is_data", False)
    era = kwargs.get("era", None)
    subsamples = kwargs.get("subsamples", {})
    eft_reweighting = kwargs.get("eft_reweighting", None)
    max_weight = kwargs.get("max_weight", None)
    genmatching_nlep = kwargs.get("genmatching_nlep", 2)
    ho_corrections = kwargs.get("ho_corrections", [])
    do_h2erratum_rwgt = kwargs.get("h2erratum_rwgt", False)
    do_top_pt_rwgt = kwargs.get("top_pt_rwgt", False)

    variations = variation_module.Variation()
    variations.register_variation([], "nom")

    if isData:
        events["weight"] = ak.ones_like(events.run)
        lumimask = LumiMask(cfg["lumiMask"])
        events = lumi_mask(events, lumimask)
    else:
        events["weight"] = events.genWeight
        events = pass_weightfilter(events, max_weight)
        events = events[events.pass_weightfilter]

    sumw = ak.sum(events.weight)
    nevents = ak.num(events.weight, axis=0)

    # LHE level selections
    if "DY" in dataset:
        outgoing_mask = (events.LHEPart.status == 1)
        lepton_mask = (abs(events.LHEPart.pdgId) == 13)
        lhe_leptons = events.LHEPart[outgoing_mask & lepton_mask]

        if ak.all(ak.num(lhe_leptons) == 2):
            lhe_mll = (lhe_leptons[:, 0] + lhe_leptons[:, 1]).mass
            if "M-50to100" in dataset:
                events = events[(50 < lhe_mll) & (lhe_mll < 100)]
            if "50_120" in dataset:
                events = events[(lhe_mll >= 50) & (lhe_mll <= 120)]
            if "120_200" in dataset:
                events = events[(lhe_mll > 120) & (lhe_mll <= 200)]
            if "200_400" in dataset:
                events = events[(lhe_mll > 200) & (lhe_mll <= 400)]
            if "400_600" in dataset:
                events = events[(lhe_mll > 400) & (lhe_mll <= 600)]
            if "600_800" in dataset:
                events = events[(lhe_mll > 600) & (lhe_mll <= 800)]
            if "800_1000" in dataset:
                events = events[(lhe_mll > 800) & (lhe_mll <= 1000)]
            if "1000_1500" in dataset:
                events = events[(lhe_mll > 1000) & (lhe_mll <= 1500)]
            if "1500_inf" in dataset:
                events = events[(lhe_mll > 1500)]
            if "1000_3000" in dataset:
                events = events[(lhe_mll > 1000)]

    # pass trigger and flags
    events = assign_run_period(events, isData, cfg, ceval_assign_run)
    events = pass_trigger(events, cfg["era"])
    events = pass_flags(events, cfg["flags"])
    events = events[events.pass_flags & events.pass_trigger]

    if isData: # each data DataSet has its own trigger_sel
        events = events[eval(trigger_sel)]

    # Require at least one good PV
    events = events[events.PV.npvsGood > 0]

    # Lepton preselection
    events = create_lepton(events)
    events = lepton_sel(events, cfg)
    events["Lepton"] = events.Lepton[events.Lepton.isLoose]

    # Jet preselection
    events = jet_sel(events, cfg) # tight ID, eta < 2.5 (2017,2018) or 2.4 (2016)
    events = clean_jet(events)
    events = remove_jets_HEM_issue(events, cfg)
    events = jet_veto(events, cfg)

    # Gen matching
    if not isData:
        events = prompt_gen_match_leptons(events)

    # Trigger matching
    events = match_trigger_object(events, cfg)

    # Muon Rochester corrections
    events, variations = vary_rochester(events, variations, isData, rochester, do_rochester_stat_variations)
    events, variations = correct_rochester(events, variations, isData, rochester)

    # Jet energy scale and resolution corrections
    if not isData:
        events, variations = correct_jets_mc(events, variations, cfg, run_variations=do_jet_variations)
    else:
        events, variations = correct_jets_data(events, variations, cfg, era)

    # Apply a skim!
    lepton_sort = ak.argsort(events.Lepton.pt, ascending=False, axis=1)
    events["Lepton"] = events.Lepton[lepton_sort]
    events = events[ak.num(events.Lepton) >= 2]
    events = events[events.Lepton[:, 0].pt >= 24]
    events = events[events.Lepton[:, 1].pt >= 10]

    if len(events) == 0:
        return {}

    # Fake lepton reweighting
    if reweight_fakes:
        variations, fakes_param = get_fake_weights(variations, cfg)

    ##################################################

    # Scale factors
    if not isData:
        # puWeight SF
        events, variations = puweight_sf(events, variations, ceval_puWeight, cfg)

        # prefire weight
        events, variations = prefireweight(events, variations)

        # lepton SF
        events, variations = lepton_sf(events, variations, ceval_lepton_sf, cfg)

        # trigger SF
        events, variations = trigger_sf(events, variations, ceval_lepton_sf, cfg)

        # puId SF
        events, variations = puid_sf(events, variations, ceval_puid, cfg)

        # btag SF
        events, variations = btag_sf(events, variations, ceval_btag, ceval_btageff, cfg, dataset)

        # Higher-order corrections
        for ho_corr in ho_corrections:
            events, variations = HO_reweight(events, variations, ho_corr)

        # H2ErratumFix
        if do_h2erratum_rwgt:
            events, variations = h2erratum_reweight(events, variations, cfg, dataset)

        # Top pT reweighting
        if do_top_pt_rwgt:
            events, variations = tt_reweight(events, variations)

        # Theory unc.
        if do_theory_variations:
            events, variations = theory_unc(events, variations)

    ##################################################

    # Set up results
    if not do_variations:
        variations.variations_dict = {
            k: v for k, v in variations.variations_dict.items() if k == "nom"
        }

    region_axis = hist.axis.StrCategory([region for region in regions], name="category")
    variation_axis = hist.axis.StrCategory(
        sorted(list(variations.get_variations_all())), name="syst")
    default_axis = [region_axis, variation_axis]

    # eft_reweighting names: points + covariance terms
    results = {}
    if subsamples != {} or eft_reweighting is not None:
        for subsample in subsamples:
            results[f"{dataset}_{subsample}"] = {
                "sumw": sumw, "nevents": nevents, "events": 0, "histos": 0,
            }

        if eft_reweighting is not None:
            results[eft_points_key(dataset)] = {
                "sumw": sumw, "nevents": nevents, "events": 0, "histos": 0,
                "eft_batch_size": EFT_BATCH_SIZE, 
                "eft_names": [name for name in eft_reweighting["points"]]
            }
            results[eft_cov_key(dataset)] = {
                "sumw": sumw, "nevents": nevents, "events": 0, "histos": 0,
                "eft_batch_size": EFT_BATCH_SIZE, 
                "eft_names": [cov_name(a,b) for a,b in eft_reweighting.get("covariance_pairs", [])]
            }
    else:
        results[dataset] = {
            "sumw": sumw, "nevents": nevents, "events": 0, "histos": 0
        }

    for dataset_name in results:
        _events, histos = {}, {}
        is_eft_combined = dataset_name == eft_points_key(dataset)
        is_eft_covariance = dataset_name == eft_cov_key(dataset)
        for variable in variables:
            _events[variable] = ak.Array([])

            if "axis" in variables[variable]:
                axis_def = variables[variable]["axis"]
                axes = axis_def if isinstance(axis_def, list) else [axis_def]

                if dataset_name in [eft_points_key(dataset), eft_cov_key(dataset)]:
                    # a list of ~eft_names/EFT_BATCH_SIZE histograms
                    batch_histos = []
                    eft_names = results[dataset_name]["eft_names"]
                    if dataset_name == eft_points_key(dataset):
                        variation_axis_ = variation_axis
                    else:
                        variation_axis_ = hist.axis.StrCategory(["nom"], name="syst")

                    for b in range(1 + len(eft_names)//EFT_BATCH_SIZE):
                        batch_len = min(EFT_BATCH_SIZE, len(eft_names) - b*EFT_BATCH_SIZE)
                        batch_axis = hist.axis.IntCategory(list(range(batch_len)), name="subsample")
                        batch_histos.append(
                            hist.Hist(
                                *axes, *[batch_axis, region_axis, variation_axis_], hist.storage.Weight())
                        )
                    histos[variable] = batch_histos
                else:
                    histos[variable] = hist.Hist(*axes, *default_axis, hist.storage.Weight())

        results[dataset_name]["histos"] = histos
        results[dataset_name]["events"] = _events

    ##################################################

    # Loop over variations
    print("Doing variations")
    originalEvents = ak.copy(events)

    for variation in sorted(variations.get_variations_all()):
        print(variation)
        events = ak.copy(originalEvents)

        for switch in variations.get_variation_subs(variation):
            if len(switch) == 2:
                variation_dest, variation_source = switch
                events[variation_dest] = events[variation_source]

        # resort Leptons
        lepton_sort = ak.argsort(events.Lepton.pt, ascending=False, axis=1)
        events["Lepton"] = events.Lepton[lepton_sort]

        # Define categories
        events["mm"] = (events.Lepton[:, 0].pdgId * events.Lepton[:, 1].pdgId) == -13 * 13
        events["mm_ss"] = (events.Lepton[:, 0].pdgId * events.Lepton[:, 1].pdgId) == 13 * 13
        events = events[events.mm | events.mm_ss]

        # Cut on pt of two leading leptons
        ptcut = (events.Lepton[:, 0].pt > 29) & (events.Lepton[:, 1].pt > 15)
        events = events[ptcut]

        # tight ID requirement
        muWP = cfg["leptonsWP"]["muWP"]
        lTight = events.Lepton[:, 0][f"isTightMuon_{muWP}"] & events.Lepton[:, 1][f"isTightMuon_{muWP}"]
        events = events[lTight]

        # isolation requirement
        l1Iso = events.Lepton[:, 0]["isTightMuon_RelIso"]
        l1IsoLoose = events.Lepton[:, 0]["isTightMuon_RelIso_loose"]
        l2Iso = events.Lepton[:, 1]["isTightMuon_RelIso"]
        l2IsoLoose = events.Lepton[:, 1]["isTightMuon_RelIso_loose"]

        if invert_one_isolation:
            lIso = (l1Iso & ~l2Iso) | (~l1Iso & l2Iso)
        elif invert_one_isolation_loose:
            lIso = (l1Iso & ~l2IsoLoose) | (~l1IsoLoose & l2Iso)
        elif invert_one_isolation_control:
            lIso = (l1Iso & l2IsoLoose & ~l2Iso) | (l1IsoLoose & ~l1Iso & l2Iso)
        elif invert_both_isolation:
            lIso = ~l1Iso & ~l2Iso
        else:
            lIso = l1Iso & l2Iso
        events = events[lIso]

        # third lepton veto
        events["Lepton"] = events.Lepton[events.Lepton.pt >= 10]
        l3Veto = ak.num(events.Lepton) < 3
        events = events[l3Veto]

        # prompt gen matching
        if not isData:
            events["prompt_gen_match_1l"] = (
                events.Lepton[:, 0].promptgenmatched | events.Lepton[:, 1].promptgenmatched
            )
            events["prompt_gen_match_2l"] = (
                events.Lepton[:, 0].promptgenmatched & events.Lepton[:, 1].promptgenmatched
            )
            if genmatching_nlep == 1:
                events = events[events.prompt_gen_match_1l]
            elif genmatching_nlep > 1:
                events = events[events.prompt_gen_match_2l]

        if len(events) == 0:
            continue

        # Jet selection and b-tag veto
        bveto_pt = cfg["bVeto"]["pt"]
        bveto_wp = cfg["bTag"][f"btag{cfg["bVeto"]["wp"]}"]

        events["Jet"] = events.Jet[events.Jet.pt >= bveto_pt]
        events["LowPtJet"] = events.Jet[~events.Jet.pass_highPt]
        events["Jet"] = events.Jet[events.Jet.pass_puId | events.Jet.pass_highPt]

        btagged = (events.Jet.btagDeepFlavB >= bveto_wp)
        events["BJet"] = events.Jet[btagged]
        events["bveto"] = ak.num(events.BJet) == 0
        events["btag"] = ak.num(events.BJet) >= 1

        ##################################################

        # Load all SFs
        if not isData:
            events["RecoSF"] = events.Lepton[:, 0].RecoSF * events.Lepton[:, 1].RecoSF
            events["IdSF"] = events.Lepton[:, 0].IdSF * events.Lepton[:, 1].IdSF
            events["IsoSF"] = events.Lepton[:, 0].IsoSF * events.Lepton[:, 1].IsoSF
            events["btagSF"] = ak.prod(events.Jet.btagSF, axis=-1)
            events["puidSF"] = ak.prod(events.LowPtJet.puidSF, axis=-1)

            events["weight"] = (
                events.weight
                * events.puWeight
                * events.prefireWeight
                * events.RecoSF
                * events.IdSF
                * events.IsoSF
                * events.TriggerSF
                * events.puidSF
                * events.btagSF
            )

            for ho_corr in ho_corrections:
                events["weight"] = events.weight * events[ho_corr["name"]]
            if do_h2erratum_rwgt:
                events["weight"] = events.weight * events.H2ErratumWeight
            if do_top_pt_rwgt:
                events["weight"] = events.weight * events.topPtWeight

        # Fake lepton reweighting (only in the same-sign region)
        if reweight_fakes:
            events["fakesRW"] = fakes_reweight(events, variation, fakes_param)
            events["fakesRW"] = ak.where(events.mm_ss, events.fakesRW, ak.ones_like(events.weight))
            events["weight"] = events.weight * events.fakesRW

        ##################################################

        # Variable and Regions definitions
        for variable in variables:
            if "func" in variables[variable]:
                events[variable] = variables[variable]["func"](events)

        for region in regions:
            regions[region]["mask"] = regions[region]["func"](events)

        # Fill histograms
        # datasets without subsamples or eft_reweighting
        if dataset in results:
            fill_histos(
                results[dataset]["histos"], events.weight, events, regions, variables, variation
            )

        # datasets with subsamples
        if subsamples != {}:
            n_subsamples = len(subsamples)
            for i, subsample in enumerate(subsamples):
                subsample_val = subsamples[subsample]
                if isinstance(subsample_val, str):
                    mask_expr, weight_expr = subsample_val, None
                elif isinstance(subsample_val, (tuple, list)) and len(subsample_val) == 2:
                    mask_expr, weight_expr = subsample_val
                else:
                    raise Exception("subsample value must be either a str (mask) or a "
                        "(mask, weight) tuple/list of length 2"
                    )

                subsample_mask = eval(mask_expr)
                if weight_expr is None:
                    subsample_weight = events.weight
                else:
                    subsample_weight = events.weight * eval(weight_expr)

                fill_histos(
                    results[f"{dataset}_{subsample}"]["histos"], subsample_weight, 
                    events, regions, variables, variation, mask=subsample_mask
                )

        # datasets with eft_reweighting
        if eft_reweighting is not None:
            points = eft_reweighting["points"]
            covariance_pairs = eft_reweighting.get("covariance_pairs", [])
            rw_weight = events[eft_reweighting["weight_branch"]]

            point_names = list(points.keys())
            idx = np.array(list(points.values()), dtype=np.intp)
            point_weights = events.weight[:,None] * rw_weight[:,idx]

            fill_histos(
                results[eft_points_key(dataset)]["histos"], point_weights, 
                events, regions, variables, variation, names=point_names, batch_size=EFT_BATCH_SIZE
            )

            # Covariance terms are only computed once (nominal)
            if variation == "nom" and covariance_pairs:
                cov_names = [cov_name(a, b) for a, b in covariance_pairs]
                idx_i = np.array([points[a] for a, b in covariance_pairs], dtype=np.intp)
                idx_j = np.array([points[b] for a, b in covariance_pairs], dtype=np.intp)
                cov_weights = (events.weight**2)[:,None] * rw_weight[:,idx_i] * rw_weight[:,idx_j]

                fill_histos(
                    results[eft_cov_key(dataset)]["histos"], cov_weights, 
                    events, regions, variables, variation, names=cov_names, batch_size=EFT_BATCH_SIZE
                )


    gc.collect()
    return results

##################################################

def chunk_str(chunk):
    drop_keys = ["read_form", "ho_corrections", "eft_reweighting"]
    return str({k:v for k,v in chunk["data"].items() if not k in drop_keys})

if __name__ == "__main__":
    new_chunks = read_chunks("chunks_job.pkl")
    print("N chunks to process", len(new_chunks))

    for i in range(len(new_chunks)):
        new_chunk = new_chunks[i]

        if new_chunk["result"] != {}:
            print(f"Skip chunk {chunk_str(new_chunk)}, was already processed")
            continue

        print(f"chunk {i+1}/{len(new_chunks)}: {new_chunk['data']['dataset']}")

        try:
            new_chunks[i]["result"] = big_process(process=process, **new_chunk["data"])
            new_chunks[i]["error"] = ""
        except Exception as e:
            nice_exception = "".join(tb.format_exception(None, e, e.__traceback__))
            print(f"\n\nError for chunk {chunk_str(new_chunk)}", file=sys.stderr)
            print(nice_exception, file=sys.stderr)
            new_chunks[i]["result"] = {}
            new_chunks[i]["error"] = nice_exception

        print(f"Done {i+1}/{len(new_chunks)}\n")

    write_chunks(new_chunks, "results.pkl")

