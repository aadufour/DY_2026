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
from fnmatch import fnmatch
from spritz.framework.framework import big_process, get_analysis_dict, get_fw_path, read_chunks, write_chunks, add_dict, EFTHistogram
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


def cov_name(name_i, name_j):
    return f"cov_{name_i}_{name_j}"

def eft_points_key(dataset):
    return f"{dataset}__eft_points"

def eft_cov_key(dataset):
    return f"{dataset}__eft_covariances"

def fill_histos(histos, weights, events, regions, variables, variation, mask=None):
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

            histos[variable].fill(
                **vals, category=region, syst=variation, weight=weights_sel,
            )

def to_numpy(array, fill_value):
    """Plain float64 numpy array: weights and variables can be option-typed
    (e.g. after ak.where on ak.mask-ed SFs), which ak.to_numpy turns into a
    masked array that numpy matmul/reduceat don't handle"""
    return np.ma.filled(ak.to_numpy(array), fill_value).astype(np.float64)

def eft_point_weights(events, eft_reweighting):
    """(n_events, n_points) numpy matrix of event weight x reweighting
    weight, one column per eft_reweighting["points"] entry."""
    branch = eft_reweighting["weight_branch"]
    rw_weight = to_numpy(ak.to_regular(events[branch]), 0.0)
    n_weights = eft_reweighting.get("n_weights")
    if n_weights is not None and rw_weight.shape[1] != n_weights:
        raise ValueError(
            f"{branch} has {rw_weight.shape[1]} columns, expected {n_weights}: this "
            "file's reweight card does not match the point -> column mapping"
        )
    idx = np.array(list(eft_reweighting["points"].values()), dtype=np.intp)
    weight = to_numpy(events.weight, 0.0)
    return weight[:, None] * rw_weight[:, idx]

def fill_eft(histos, point_weights, events, regions, variables, variation, cov_pairs=None):
    """Fill EFTHistograms: points (histos[variable] has one name per
    column of point_weights) or, with cov_pairs = (idx_i, idx_j), the
    covariance terms sum_e w_ei * w_ej between those columns."""
    for region in regions:
        mask = np.ma.filled(ak.to_numpy(regions[region]["mask"]), False)
        if not mask.any():
            continue
        weights_sel = point_weights[mask]

        for variable, h in histos.items():
            if region not in h.regions:
                continue
            vals = {ax.name: to_numpy(events[ax.name], np.nan)[mask] for ax in h.axes}
            if cov_pairs is None:
                h.fill(vals, weights_sel, region, variation)
            else:
                h.fill_covariance(vals, weights_sel, *cov_pairs, region, variation)

# Columns that only enter the event weight (on top of each dataset's
# ho_corrections columns): a variation replacing nothing else selects the
# same events as the nominal one, so it reuses the nominal selection and
# only recomputes the weight
WEIGHT_COLUMNS = {
    "weight", "puWeight", "prefireWeight", "TriggerSF", "topPtWeight", "H2ErratumWeight",
    ("Lepton", "RecoSF"), ("Lepton", "IdSF"), ("Lepton", "IsoSF"),
    ("Jet", "btagSF"), ("Jet", "puidSF"),
}

def declared_kinematic(analysis_cfg):
    """variation name -> "kinematic" flag of the nuisance (or the "_before"
    variation of the correction) it belongs to, for the entries that set one"""
    flags = {}
    for nuis in analysis_cfg.get("nuisances", {}).values():
        if "kinematic" not in nuis:
            continue
        name = nuis.get("name")
        tags = [f"{name}_up", f"{name}_down"] if nuis.get("kind") in ("weight", "suffix") else []
        tags += [nuis[k] for k in ("tag", "tag_up", "tag_down") if k in nuis]
        tags += [v["tag"] for v in nuis.get("variations", [])]
        for tag in tags:
            flags[tag] = flags.get(tag, False) or nuis["kinematic"]
    for corr in analysis_cfg.get("corrections", {}).values():
        if "kinematic" in corr:
            tag = f"{corr.get('name')}_before"
            flags[tag] = flags.get(tag, False) or corr["kinematic"]
    return flags

_warned = set()

def split_variations(variations, weight_columns):
    """(variations that need the full selection, weight-only variations).
    A variation reuses the nominal selection only if it replaces nothing but
    weight columns AND the config does not declare it kinematic; any
    disagreement between the two is reported once and resolved by running
    the full selection."""
    full, reweight = [], []
    for variation in sorted(variations.get_variations_all()):
        replaced = [sw[0] for sw in variations.get_variation_subs(variation) if len(sw) == 2]
        detected = any(column not in weight_columns for column in replaced)
        declared = declared_flags.get(variation)
        if declared is not None and declared != detected and variation not in _warned:
            _warned.add(variation)
            if declared:
                print(f"WARNING: {variation} is declared kinematic but only replaces weight "
                      f"columns {replaced}: running the full selection anyway")
            else:
                print(f"WARNING: {variation} is declared not kinematic but replaces "
                      f"{[c for c in replaced if c not in weight_columns]}: running the full selection")
        if variation == "nom" or detected or declared:
            full.append(variation)
        else:
            reweight.append(variation)
    return full, reweight

def nested_fields(events):
    """field name -> collections holding it (e.g. Jet and the LowPtJet/BJet
    built from it during the selection)"""
    nested = {}
    for collection in ak.fields(events):
        for field in ak.fields(events[collection]):
            nested.setdefault(field, []).append(collection)
    return nested

def replace_weight_columns(events, switches, nested):
    """Apply the substitutions of a weight-only variation to already selected
    events: a collection field is replaced in every collection holding both
    the nominal and the varied field"""
    for switch in switches:
        if len(switch) != 2:
            continue
        dest, source = switch
        if isinstance(dest, tuple):
            for collection in nested.get(dest[-1], []):
                if collection in nested.get(source[-1], []):
                    events[(collection, dest[-1])] = events[(collection, source[-1])]
        else:
            events[dest] = events[source]
    return events

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
keep_variations = special_analysis_cfg.get("variations", None) # fnmatch patterns, None = all
do_rochester_stat_variations = special_analysis_cfg.get("do_rochester_stat_variations", False)
do_jet_variations = special_analysis_cfg.get("do_jet_variations", False)
do_theory_variations = special_analysis_cfg.get("do_theory_variations", False)
invert_one_isolation = special_analysis_cfg.get("invert_one_isolation", False)
invert_one_isolation_loose = special_analysis_cfg.get("invert_one_isolation_loose", False)
invert_one_isolation_control = special_analysis_cfg.get("invert_one_isolation_control", False)
invert_both_isolation = special_analysis_cfg.get("invert_both_isolation", False)
declared_flags = declared_kinematic(analysis_cfg)

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

    # genWeight is float32 in NanoAOD: sum in float64, a float32 accumulator
    # loses ~1e-6 relative per 1e5 events, and more when summing all chunks
    sumw = ak.sum(ak.values_astype(events.weight, np.float64))
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

    # Compact the selected events: the selections above only index into the
    # full chunk, and everything below would keep working through those indices
    events = ak.to_packed(events)
    gc.collect()

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

        # ...and all of them switched off at once
        if ho_corrections:
            for ho_corr in ho_corrections:
                events[f"{ho_corr['name']}_ho_before"] = ak.ones_like(events.weight)
            variations.register_variation([ho_corr["name"] for ho_corr in ho_corrections], "ho_before")

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
    elif keep_variations is not None:
        variations.variations_dict = {
            k: v for k, v in variations.variations_dict.items()
            if k == "nom" or any(fnmatch(k, pattern) for pattern in keep_variations)
        }

    region_axis = hist.axis.StrCategory([region for region in regions], name="category")
    variation_axis_names = sorted(list(variations.get_variations_all()))
    variation_axis = hist.axis.StrCategory(variation_axis_names, name="syst")
    default_axis = [region_axis, variation_axis]

    results = {}
    if subsamples != {}:
        for subsample in subsamples:
            results[f"{dataset}_{subsample}"] = {
                "sumw": sumw, "nevents": nevents, "events": 0, "histos": 0,
            }
    elif eft_reweighting is None:
        results[dataset] = {
            "sumw": sumw, "nevents": nevents, "events": 0, "histos": 0
        }

    for dataset_name in results:
        _events, histos = {}, {}
        for variable in variables:
            _events[variable] = ak.Array([])

            if "axis" in variables[variable]:
                axis_def = variables[variable]["axis"]
                axes = axis_def if isinstance(axis_def, list) else [axis_def]
                histos[variable] = hist.Hist(*axes, *default_axis, hist.storage.Weight())

        results[dataset_name]["histos"] = histos
        results[dataset_name]["events"] = _events

    # eft_reweighting: one EFTHistogram per variable holding all points, and
    # one holding all covariance terms (nominal only, and by default only
    # for the cards regions/variables, the only ones spritz-cov-matrix reads)
    if eft_reweighting is not None:
        eft_points = list(eft_reweighting["points"])
        cov_pairs = eft_reweighting.get("covariance_pairs", [])
        hist_variables = [v for v in variables if "axis" in variables[v]]
        # config's covariance_regions/variables, default cards_regions/variables
        cov_regions = analysis_cfg.get("covariance_regions", analysis_cfg.get("cards_regions", list(regions)))
        cov_regions = [r for r in regions if r in cov_regions]
        cov_variables = analysis_cfg.get("covariance_variables", analysis_cfg.get("cards_variables", hist_variables))
        cov_variables = [v for v in hist_variables if v in cov_variables]

        def get_axes(variable):
            axis_def = variables[variable]["axis"]
            return axis_def if isinstance(axis_def, list) else [axis_def]

        results[eft_points_key(dataset)] = {
            "sumw": sumw, "nevents": nevents, "events": 0,
            "eft_names": eft_points,
            "histos": {
                v: EFTHistogram(get_axes(v), list(regions), variation_axis_names, len(eft_points), variances=True)
                for v in hist_variables
            },
        }
        if cov_pairs:
            point_pos = {p: i for i, p in enumerate(eft_points)}
            cov_idx = (
                np.array([point_pos[a] for a, b in cov_pairs], dtype=np.intp),
                np.array([point_pos[b] for a, b in cov_pairs], dtype=np.intp),
            )
            results[eft_cov_key(dataset)] = {
                "sumw": sumw, "nevents": nevents, "events": 0,
                "eft_names": [cov_name(a, b) for a, b in cov_pairs],
                "histos": {
                    v: EFTHistogram(get_axes(v), cov_regions, ["nom"], len(cov_pairs))
                    for v in cov_variables
                },
            }

    ##################################################

    # Loop over variations
    print("Doing variations")
    # No deep copies: awkward arrays are immutable and field assignments only
    # change the ak.Array wrapper they are made on, so each variation can start
    # from a new wrapper sharing originalEvents' buffers (ak.copy duplicated
    # every buffer, including those of the full chunk behind the selections)
    originalEvents = events

    def select(events):
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
            return None

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

        return events

    def add_weight(events, variation):
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

        return events

    def fill(events, variation):
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
            point_weights = eft_point_weights(events, eft_reweighting)
            fill_eft(
                results[eft_points_key(dataset)]["histos"], point_weights,
                events, regions, variables, variation,
            )

            # Covariance terms are only computed once (nominal)
            if variation == "nom" and eft_cov_key(dataset) in results:
                fill_eft(
                    results[eft_cov_key(dataset)]["histos"], point_weights,
                    events, regions, variables, variation, cov_pairs=cov_idx,
                )

    # Weight-only variations reuse the nominal selection (see split_variations)
    weight_columns = WEIGHT_COLUMNS | {ho_corr["name"] for ho_corr in ho_corrections}
    full_variations, reweight_variations = split_variations(variations, weight_columns)

    for variation in full_variations:
        print(variation)
        events = ak.Array(originalEvents)

        for switch in variations.get_variation_subs(variation):
            if len(switch) == 2:
                variation_dest, variation_source = switch
                events[variation_dest] = events[variation_source]

        events = select(events)
        if events is None:
            continue

        # Variable and Regions definitions
        for variable in variables:
            if "func" in variables[variable]:
                events[variable] = variables[variable]["func"](events)

        for region in regions:
            regions[region]["mask"] = regions[region]["func"](events)

        selected = events
        nested = nested_fields(selected)
        for weight_variation in [variation] + (reweight_variations if variation == "nom" else []):
            if weight_variation != variation:
                print(weight_variation)
                events = replace_weight_columns(
                    ak.Array(selected), variations.get_variation_subs(weight_variation), nested
                )
            else:
                events = ak.Array(selected)
            events = add_weight(events, weight_variation)
            fill(events, weight_variation)

    gc.collect()
    return results

##################################################

def prune_read_form(data):
    """read_form without the per-event weight vectors this chunk never uses:
    the theory weights (LHEPdfWeight alone is 103 floats per event) when
    theory variations are off, the EFT reweighting weights for non-EFT datasets"""
    unused = set()
    if not do_theory_variations:
        unused |= {"LHEPdfWeight", "LHEScaleWeight", "PSWeight"}
    eft_reweighting = data.get("eft_reweighting")
    if eft_reweighting is None or eft_reweighting["weight_branch"] != "LHEReweightingWeight":
        unused.add("LHEReweightingWeight")
    return {coll: form for coll, form in data["read_form"].items() if coll not in unused}

def chunk_str(chunk):
    drop_keys = ["read_form", "ho_corrections", "eft_reweighting"]
    return str({k:v for k,v in chunk["data"].items() if not k in drop_keys})

if __name__ == "__main__":
    new_chunks = read_chunks("chunks_job.pkl")
    print("N chunks to process", len(new_chunks))

    # Results of chunks of the same dataset are summed into the first chunk
    # of that dataset that succeeded, to keep the output (and memory) small.
    # Never into a failed one: check_errors only flags chunks whose result
    # is {}, so merging into a failed chunk would hide its failure.
    first_ok = {}

    for i in range(len(new_chunks)):
        new_chunk = new_chunks[i]
        dataset = new_chunk["data"]["dataset"]

        if new_chunk["result"] != {}:
            print(f"Skip chunk {chunk_str(new_chunk)}, was already processed")
            first_ok.setdefault(dataset, i)
            continue

        print(f"chunk {i+1}/{len(new_chunks)}: {dataset}")

        try:
            data = {**new_chunk["data"], "read_form": prune_read_form(new_chunk["data"])}
            result = big_process(process=process, **data)
            new_chunks[i]["error"] = ""
            j = first_ok.get(dataset)
            if j is not None and result != {}:
                new_chunks[j]["result"] = add_dict(new_chunks[j]["result"], result)
                new_chunks[i]["result"] = { "real_results": {}, 
                    "performance": {k:v for k,v in result["performance"].items()} }
                del result
            else:
                new_chunks[i]["result"] = result
                if result != {}:
                    first_ok[dataset] = i

        except Exception as e:
            nice_exception = "".join(tb.format_exception(None, e, e.__traceback__))
            print(f"\n\nError for chunk {chunk_str(new_chunk)}", file=sys.stderr)
            print(nice_exception, file=sys.stderr)
            new_chunks[i]["result"] = {}
            new_chunks[i]["error"] = nice_exception

        print(f"Done {i+1}/{len(new_chunks)}\n")

    write_chunks(new_chunks, "results.pkl")

