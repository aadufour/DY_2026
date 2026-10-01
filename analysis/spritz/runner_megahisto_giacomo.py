import gc
import json
import sys
import time
_t_module_start = time.time()
import traceback as tb
import awkward as ak
import correctionlib
import hist
import numpy as np
import onnxruntime as ort
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
from spritz.modules.fake_leptons import get_fake_weights, fakes_reweight
from spritz.modules.lepton_sel import create_lepton, lepton_sel
from spritz.modules.lepton_sf import lepton_sf
from spritz.modules.prompt_gen import prompt_gen_match_leptons
from spritz.modules.puweight import puweight_sf
from spritz.modules.rochester import (
    correct_rochester,
    get_rochester,
    vary_rochester,
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
print(f"  [timing] imports done +{time.time()-_t_module_start:.2f}s", flush=True)

# Matches configs/dy-eft-full-smeftsim-propcorr-2018/config.py's own
# covariance_name() exactly -- must, since post_process.py's `samples` dict
# references these same f"cov_{a}_{b}" keys to find the covariance terms.
def covariance_name(name_i, name_j):
    return f"cov_{name_i}_{name_j}"


def EFT_COMBINED_KEY(dataset):
    """Key under which the ONE combined (name-axis) eft_reweighting result
    lives in `results`, in place of 83,027 separate f"{dataset}_{name}"
    entries. Kept as a function (not a plain f-string inline) so every use
    site stays byte-for-byte consistent."""
    return f"{dataset}__eft_combined"


# `eft_reweighting` (an alternative to `subsamples` for EFT reweight points
# specifically -- both can be used on the same dataset) is a dict:
#   {
#     "weight_branch": "LHEReweightingWeight",
#     "points": {"sm": 0, "w1_cHDD": 1, "wm1_cHDD": 28, ...},  # name -> column index
#     "covariance_pairs": [("sm", "sm"), ("sm", "w1_cHDD"), ...],  # optional
#     "n_weights": 406,  # optional: exact width the weight branch must have
#   }
# `n_weights` guards against files whose reweight card differs from the one
# `points` was written for (e.g. an incomplete reweighting in one production
# job): the column indices above are only meaningful for that card layout, so
# a chunk whose branch has a different width fails with a chunk error instead
# of being silently mapped onto the wrong points.
# Unlike `subsamples` (a Python-eval'd (mask_expr, weight_expr) string pair
# per name -- the right tool for genuine per-name event *selections*, e.g.
# splitting a sample into Z->ee vs Z->mumu), every "point" here shares the
# exact same event selection (all events) and differs only in which column
# of the same 2D weight branch it multiplies in. That structure means the
# weight for every point -- and every requested covariance term -- can be
# computed as a single vectorized matrix operation instead of one eval()
# call per name, and histograms can be filled in batches instead of one
# fill() call per name. Both were measured to dominate the per-chunk cost
# once the number of names reaches the tens of thousands (406 real
# templates + one covariance term per unordered pair of them, for a full
# 27-operator EFT morphing fit, is 83,027 names).
#
# covariance_pairs is only meant to be passed for the "nom" variation --
# the covariance terms are MC-stat bookkeeping between templates built from
# the same nominal event weights, not a quantity with a meaningful
# systematic-varied version -- callers should omit it (or pass an empty
# list) for non-nominal variations so that cost (by far the largest part,
# e.g. 82,621 of the 83,027 names above) is only ever paid once per chunk.
FILL_BATCH_SIZE = 2000


def eft_reweighting_arrays(events, eft_reweighting):
    """Vectorized replacement for evaluating each point's/covariance term's
    weight_expr individually: pulls the weight branch out as one plain
    numpy 2D array, then computes every requested weight column as a single
    broadcasted operation. Returns (names, weights) where weights has shape
    (n_events, len(names))."""
    points = eft_reweighting["points"]
    covariance_pairs = eft_reweighting.get("covariance_pairs", [])

    W = ak.to_numpy(ak.to_regular(events[eft_reweighting["weight_branch"]]))
    n_weights = eft_reweighting.get("n_weights")
    if n_weights is not None and W.shape[1] != n_weights:
        raise ValueError(
            f"{eft_reweighting['weight_branch']} has {W.shape[1]} columns, expected "
            f"{n_weights}: this file's reweight card does not match the one the "
            "point -> column mapping was written for"
        )
    base_weight = ak.to_numpy(events.weight)

    point_names = list(points.keys())
    point_idx = np.array([points[n] for n in point_names], dtype=np.intp)
    point_weights = base_weight[:, None] * W[:, point_idx]

    if covariance_pairs:
        cov_names = [covariance_name(a, b) for a, b in covariance_pairs]
        idx_i = np.array([points[a] for a, b in covariance_pairs], dtype=np.intp)
        idx_j = np.array([points[b] for a, b in covariance_pairs], dtype=np.intp)
        cov_weights = (base_weight**2)[:, None] * W[:, idx_i] * W[:, idx_j]
        names = point_names + cov_names
        weights = np.concatenate([point_weights, cov_weights], axis=1)
    else:
        names = point_names
        weights = point_weights

    return names, weights

path_fw = get_fw_path()
with open("cfg.json") as file:
    txt = file.read()
    txt = txt.replace("RPLME_PATH_FW", path_fw)
    cfg = json.loads(txt)

ceval_puWeight = correctionlib.CorrectionSet.from_file(cfg["puWeights"])
ceval_lepton_sf = correctionlib.CorrectionSet.from_file(cfg["leptonSF"])
ceval_assign_run = correctionlib.CorrectionSet.from_file(cfg["run_to_era"])
print(f"  [timing] correctionlib CorrectionSets loaded +{time.time()-_t_module_start:.2f}s", flush=True)

rochester = get_rochester(cfg)
print(f"  [timing] rochester loaded +{time.time()-_t_module_start:.2f}s", flush=True)

analysis_path = sys.argv[1]
analysis_cfg = get_analysis_dict(analysis_path)
print(f"  [timing] config.py loaded (get_analysis_dict) +{time.time()-_t_module_start:.2f}s", flush=True)
special_analysis_cfg = analysis_cfg["special_analysis_cfg"]
sess_opt = ort.SessionOptions()
sess_opt.intra_op_num_threads = 1
sess_opt.inter_op_num_threads = 1


def process(events, **kwargs):

    dataset = kwargs["dataset"]
    trigger_sel = kwargs.get("trigger_sel", "")
    isData = kwargs.get("is_data", False)
    era = kwargs.get("era", None)
    subsamples = kwargs.get("subsamples", {})
    special_weight = eval(kwargs.get("weight", "1.0"))
    # Structured, vectorized alternative to `subsamples` for EFT reweighting
    # specifically -- see module docstring below for the shape and why this
    # exists alongside (not instead of) the generic subsamples mechanism.
    eft_reweighting = kwargs.get("eft_reweighting", None)

    # Diagnostic instrumentation: localizing a ~25x baseline-phase slowdown
    # observed on real cluster jobs vs local runs (isolated XRootD read
    # speed tests ruled out network contention as the cause -- 20 concurrent
    # reads of different files were just as fast as one alone). Temporary,
    # remove once the actual bottleneck step is identified.
    _t_checkpoint = [time.time()]

    def checkpoint(label):
        now = time.time()
        print(f"  [{dataset}] timing: {label} +{now - _t_checkpoint[0]:.2f}s", flush=True)
        _t_checkpoint[0] = now

    variations = variation_module.Variation()
    variations.register_variation([], "nom")

    if isData:
        events["weight"] = ak.ones_like(events.run)
    else:
        events["weight"] = events.genWeight

    if isData:
        lumimask = LumiMask(cfg["lumiMask"])
        events = lumi_mask(events, lumimask)
    else:
        events = pass_weightfilter(events, kwargs.get("max_weight", None))
        events = events[events.pass_weightfilter]

    sumw = ak.sum(events.weight)
    nevents = ak.num(events.weight, axis=0)

    # Add special weight for each dataset (not subsamples)
    if special_weight != 1.0:
        print(f"Using special weight for {dataset}: {special_weight}")

    events["weight"] = events.weight * special_weight

    # pass trigger and flags
    events = assign_run_period(events, isData, cfg, ceval_assign_run)
    events = pass_trigger(events, cfg["era"])
    events = pass_flags(events, cfg["flags"])

    events = events[events.pass_flags & events.pass_trigger]

    if isData: # each data DataSet has its own trigger_sel
        events = events[eval(trigger_sel)]

    checkpoint("read+trigger+flags")

    events = create_lepton(events)
    events = lepton_sel(events, cfg)
    events["Lepton"] = events.Lepton[events.Lepton.isLoose]

    # Apply a skim!
    events = events[ak.num(events.Lepton) >= 2]
    events = events[events.Lepton[:, 0].pt >= 24]
    events = events[events.Lepton[:, 1].pt >= 10]

    checkpoint("lepton_sel+skim")

    if dataset in ["DYmm_M-50to100", "DYmm_M-50"]: # for mll > 100 GeV we have separate DY samples
        outgoing_mask = (events.LHEPart.status == 1)
        lepton_mask = (abs(events.LHEPart.pdgId) == 13)
        lhe_leptons = events.LHEPart[outgoing_mask & lepton_mask]

        assert ak.all(ak.num(lhe_leptons) == 2)
        lhe_mll = (lhe_leptons[:, 0] + lhe_leptons[:, 1]).mass
        events = events[(50 < lhe_mll) & (lhe_mll < 100)]

    # The mll-binned EFT productions are generated in a *wider* truth-level
    # window than their name states (e.g. the "mll50_120" gridpack actually
    # contains events from ~30 to ~150 GeV, not a hard 50-120 cut) -- the
    # declared xsec in samples.json corresponds to the *stated* window only,
    # so without this cut each sample leaks extra, wrongly-normalized events
    # into its neighbors' territory, producing jumps/double-counting right at
    # the mll50_120/120_200/... boundaries. Same substring-matching pattern as
    # spritz's runner_3DY_trees_singleTriggers.py.
    if not isData:
        outgoing_mask = (events.LHEPart.status == 1)
        lepton_mask = (
            (abs(events.LHEPart.pdgId) == 11)
            | (abs(events.LHEPart.pdgId) == 13)
            | (abs(events.LHEPart.pdgId) == 15)
        )
        lhe_leptons = events.LHEPart[outgoing_mask & lepton_mask]
        if ak.all(ak.num(lhe_leptons) == 2):
            lhe_mll = (lhe_leptons[:, 0] + lhe_leptons[:, 1]).mass

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
            # SMEFTsim-only window (replaces SMEFTatNLO's 1000_1500 + 1500_inf
            # with one bin extending to 3000 GeV)
            if "1000_3000" in dataset:
                events = events[(lhe_mll > 1000) & (lhe_mll <= 3000)]

    checkpoint("LHE truth-window cut")

    if not isData:
        events = prompt_gen_match_leptons(events)

    checkpoint("prompt_gen_match")

    # Require at least one good PV
    events = events[events.PV.npvsGood > 0]

    # Top pT reweighting
    if kwargs.get("top_pt_rwgt", False):
        events, variations = tt_reweight(events, variations)
    else:
        events['topPtWeight'] = ak.ones_like(events.weight)

    # Correct Muons with rochester
    if special_analysis_cfg.get("do_rochester_variations", False):
        events, variations = vary_rochester(events, variations, isData, rochester)
    events, variations = correct_rochester(events, variations, isData, rochester)

    checkpoint("rochester")

    # Trigger matching
    events = match_trigger_object(events, cfg)

    checkpoint("trigger matching")

    # Fake lepton reweighting
    if special_analysis_cfg.get("reweight_fakes", False):
        variations, fakes_param = get_fake_weights(variations, cfg)

    if not isData:
        # puWeight
        events, variations = puweight_sf(events, variations, ceval_puWeight, cfg)

        checkpoint("puweight_sf")

        # add trigger SF
        events, variations = trigger_sf(events, variations, ceval_lepton_sf, cfg)

        checkpoint("trigger_sf")

        # add LeptonSF
        events, variations = lepton_sf(events, variations, ceval_lepton_sf, cfg)

        checkpoint("lepton_sf")

        # prefire
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

        # Theory unc.
        if special_analysis_cfg.get("do_theory_variations", True):
            events, variations = theory_unc(events, variations)

    checkpoint("prefire+theory_unc (end of baseline phase)")

    # Regions definitions
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

    # eft_reweighting names are declared once, up front, unconditionally --
    # this is the FULL (nominal-scope) name list, points + every requested
    # covariance term, regardless of which variation is currently being
    # processed. The combined histogram's "subsample" axis is built once
    # with this fixed category list; a non-nominal variation just never
    # fills the covariance slots (they stay at zero for that syst), same
    # semantics as before, just without needing a separately-sized object.
    eft_names_full = []
    eft_n_batches = 0
    if eft_reweighting is not None:
        eft_names_full = list(eft_reweighting["points"]) + [
            covariance_name(a, b) for a, b in eft_reweighting.get("covariance_pairs", [])
        ]
        eft_n_batches = -(-len(eft_names_full) // FILL_BATCH_SIZE)  # ceil div

    results = {dataset: {"sumw": sumw, "nevents": nevents, "events": 0, "histos": 0}}
    if subsamples != {} or eft_reweighting is not None:
        results = {}
        for subsample in subsamples:
            results[f"{dataset}_{subsample}"] = {
                "sumw": sumw,
                "nevents": nevents,
                "events": 0,
                "histos": 0,
            }
        if eft_reweighting is not None:
            # ONE combined result entry for the whole eft_reweighting name
            # list, instead of one per name -- see EFT_COMBINED_KEY's
            # histos below for the actual (batch-)histograms.
            results[EFT_COMBINED_KEY(dataset)] = {
                "sumw": sumw,
                "nevents": nevents,
                "events": 0,
                "histos": 0,
                "eft_names": eft_names_full,
                "eft_batch_size": FILL_BATCH_SIZE,
            }

    _t_histo_create = time.time()
    for dataset_name in results:
        _events = {}
        histos = {}
        is_eft_combined = eft_reweighting is not None and dataset_name == EFT_COMBINED_KEY(dataset)
        for variable in variables:
            _events[variable] = ak.Array([])

            if "axis" in variables[variable]:
                axis_def = variables[variable]["axis"]
                axes = list(axis_def) if isinstance(axis_def, list) else [axis_def]
                if is_eft_combined:
                    # A LIST of ~n_names/FILL_BATCH_SIZE (e.g. 42 for 83,027
                    # names) small-category-axis histograms, not one huge
                    # one and not 83,027 separate ones. A boost-histogram
                    # category axis's per-fill lookup cost was measured to
                    # scale with how many categories the axis was DECLARED
                    # with (not how many are actually touched per fill) --
                    # one histogram with an 83,027-category axis showed
                    # severe superlinear fill cost (0.4s for the first 2000
                    # names, but +26s just for names 10,000-20,000), while
                    # 42 independent small (<=2000-category) histograms
                    # filled the same 83,027 names in 15.8s flat, linear.
                    # This keeps that fill-time property while still
                    # cutting object count (hence creation+pickle cost) by
                    # ~2000x versus one histogram per name.
                    batch_histos = []
                    for b in range(eft_n_batches):
                        batch_len = min(FILL_BATCH_SIZE, len(eft_names_full) - b * FILL_BATCH_SIZE)
                        batch_axis = hist.axis.IntCategory(list(range(batch_len)), name="subsample")
                        batch_histos.append(hist.Hist(*axes, batch_axis, *default_axis, hist.storage.Weight()))
                    histos[variable] = batch_histos
                else:
                    histos[variable] = hist.Hist(*axes, *default_axis, hist.storage.Weight())

        results[dataset_name]["histos"] = histos
        results[dataset_name]["events"] = _events

    print(
        f"  [{dataset}] timing: created {len(results)} histogram-dict entries in "
        f"{time.time() - _t_histo_create:.2f}s",
        flush=True,
    )

    originalEvents = ak.copy(events)

    print("Doing variations", flush=True)
    for variation in sorted(variations.get_variations_all()):
        events = ak.copy(originalEvents)

        print(variation, flush=True)
        for switch in variations.get_variation_subs(variation):
            if len(switch) == 2:
                variation_dest, variation_source = switch
                events[variation_dest] = events[variation_source]

        # resort Leptons
        lepton_sort = ak.argsort(events[("Lepton", "pt")], ascending=False, axis=1)
        events["Lepton"] = events.Lepton[lepton_sort]

        # l2tight
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


        # Define categories
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

        # Analysis level cuts
        leptoncut = (events.mm | events.mm_ss)

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

        # Cut on pt of two leading leptons
        leptoncut = leptoncut & (events.Lepton[:, 0].pt > 29) & (events.Lepton[:, 1].pt > 15)
        events = events[leptoncut]

        ##################################################

        if len(events) == 0:
            continue

        if not isData:
            # Load all SFs
            events["RecoSF"] = events.Lepton[:, 0].RecoSF * events.Lepton[:, 1].RecoSF
            events["IdSF"] = events.Lepton[:, 0].IdSF * events.Lepton[:, 1].IdSF
            events["IsoSF"] = events.Lepton[:, 0].IsoSF * events.Lepton[:, 1].IsoSF

            events["weight"] = (
                events.weight
                * events.puWeight
                * events.topPtWeight
                * events.RecoSF
                * events.IdSF
                * events.IsoSF
                * events.prefireWeight
                * events.TriggerSF
            )

        # Fake lepton reweighting is only applied in the same-sign region
        if special_analysis_cfg.get("reweight_fakes", False):
            events["fakesRW"] = fakes_reweight(events, variation, fakes_param)
            events["fakesRW"] = ak.where(events.mm_ss, events.fakesRW, ak.ones_like(events.weight))

            events["weight"] = events.weight * events.fakesRW
        
        ##################################################
        # Variable definitions

        for variable in variables:
            if "func" in variables[variable]:
                events[variable] = variables[variable]["func"](events)
                
        ##################################################
        
        # Compute HO corrections if any 
        
        ho_corrections = kwargs.get("ho_corrections", False)
        if ho_corrections:
            print("Sono qui")
            import importlib
            for h__ in ho_corrections:
                ho_corr_module = importlib.import_module(f"spritz.modules.{h__['module']}")
                print(dir(ho_corr_module))
                events, variations = ho_corr_module.HO_reweight(events, variations, h__["weight"],  h__["weight_err"],  h__["edges"], name=h__["name"], observable=h__.get("observable"))
                
                print(f"Applying additional weight to {dataset}")
                print(f'{h__["name"]}: {events[h__["name"]]}')
                events["weight"] = events.weight * events[h__["name"]] # multiplicative correction

        for region in regions:
            regions[region]["mask"] = regions[region]["func"](events)

        def fill_histos(dataset_name, subsample_mask, subsample_weight):
            for region in regions:
                mask = regions[region]["mask"] & subsample_mask

                # events[mask] boolean-filters every field of the whole
                # events record just to check emptiness -- expensive, and
                # called once per subsample (83,027 times for a full EFT
                # morphing fit). ak.sum(mask) only touches the mask array
                # itself; same emptiness check, none of the full-record
                # filtering cost.
                if ak.sum(mask) == 0:
                    continue

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
                            weight=subsample_weight[mask],
                        )
                    else:
                        var_name = variables[variable]["axis"].name
                        results[dataset_name]["histos"][variable].fill(
                            events[var_name][mask],
                            category=region,
                            syst=variation,
                            weight=subsample_weight[mask],
                        )

        # Bare dataset (no subsample split): only present in `results` when
        # this dataset has no subsamples at all (see the results = {...}
        # setup above -- results only ever holds f"{dataset}_{subsample}"
        # keys once subsamples are used, never the bare dataset name too).
        if dataset in results:
            fill_histos(dataset, ak.ones_like(events.run) == 1.0, events.weight)

        # Subsample values can be either:
        #   - a str: a boolean mask expression (existing behavior, e.g. jet
        #     flavour splits) -- all subsamples share the same events.weight.
        #   - a (mask_expr, weight_expr) tuple/list: mask_expr selects events
        #     as before, and weight_expr replaces events.weight for that
        #     subsample only (e.g. an EFT LHEReweightingWeight point). This is
        #     how per-dataset EFT reweight predictions are built directly as
        #     histograms instead of via saved per-event trees.
        #
        # Each subsample's mask/weight is computed and immediately consumed
        # by fill_histos(), then goes out of scope -- deliberately a single
        # pass, not "compute all subsamples into a dict, then fill from the
        # dict afterward". A full EFT morphing fit can have 80k+ subsamples
        # (406 raw reweight points + one covariance term per unordered pair
        # of them); materializing all of their mask+weight arrays at once
        # (each the length of the whole chunk) pushed RSS over a 2GB condor
        # request_memory limit and put jobs on hold, even though the actual
        # per-subsample reweighting arithmetic is cheap. One pass keeps at
        # most one subsample's arrays alive at a time, independent of how
        # many subsamples there are.
        if subsamples != {}:
            n_subsamples = len(subsamples)
            print_every = max(1, n_subsamples // 20)  # ~20 progress lines regardless of scale
            t_sub_start = time.time()
            # Many subsamples can share the exact same mask expression (e.g.
            # every one of this config's 83,027 subsamples uses the same
            # "select everything" mask, only the weight differs) -- caching
            # by the expression string avoids re-eval'ing (and, for
            # ak.ones_like-style masks, re-allocating) an identical array
            # tens of thousands of times per chunk. Falls back to a plain
            # eval() for any expression not seen before, so this stays
            # correct and general for configs whose masks genuinely differ
            # per subsample too.
            mask_cache = {}
            for i, subsample in enumerate(subsamples):
                subsample_val = subsamples[subsample]
                if isinstance(subsample_val, str):
                    mask_expr, weight_expr = subsample_val, None
                elif isinstance(subsample_val, (tuple, list)) and len(subsample_val) == 2:
                    mask_expr, weight_expr = subsample_val
                else:
                    raise Exception(
                        "subsample value must be either a str (mask) or a "
                        "(mask, weight) tuple/list of length 2"
                    )

                if mask_expr not in mask_cache:
                    mask_cache[mask_expr] = eval(mask_expr)
                subsample_mask = mask_cache[mask_expr]

                if weight_expr is None:
                    subsample_weight = events.weight
                else:
                    subsample_weight = events.weight * eval(weight_expr)

                fill_histos(f"{dataset}_{subsample}", subsample_mask, subsample_weight)

                if (i + 1) % print_every == 0 or (i + 1) == n_subsamples:
                    elapsed = time.time() - t_sub_start
                    rate = (i + 1) / elapsed if elapsed > 0 else 0
                    print(
                        f"  [{dataset}] subsamples: {i + 1}/{n_subsamples} "
                        f"({100 * (i + 1) / n_subsamples:.1f}%) elapsed={elapsed:.1f}s "
                        f"rate={rate:.0f}/s",
                        flush=True,
                    )

        if eft_reweighting is not None:
            t_eft_start = time.time()
            # Covariance terms are MC-stat bookkeeping between templates
            # built from the same nominal event weights -- only computed
            # once, at nominal, regardless of what's in the config (see
            # module docstring above). Filling fewer names for a
            # non-nominal variation just leaves those "subsample" axis
            # categories at zero for this syst -- the combined histogram's
            # category list (eft_names_full) always covers the nominal
            # (full) scope, built once outside the variation loop.
            cov_pairs = eft_reweighting.get("covariance_pairs", []) if variation == "nom" else []
            names, weights = eft_reweighting_arrays(
                events, {**eft_reweighting, "covariance_pairs": cov_pairs}
            )
            eft_histos = results[EFT_COMBINED_KEY(dataset)]["histos"]

            for region in regions:
                mask_np = ak.to_numpy(regions[region]["mask"])
                n_sel = int(mask_np.sum())
                if n_sel == 0:
                    continue
                sel_weights_full = weights[mask_np]

                for variable, batch_histos in eft_histos.items():
                    axis_def = variables[variable]["axis"]
                    if isinstance(axis_def, list):
                        var_names = [k.name for k in axis_def]
                        sel_vals = {vn: ak.to_numpy(events[vn])[mask_np] for vn in var_names}
                    else:
                        var_name = axis_def.name
                        sel_vals = {var_name: ak.to_numpy(events[var_name])[mask_np]}

                    # `names` is a prefix of eft_names_full for a
                    # non-nominal variation (points only, no covariance),
                    # so `start // FILL_BATCH_SIZE` always lands on the
                    # correct persistent batch histogram -- batch
                    # boundaries are defined identically (same
                    # FILL_BATCH_SIZE) on both the histogram-creation side
                    # and here.
                    for start in range(0, len(names), FILL_BATCH_SIZE):
                        batch_idx = start // FILL_BATCH_SIZE
                        batch_names = names[start:start + FILL_BATCH_SIZE]
                        b = len(batch_names)
                        batch_weights = sel_weights_full[:, start:start + b]

                        tiled_vals = {vn: np.tile(v, b) for vn, v in sel_vals.items()}
                        local_idx_tiled = np.repeat(np.arange(b, dtype=np.intp), n_sel)
                        weight_tiled = batch_weights.T.reshape(-1)

                        batch_histos[batch_idx].fill(
                            **tiled_vals,
                            subsample=local_idx_tiled,
                            category=region,
                            syst=variation,
                            weight=weight_tiled,
                        )

            elapsed = time.time() - t_eft_start
            rate = len(names) / elapsed if elapsed > 0 else 0
            print(
                f"  [{dataset}] eft_reweighting ({variation}): {len(names)} names "
                f"elapsed={elapsed:.1f}s rate={rate:.0f}/s",
                flush=True,
            )

    gc.collect()
    return results


if __name__ == "__main__":
    chunks_readable = False
    new_chunks = read_chunks("chunks_job.pkl", readable=chunks_readable)
    print("N chunks to process", len(new_chunks), flush=True)

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

        print(f"chunk {i+1}/{len(new_chunks)}: {new_chunk['data']['dataset']}", flush=True)

        try:
            new_chunks[i]["result"] = big_process(process=process, **new_chunk["data"])
            new_chunks[i]["error"] = ""
        except Exception as e:
            print("\n\nError for chunk", new_chunk, file=sys.stderr)
            nice_exception = "".join(tb.format_exception(None, e, e.__traceback__))
            print(nice_exception, file=sys.stderr)
            new_chunks[i]["result"] = {}
            new_chunks[i]["error"] = nice_exception

        print(f"Done {i+1}/{len(new_chunks)}", flush=True)

    datasets = list(filter(lambda k: "root:/" not in k, results.keys()))

    _t_write = time.time()
    write_chunks(new_chunks, "results.pkl", readable=chunks_readable)
    print(f"  [timing] write_chunks(results.pkl) took {time.time()-_t_write:.2f}s", flush=True)
