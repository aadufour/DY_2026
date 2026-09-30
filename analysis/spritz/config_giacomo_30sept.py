# ruff: noqa: E501
#
# Full-scale local test bed for the "megahisto" batch-histogram
# eft_reweighting runner (runner_3DY_eft_full_morphing_megahisto.py).
# Identical scope, datasets, samples, regions, variables, njobs and
# data/chunks.pkl as configs/dy-eft-full-smeftsim-propcorr-2018-vectorized-demo
# (read that config's docstring, and dy-eft-full-smeftsim-propcorr-2018's
# before it) -- the only difference from -vectorized-demo is the runner:
#
#   - -vectorized-demo: one hist.Hist per eft_reweighting name (406
#     templates + 82,621 covariance terms = 83,027 objects), created/filled/
#     pickled individually every chunk. Correct, but ~62s of a ~93s
#     single-chunk runtime was pure object creation (40.5s) +
#     serialization (21.5s) overhead, not weight computation.
#   - here: eft_names is split into ~42 batches of <=2000 names, each
#     batch sharing ONE small hist.Hist with its own "subsample"
#     IntCategory axis, cutting object count ~2000x. Validated
#     bit-for-bit identical to the -vectorized-demo output (0/83,027
#     mismatches, values and variances) and ~3x faster end-to-end
#     (93.65s -> 30.47s per chunk; 14x faster than the original
#     eval()-based runner_3DY_eft_full_morphing.py's 428.81s).
#
# Downstream compatibility: post_process.py transparently unpacks the
# batch-histogram layout back into ordinary per-name histograms via
# spritz.framework.framework.expand_eft_combined() before anything else
# touches `results` -- build_covariance_matrix.py needs no changes, and
# older per-name results (this config's own -vectorized-demo, or the
# original subsamples-based dy-eft-full-smeftsim-propcorr-2018) pass
# through expand_eft_combined() untouched. merge.py's add_dict() was
# extended to sum the batch-histogram lists elementwise (not concatenate
# them) and to treat eft_names/eft_batch_size as invariant metadata rather
# than additive quantities.
#
# Memory: real per-job chunks are exactly 1000 events (matching njobs=6100
# for ~6059 real chunks); measured peak RSS is ~2.9GB/job regardless of
# which of the 7 mass-bin datasets, comfortably under a 4608MB condor
# request_memory and trivial for local runs on a machine with tens of GB
# free. Memory scales with events-per-chunk (the vectorized weight matrix
# is n_events x n_names), NOT with total dataset size or chunk count --
# don't lower njobs much below ~6100 or chunks get bigger and this budget
# no longer holds.

import itertools
import json

import hist
import numpy as np
from spritz.framework.framework import cmap_petroff, get_fw_path

fw_path = get_fw_path()
with open(f"{fw_path}/data/common/lumi.json") as file:
    lumis = json.load(file)

year = "Full2018v9"
lumi = lumis[year]["tot"] / 1000
lumi_unc = lumis[year]["rel_unc"]
plot_label = "DY EFT full fit -- SMEFTsim propcorr, 27 operators (megahisto runner)"
year_label = "2018"
njobs = 6100

runner = f"{fw_path}/src/spritz/runners/runner_3DY_eft_full_morphing_megahisto.py"

special_analysis_cfg = {
    "do_theory_variations": False,
}

VARIANT = "propcorr"

# -----------------------------
# Same 27 operators, same reweight-card index mapping, same naming
# convention as configs/dy-eft-full-smeftsim-propcorr-2018 -- see that
# config's docstring for how this was derived/verified.
# -----------------------------
OPERATORS = [
    "cHDD", "cHWB", "cbWRe", "cbBRe", "cHj1", "cHQ1", "cHj3", "cHQ3", "cHu",
    "cHd", "cHbq", "cHl1", "cHl3", "cHe", "cll1", "clj1", "clj3", "cQl1",
    "cQl3", "ceu", "ced", "cbe", "cje", "cQe", "clu", "cld", "cbl",
]
assert len(OPERATORS) == 27

_eft_points = (
    ["sm"]
    + [f"w1_{op}" for op in OPERATORS]
    + [f"wm1_{op}" for op in OPERATORS]
    + [f"w11_{a}_{b}" for a, b in itertools.combinations(OPERATORS, 2)]
)
assert len(_eft_points) == 406

_index_of = {"sm": 0}
for i, op in enumerate(OPERATORS):
    _index_of[f"w1_{op}"] = 1 + i
    _index_of[f"wm1_{op}"] = 28 + i
for k, (a, b) in enumerate(itertools.combinations(OPERATORS, 2)):
    _index_of[f"w11_{a}_{b}"] = 55 + k
assert set(_index_of) == set(_eft_points)


def covariance_name(name_i, name_j):
    return f"cov_{name_i}_{name_j}"


# -----------------------------
# The structured eft_reweighting spec, in place of the giant `subsamples`
# eval-string dict: "points" is directly {name: column index} (no need for
# the "events.LHEReweightingWeight[:, N]" string wrapper the eval-based
# mechanism needed -- the vectorized runner indexes the extracted weight
# matrix by column number directly), and "covariance_pairs" is the same
# set of unordered pairs as before, just as plain (name_i, name_j) tuples
# instead of pre-built weight-expression strings.
# -----------------------------
_covariance_pairs = list(itertools.combinations_with_replacement(_eft_points, 2))
assert len(_covariance_pairs) == 406 * 407 // 2  # 82,621

_eft_reweighting = {
    "weight_branch": "LHEReweightingWeight",
    "points": _index_of,
    "covariance_pairs": _covariance_pairs,
}

# -----------------------------
# Datasets -- all 7 SMEFTsim mll-binned 2018 datasets, propcorr variant.
# -----------------------------
_mll_bins = ["50_120", "120_200", "200_400", "400_600", "600_800", "800_1000", "1000_3000"]
datasets = {}
for b in _mll_bins:
    name = f"DYMuMu_LO_EFT_SMEFTsim_{VARIANT}_mll{b}_Photos_startingOne"
    datasets[name] = {
        "files": name,
        "task_weight": 8,
        "eft_reweighting": _eft_reweighting,
    }

for dataset in datasets:
    datasets[dataset]["read_form"] = "mc"

# -----------------------------
# Samples -- identical to dy-eft-full-smeftsim-propcorr-2018: same names,
# same covariance_of/is_variance/exclude_from_datacard flags. The runner
# that produced the underlying per-dataset histograms is invisible here.
# -----------------------------
samples = {
    point: {
        "samples": [f"{dataset}_{point}" for dataset in datasets],
        **({"is_signal": True} if point != "sm" else {}),
    }
    for point in _eft_points
}

samples.update({
    covariance_name(name_i, name_j): {
        "samples": [f"{dataset}_{covariance_name(name_i, name_j)}" for dataset in datasets],
        "is_variance": True,
        "exclude_from_datacard": True,
        "covariance_of": (name_i, name_j),
    }
    for name_i, name_j in _covariance_pairs
})
assert len(samples) == 406 + len(_covariance_pairs)

colors = {name: cmap_petroff[i % len(cmap_petroff)] for i, name in enumerate(_eft_points)}

# -----------------------------
# Regions -- identical to dy-eft-full-smeftsim-propcorr-2018.
# -----------------------------
preselections = lambda events: (events.mll > 50)  # noqa: E731

regions = {
    "inc_mm": {
        "func": lambda events: preselections(events) & events["mm"],
        "mask": 0,
    },
}

# -----------------------------
# Variables -- identical to dy-eft-full-smeftsim-propcorr-2018 (same
# triple_diff binning, same reasoning for why it's this coarse).
# -----------------------------
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
    },
    "costhetastar": {
        "func": lambda events: cos_theta_star(events.Lepton[:, 0], events.Lepton[:, 1]),
    },
    "rapll_abs": {
        "func": lambda events: abs((events.Lepton[:, 0] + events.Lepton[:, 1]).rapidity),
    },
    "triple_diff": {
        "axis": [
            hist.axis.Variable([50, 120, 200, 400, 600, 800, 1000, 3000], name="mll"),
            hist.axis.Variable([-1.0, 0.0, 1.0], name="costhetastar"),
            hist.axis.Variable([0.0, 1.2, 2.4], name="rapll_abs"),
        ],
        "label": ["$m_{\\ell\\ell}$", "$cos\\,\\theta^{\\ast}$", "$|y_{\\ell\\ell}|$"],
        "unit": ["GeV", "", ""],
        "xlog": True,
    },
}

cards_regions = ["inc_mm"]
cards_variables = ["triple_diff"]
covariance_file = "covariance.root"

nuisances = {}

nuisances["lumi"] = {
    "name": "lumi",
    "type": "lnN",
    "samples": dict((skey, "1.02") for skey in samples),
}

nuisances["stat"] = {
    "type": "auto",
    "maxPoiss": "10",
    "includeSignal": "0",
    "samples": {},
}

check_weights = {}
