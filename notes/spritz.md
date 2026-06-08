# Spritz Analysis on LLR T3

## Overview

Spritz is a CMS NanoAOD analysis framework. We use it to process NanoAOD samples and produce histograms for offline SMEFT analysis. There are two parallel tracks:

- **EFT signal track** (Giacomo's spritz, `analysis/spritz/`): processes private DY SMEFTsim LO NanoAOD, produces EFT-weighted histograms for morphing combine workflow.
- **Background track** (Fabian's spritz, `spritz_fabian/`): processes standard CMS MC backgrounds (DY MiNNLO, tt, WW, WZ, ZZ, single top, GG→ll) + data with fake lepton estimation. Target: reproduce Fabian's pag.5 plot (29 May) at LO, no systematics, lumi only.

The full pipeline (both tracks):
```
spritz-fileset → spritz-chunks → spritz-batch-llr → condor_submit → spritz-merge → spritz-postproc → spritz-plot
```

---

## Physics strategy (from Giacomo, June 2026)

### EFT + backgrounds

The goal is to build a complete EFT analysis on top of Fabian's background framework:
- If Wilson coefficients → 0, recover Fabian's SM plot
- EFT signal = our DY SMEFTsim LO samples, rescaled by a k-factor

### K-factor rescaling (bin-by-bin)

Our LO MadGraph SM prediction is less precise than MiNNLO (Fabian's DY). To match precision:

```
k = MiNNLO / SM_MG        (bin-by-bin ratio)
quad       → quad       × k
sm_lin_quad → sm_lin_quad × k
```

**Important**: cannot rescale SM alone inside the EFT parametrisation — the decomposition would no longer close. Must rescale all EFT components consistently.

### Fake leptons (reducible background)

- **Irreducible backgrounds** (same final state, e.g. ZZ→4l with 2 lost): cut away with kinematic selections.
- **Reducible backgrounds** (e.g. W+jet where a jet fakes a lepton, or Z→ll with one lost lepton faking a W): estimated **data-driven**, not from MC (MC QCD/PDF uncertainties are too large).
- Method: **same-sign / opposite-sign (SS/OS) ratio**. In the SS region, prompt dileptons are impossible (Z→l⁺l⁺ doesn't exist), so any SS events are fakes. After subtracting prompt MC from SS, the remainder gives the fake rate.
- Muon fake criterion: isolation — require no hadronic activity in ΔR < 0.4 cone. Non-isolated muons are likely non-prompt (e.g. from b-decays).
- Use loose isolation WP, measure SS/OS ratio in data after prompt MC subtraction.
- **Does not work for all topologies** — e.g. fails for ttbar (see Fabian pag.9, 29 May).

### Step 1: reproduce Fabian's background plot ✓ DONE

Run Fabian's config on standard CMS backgrounds (DY MiNNLO, TT, WW, WZ, ZZ, single top, GGToLL, data):
- No systematics (only lumi uncertainty for now)
- LO, no NLO corrections
- Config: `analysis/spritz/config_fabian_test.py` → deployed as `spritz_fabian/configs/test_v2/config.py`
- Reproduced Fabian's pag.5 (29 May) background stack plot ✓

### Step 2: add EFT signal on top ✓ DONE

Built EFT signal on top of Fabian's background framework using our DY SMEFTsim LO NanoAOD.

**Config:** `analysis/spritz/config.py` → deployed as `spritz_fabian/configs/test_v3/config.py`

**Runner:** `analysis/spritz/runner.py` — Fabian's `runner_3DY.py` with 4 EFT additions (marked `### EFT ###`):
1. Filter events to exactly 406 LHEReweightingWeights when `EFT` kwarg is set
2. LHE mll filter to avoid double-counting between the 7 mll-binned EFT samples
3. Subsample handling: supports `(mask_expr, weight_expr)` tuples for EFT weights
4. Per-subsample weight in histogram fill

`do_theory_variations=False` (lumi only). No systematics in this config — see Giacomo's v8 for theory systs.

**K-factor:** k(bin) = MiNNLO(bin) / SM_MG(bin), applied bin-by-bin at plotting stage to all EFT components (sm, w1_{op}, wm1_{op}) consistently. See notes on why all three must be rescaled.

**3000 condor jobs ran successfully** (batch 807506). Result in `spritz_fabian/configs/test_v3/`.

---

## Setup

### Apptainer image
The Spritz environment lives in a frozen Apptainer (Singularity) image:
```
/grid_mnt/data__data.polcms/cms/adufour/spritz-env.sif
```
All spritz commands (`spritz-fileset`, `spritz-chunks`, etc.) must be run **inside** the apptainer.

### Activate apptainer shell
Use the `spritz-shell` alias (defined in `~/.bashrc`):
```bash
spritz-shell
```

Which expands to:
```bash
apptainer exec \
  -B /etc/grid-security/certificates:/etc/grid-security/certificates \
  -B /cvmfs \
  -B /grid_mnt \
  /grid_mnt/data__data.polcms/cms/adufour/spritz-env.sif bash --rcfile ~/.bashrc
```

Using `exec bash --rcfile` instead of `shell` forces `~/.bashrc` to be sourced, which puts `spritz-batch-llr` in PATH.

### Stacking apptainer + analysis_venv
For plotting (needs mplhep), activate the analysis venv **inside** the apptainer:
```bash
spritz-shell
dy_analysis   # activates /grid_mnt/data__data.polcms/cms/adufour/analysis_venv
```
This gives access to both `spritz.framework` (from apptainer) and `mplhep` (from venv).

### Exit apptainer
```bash
exit
```

**Important**: `condor_submit` is NOT available inside the apptainer. Always exit before submitting jobs.

---

## Fabian's spritz setup (`spritz_fabian`)

Cloned from `https://github.com/fstaeg/spritz.git` into `/grid_mnt/data__data.polcms/cms/adufour/spritz_fabian`.

### Required fixes vs Fabian's upstream

| File | Fix |
|------|-----|
| `src/spritz/utils/rucio_utils.py` line 94 | `site.get("rse")` instead of `site["rse"]` (KeyError fix) |
| `src/spritz/scripts/batch.py` line ~160 | Add `year = an_dict["year"]` before command block (NameError fix) |
| `src/spritz/modules/lepton_sf.py` line 47 | `"isTightMuon_" + cfg["leptonsWP"]["muWP"]` instead of hardcoded `"isTightMuon_RelIso"` |
| `data/Full2018v9/samples/samples.json` | Has 6 colors in `cmap_petroff` — hardcode colors beyond index 5 in config |
| `batch_config.json` (root of repo) | Must point to LLR setup: condor, LLR proxy, LLR sif |
| `src/spritz/scripts/post_process.py` line ~331 | f-string fix: `read_chunks(get_batch_cfg()["BATCH_SYSTEM"] + "/results_merged_new.pkl")` |

### `~/.bashrc` additions for Fabian's spritz

```bash
export SPRITZ_PATH=/grid_mnt/data__data.polcms/cms/adufour/spritz_fabian
export PYTHONPATH=/grid_mnt/data__data.polcms/cms/adufour/spritz_fabian/src:/grid_mnt/data__data.polcms/cms/adufour/spritz_fabian:$PYTHONPATH
```

### `batch_config.json` (at `spritz_fabian/batch_config.json`)

```json
{
  "X509_USER_PROXY": "/grid_mnt/data__data.polcms/cms/adufour/proxy.pem",
  "SINGULARITY_IMAGE": "/grid_mnt/data__data.polcms/cms/adufour/spritz-env.sif",
  "BATCH_SYSTEM": "condor"
}
```

### Fresh start procedure (always do this)

When restarting a config from scratch, keep only:
```
config.py
```
Delete everything else (`condor/`, `slurm/`, `data/`, `__pycache__/`, `cfg.json`). Then:

```bash
# 0. Refresh proxy BEFORE submitting (critical for EOS/xrootd access)
cp /tmp/x509up_u1279 /grid_mnt/data__data.polcms/cms/adufour/proxy.pem

spritz-shell
cd /grid_mnt/.../configs/myconfig

# Option A: reuse existing fileset (saves ~10 min, skip if samples changed)
mkdir -p data
cp /grid_mnt/.../configs/test_v2/data/fileset.json data/

# Option B: regenerate fileset from scratch
spritz-fileset

spritz-chunks
spritz-batch-llr     # generates condor/ with correct run.sh and submit.jdl
exit
cd condor
condor_submit submit.jdl
```

**Key settings in config:**
- `njobs = 2000` — keeps each job short enough for the `short` queue (~1-2h each)
- Proxy must be fresh: `timeleft > 0` in `proxy.pem`

**After jobs finish:**
```bash
spritz-shell
cd /grid_mnt/.../configs/myconfig
spritz-merge       # run from config dir, NOT from condor/
spritz-postproc
spritz-plot        # or dy_analysis + custom plot script
```

---

## Full Pipeline (Fabian's spritz — combined EFT+bkg, test_v3)

### Active files
- Config: `analysis/spritz/config.py` → `spritz_fabian/configs/test_v3/config.py`
- Runner: `analysis/spritz/runner.py` → `spritz_fabian/configs/test_v3/condor/runner.py`
- 3000 jobs, lumi-only uncertainty, no theory systematics

### Required patches to Fabian's spritz (one-time, already applied)

| File | Patch |
|------|-------|
| `src/spritz/scripts/fileset.py` | Added `xrdfs ls -R` for xrootd paths + parallel uproot opens with tqdm + `-r` recycle flag |
| `data/common/forms.json` | Added `LHEReweightingWeight` to the mc form |
| `data/Full2018v9/samples/samples.json` | Added 7 DYSMEFTsim_LO_mll_* entries with xsec, kfact, ref, path |
| `src/spritz/scripts/make_cards.py` | Replaced hardcoded `good_regions`/`good_variables` with `list(analysis_dict["regions/variables"].keys())` |
| `src/spritz/scripts/make_cards.py` | Added `.replace(" ", "_")` to process name (fixes "Single Top" space → combine parse error) |

### Workflow

```bash
# Inside apptainer (spritz-shell)
cd /grid_mnt/data__data.polcms/cms/adufour/spritz_fabian/configs/test_v3

spritz-fileset     # discovers files via xrdfs ls -R for EFT samples
spritz-chunks
spritz-batch-llr
exit               # condor_submit must be outside apptainer
cd condor && condor_submit submit.jdl

# After jobs finish — back inside apptainer
spritz-merge
# Patch fileset.json for EFT subsamples (see snippet below)
spritz-postproc    # produces histos.root
spritz-cards       # produces datacards/inc_mm/mll/datacard.txt + shapes.root
```

### fileset.json patch (before spritz-postproc)

```python
import json
path = "data/fileset.json"
with open(path) as f: s = json.load(f)
MLL_BINS = ["50_120","120_200","200_400","400_600","600_800","800_1000","1000_3000"]
OPERATORS = ["cHDD","cHWB","cbWRe","cbBRe","cHj1","cHQ1","cHj3","cHQ3",
             "cHu","cHd","cHbq","cHl1","cHl3","cHe","cll1","clj1","clj3",
             "cQl1","cQl3","ceu","ced","cbe","cje","cQe","clu","cld","cbl"]
for b in MLL_BINS:
    base = f"DYSMEFTsim_LO_mll_{b}"
    parent = s[base].copy()
    s[f"{base}_sm"] = parent.copy()
    for op in OPERATORS:
        s[f"{base}_w1_{op}"]  = parent.copy()
        s[f"{base}_wm1_{op}"] = parent.copy()
with open(path, "w") as f: json.dump(s, f, indent=2)
```

### spritz-cards → combine

```bash
# Still inside apptainer, from test_v3/
dy_analysis
spritz-cards     # → datacards/inc_mm/mll/datacard.txt + shapes.root

dy_combine_morphing
cd datacards/inc_mm/mll
createCombineJson.py --datacard datacard.txt --binname w1_ --output jsonComb.json
# Create metadata.json (scan ranges per operator)
createWS.py 1
runScans.py 1 initial
runScans.py 1 scan
runScans.py 1 initial --stat
runScans.py 1 scan --stat
runPlots_compare.py 1

# Summary plot (analysis_venv)
dy_analysis
makeSummary.py --indir .
xrdcp eft_summary_two_panel.p* root://eosuser.cern.ch//eos/user/a/aldufour/www/
```

---

## Config versions

| Version | Config | Runner | Subsamples | Notes |
|---------|--------|--------|------------|-------|
| v4 | `old/config_dy_smeft_eft.py` | `runner_dy_smeft.py` | SM + op01_lin + op01_quad (×27) | Lin/quad precomputed — normalisation bug |
| v5 | `config_dy_smeft_v5.py` | `runner_dy_smeft.py` | SM + cHDD + cHDD_m1 (×27) | Raw weights — correct normalisation ✓ |
| syst | `config_dy_smeft_syst.py` | `runner_dy_smeft_syst.py` | SM + op01_lin + op01_quad (×27) | Prototype syst config — superseded by v7 |
| v7 | `config_dy_smeftsim_v7.py` | `runner_dy_smeft_v7.py` | sm + w1_{op} + wm1_{op} (×27) | Morphing names + PDF/QCD scale systematics. **Requires manual pkl 3-axis fix before postproc.** |
| **v8** | **`config_dy_smeftsim_v8.py`** | **`runner_dy_smeft_v8.py`** | sm + w1_{op} + wm1_{op} (×27) | **Same as v7 but histograms are 3-axis natively → no pkl preprocessing needed ← use for new runs** |

**Active runner: v8** (`runner_dy_smeft_v8.py`) — use for any new condor submission.
v7 results (already processed) remain valid.

---

## EFT subsample approach (v7 — active)

### Morphing-convention names

v7 renames subsamples and postproc samples to match `AnomalousCouplingMorphing.py` directly:

| v5 name | v7 name | Meaning |
|---------|---------|---------|
| `DYSMEFTsim_SM` | `sm` | SM template w(0) |
| `DYSMEFTsim_cHDD` | `w1_cHDD` | c=+1 template w(+1) |
| `DYSMEFTsim_cHDD_m1` | `wm1_cHDD` | c=-1 template w(-1) |

This means `histos.root` from `spritz-postproc` is already combine-ready. `build_shapes_morphing.py` only needs to add `histo_Data` (Asimov = SM) and write the datacard.

### Theory systematics (v7)

`do_theory_variations: True` in `special_analysis_cfg` enables:
- QCD scale variations: 6-point envelope via `LHEScaleWeight` (type `lheScaleWeight` in nuisances)
- PDF variations: replicas via `LHEPdfWeight` (type `lhePdfWeight` in nuisances)

`runner_dy_smeft_v7.py` changes vs v6:
1. `doTheoryVariations = special_analysis_cfg.get("do_theory_variations", False)` — no longer locked to `dataset == "Zjj"`
2. Nom-only filter is conditional: `if not doTheoryVariations:` — when theory vars are on, all variation slices are filled in the `syst` axis of each histogram

To inspect the exact variation names registered by `theory_unc`:
```bash
python3 -c "
from spritz.modules.theory_unc import theory_unc
import inspect
print(inspect.getsource(theory_unc))
"
```

### Subsample structure (v7)

55 subsamples per dataset:
```python
subsamples = {
    "sm":         (mask, "LHEReweightingWeight[:, 0]"),
    "w1_cHDD":    (mask, "LHEReweightingWeight[:, 1]"),   # c=+1
    "wm1_cHDD":   (mask, "LHEReweightingWeight[:, 28]"),  # c=-1
    "w1_cHWB":    (mask, "LHEReweightingWeight[:, 2]"),
    "wm1_cHWB":   (mask, "LHEReweightingWeight[:, 29]"),
    ...
}
```

### fileset.json patch for v7

Before running `spritz-postproc`, patch `data/fileset.json` to add all `dataset_subsample` keys:

```python
import json
path = "data/fileset.json"
with open(path) as f: s = json.load(f)
MLL_BINS = ["50_120","120_200","200_400","400_600","600_800","800_1000","1000_3000"]
OPERATORS = ["cHDD","cHWB","cbWRe","cbBRe","cHj1","cHQ1","cHj3","cHQ3",
             "cHu","cHd","cHbq","cHl1","cHl3","cHe","cll1","clj1","clj3",
             "cQl1","cQl3","ceu","ced","cbe","cje","cQe","clu","cld","cbl"]
for b in MLL_BINS:
    base = f"DYSMEFTsim_LO_mll_{b}"
    parent = s[base].copy()
    s[f"{base}_sm"] = parent.copy()
    for op in OPERATORS:
        s[f"{base}_w1_{op}"]  = parent.copy()
        s[f"{base}_wm1_{op}"] = parent.copy()
with open(path, "w") as f: json.dump(s, f, indent=2)
```

---

## EFT subsample approach (v5 — kept for reference)

### Why raw weights instead of lin/quad

Storing lin/quad precomputed in the runner causes a normalisation bug at postproc time:
- `sumw_lin = sum(genWeight × 0.5×(w_k+ − w_k−))` → a **difference**, can be tiny or near-zero
- postproc normalises each histogram by its own sumw → lin histogram gets blown up
- Result: visible bumps at mll bin boundaries in the EFT/SM ratio

**v5 fix (Giacomo's approach):** store raw weights at c=+1 and c=−1:
- `sumw_p1 = sum(genWeight × LHEReweightingWeight[:,k])` ≈ `sumw_SM` ≈ `sum(genWeight)`
- postproc normalisation is consistent for SM and all operator histograms
- Lin/quad decomposition is done **at plotting time**

### Operator names

27 operators extracted from LHE reweighting block (`<weight id='...'>` in header):

| idx | name | idx | name | idx | name |
|-----|------|-----|------|-----|------|
| 1 | cHDD | 10 | cHd | 19 | cQl3 |
| 2 | cHWB | 11 | cHbq | 20 | ceu |
| 3 | cbWRe | 12 | cHl1 | 21 | ced |
| 4 | cbBRe | 13 | cHl3 | 22 | cbe |
| 5 | cHj1 | 14 | cHe | 23 | cje |
| 6 | cHQ1 | 15 | cll1 | 24 | cQe |
| 7 | cHj3 | 16 | clj1 | 25 | clu |
| 8 | cHQ3 | 17 | clj3 | 26 | cld |
| 9 | cHu | 18 | cQl1 | 27 | cbl |

LHEReweightingWeight indexing (406 total):
- Index 0 = SM
- Indices 1..27 = op_k at c=+1
- Indices 28..54 = op_k at c=−1
- Indices 55..405 = pairs (op_i, op_j) — not used

### Subsample structure (v5)

55 subsamples per dataset:
```python
subsamples = {
    "SM":        (mask, "LHEReweightingWeight[:, 0]"),
    "cHDD":      (mask, "LHEReweightingWeight[:, 1]"),   # c=+1
    "cHDD_m1":   (mask, "LHEReweightingWeight[:, 28]"),  # c=-1
    "cHWB":      (mask, "LHEReweightingWeight[:, 2]"),
    "cHWB_m1":   (mask, "LHEReweightingWeight[:, 29]"),
    ...
}
```

### EFT reconstruction at plotting time

```python
vals_lin  = 0.5 * (vals_p1 - vals_m1)
vals_quad = 0.5 * (vals_p1 + vals_m1) - vals_sm
# EFT distribution at coupling c:
H_eft = H_sm + c * H_lin + c**2 * H_quad
```

---

## Cross sections (from LHE `<init>` block)

```
DYSMEFTsim_LO_mll_50_120    → 1188.67 pb
DYSMEFTsim_LO_mll_120_200   → 13.187  pb
DYSMEFTsim_LO_mll_200_400   → 4.287   pb
DYSMEFTsim_LO_mll_400_600   → 2.076   pb
DYSMEFTsim_LO_mll_600_800   → 1.565   pb
DYSMEFTsim_LO_mll_800_1000  → 1.240   pb
DYSMEFTsim_LO_mll_1000_3000 → 4.197   pb
```

Extract from LHE file:
```bash
for d in /grid_mnt/data__data.polcms/cms/adufour/LHE/SYST_slc7/DYSM*/; do
    echo -n "$(basename $d): "
    grep -A 4 "<init>" $d/unweighted_events.lhe | grep -v "[<>]" | sed -n '2p'
done
```

---

## samples.json

Location: `/grid_mnt/data__data.polcms/cms/adufour/spritz/data/Full2018v9/samples/samples.json`

Structure: `{"headers": {...}, "samples": {"dataset_name": {"path": ..., "xsec": ..., "kfact": ..., "ref": ...}}}`

For subsamples, postproc reads `xss` from `data/fileset.json` (local to each config dir), NOT from the global samples.json. When running postproc with subsamples, all `dataset_subsample` keys must be present in `data/fileset.json`. Add them with:

```python
import json
path = "data/fileset.json"
with open(path) as f: s = json.load(f)
MLL_BINS = ["50_120","120_200","200_400","400_600","600_800","800_1000","1000_3000"]
OPERATORS = ["cHDD","cHWB","cbWRe","cbBRe","cHj1","cHQ1","cHj3","cHQ3",
             "cHu","cHd","cHbq","cHl1","cHl3","cHe","cll1","clj1","clj3",
             "cQl1","cQl3","ceu","ced","cbe","cje","cQe","clu","cld","cbl"]
for b in MLL_BINS:
    base = f"DYSMEFTsim_LO_mll_{b}"
    parent = s[base].copy()
    s[f"{base}_SM"] = parent.copy()
    for op in OPERATORS:
        s[f"{base}_{op}"] = parent.copy()
        s[f"{base}_{op}_m1"] = parent.copy()
with open(path, "w") as f: json.dump(s, f, indent=2)
```

---

## Changes vs Giacomo's original repo

### 1. `runner_dy_smeft.py`
Based on Giacomo's `runner_3DY_trees_singleTriggers_EFT.py` with:
- Removed `onnxruntime`/DNN imports (not needed for DY EFT)
- Added try/except guard for `trigger_sf_latinos` import
- Added `1000_3000` mll LHE filter (was missing in original)
- **The 2-line fix**: histogram fill uses `events[f"weight_{dataset_name}"][mask]` instead of `events[cwgt][mask]` — makes the subsample EFT weight actually used in filling

### 2. `runner_dy_smeft_v7.py` ← active
Same as v6 with two additional changes:
- Removes the `dataset == "Zjj"` gate on theory variations:
  ```python
  doTheoryVariations = special_analysis_cfg.get("do_theory_variations", False)
  ```
- Makes the nom-only filter conditional so theory variation histograms are actually filled:
  ```python
  if not doTheoryVariations:
      variations.variations_dict = {k: v for k, v in variations.variations_dict.items() if k == "nom"}
  ```

### `runner_dy_smeft_syst.py` (superseded by v7)
Intermediate prototype — had `dataset == "Zjj"` removed but still stripped all non-nom variations unconditionally.

### 3. `condor/run.sh`
```bash
#!/bin/bash
export X509_USER_PROXY=/grid_mnt/data__data.polcms/cms/adufour/proxy.pem
time apptainer exec \
    -B /etc/grid-security/certificates:/etc/grid-security/certificates \
    -B /cvmfs \
    -B /grid_mnt \
    /grid_mnt/data__data.polcms/cms/adufour/spritz-env.sif \
    python runner.py .
```

### 4. `condor/submit.jdl`
Fixes for LLR T3:
- Remove `MY.SingularityImage`, `Requirements`, `+JobFlavour`, `/tmp/x509up` from inputs
- Fix `cfg.json` path from `/opt/spritz/` to gridmount path
- Add T3Queue/WNTag/include lines

All automated by `spritz-batch-llr`.

---

## Full Pipeline (v7)

### 0. Copy config to gridmount
```bash
# create config dir
mkdir -p /grid_mnt/data__data.polcms/cms/adufour/spritz/configs/dy_smeftsim_v7
cp /grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/spritz/config_dy_smeftsim_v7.py \
   /grid_mnt/data__data.polcms/cms/adufour/spritz/configs/dy_smeftsim_v7/config.py
```

### 1. Register fileset (inside apptainer, from config dir)
```bash
cd /grid_mnt/data__data.polcms/cms/adufour/spritz/configs/dy_smeftsim_v7
spritz-fileset
```

### 2. Build chunks
```bash
spritz-chunks
```

### 3. Build condor jobs
```bash
spritz-batch-llr
```

### 4. Submit (OUTSIDE apptainer)
```bash
exit
cd condor && condor_submit submit.jdl
```

### 5. Monitor
```bash
condor_q | grep adufour
find job_* -name chunks_job.pkl -size +14k | wc -l   # successful jobs
```

### 6. Merge (inside apptainer)
```bash
spritz-merge
```
Output: `condor/results_merged_new.pkl` (spritz compressed format)

### 7. Patch fileset.json for subsamples (before postproc)

See the v7 fileset.json patch snippet in the subsample section above (uses `sm`, `w1_{op}`, `wm1_{op}` keys).

### 8. Post-process (inside apptainer)
```bash
spritz-postproc
```
Output: `histos.root` with morphing-convention histogram names

### 9. Plot EFT shapes (analysis_venv, from config dir)
```bash
dy_analysis
python3 /grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/spritz/plot_eft_shapes.py \
    --root histos.root --region inc_mm --outdir check
```

Produces one PNG+PDF per operator in `check/`: SM (black), c=+1 (orange), c=−1 (blue), ratio panel.

### 10. Build combine shapes + datacard
```bash
dy_analysis
python3 /grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/spritz/build_shapes_morphing.py \
    --input histos.root \
    --outdir /grid_mnt/data__data.polcms/cms/adufour/spritz/configs/dy_smeftsim_v7/datacards_morphing \
    --region inc_mm --variable mll
```

Output: `datacards_morphing/inc_mm/mll/shapes.root` + `datacard.txt`

### 11. Run combine (morphing workflow)
```bash
dy_combine_morphing
cd /grid_mnt/data__data.polcms/cms/adufour/spritz/configs/dy_smeftsim_v7/datacards_morphing/inc_mm/mll
createJson.py --datacard datacard.txt --binname w1_
createCombineJson.py --datacard datacard.txt
createWS.py 1
runScans.py 1 initial
runScans.py 1 scan
runPlots.py 1
# or for stat vs syst comparison:
runScans.py 1 initial --stat && runScans.py 1 scan --stat
runPlots_compare.py 1 --label "Stat + Syst" --compare-stat
```

---

## histos.root structure

```
region/variable/nominal/histo_SAMPLENAME
```

Example (v7 — morphing convention):
```
inc_mm/mll/nominal/histo_sm
inc_mm/mll/nominal/histo_w1_cHDD
inc_mm/mll/nominal/histo_wm1_cHDD
...
```

Example (v5 — old naming, kept for reference):
```
inc_mm/mll/nominal/histo_DYSMEFTsim_SM
inc_mm/mll/nominal/histo_DYSMEFTsim_cHDD
inc_mm/mll/nominal/histo_DYSMEFTsim_cHDD_m1
```

Read with uproot (in analysis_venv, no apptainer needed):
```python
import uproot
f = uproot.open("histos.root")
h = f["inc_mm/mll/nominal/histo_sm"]
vals, edges = h.values(), h.axes[0].edges()
```

When theory variations are enabled (v7), histograms have a `syst` axis with entries for each
variation (`nom`, `QCDscale_*`, `PDF_*`, etc.). `spritz-postproc` collapses this to the
`nominal` directory for the nominal slice and separate directories for up/down variations.

---

## Useful one-liners

### Check job success rate
```bash
find job_* -name chunks_job.pkl -size +14k | wc -l   # good jobs (>14k)
find job_* -name chunks_job.pkl -size -14k | wc -l   # failed jobs
```

### Validate samples.json
```bash
python3 -c "import json; json.load(open('samples.json')); print('JSON OK')"
```

### Extract xsec from LHE init block
```bash
grep -A 4 "<init>" unweighted_events.lhe | grep -v "[<>]" | sed -n '2p'
```

### Get operator names from LHE reweighting block
```bash
grep "<weight " unweighted_events.lhe | head -60
```
