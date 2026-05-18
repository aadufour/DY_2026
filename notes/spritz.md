# Spritz Analysis on LLR T3

## Overview

Spritz is a CMS NanoAOD analysis framework. We use it to process private DY SMEFTsim LO NanoAOD samples and produce EFT-weighted histograms for offline SMEFT analysis.

The full pipeline is:
```
spritz-fileset → spritz-chunks → spritz-batch-llr → condor_submit → spritz-merge → spritz-postproc → spritz-plot / plot_eft_shapes.py
```

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

## Config versions

| Version | Config | Runner | Subsamples | Notes |
|---------|--------|--------|------------|-------|
| v4 | `config_dy_smeft_eft.py` | `runner_dy_smeft.py` | SM + op01_lin + op01_quad (×27) | Lin/quad precomputed — normalisation bug |
| v5 | `config_dy_smeft_v5.py` | `runner_dy_smeft.py` | SM + cHDD + cHDD_m1 (×27) | Raw weights — correct normalisation ✓ |
| syst | `config_dy_smeft_syst.py` | `runner_dy_smeft_syst.py` | SM + op01_lin + op01_quad (×27) | Theory variations enabled (future) |

**Active config: v5** (`config_dy_smeft_v5.py`)

---

## EFT subsample approach (v5)

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

### 2. `runner_dy_smeft_syst.py`
Same as above but removes the `dataset == "Zjj"` gate on theory variations:
```python
doTheoryVariations = special_analysis_cfg.get("do_theory_variations", False)
```

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

## Full Pipeline

### 0. Copy config to gridmount
```bash
cp /grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/spritz/config_dy_smeft_v5.py \
   /grid_mnt/data__data.polcms/cms/adufour/spritz/configs/dy_smeftsim_v5/config.py
```

### 1. Register fileset (inside apptainer, from config dir)
```bash
cd /grid_mnt/data__data.polcms/cms/adufour/spritz/configs/dy_smeftsim_v5
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

### 7. Post-process (inside apptainer)
```bash
spritz-postproc
```
Output: `histos.root`

**Note**: before running postproc with subsamples, patch `data/fileset.json` (see samples.json section above).

### 8. Plot EFT shapes (analysis_venv, from config dir)
```bash
dy_analysis
python3 /grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/spritz/plot_eft_shapes.py \
    --root histos.root --region inc_mm --outdir check
```

Produces one PNG+PDF per operator in `check/`: SM (black), c=+1 (orange), c=−1 (blue), ratio panel.

---

## histos.root structure

```
region/variable/nominal/histo_SAMPLENAME
```

Example (v5):
```
inc_mm/mll/nominal/histo_DYSMEFTsim_SM
inc_mm/mll/nominal/histo_DYSMEFTsim_cHDD
inc_mm/mll/nominal/histo_DYSMEFTsim_cHDD_m1
...
```

Read with uproot (in analysis_venv, no apptainer needed):
```python
import uproot
f = uproot.open("histos.root")
h = f["inc_mm/mll/nominal/histo_DYSMEFTsim_SM"]
vals, edges = h.values(), h.axes[0].edges()
```

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
