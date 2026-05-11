# Spritz Analysis on LLR T3

## Overview

Spritz is a CMS NanoAOD analysis framework. We use it to process private DY SMEFTsim LO NanoAOD samples and save per-event arrays (mll, EFT weights, PDF/scale weights) for offline SMEFT analysis.

The full pipeline is:
```
spritz-fileset → spritz-chunks → spritz-batch-llr → condor_submit → spritz-merge → spritz-postproc → spritz-plot
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
  -B /grid_mnt/data__data.polcms/cms/adufour/spritz/data/Full2018v9/samples/samples.json:/opt/spritz/data/Full2018v9/samples/samples.json \
  /grid_mnt/data__data.polcms/cms/adufour/spritz-env.sif bash --rcfile ~/.bashrc
```

Using `exec bash --rcfile` instead of `shell` forces `~/.bashrc` to be sourced, which puts `spritz-batch-llr` in PATH.

The third `-B` bind-mount is critical: it overlays our editable `samples.json` (with private sample paths) over the read-only frozen one inside the image.

### Exit apptainer
```bash
exit
```

**Important**: `condor_submit` is NOT available inside the apptainer. Always exit before submitting jobs.

---

## Changes vs Giacomo's original repo

### 1. `samples.json` (bind-mounted over `/opt/spritz/data/Full2018v9/samples/samples.json`)
Added 7 private DY SMEFTsim LO samples. Private samples use the `"path"` key (xrootd path to the directory) instead of a DAS query. Example entry:
```json
"DYSMEFTsim_LO_mll_50_120": {
    "path": "root://eos.grif.fr:1094//eos/grif/cms/llr/store/user/aldufour/3DY_SMEFTsim_LO/DYSMEFTMll-nanoaod18_SMEFTsim_mll_50_120/DYSMEFTMll-nanoaod18_SMEFTsim_mll_50_120/260504_081708",
    "xsec": "1.0",
    "kfact": "1.000",
    "ref": "A1"
}
```
Note: `xsec=1.0` is a placeholder. Use the LO cross section from the MadGraph banner file for the restriction card hypothesis (no reweighting).

### 2. Config (`analysis/spritz/config_dy_smeft_lo.py`)
Key differences from a standard Spritz config:
- `save_events=True` on all variables → saves per-event arrays instead of histograms
- Fine mll binning around Z peak (2 GeV steps from 76–105 GeV)
- Saves all 406 `LHEReweightingWeight` entries as `w_0`…`w_405`
  - Index 0 = SM, 1–27 = op_i(+1), 28–54 = op_i(−1), 55–405 = pairs
- Saves 8 QCD scale weights as `scale_w_0`…`scale_w_7` (from `LHEScaleWeight`)
- Saves 103 PDF weights as `pdf_w_0`…`pdf_w_102` (from `LHEPdfWeight`, NNPDF3.1)
  - Index 0 = central, 1–100 = replicas, 101–102 = alphaS up/down
- `do_theory_variations: False` (we handle theory weights ourselves)
- Runner: `runner_3DY_trees_singleTriggers.py`

### 3. `condor/run.sh`
The default spritz-batch generates a broken `run.sh` that doesn't call apptainer. Our working version:
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
Note: proxy must be on the gridmount (not `/tmp`) because the condor schedd runs on a different machine than the login node.

### 4. `condor/submit.jdl`
The default spritz-batch generates a jdl for CERN/Pisa clusters. For LLR T3:
- Remove `MY.SingularityImage` (not supported at LLR)
- Remove `Requirements = (machine == "pccms...")`
- Remove `+JobFlavour = "workday"`
- Remove `/tmp/x509up_u...` from `transfer_input_files`
- Fix `cfg.json` path: `/opt/spritz/...` → `/grid_mnt/data__data.polcms/cms/adufour/spritz/...`
- Add after the `log` line:
  ```
  T3Queue = short
  WNTag   = el9
  include : /opt/exp_soft/cms/t3/t3queue |
  ```

All of this is automated by `spritz-batch-llr` (see below).

---

## Full Pipeline

### 0. Copy config to gridmount
```bash
cp /grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/spritz/config_dy_smeft_lo.py \
   /grid_mnt/data__data.polcms/cms/adufour/spritz/configs/dy_smeftsim_vX/config.py
```

### 1. Register fileset (inside apptainer, from config dir)
```bash
cd /grid_mnt/data__data.polcms/cms/adufour/spritz/configs/dy_smeftsim_vX
spritz-fileset
```
Scans the xrootd paths in `samples.json` and builds the fileset.

### 2. Build chunks (inside apptainer, from config dir)
```bash
spritz-chunks
```
Splits the fileset into job chunks.

### 3. Build condor jobs — use our wrapper (inside apptainer, from config dir)
```bash
spritz-batch-llr
```
This runs `spritz-batch` and automatically fixes both `submit.jdl` and `run.sh` for LLR T3.

Note: `spritz-chunks` must be run before `spritz-batch-llr` (it needs `data/chunks.pkl`).

### 4. Submit to condor (OUTSIDE apptainer)
```bash
exit   # exit apptainer first!
cd /grid_mnt/data__data.polcms/cms/adufour/spritz/configs/dy_smeftsim_vX/condor
condor_submit submit.jdl
```

### 5. Monitor jobs
```bash
condor_q | grep adufour
```

### 6. Merge results (inside apptainer, from config dir)
```bash
cd /grid_mnt/data__data.polcms/cms/adufour/spritz/configs/dy_smeftsim_vX
spritz-merge
```
Output: `condor/results_merged_new.pkl`

### 7. Post-process (inside apptainer, from config dir)
```bash
spritz-postproc
```
Output: `histos.root`

### 8. Plot (inside apptainer, from config dir)
```bash
mkdir -p plots
spritz-plot
```

---

## Output data structure

The merged pkl uses spritz's compressed format. Read it with `read_chunks`, not plain `pickle.load`:

```python
from spritz.framework.framework import read_chunks
r = read_chunks('condor/results_merged_new.pkl')
# r is a list of job results
# r[i]['result']['real_results'][dataset][region][variable] → awkward array
```

Structure:
```
r[i]['result']['real_results']
  └── 'DYSMEFTsim_LO_mll_50_120'
        ├── 'sumw'      # float32, sum of gen weights
        ├── 'nevents'   # int
        ├── 'events'
        │     ├── 'inc_ee'   # dielectron region
        │     │     ├── 'mll'        # awkward array, per-event reco mll [GeV]
        │     │     ├── 'weight'     # per-event nominal weight
        │     │     ├── 'w_0'        # SM LHEReweightingWeight
        │     │     ├── 'w_1'…'w_405' # EFT weights
        │     │     ├── 'scale_w_0'…'scale_w_7'
        │     │     └── 'pdf_w_0'…'pdf_w_102'
        │     ├── 'inc_mm'   # dimuon region
        │     └── 'inc_em'   # mixed region
        └── 'histos'    # pre-binned histograms for mll, costhetastar_bins, yZ_bins
```

---

## EFT reweighting formula

For operator k at Wilson coefficient value c:
```python
w_kp = ev['w_{k}']        # weight with op k at c=+1  (index 1..27)
w_km = ev['w_{k+27}']     # weight with op k at c=-1  (index 28..54)
w_sm = ev['w_0']

w_lin  = 0.5 * (w_kp - w_km)           # linear term
w_quad = 0.5 * (w_kp + w_km) - w_sm    # quadratic term
w_eft  = w_sm + c * w_lin + c**2 * w_quad  # full weight at coupling c
```

---

## Useful one-liners

### Inspect a job pkl (outside apptainer, using apptainer exec)
```bash
apptainer exec -B /grid_mnt /grid_mnt/data__data.polcms/cms/adufour/spritz-env.sif python -c "
from spritz.framework.framework import read_chunks
import numpy as np
r = read_chunks('job_5/chunks_job.pkl')
ev = r[0]['result']['real_results']['DYSMEFTsim_LO_mll_50_120']['events']['inc_mm']
print('nevents:', len(ev['mll']))
print('mll[:5]:', np.array(ev['mll'][:5]))
print('w_0[:5]:', np.array(ev['w_0'][:5]))
"
```

### Check LHEPdfWeight and LHEScaleWeight counts in a NanoAOD file
```bash
apptainer exec -B /grid_mnt /grid_mnt/data__data.polcms/cms/adufour/spritz-env.sif python -c "
import uproot
f = uproot.open('root://eos.grif.fr:1094//eos/grif/cms/llr/store/user/aldufour/3DY_SMEFTsim_LO/DYSMEFTMll-nanoaod18_SMEFTsim_mll_50_120/DYSMEFTMll-nanoaod18_SMEFTsim_mll_50_120/260504_081708/0000/SMP-RunIISummer20UL18NanoAODv9-00051_1.root')
t = f['Events']
print('LHEPdfWeight length:', len(t['LHEPdfWeight'].array()[0]))
print('LHEScaleWeight length:', len(t['LHEScaleWeight'].array()[0]))
"
```

### Check how many jobs failed
```bash
apptainer exec -B /grid_mnt /grid_mnt/data__data.polcms/cms/adufour/spritz-env.sif python -c "
import pickle, glob
jobs = sorted(glob.glob('job_*/chunks_job.pkl'))
bad = []
for j in jobs:
    try:
        with open(j,'rb') as f: r = pickle.load(f)
        if r is None: bad.append(('none', j))
    except Exception as e: bad.append(('error', j, str(e)))
print(f'{len(bad)} failed / {len(jobs)} total')
for b in bad[:5]: print(b)
"
```
Note: job pkls use spritz compression, so plain `pickle.load` will fail even for good jobs. Use `read_chunks` instead to validate results.
