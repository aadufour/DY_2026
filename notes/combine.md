# Combine — DY SMEFT Analysis Notes

---

## Conceptual Overview

*(from Giacomo's lectures, March 2026)*

**Likelihood** = ∏(channels, phase spaces) × ∏(bins) { L(μ, θ | data) × prior(θ) }

Phase-space regions must be statistically independent blocks (e.g. mll < 50 vs mll > 50, years 2016/2017/2018, orthogonal jet multiplicities, …).

**Nuisance parameters θ:**
- *Theory*: QCD scale (renorm. + fact. scales), αs, PDFs, …
- *Experimental*: trigger, muon/electron reco, JES, JER, MET, pile-up reweighting, jet pile-up ID, …

**Rate parameters** (flat prior, no Gaussian constraint): used for data-driven background normalisation.

**autoMCStats**: MC statistical uncertainty (Barlow-Beeston lite). Added to the datacard as `<bin> autoMCStats 10 0 1`.

---

## EFT Specifics

The EFT cross section expands as:

```
σ(c) = σ_SM  +  c · σ_lin  +  c² · σ_quad
```

The linear term `σ_lin` can be **negative** (interference), which breaks combine (PDFs must be non-negative).

**Solution (morphing convention):** store raw templates at c = 0 (SM), c = +1, c = −1, and reconstruct lin/quad analytically. The `AnomalousCouplingMorphing` physics model handles this internally.

---

## Environment Setup (LLR T3)

```bash
# Analysis (spritz, plotting, build scripts)
dy_analysis        # activates analysis_venv inside apptainer

# Combine — old lin/quad model (kept for reference)
dy_combine

# Combine — morphing model  ← ACTIVE for v7
dy_combine_morphing
```

`dy_combine_morphing` sources `analysis/combine_tools/env_llr_morphing.sh`, which:
1. Runs `cmsenv` in `CMSSW_spritz/CMSSW_14_1_0_pre4` (has `AnalyticAnomalousCoupling` with `AnomalousCouplingMorphing`)
2. Prepends `tools/combine_helpers`, `tools/combination`, `tools/plotters` to PATH so the morphing-aware scripts are picked up

---

## Script Inventory (`analysis/combine_tools/`)

| Script | Purpose | Notes |
|--------|---------|-------|
| `createJson.py` | Interactive: create `metadata.json` with operator scan ranges | Copy from previous version; edit ranges manually |
| `createCombineJson.py` | Parse datacard → `jsonComb.json` (bin → operator list) | Updated with `--binname` option (default `w1_`) |
| `createWS.py` | Wrapper for `text2workspace.py` → one `model_{op}.root` per operator | Updated to use `AnomalousCouplingMorphing_comb` |
| `runScans.py` | Runs `combine -M MultiDimFit` (initial fit + grid scan) | Reads ranges from `metadata.json` |
| `runPlots.py` | Makes likelihood scan plots | From morphing tools |
| `env_llr.sh` | Old combine env setup | |
| `env_llr_morphing.sh` | Morphing combine env setup ← use this | |

---

## Full Morphing Workflow (v7 — active)

### Prerequisites
- `histos.root` produced by spritz v7 (see `notes/spritz.md`)
- `dy_combine_morphing` environment active

### Step 1 — Build shapes.root + datacard.txt

```bash
dy_analysis
cd /grid_mnt/data__data.polcms/cms/adufour/spritz/configs/dy_smeftsim_v7

python3 /grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/spritz/build_shapes_morphing.py \
    --input  histos.root \
    --outdir datacards_morphing \
    --region inc_mm --variable mll
```

Output:
- `datacards_morphing/inc_mm/mll/shapes.root` — `histo_sm`, `histo_w1_{op}`, `histo_wm1_{op}`, `histo_Data`, plus `histo_{proc}_QCDscaleUp/Down`, `histo_{proc}_PDFUp/Down`
- `datacards_morphing/inc_mm/mll/datacard.txt` — process indices: sm=1 (background/reference), w1_op1=0, wm1_op1=−1, …

### Step 2 — Prepare metadata.json

Copy from a previous version and adapt:

```bash
cp .../dy_smeftsim_v6/datacards_morphing/inc_mm/mll/metadata.json \
   datacards_morphing/inc_mm/mll/metadata.json

# Fix analysis tag and nuisances
python3 -c "
import json
path = 'datacards_morphing/inc_mm/mll/metadata.json'
with open(path) as f: m = json.load(f)
m['analysis'] = 'dy_smeft_lo'
m['nuisances'] = ['QCDscale', 'PDF']
with open(path, 'w') as f: json.dump(m, f, indent=4)
"
```

`metadata.json` structure:
```json
{
    "analysis": "dy_smeft_lo",
    "card": "datacard.txt",
    "operators": {
        "cHDD": [-0.03, 0.03],
        "cHWB": [-0.01, 0.01],
        ...
    },
    "nuisances": ["QCDscale", "PDF"]
}
```

The operator ranges define the x-axis of the likelihood scan.

### Step 3 — Switch to morphing combine env

```bash
dy_combine_morphing
cd /grid_mnt/data__data.polcms/cms/adufour/spritz/configs/dy_smeftsim_v7/datacards_morphing/inc_mm/mll
```

### Step 4 — Create jsonComb.json

Parses the datacard for processes starting with `w1_` and maps them to the bin:

```bash
createCombineJson.py --datacard datacard.txt --binname w1_ --output jsonComb.json
```

Output: `{"inc_mm_mll": ["cHDD", "cHWB", ..., "cbl"]}` (27 operators)

### Step 5 — Build per-operator workspaces

```bash
createWS.py 1
```

Runs `text2workspace.py` with `AnomalousCouplingMorphing_comb` for each of the 27 operators.
Produces `model_cHDD.root`, `model_cHWB.root`, … in the current directory.

To run for 2D operator combinations: `createWS.py 2`.

### Step 6 — Initial fit

```bash
runScans.py 1 initial
```

Finds the global minimum (best fit) for each operator workspace.

### Step 7 — Likelihood scan

```bash
runScans.py 1 scan
```

Runs `combineTool.py` grid scan (50 points per operator by default).

Stat-only version:
```bash
runScans.py 1 initial --stat
runScans.py 1 scan    --stat
```

### Step 8 — Plots

Full-syst as main curve, stat-only overlaid (dashed red):
```bash
runPlots_compare.py 1 --label "Stat + Syst" --compare-stat
```

Or the inverse:
```bash
runPlots_compare.py 1 --stat --label "Stat only" --compare-syst
```

Plain plots (no comparison):
```bash
runPlots.py 1
```

To redo a single operator after fixing its scan range in `metadata.json`:
```bash
runScans.py 1 initial --doOnly cll1
runScans.py 1 scan    --doOnly cll1
runScans.py 1 initial --doOnly cll1 --stat
runScans.py 1 scan    --doOnly cll1 --stat
runPlots_compare.py 1 --label "Stat + Syst" --compare-stat
```

---

## Process Index Convention (AnomalousCouplingMorphing)

| Process | Index | Role in combine |
|---------|-------|----------------|
| `sm` | 1 | Background (reference/SM template) |
| `w1_op1` | 0 | Signal (c=+1 template) |
| `wm1_op1` | −1 | Signal (c=−1 template) |
| `w1_op2` | −2 | Signal |
| `wm1_op2` | −3 | Signal |
| … | … | … |

Combine requires ≥ 1 positive process index (background). `sm=1` fills that role.

---

## Old LHE-based Workflow (kept for reference)

Used before spritz v7. Built datacards directly from LHE files via `build_cache.py` + `build_datacard.py`.

```bash
# Steps 1-2: build histograms + datacard from LHE cache
dy_analysis
cd /grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/combine
python3 build_cache.py
python3 build_datacard.py --op cHDD

# Steps 3-7: combine (old EFTNegative model)
dy_combine
createJson.py --datacard datacard.txt
createCombineJson.py --datacard datacard.txt   # old version: looks for quad_ processes
createWS.py 1                                  # old version: uses AnomalousCouplingEFTNegative_comb
runScans.py 1 initial
runScans.py 1 scan
runPlots.py 1
```

---

## Key Differences: Old vs Morphing Workflow

| | Old (lin/quad) | Morphing (v7) |
|-|---------------|--------------|
| Process names | `sm_lin_quad_cHDD`, `quad_cHDD`, `lin_cHDD` | `sm`, `w1_cHDD`, `wm1_cHDD` |
| Physics model | `AnomalousCouplingEFTNegative_comb` | `AnomalousCouplingMorphing_comb` |
| `createCombineJson.py` | `--binname quad_` (default) | `--binname w1_` |
| Workspace | One per operator, separate lin/quad templates | One per operator, morphing from 3 templates |
| Theory systs | Not included | QCDscale + PDF shape nuisances |
| Datacard source | `build_datacard.py` (LHE-based) | `build_shapes_morphing.py` (spritz histos.root) |
