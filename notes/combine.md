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

## Two Pipelines — Critical Distinction

| | **RECO pipeline** | **LHE pipeline** |
|-|-------------------|-----------------|
| Input | spritz `histos.root` | LHE cache (`lhe_cache_syst.pkl`) |
| Process names | `sm`, `w1_{op}`, `wm1_{op}` | `sm_lin_quad_{op}`, `quad_{op}` |
| Physics model | `AnomalousCouplingMorphing_comb` | `AnomalousCouplingEFTNegative_comb` |
| Combine env | `dy_combine_morphing` | `dy_combine` |
| `createWS` script | `createWS.py` | `createWS_lhe.py` |
| `createCombineJson` | `--binname w1_` | `--binname quad_` |
| Datacard builder | `build_shapes_morphing.py` | `build_datacard_reco_bins.py` |
| Theory systs | QCDscale + PDF shape nuisances | QCDscale + PDF shape nuisances |
| Normalization | spritz-postproc handles it correctly | Weights must be divided by N_gen (fix in `build_datacard_reco_bins.py`) |
| Binning | 34 bins, 50–3000 GeV (RECO binning) | 34 bins, 50–3000 GeV (matched to RECO for comparison) |

**Never mix these two pipelines.** `AnomalousCouplingMorphing_comb` cannot read `quad_` process names and vice versa.

---

## Environment Setup (LLR T3)

```bash
# Analysis (spritz, plotting, build scripts)
dy_analysis        # activates analysis_venv inside apptainer

# Combine — RECO morphing pipeline  ← use for spritz v7/v8 datacards
dy_combine_morphing

# Combine — LHE pipeline (EFTNegative model)
dy_combine
```

`dy_combine_morphing` sources `analysis/combine_tools/env_llr_morphing.sh`, which:
1. Runs `cmsenv` in `CMSSW_spritz/CMSSW_14_1_0_pre4` (has `AnalyticAnomalousCoupling` with `AnomalousCouplingMorphing`)
2. Prepends `tools/combine_helpers`, `tools/combination`, `tools/plotters` to PATH

`dy_combine` sources `analysis/combine_tools/env_llr.sh`:
1. Runs `cmsenv` in the old CMSSW (has `CombinedLimit` with `AnomalousCouplingEFTNegative_comb`)

---

## Script Inventory (`analysis/combine_tools/`)

| Script | Pipeline | Purpose |
|--------|----------|---------|
| `createJson.py` | both | Interactive: create `metadata.json` with operator scan ranges |
| `createCombineJson.py` | both | Parse datacard → `jsonComb.json`. Use `--binname w1_` (RECO) or `--binname quad_` (LHE) |
| `createWS.py` | **RECO only** | `text2workspace.py` with `AnomalousCouplingMorphing_comb` |
| `createWS_lhe.py` | **LHE only** | `text2workspace.py` with `AnomalousCouplingEFTNegative_comb` |
| `runScans.py` | both | Runs `combine -M MultiDimFit` (initial fit + grid scan) |
| `runPlots.py` | both | Makes likelihood scan plots |
| `build_shapes_morphing.py` | RECO | Reads spritz `histos.root` → `shapes.root` + `datacard.txt` |
| `build_datacard_reco_bins.py` | LHE | Reads LHE cache → `histograms.root` + `datacard.txt` (RECO binning, N_gen fix applied) |
| `build_datacard_syst.py` | LHE | Same as above but with coarse 7-bin LHE binning |
| `rank_operators.py` | both | Ranked sensitivity plot from scan ROOT files |
| `env_llr_morphing.sh` | RECO | Morphing combine env setup |
| `env_llr.sh` | LHE | Old combine env setup |

---

## RECO Morphing Workflow (active — spritz v7/v8)

### Prerequisites
- `histos.root` produced by spritz v7/v8 (see `notes/spritz.md`)
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

### Step 3 — Switch to morphing combine env

```bash
dy_combine_morphing
cd /grid_mnt/data__data.polcms/cms/adufour/spritz/configs/dy_smeftsim_v7/datacards_morphing/inc_mm/mll
```

### Step 4 — Create jsonComb.json

```bash
createCombineJson.py --datacard datacard.txt --binname w1_ --output jsonComb.json
```

Output: `{"inc_mm_mll": ["cHDD", "cHWB", ..., "cbl"]}` (27 operators)

### Step 5 — Build per-operator workspaces

```bash
createWS.py 1
```

Produces `model_cHDD.root`, `model_cHWB.root`, … using `AnomalousCouplingMorphing_comb`.

### Step 6 — Initial fit

```bash
runScans.py 1 initial
```

### Step 7 — Likelihood scan

```bash
runScans.py 1 scan
```

Stat-only version:
```bash
runScans.py 1 initial --stat
runScans.py 1 scan    --stat
```

### Step 8 — Plots

```bash
# Full-syst as main curve, stat-only overlaid
runPlots_compare.py 1 --label "Stat + Syst" --compare-stat

# Stat-only as main curve, full-syst overlaid
runPlots_compare.py 1 --stat --label "Stat only" --compare-syst

# Plain (no comparison)
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

## Process Index Convention (AnomalousCouplingMorphing — RECO only)

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

## LHE Workflow (parton-level validation)

Used to validate EFT templates at parton level and compare constraints with the RECO analysis.
Uses the old `AnomalousCouplingEFTNegative_comb` model with `quad_` / `sm_lin_quad_` process names.

### Normalization fix

LHE cache weights are **not divided by N_gen**, so raw `sum(w) * LUMI` gives ~10¹² events instead of ~6×10⁷. `build_datacard_reco_bins.py` applies the fix automatically:

```python
N_gen = len(w_SM)
w_SM = w_SM / N_gen   # and same for all other weight arrays
```

This makes `sum(w_SM) * LUMI ≈ σ_DY * LUMI ≈ 6×10⁷` — the correct expected yield.
Without this fix, MINUIT cannot converge with `--lumi 59740`.

### Step 1 — Build histograms + datacard

```bash
dy_analysis
cd /grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/combine/5_flav_bins

python3 /grid_mnt/.../DY_2026/analysis/combine_tools/build_datacard_reco_bins.py \
    --all_op --lumi 59740 --pdf-flavour 5 \
    --output histograms_reco_bins.root \
    --datacard datacard_reco_bins.txt
```

Rename to `datacard.txt` if `createWS_lhe.py` does not find it:
```bash
cp datacard_reco_bins.txt datacard.txt
```

### Step 2 — Prepare metadata.json

```bash
cp ../5_flav/metadata.json .
# Edit operator scan ranges manually if needed
```

### Step 3 — Switch to old combine env

```bash
dy_combine
cd /grid_mnt/.../combine/5_flav_bins
```

### Step 4 — Create jsonComb.json

```bash
createCombineJson.py --datacard datacard.txt --binname quad_ --output jsonComb.json
```

### Step 5 — Build per-operator workspaces

```bash
createWS_lhe.py 1
```

Uses `AnomalousCouplingEFTNegative_comb`. **Never use `createWS.py` here** — that script uses the morphing model and will fail on `quad_` process names.

### Steps 6–8 — Scans and plots

Same commands as RECO workflow:
```bash
runScans.py 1 initial
runScans.py 1 scan
runScans.py 1 initial --stat
runScans.py 1 scan    --stat
runPlots_compare.py 1 --label "Stat + Syst" --compare-stat
```

### Ranking plot (both pipelines)

```bash
# RECO (combine TTree output)
rank_operators.py --indir . --outdir ranking --stat --wide \
    --pattern-stat "higgsCombine.{op}_stat.individual.MultiDimFit.mH125.root"

# LHE (mkEFTScan TGraph output)
rank_operators.py --indir . --outdir ranking --stat --wide \
    --pattern "scan_{op}.root" \
    --pattern-stat "scan_{op}.root" \
    --tgraph-key-syst "Stat + Syst" \
    --tgraph-key-stat "Stat only"
```

---

## Key Differences: LHE vs RECO Workflow

| | **LHE** | **RECO (morphing)** |
|-|---------|---------------------|
| Process names | `sm_lin_quad_{op}`, `quad_{op}` | `sm`, `w1_{op}`, `wm1_{op}` |
| Physics model | `AnomalousCouplingEFTNegative_comb` | `AnomalousCouplingMorphing_comb` |
| Combine env | `dy_combine` | `dy_combine_morphing` |
| `createWS` script | `createWS_lhe.py` | `createWS.py` |
| `createCombineJson` flag | `--binname quad_` | `--binname w1_` |
| Datacard builder | `build_datacard_reco_bins.py` | `build_shapes_morphing.py` |
| Normalization | Manual N_gen fix required | Handled by spritz-postproc |
| Scan output format | TGraph in `scan_{op}.root` | TTree in `higgsCombine.{op}.*.root` |
| Theory systs | QCDscale + PDF | QCDscale + PDF |
| Binning | 34 bins, 50–3000 GeV (matched to RECO) | 34 bins, 50–3000 GeV |
