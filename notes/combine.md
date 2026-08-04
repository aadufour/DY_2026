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
| `readapt_double_boundaries.py` | both | 2D scan boundary sizing from likelihood grids — see [2D (Double) EFT Scans](#2d-double-eft-scans) |
| `check_condor_scan_failures.py` | both | Flags condor scan jobs that hit the EFT negative-yield failure — see below |

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

### Step 9 — Summary plot (all operators, two-panel)

Requires both full and stat-only scans to be complete (`higgsCombine.{op}.individual.*.root` AND `higgsCombine.{op}_stat.individual.*.root`).

Switch back to analysis env (NOT combine env — uses matplotlib):
```bash
dy_analysis
cd /path/to/datacards/inc_mm/mll
makeSummary.py --indir .
```

Output: `eft_summary_two_panel.pdf` and `eft_summary_two_panel.png`
- Left panel: 1σ/2σ Wilson coefficient intervals per operator (blue = full, red = stat-only)
- Right panel: Λ reach at 95% CL in TeV (stacked bars, log scale)

Copy to EOS web area:
```bash
# Create dir if it doesn't exist:
xrdfs root://eosuser.cern.ch mkdir -p /eos/user/a/aldufour/www/my_output_dir
xrdcp eft_summary_two_panel.p* root://eosuser.cern.ch//eos/user/a/aldufour/www/my_output_dir/
```

**Note:** `makeSummary.py` crashes with `Singular matrix` if `results` is empty — this means no scan files were found (wrong `--indir`) or stat-only files are missing (skip operator). Always run stat scans before calling it.

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

---

## 2D (Double) EFT Scans

`runScans.py`/`runPlots.py 2` build every operator **pair** (`C(N,2)` combinations, e.g. 351 for 27 operators) instead of one scan per operator. Same environment/pipeline rules as the 1D workflow above apply — this section only covers what's different for `mode=2`.

### CL thresholds are different from 1D

A 1D scan uses `2*deltaNLL = 1.0` (68%) / `3.84` (95%) — the chi2(1 dof) quantiles.
A 2D **joint** region (both operators simultaneously) uses chi2(2 dof) instead:

| CL | 1D threshold | 2D joint threshold |
|----|-------------|---------------------|
| 68% | 1.00 | **2.30** |
| 95% | 3.84 | **5.99** |

`readapt_double_boundaries.py` and the `scanEFT` class inside `mkEFTScan.py` (`HiggsAnalysis.AnalyticAnomalousCoupling.utils.scan`) both use 2.30/5.99 for the 2D contours.

### `readapt_double_boundaries.py` — sizing the scan boundaries

Reads every `higgsCombine.<op1>_<op2>.individual.MultiDimFit.mH125.root` grid, finds the 68%/95% contours via ROOT's `TH2::SetContour` + `CONT Z LIST` (same mechanism `mkEFTScan.py` uses), and proposes a per-operator boundary sized so the 95% CL comfortably fits inside the box rather than touching the edge or rattling around in a needlessly huge one.

Key options:
- `--target-fill` (default `0.5`) — when a boundary needs resizing, the new box is sized so the 95% CL reach fills this fraction of it.
- `--min-fill` / `--max-fill` (default `0.25` / `0.90`) — an operator is left untouched ("good as is") only if, for **every** pair it appears in, the 95% CL fill fraction is inside this band. Outside it, the box gets resized toward `--target-fill`. A single unsatisfied pair marks the whole operator for redo, since a boundary is shared across all 26 pairs an operator appears in.
- `--fallback-factor` (default `10.0`) — multiplier used only when a pair's likelihood is essentially flat within the current box (no usable curvature to extrapolate from at all).
- `--fit-zcap` (default `50.0`) — grid points above this `2*deltaNLL` are excluded from the local quadratic extrapolation fit used when a contour hasn't closed yet, so saturated/non-converged points don't skew it.
- `--pairs op1_op2,op3_op4,...` — restrict to specific pairs, for testing before a full run.
- `--metadata` / `--out-metadata` / `--scan-dir` / `--report` — as named.

```bash
# quick check on one pair
readapt_double_boundaries.py --metadata metadata.json --scan-dir . \
    --out-metadata metadata.json --report report.json --pairs=cbBRe_cbe

# full run, all pairs
readapt_double_boundaries.py --metadata metadata.json --scan-dir . \
    --out-metadata metadata_new.json --report report_full.json
```

Output includes a `GOOD AS IS` / `NEEDS REDO` list per operator, with which partner pair drove the lo/hi bound.

**Two things it corrects for that aren't obvious from the algorithm alone:**
- `--doSplitPoints` jobs each re-write their own copy of the best-fit reference point, so `hadd`-ing them back together stacks N duplicate copies of the origin point on top of the true minimum — this degrades both the `TGraph2D` Delaunay triangulation used for contour extraction and the local quadratic fit. `load_grid()` dedupes on rounded `(x,y)`, keeping the lowest `deltaNLL` per coordinate.
- If a pair's 95% contour hasn't closed within the current box, the reach estimate (extrapolated or fallback) is clamped to never propose a box *smaller* than what was already scanned without finding convergence — logically, if the box wasn't big enough, the true reach can't be smaller than it.

**Resolution matters more than you'd think.** A ~10×10-point grid (`--points=100` on a 2-POI scan) can fail to resolve a contour that's actually there, producing a spuriously narrow or unstable result — the fix in one case here was simply using far more points (thousands), not a smaller box.

### Running the scan at scale — HTCondor support

`runScans.py 2 scan --points=N --doSplitPoints=M` builds `pairs × M` individual `combineTool.py` commands. By default these run locally via Python `multiprocessing` on whatever machine invokes the script — fine for a handful of pairs, not for the full 351 (e.g. `351 × 50 = 17,550` jobs). Add `--condor` to generate (never auto-submit) an LLR T3 HTCondor submission instead, matching the `T3Queue`/`WNTag`/`include: /opt/exp_soft/cms/t3/t3queue` conventions used in `analysis/gridpack/*.sub`.

Options:
- `--condor` — write per-job scripts + a `.sub` file under `--condor-dir` (default `condor_jobs/`), print the `condor_submit` command, and exit without submitting or running anything.
- `--condor-dir` (default `condor_jobs`)
- `--cmssw-base` (default `/grid_mnt/data__data.polcms/cms/adufour/CMSSW_spritz/CMSSW_14_1_0_pre4`) — `cmsenv`'d into on the worker node; must be visible from LLR worker nodes.
- `--proxy` (default `/home/llr/cms/adufour/.t3/proxy.cert`) — exported as `X509_USER_PROXY` in each job.
- `--condor-queue` (default `long`) — LLR `T3Queue` value.
- `--request-memory` (default `2G`), `--request-cpus` (default `1`) — a combine grid-point fit is much lighter than gridpack generation, so these default well below the gridpack templates' `20G`/`8`.
- `--hadd-only` — skip building/running scan commands entirely and just hadd whatever `higgsCombine.<pair>.individual.POINTS.*.root` files are already on disk. Use this after the condor jobs finish.

```bash
runScans.py 2 initial --points=10000

runScans.py 2 scan --points=10000 --doSplitPoints=50 --condor
# review condor_jobs/scan.sub and a job_*.sh before submitting
condor_submit condor_jobs/scan.sub

# ... wait for all jobs to finish ...

python3 check_condor_scan_failures.py --condor-dir condor_jobs   # see below, check before hadd'ing

runScans.py 2 scan --points=10000 --doSplitPoints=50 --hadd-only
```

Each generated `job_*.sh` is self-contained: exports the proxy, sources `cvmfs cmsset_default.sh`, `cmsenv`'s into `--cmssw-base`, `cd`s into the working directory the scan was launched from, then runs its one `combineTool.py` invocation. Always sanity-check one script's content, and test on a small `--doOnly` slice, before submitting the full batch.

### The EFT negative-yield failure mode

A condor scan job can finish with `return value 0` (no crash, no eviction) but still come back with far fewer points than requested. Cause: at large enough Wilson coefficient values, the quadratic SMEFT parametrization (`σ_SM + c·σ_lin + c²·σ_quad`, see [EFT Specifics](#eft-specifics) above) can predict a **negative** event yield in some bin — a genuine EFT-validity boundary, not a bug. combine's grid loop appears to stop processing the rest of that job's point range once it hits this, rather than skipping the bad point and continuing, so a chunk silently comes back partial with a clean exit code. Signature in the job's `.out` log:

```
FASTEXIT from pdf_bininc_mm_mll
RooAbsMinimizerFcn: Minimized function has error status.
...
Number of events is negative or error @ params=(... k_cbBRe = -14.55 ...)
```

`check_condor_scan_failures.py --condor-dir condor_jobs` scans every job's `.out` log for this signature (cross-referenced against `condor_jobs/joblist.txt` to know which operator pair each job belongs to) and reports which pairs hit it — run this **before** `--hadd-only`, since the per-chunk logs are the only place this is visible; after merging, an affected pair just looks like it has somewhat fewer total points than requested, not an obvious failure.

Operators whose scan box already reaches a large absolute Wilson coefficient value (either because the operator itself has weak sensitivity, or because a boundary-sizing heuristic scaled up an already-wide 1D range) are the ones most likely to hit this — there's no way to predict it exactly without evaluating the per-bin quadratic yield formula directly, so the practical approach is: submit the full batch, run the checker afterward, and narrow just the flagged operators' ranges before re-running those specific pairs.

### `runPlots.py` additions

`runPlots.py` (all modes) now also supports, mirroring `runPlots_compare.py`:
- `--doOnly op1,op2,...` — comma-separated operator names (not pair names) to restrict which combinations get (re)plotted.
- `-j` / `--cores N` — parallel `mkEFTScan.py` jobs via `ThreadPoolExecutor`, capped at half the machine's CPU count by default.

```bash
runPlots.py 2 --doOnly=cbBRe,cbe -j 4
```
