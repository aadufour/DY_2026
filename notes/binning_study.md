# Same-binning study ("agentic analysis", Govoni)

*(started 2026-11-20)*

## Goal
How much of the EFT sensitivity comes from the binning? Re-run the whole framework with the **same number of bins N**
in `mll`, `costhetastar`, `rapll_abs` and `triple_diff`, then compare the limits with v3
(`propcorr_v3.md`: mll 14, costhetastar 5, rapll_abs 4, triple_diff 14x5x4 = 280 bins).
Earlier idea: deferred memory "binning-vs-limits study, two extreme binnings".

## Step 1: LHE-level check of the candidate binnings (2026-11-20)
Tool: `analysis/combine_tools/lhe_binning_check.py` (binnings defined in `BINNINGS` inside the script).
Input: propcorr LHE cache `~/work/LHE/propcorr/CACHE/lhe_cache_propcorr_parallel.pkl` (660k events, fiducial
|eta|<2.4, pT>25/10 GeV, 50<mll<3000), weights / N_gen (1e5) x 59.74 fb^-1. Normalisation checked continuous across
the 7 generator mll slices.
Output: LLR `~/work/binning_study/lhe_N40v5/` (`<binning>_<var>.png`, `sigma_table.md`, `binnings.json`).

Per binning and variable the plot shows: SM events per bin, LHE MC stat of the SM (sqrt(sum w^2)/sum w) and raw LHE
events, the linear EFT shape (lin/SM scaled to max 1) for cHj3, cHe, clj1, ced, clu.
Figure of merit: stat-only linear expected 1 sigma per operator, 1/sqrt(sum_i lin_i^2/SM_i), with lin_i^2 corrected
for the LHE MC variance of lin_i (without it fine bins look up to 2x better: cHd, cQe, cbl, cHDD).
Limits of the method: parton level (no reco efficiency or resolution, no backgrounds, no systematics). Only for
comparing binnings.

### Candidate N = 40
| variable | binning |
|---|---|
| mll | in the spirit of `eft_bkg_fullsyst_v5` (34 bins): 50-75 by 5, 75, 78, 80-106 by 2 (Z peak), 110, 115, 120, 130, 140, 150, 175, 200, 250, 300, 350, 400, 500, 600, 700, 800, 1000, 1200, 1500, 3000 |
| costhetastar | 40 uniform in [-1, 1] |
| rapll_abs | 40 uniform in [0, 2.4] |
| triple_diff A | mll [50,70,80,86,90,94,100,120,200,500,3000] x cos [-1,0,1] x \|y\| [0,0.9,2.4] (10x2x2) |
| triple_diff B | mll [50,80,100,200,500,3000] x cos [-1,-0.5,0,0.5,1] x \|y\| [0,0.9,2.4] (5x4x2) |
| triple_diff C | mll [50,76,86,96,106,150,250,400,600,1000,3000] x cos [-1,0,1] x \|y\| [0,0.9,2.4] (10x2x2) |

0.9 = SM-weighted median of |y_ll|. Unrolled triple_diff plots use the same order as `hist_unroll` in
`post_process_minnlo.py` (`view().T.flatten()`): mll fastest, then costhetastar, |y| outermost, i.e. one block per
(cos, |y|) cell with the full mll spectrum inside (fixed 2026-11-20; the order does not change the 1 sigma numbers). A first try with 40 log-spaced mll bins was bad (one bin over the whole Z peak,
last bins with ~3 SM events and 40-90% MC stat): replaced by the v5-style binning.

### Statistics
- mll: MC stat < 3% from 75 GeV to 1 TeV; 6-11% at 50-75 GeV (few LHE events below the peak); 9-17% above 1 TeV
  (17% in 1500-3000, same as the v3 last bin).
- costhetastar: ~3% except the 3 outer bins per side (|cos| > 0.85: 7-11%).
- rapll_abs: ~3% up to |y| = 2, last bin (2.34-2.4) 15%. The four-fermion lin terms are pure MC noise in |y| (lin/SM ~ 0.002).
- triple_diff: A max 10%, B 15%, C 18% (top mll slice). v3 triple_diff (280 bins): up to 63%, some bins empty.

### Expected 1 sigma, ratio N40 / v3 (< 1 = better than v3)
triple_diff A/B/C are all compared with the v3 triple_diff (280 bins).

| op | mll | costhetastar | rapll_abs | td A | td B | td C |
|---|---|---|---|---|---|---|
| cHDD | 0.94 | 0.90 | 0.93 | 0.93 | 1.01 | 0.96 |
| cHWB | 0.91 | 0.90 | 0.91 | 0.93 | 1.05 | 0.98 |
| cbWRe | 0.76 | 1.00 | 0.99 | 0.83 | 1.00 | 0.93 |
| cbBRe | 0.76 | 1.00 | 1.01 | 0.82 | 1.00 | 0.93 |
| cHj1 | 0.94 | 0.99 | 0.99 | 0.96 | 1.01 | 1.01 |
| cHQ1 (=cHQ3) | 0.78 | 1.00 | 1.00 | 0.83 | 1.00 | 0.93 |
| cHj3 | 0.78 | 1.00 | 1.00 | 0.83 | 1.00 | 0.94 |
| cHu | 0.73 | 1.00 | 1.00 | 0.79 | 1.00 | 0.90 |
| cHd | 0.73 | 1.11 | 0.93 | 1.11 | 1.37 | 1.35 |
| cHbq | 0.82 | 1.00 | 1.01 | 0.85 | 1.00 | 0.94 |
| cHl1 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| cHl3 | 0.71 | 0.99 | 1.00 | 0.75 | 0.99 | 0.88 |
| cHe | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| cll1 | 0.87 | 1.00 | 1.00 | 0.90 | 1.00 | 0.96 |
| clj1 | 0.95 | 0.94 | 0.98 | 1.46 | 1.49 | 1.03 |
| clj3 | 0.98 | 0.48 | 1.00 | 1.14 | 1.16 | 0.98 |
| cQl1 (=cQl3) | 0.96 | 0.48 | 0.97 | 1.01 | 1.08 | 1.01 |
| ceu | 0.99 | 0.65 | 1.01 | 1.31 | 1.34 | 1.10 |
| ced | 0.99 | 0.46 | 0.99 | 1.18 | 1.21 | 1.06 |
| cbe | 0.98 | 0.48 | 0.99 | 1.06 | 1.17 | 1.00 |
| cje | 0.99 | 0.45 | 1.02 | 1.27 | 1.29 | 1.11 |
| cQe | 0.90 | 0.30 | 0.43 | 0.82 | 0.95 | 0.87 |
| clu | 0.99 | 0.74 | 1.00 | 1.32 | 1.34 | 1.11 |
| cld | 0.98 | 0.80 | 1.01 | 1.18 | 1.21 | 1.05 |
| cbl | 1.01 | 0.90 | 1.05 | 1.07 | 0.92 | 1.04 |

### Findings
- mll: the fine Z-peak bins gain 20-30% for the Z-coupling operators (cHu, cHd, cHl3, cHj3, cHQ1, cbWRe). Partly a
  parton-level effect: at reco the peak is smeared by the muon resolution (~1-2 GeV), so expect less.
- costhetastar 40 bins: big gains (x0.3-0.5) for the four-fermion operators. Not angular information: with
  |eta| < 2.4, large |cos theta*| is only reachable at high mll, so the outer cos bins act as a hidden high-mll selection.
- triple_diff with 40 bins: the split decides who wins. Peak-heavy (A) helps Z couplings, loses 15-45% for the
  four-fermion operators (all of mll > 500 in one bin). C (fewer peak bins, more at high mll) is the balanced one:
  Z couplings -5 to -12%, four-fermion within +-11%; only cHd loses (x1.35), it needs the angular bins.
- cHl1, cHe: binning-independent (pure normalisation effects).

## Step 2: how many bins can a single variable use? (2026-11-20)
Question (user): why not 280 bins for each 1D variable too, as in the v3 triple_diff? Tool:
`analysis/combine_tools/lhe_binning_saturation.py` (LLR `~/work/binning_study/lhe_saturation/`): stat-only linear
1 sigma vs N (5 ... 280) per 1D variable; mll with equal LHE events per bin and a 1 GeV minimum width, cos and |y| uniform.

| var | N=5 | 10 | 20 | 40 | 80 | 140 | 280 |
|---|---|---|---|---|---|---|---|
| mll: median sigma(N)/sigma(280) | 1.25 | 1.13 | 1.08 | 1.04 | 1.01 | 1.01 | 1.00 |
| mll: worst MC stat | 5% | 12% | 27% | 51% | 75% | 86% | 93% |
| mll: min SM events/bin | 399 | 89 | 20 | 5 | 1.6 | 0.6 | 0.2 |
| costhetastar: median sigma(N)/sigma(280) | 1.32 | 1.31 | 1.13 | 1.09 | 1.04 | 1.02 | 1.00 |
| costhetastar: worst MC stat | 2% | 3% | 8% | 11% | 26% | 27% | 44% |
| costhetastar: min SM events/bin | 6e6 | 8e5 | 8e3 | 500 | 18 | 0 | 0 |
| rapll_abs: median sigma(N)/sigma(280) | 1.03 | 1.03 | 1.05 | 1.02 | 1.03 | 1.04 | 1.00 |
| rapll_abs: worst MC stat | 2% | 4% | 8% | 15% | 34% | 65% | 99% |

- The 1D information saturates: mll and costhetastar within ~4% of their 280-bin value at N = 40-80; |y| alone carries
  no shape information whatever the binning. Beyond that, more bins only add MC noise (EFT templates are `noStat`, so it
  goes untreated into the fit) and empty bins (mll: < 1 SM event per bin at N >= 140; cos: empty bins beyond the
  acceptance edge at |cos| ~ 0.98). At reco, mll bins below ~1-2 GeV are also below the muon resolution.
- 280 1D bins are therefore not useful; 40 is on the low side for triple_diff (four-fermion operators lose 10-45% vs v3
  at 40 bins, step 1). Natural compromise: N ~ 80 (1D saturated, triple_diff e.g. 10x4x2 or 20x2x2).
- Caveat: one operator at a time, linear, stat only. A finely binned 1D mll is already as good as the 280-bin triple_diff
  for most single operators (exceptions: cHd 1.30x); the strength of triple_diff is in multi-operator (profiled) fits,
  where the angular/rapidity bins break degeneracies. This Fisher does not see that.

## Step 3: candidate N = 80 (2026-11-20)
`lhe_binning_check.py --set N80`, output LLR `~/work/binning_study/lhe_N80/`.

| variable | binning |
|---|---|
| mll | 50-86 by 2, 86-96 by 1 (peak core), 96-106 by 2, 106, 110, 115, 120-180 by 5, 180-200, 200-300 by 10, 300-500 by 20, 500, 550, 600, 650, 700, 750, 800, 900, 1000, 1200, 1500, 3000 (final, after step 4) |
| costhetastar | 80 uniform in [-1, 1] |
| rapll_abs | 80 uniform in [0, 2.4] |
| triple_diff A | mll [50,76,86,96,106,150,250,400,600,1000,3000] x cos [-1,-0.5,0,0.5,1] x \|y\| [0,0.9,2.4] (10x4x2) |
| triple_diff B | mll [50,60,70,76,82,86,89,92,95,100,106,120,150,200,250,300,400,500,700,1000,3000] x cos [-1,0,1] x \|y\| [0,0.9,2.4] (20x2x2) |

mll tail: a first version with 2000-2500-3000 had 77% MC stat in the last bin (about 1 SM event), merging 2-3 TeV
still 39% (few large-weight events): last bin 1500-3000 as in v3. Worst MC stat now: mll 23% (50-52 GeV),
costhetastar 26% and rapll_abs 34% (outermost bins), triple_diff A 26%, B 18% (top mll bins).
1 GeV peak bins are at the muon resolution at reco: expect part of the mll gain to disappear there.

Expected 1 sigma ratio N80 / v3 (triple_diff A and B vs v3 triple_diff, 280 bins):

| op | mll | costhetastar | rapll_abs | td A | td B |
|---|---|---|---|---|---|
| cHDD | 0.94 | 0.83 | 0.77 | 0.96 | 0.93 |
| cHWB | 0.91 | 0.87 | 0.84 | 0.98 | 0.91 |
| cbWRe | 0.75 | 0.98 | 0.92 | 0.93 | 0.79 |
| cbBRe | 0.75 | 1.00 | 0.97 | 0.93 | 0.79 |
| cHj1 | 0.94 | 1.00 | 1.00 | 1.01 | 0.96 |
| cHQ1 (=cHQ3) | 0.78 | 1.00 | 1.00 | 0.93 | 0.82 |
| cHj3 | 0.76 | 1.00 | 1.00 | 0.94 | 0.81 |
| cHu | 0.72 | 1.00 | 1.00 | 0.90 | 0.75 |
| cHd | 0.73 | 1.18 | 0.81 | 1.20 | 1.15 |
| cHbq | 0.81 | 0.98 | 0.97 | 0.94 | 0.86 |
| cHl1 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| cHl3 | 0.70 | 1.00 | 1.01 | 0.88 | 0.74 |
| cHe | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| cll1 | 0.86 | 1.00 | 1.00 | 0.96 | 0.89 |
| clj1 | 0.95 | 0.87 | 0.95 | 1.14 | 1.03 |
| clj3 | 0.98 | 0.47 | 1.00 | 0.99 | 0.96 |
| cQl1 (=cQl3) | 0.96 | 0.47 | 0.96 | 1.01 | 0.98 |
| ceu | 0.99 | 0.65 | 1.02 | 1.10 | 1.08 |
| ced | 0.99 | 0.44 | 0.94 | 1.06 | 1.04 |
| cbe | 0.98 | 0.46 | 0.99 | 1.00 | 0.98 |
| cje | 0.98 | 0.42 | 0.94 | 1.11 | 1.08 |
| cQe | 0.84 | 0.26 | 0.33 | 0.91 | 0.81 |
| clu | 0.98 | 0.73 | 1.00 | 1.10 | 1.08 |
| cld | 0.97 | 0.79 | 1.01 | 1.04 | 1.03 |
| cbl | 1.06 | 1.09 | 1.09 | 1.00 | 1.13 |

- 1D at 80 vs 40: mll, cos, |y| essentially unchanged (saturated, step 2).
- triple_diff B (20x2x2) is within -25% / +15% of the 280-bin v3 for every operator: better for Z couplings, at most
  8% worse for four-fermion operators (except cHd 1.15, cbl 1.13). A (10x4x2) gains less on Z couplings, similar elsewhere.

## Step 4: generator mll slices, normalisation and MC stat at their boundaries (2026-11-20)
Concern (user): steps at the gridpack mll slice boundaries (120, 200 GeV) in the N80 plots, i.e. events not
normalised to the right cross section of their slice.
- Where the per-slice cross sections come from: `spritz_fabian/data/Full2018v9/samples/samples.json`, entries
  `DYSMEFTsim_LO_propcorr_mll_<slice>` (xsec 1088.54, 13.138, 4.2912, 2.0757, 1.5693, 1.2393, 4.2038 pb for 50_120 ...
  1000_3000), read from the `<init>` block (XSECUP) of the regenerated propcorr LHE files (`propagator_correction.md`
  section 17). This is the cross section at the generation (reweighting base) point; the sm/w1/wm1 templates are that
  times the LHEReweightingWeight ratio, as usual.
- LHE cache check: every slice has 10 files x 10000 events = 1e5, every event weight = XSECUP (event_norm average), so
  weight / 1e5 is the right per-slice normalisation, the same xsec as samples.json.
- Continuity of dsigma/dm (SM, LHE) across the slice boundaries, windows of +-4% of the boundary:

  | boundary [GeV] | ratio above/below | smooth power-law expectation |
  |---|---|---|
  | 120 | 0.78 +- 0.06 | 0.70 |
  | 200 | 0.86 +- 0.08 | 0.83 |
  | 400 | 0.85 +- 0.04 | 0.85 |
  | 600 | 0.89 +- 0.05 | 0.84 |
  | 800 | 0.80 +- 0.07 | 0.83 |
  | 1000 | 0.72 +- 0.13 | 0.82 |

  -> no normalisation step; lin/SM (cHe) also continuous. The steps in the plots were the bin-width changes (events per
  bin, 2 -> 5 GeV at 120, 5 -> 20 at 200): the 1D top panel now shows events / GeV.
- Real boundary effect: MC stat. Each slice has 1e5 events, but the spectrum falls inside the slice, so the top of each
  slice is thin (106-120: 5-13%, the Z peak takes the 50-120 events; 180-200: ~10%). N80 mll rebinned: wide bins at
  106-120 and 180-200, finer bins at 200-500 where the samples are dense. Now 106-120 6-8%, 180-200 5%.
- Found on the way: the SM template is reweighted from the generation point (all operators = 1), and the SM weights have
  outliers. 200-400 slice: one event (mll 230.5) with SM/generation weight ratio 782 (median 0.17) -> SM n_eff = 990 of
  95k events (1%), a 23% MC-stat spike in the 230-240 bin. SM n_eff per slice: 50-120 48%, 120-200 29%, 200-400 1%,
  400-600 11%, 600-800 3.4%, 800-1000 1.3%, 1000-3000 0.2%. Not a binning issue, but it affects the sm template
  (noStat) in every binning, v3 included: to check at reco level (open, `todo.md`).

## Config (2026-11-20)
N = 80, triple_diff B. `analysis/spritz/config_same_binning_N80.py` (= config_propcorr_v1.py + runner_pdfcut, new axes
only), deployed to LLR `spritz_fabian/configs/same_binning_N80/config.py`; loads in the spritz env (64 samples,
47 datasets, 14 nuisances; 80 bins in every variable).
Caveat for the postproc: the MiNNLO k = MiNNLO_V / LO_sm_V is per bin, so with 80 bins it also picks up the MC noise of
both samples (and the sm weight outliers above).

## Launch (2026-11-20 15:48)
- Condor cluster **823662**: 3000 jobs (same fileset/chunks as v3_pdfcut), long queue, request_memory 4500
  (v3: median 93 min, max 155 min, up to 4.2 GB). Proxy renewed for 7 days before submitting.
- Unattended chain in tmux `n80_chain` on llrcms: `analysis/spritz/run_chain_llr.sh same_binning_N80 823662
  propcorr_v3_pdfcut N80 v3`, log `same_binning_N80/chain.log`. Steps: wait, retry once bad jobs, merge, MiNNLO
  postproc, cards, 1D fits (100 points, stat+syst and stat-only) for mll / rapll_abs / costhetastar / triple_diff,
  summaries (`--logscale`) in `summary_plots/`, comparison with v3 in `compare_ref/`.
- Condor done 11-20 23:05 (all 3000 ok, no retry), merge 00:10, postproc + cards 00:16. The 1D fits then failed: no
  `metadata.json` (scan ranges) in the new cards. Fixed the chain to copy v3's per-variable `metadata.json` (same
  ranges -> like-for-like comparison) and to run the summaries in `dy_analysis`; relaunched 11-21 10:13 (first log
  kept as `chain_try1.log`).

## First results and a scan-range artefact (2026-11-21)
- Chain done 11-21 10:26 (1D fits 100 points, v3 per-variable ranges from `metadata.json`). Summaries had to be rerun in
  `dy_combine_morphing` (the `dy_analysis` venv has no ROOT); `run_chain_llr.sh` fixed.
- v3 1D binnings were coarse: mll 14 bins (whole Z peak in 80-100), cos 5, |y| 4. N80 vs v3 came out up to 50x tighter
  (stat+syst) and 4-10x (stat-only) even for |y|, where the LHE study shows no shape information beyond 4-5 bins.
- Cause: scan-grid resolution. With the fixed ranges the points next to the minimum already sit at 2dNLL >> 4
  (N80 |y| cHe stat-only: points at +-0.05 with 2dNLL ~800) and the 2dNLL = 4 crossing is a linear interpolation
  between 0 and 800 -> +-0.0002 instead of the parabolic +-0.0035 (= what mll, cos, triple_diff give). The artefact
  scales like sigma^2, so ratios get exaggerated. It also affects some v3 scans (even with 1000 stat-only points).
- Template check (`~/work/scratch_missing/noise2.py`): N80 |y| lin/sm is flat for cHe (-9%), cHl1 (+14%), cHj1 (+1.5%):
  same template information in 80 or 4 bins, so the |y| gains were not from shape or MC noise.
- Fix: `adapt_scan_ranges.py` (per-operator range = +-3.5 x the 95% half-width per side, parabolic estimate when the
  old scan does not resolve it; separate stat-only ranges) + `rescan_1d.sh` (two passes, both productions, in
  `<vardir>/rescan/`, original scans untouched; summaries, comparison and `compare_1d_widths.py` table in
  `same_binning_N80/compare_ref_rescan/`). Launched 11-21 12:22 in tmux `rescan1d`.

## Next
- When the chain is done: compare N80 vs v3 per variable (`compare_ref/`), push to EOS.
- Inspect the sm weight outliers at reco level (the 200-400 slice event with ratio 782).
