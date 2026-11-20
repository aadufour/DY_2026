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

## Next
- Choose N and the triple_diff split (current preference: N = 40, triple_diff C).
- New spritz config (copy of `propcorr_v3_pdfcut` with the new axes), then the full chain: jobs, merge, MiNNLO
  postproc, cards, 1D fits, summary plots; compare with v3.
