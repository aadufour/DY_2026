# Profiled fits on the triple differential (started 2026-11-17)

Goal (meetings 13 Oct, Govoni 23 Sep): 1D scans of each operator with the other operators floating
("profiled"), on the largest set of operators that can be constrained simultaneously. Uses the
propcorr v3 2D datacards (with all w11 cross terms), triple_diff. Asimov, all nuisances profiled.

## Which operators can be profiled together?

`analysis/combine_tools/fisher_profiling.py`: linear Fisher matrix from the templates,
F_ij = sum_b lin_i lin_j / (sm+bkg) (stat only, linear terms; optional `--free-norm`).
Per operator sigma_individual = 1/sqrt(F_ii), sigma_profiled = sqrt((F^-1)_ii); eigenvectors = weakest
combinations; `--greedy` drops the worst operator until all profiled/individual ratios < `--max-ratio`.

Result, triple_diff, 25 operators (cHQ3, cQl3 dropped as duplicates):
- 3 EXACTLY flat linear directions (eigenvalues ~1e-11, i.e. zero up to floating-point precision), only among
  cHl1, cll1, cHe, cHj3, cHl3, cHDD, cHWB (+ a little cHQ1): at linear order Drell-Yan sees fewer combinations
  of these than operators. All 25 cannot be profiled together.
- Greedy, max ratio 3 (stat only): 12 operators
  `cHDD,cHbq,cHd,cHj1,cHu,cQe,cQl1,cbe,cbl,ceu,cje,clj1` (worst ratios cHu 2.8, cHj1 2.5).
- With a free normalisation (worst case for lumi/xsec nuisances): 10 operators (cHj1, cHu also dropped).
- mll alone: condition number ~1e301, hopeless (as expected from the 2D scans).
Caveats: linear and stat-only, a guide for choosing the subset, not a result. Quadratic terms can close
linear flat directions far out; the real covariance comes from the combine fit.

### With systematics (2026-11-21, `--systs`, triple_diff_w11export)
The Fisher input (sm, bkg, w1/wm1) does not involve w11, so the subsets above were not affected by the missing
w11 step. New option `--systs`: the 12 shape nuisances + lumi lnN are parameters with unit Gaussian constraint and
profiled (Schur complement), MC stat added to the per-bin variance. Outputs: LLR
`propcorr_v3_2d/datacards/inc_mm/fisher_triple_diff/{stat,systs,stat_test2,systs_test2}.txt`.

Validation against combine (test2 subset, individual 68% from the 1D scans, profiled/individual from test2):

| op | Fisher+systs sigma_ind | combine 68% | Fisher+systs ratio | combine ratio | stat-only ratio |
|---|---|---|---|---|---|
| cHDD | 0.0087 | 0.0085 | 1.1 | 1.13 | 1.1 |
| cHj1 | 0.052 | 0.044 | 1.0 | 1.14 | 2.2 |
| cHu | 0.067 | 0.065 | 1.1 | 1.12 | 2.6 |
| cHd | 0.048 | 0.046 | 1.0 | 1.02 | 1.6 |
| clj1 | 0.0063 | ~0.006 | 1.0 | 1.04 | 1.0 |

So with systematics the linear Fisher reproduces combine; the stat-only one overestimates the profiling cost
(systematics dominate, and they absorb the operator correlations).

Result with systematics, 25 operators:
- The 3 exactly flat directions are still there (eigenvalues ~1e-10): cHDD, cHWB, cHl3, cHe, cll1, cHl1, cHj3.
  4th weakest: clj3 - ced (+clj1), x16.
- Greedy, max ratio 3: 16 operators
  `cHWB,cHbq,cHd,cHe,cHj1,cHj3,cHu,cQe,cQl1,cbBRe,cbWRe,cbe,cbl,ceu,cje,clj1` (worst ratios cbWRe, cbBRe 2.3).
  Dropped, in order: cHDD, cHl3, cll1 (the 3 flat directions), then clj3, clu, cHl1, ced, cHQ1, cld (ratio > 3).
  Which 3 of the flat-direction operators to drop is a choice: greedy drops the worst one each time, but any
  3 that break the 3 flat combinations work (e.g. keep cHDD and drop cHWB instead).

## How it is run (machinery already in runScans.py / createWS.py)

Folder: LLR `propcorr_v3_2d/datacards/inc_mm/triple_diff/profiled/` (symlinks to the 2D datacard, shapes,
jsonComb.json, metadata.json).

1. Workspace with only the chosen operators as POIs (all 27 kept in addToCompleteOperators):
   `python3 $CT/createWS.py 3 --ops cHDD,cHj1,cHu,cHd,clj1` -> `model_cHDD_cHj1_cHu_cHd_clj1.root`
   (15 s; contains the 10 mixed terms). `--ops` added to createWS.py 2026-11-17.
2. Scan ranges: `metadata_profiled.json` (the scanned operator uses metadata["operators"][op]; the profiled
   ones get +-300 in runScans.py, i.e. effectively free).
3. Per operator: `runScans.py 3 initial` and `runScans.py 3 scan` with `--doOnly <subset> --signalPOIs <op>
   --floatOtherPOIs=1` (added automatically for mode 3), `--savePOIs` (saves the profiled POIs),
   `--points N`. Wrapper: `analysis/combine_tools/run_profiled_scans.sh OPS [POINTS] [METADATA]`.
   Output: `higgsCombine.<OPS_>__<op>_profiled.individual.MultiDimFit.mH125.root`.
4. Comparison with the individual 1D scans (propcorr_v3_pdfcut, same datacard apart from w11):
   `analysis/combine_tools/compare_profiled_individual.sh OPS <v3_pdfcut>/datacards/inc_mm/triple_diff`
   -> `plots/scan_profiled_<op>.png` (mkEFTScan, individual dashed red) + interval table.

Timing: initial fit ~3.5 s; a 100-point profiled scan ~1-3 min per operator on one core (local, tmux).

Gotchas
- Use the repo runScans.py (`python3 $CT/runScans.py`, or `dy_combine_morphing` after the 2026-11-17 PATH fix).
  The stale tools copy has no --points/--condor and fails silently ("Number of processes must be at least 1").
- The option is `--savePOIs` (not --savepois).
- Do not take scan ranges from the 2D boxes: they are far wider than the 1D intervals.

## Runs

### test1 (2026-11-17, 5 operators: cHDD, cHj1, cHu, cHd, clj1) — NOT usable
Ranges 3x the 2D metadata boxes, 100 points: e.g. cHj1 +-11.5 (step 0.23) for an interval of +-0.1, so each
curve has 1-3 points near the minimum; cHd minimum apparently at +0.04, cHj1 curve meaningless.
Machinery validated (workspace, mixed terms, profiled fits, outputs). Moved to `profiled/test1_coarse/`.
Fisher prediction for this subset (profiled/individual): cHu 2.6, cHj1 2.2, cHd 1.6, cHDD 1.1, clj1 1.0.

### test2 (2026-11-17, same 5 operators) — SUPERSEDED
Used the v3_2d triple_diff card where w11 had no shape variations (missing step, see `propcorr_v3.md`,
"w11 systematics step"). Numbers below kept for the record; redo on `triple_diff_w11export`.

Ranges = 4x the individual 95% interval (cHDD +-0.066, cHj1 +-0.40, cHu +-0.50, cHd +-0.376, clj1 +-0.05),
200 points, ~1 min per operator (local, 1 core). All 200 points converged, smooth parabolas, best fit at 0.

| op | individual 68% | profiled 68% | ratio | individual 95% | profiled 95% | ratio |
|---|---|---|---|---|---|---|
| cHDD | ±0.0085 | ±0.0096 | 1.13 | ±0.0165 | ±0.0188 | 1.14 |
| cHj1 | ±0.044 | ±0.050 | 1.14 | ±0.100 | ±0.102 | 1.02 |
| cHu | ±0.065 | ±0.0725 | 1.12 | ±0.125 | ±0.1425 | 1.14 |
| cHd | ±0.046 | ±0.047 | 1.02 | [-0.090,+0.094] | [-0.092,+0.096] | 1.02 |
| clj1 | [-0.0055,+0.0065] | [-0.0058,+0.0068] | 1.04 | [-0.0085,+0.0125] | [-0.0113,+0.0143] | 1.21 |

Profiling is real (--savePOIs: the other POIs move within the 95% region, e.g. cHu up to +-0.06 along the
cHDD scan), but the widening (2-21%) is far below the stat-only Fisher estimate (up to x2.6). Likely reason:
the individual intervals are already systematics-dominated (cHj1 68% +-0.044 vs stat-only Fisher +-0.012),
and the operator correlations matter much less once the nuisances are profiled. To check: stat-only profiled
scans (`--stat`), which should reproduce the Fisher ratios. Plots `profiled/plots/scan_profiled_<op>.png`.

### test2 on triple_diff_w11export (2026-11-23, same 5 operators, ranges, 200 points)
LLR `propcorr_v3_2d/datacards/inc_mm/triple_diff_w11export/profiled/test2/` (run.sh, post.sh, plots/).
~1 min per operator, all 200 points converged, no error status. Intervals IDENTICAL to the old test2 at the grid
resolution (table above still holds); point by point the 2dNLL curves differ by < 0.11 (clj1: up to 1.5 at the
far edge, 2dNLL ~111). Expected: w11 enters as c_i*c_j, so its missing shape variations matter far from the
minimum (where the 2D scans broke), not near it. Fisher+systs reproduces these ratios (table above).

### set16 (2026-11-23, greedy Fisher+systs set, triple_diff_w11export)
`cHWB,cHbq,cHd,cHe,cHj1,cHj3,cHu,cQe,cQl1,cbBRe,cbWRe,cbe,cbl,ceu,cje,clj1` (cHDD, cHl3, cll1 dropped for the 3
flat directions). LLR `triple_diff_w11export/profiled/set16/` (metadata_profiled.json = 4x individual 95% from
v3_pdfcut 1D, condor/ = 16 short jobs, cluster 824068, post.sh). 200 points; 1-3 missing points in 6 scans, inside
95% only for cHWB (3, 2dNLL~2.5), cHe (2, ~4.6-5.1), cHj3 (1, ~0.6): no interval crossing falls in a gap.
EOS: https://aldufour.web.cern.ch/RECO3D_propcorr_v3/profiled_triple_diff/ (plots + intervals.txt).

| op | ind 68% | prof 68% | ratio 68% | ratio 95% | Fisher+systs |
|---|---|---|---|---|---|
| cHj3 | ±0.0232 | [-0.0373,+0.0591] | 2.08 | 1.84 | 1.8 |
| cbWRe | [-0.34,+0.39] | [-0.56,+0.90] | 2.02 | 1.57 | 2.3 |
| cHWB | ±0.0054 | [-0.0125,+0.0074] | 1.85 | 1.64 | 1.6 |
| ceu | ±0.0042 | [-0.0056,+0.0088] | 1.71 | 1.84 | 1.4 |
| cHe | ±0.0161 | [-0.0233,+0.0283] | 1.60 | 1.52 | 1.5 |
| cHbq | ±0.126 | [-0.170,+0.201] | 1.47 | 1.43 | 1.5 |
| cbBRe | [-0.53,+0.65] | [-0.81,+0.90] | 1.46 | 1.28 | 2.3 |
| cQe | [-0.17,+0.19] | ±0.254 | 1.41 | 1.32 | 1.7 |
| cQl1 | [-0.26,+0.19] | [-0.32,+0.29] | 1.36 | 1.30 | 1.4 |
| cbe | [-0.24,+0.19] | [-0.30,+0.28] | 1.36 | 1.29 | 1.3 |
| cHu | ±0.065 | ±0.0875 | 1.35 | 1.38 | 1.3 |
| cbl | [-0.25,+0.21] | [-0.30,+0.29] | 1.28 | 1.24 | 1.2 |
| cHj1 | ±0.044 | ±0.054 | 1.23 | 1.10 | 1.1 |
| clj1 | [-0.0055,+0.0065] | [-0.0068,+0.0072] | 1.17 | 1.31 | 1.1 |
| cje | [-0.0065,+0.0146] | [-0.0079,+0.0163] | 1.15 | 1.19 | 1.4 |
| cHd | ±0.046 | ±0.0508 | 1.10 | 1.10 | 1.1 |

All 16 constrained with others floating; widening 1.1-2.1, in line with the linear Fisher+systs (largest
differences cbBRe, cQe lower; ceu higher). Profiled intervals become asymmetric (quadratic terms + correlations).

## Next
- Check test2 against the Fisher prediction; then the 12-operator set (condor, `--condor` +
  `bundle_scan_jobs.py`, short queue), and stat-only profiled for comparison.
- Decide the treatment of the exactly flat EW-input directions (drop operators, or fit rotated combinations / PCA).
