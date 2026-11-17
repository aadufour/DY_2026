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

### test2 (2026-11-17, same 5 operators) — running
Ranges = 4x the individual 95% interval (cHDD +-0.066, cHj1 +-0.40, cHu +-0.50, cHd +-0.376, clj1 +-0.05),
200 points, tmux `profiled` on LLR, log `profiled/profiled_test.log`.

## Next
- Check test2 against the Fisher prediction; then the 12-operator set (condor, `--condor` +
  `bundle_scan_jobs.py`, short queue), and stat-only profiled for comparison.
- Decide the treatment of the exactly flat EW-input directions (drop operators, or fit rotated combinations / PCA).
