# propcorr v3 — reference version (2026-11-12 to 11-19)

v3 is the current (intended final) analysis version: propcorr_v1 setup + MiNNLO theory uncertainties
+ a soft cut on unphysical PDF weights, and a second production with the 2D cross terms.
Details of individual steps live in the topical notes linked below; this note is the map.

## What v3 is

| | |
|---|---|
| Samples | own propcorr samples (SMEFTsim topU3l, MwScheme, linearised propagator corrections, fixed gridpacks, see `propagator_correction.md`), 2018 only, 59.74 fb^-1 |
| Channel / variables | inc_mm; `mll`, `rapll_abs`, `costhetastar`, `triple_diff` (280 unrolled bins) |
| Operators | 27 (cHQ1≡cHQ3 and cQl1≡cQl3 for this process: bit-identical weights, see `combine_2Dscans.md`) |
| Theory uncertainties | MiNNLO way: per-variation k = MiNNLO_V / LO_sm_V applied to sm and EFT templates (SM at NNLO, EFT/SM at LO); `analysis/spritz/post_process_minnlo.py` |
| PDF outliers | soft cut: LHEPdfWeight members with |ratio to nominal| > 100 reset to the nominal member (event kept); `analysis/spritz/runner_pdfcut.py`. Affects 0.27% of events in 2000toInf, <= 0.02% below 1.5 TeV. PDF unc. last mll bin 24% -> 4.5%, nominal yields unchanged (`minnlo_pdf_outliers.md`) |
| Nuisances | 13 + autoMCStats (10 0 1): lumi lnN 1.0084; QCDScale, PDFweight, alphaS, PSWeight, mu_reco, mu_idiso, mu_trig, PU, prefireWeight, tt_ptrw, rochester_stat, rochester_syst (shape) |
| Model | `AnomalousCouplingMorphing_comb` (sm, w1_op, wm1_op, w11_op1_op2); Asimov fits |
| Cross terms | `w11_*` for all 351 pairs, produced nominal only; the planned step giving them the sm relative variations was forgotten until 2026-11-19, now applied with `exportNuisances.py` (see "w11 systematics step" below). Use the `*_w11export` cards for every 2D / profiled fit |

## Productions (LLR `spritz_fabian/configs/`)

| dir | what | how |
|---|---|---|
| `propcorr_v3_pdfcut/` | 1D analysis: all samples, sm/w1/wm1 with all variations | `config.py` = propcorr_v1 + runner_pdfcut; 3000 jobs long queue; merge, `pp_minnlo.py`, cards, `run_1d.sh`, `summaries.sh`, `nuis_plots.sh` |
| `propcorr_v3_eft2d/` | signal-only run with the 351 cross terms | `analysis/spritz/config_propcorr_v3_eft2d.py` + `runner_2d.py` (batched EFT fill); `combine_2Dscans.md` "2D production v3" |
| `propcorr_v3_2d/` | 2D datacards = v3_pdfcut + w11 from the eft2d run | `analysis/spritz/merge_eft2d.py --bkg propcorr_v3_pdfcut --sig propcorr_v3_eft2d --out propcorr_v3_2d`, MiNNLO postproc, spritz-cards-eft |

| `propcorr_v3_2d/datacards/inc_mm/{mll,triple_diff}_w11export/` | **the 2D cards to use**: v3_2d cards + w11 shape variations copied from sm | `exportNuisances.py --from sm --to "w11_.*" --nuis ".*"` (tools, morphing_model) + copy `jsonComb.json`, `metadata.json` by hand; shapes file is `shapes.root_0` |

Datacards: `<dir>/datacards/inc_mm/<var>/{datacard.txt, shapes.root, jsonComb.json, metadata.json}`.
v3_pdfcut and v3_2d shapes are identical except for the w11 templates.
`mll_w11sm/` (same dir) = the same fix made with `analysis/combine_tools/add_w11_syst_from_sm.py`; bit-identical to
`mll_w11export` (0 differences in datacard entries and in all 9474 histograms), kept only for the cHd_cHbq test.

## w11 systematics step (forgotten, applied 2026-11-19)
- The plan was: produce w11 nominal only (filling all variations for 351 cross terms did not merge on condor) and then
  give w11 the SM relative variations, as for the other EFT templates. That last step was forgotten, so the v3_2d cards
  had w11 with only `lumi`; all 12 shape nuisances were `-` for w11.
  sm/w1/wm1 had all variations. When the nuisances move, sm/w1/wm1 shift but w11 does not, so the morphing sum
  (coefficient of w11 = c_a*c_b, large cancellations at large |c|) is no longer a consistent EFT prediction: the implied
  mixed term shifts by -(d sm + d lin + d quad) x c_a c_b, i.e. spurious freedom for the nuisances growing with c_a c_b.
  Giacomo: "questo rompe l'algebra".
- Symptoms (mll, worst; triple_diff mildly): ridges/stripes along c_a ~ 0 / c_b ~ 0 (continuous ridges of the profiled
  NLL, e.g. cHd_cHbq at cHd=5.44: 2dNLL 164 at cHbq=0, 912 at 0.20, 56 at 0.67, pulls up to 18 sigma), inflated and
  irregular contours, FASTEXIT/error status (logged), fallbacks and failed fits.
- How it was found: stat-only likelihood (all nuisances fixed, evaluated directly from the templates; combine cannot
  "fit" with zero free parameters) is perfectly smooth: the EFT model is fine, the nuisances are the problem.
- Fix: give w11 the sm relative variations per bin and per shape nuisance (w11_V = w11 x sm_V/sm), consistent with how
  the theory variations reach w1/wm1 (MiNNLO k-factor) and a good approximation for the experimental ones (same events).
  Official tool: `tools/systematic_tools/exportNuisances.py` (Giacomo, morphing_model c5d93cd, 2026-11-19).
- Test, cHd_cHbq mll, full 10k grid, -v 3: stripes gone, smooth banana-shaped contours (68% cHd -1..1.6, cHbq +1.5..-3.5;
  95% to cHbq ~ -11), failed fits 10 -> 0, fallbacks 42 -> 13, points far above neighbours 116 -> 3 (all at 2dNLL ~
  500-4500, far outside the contours). FASTEXITs remain only at 2dNLL >> 100 (huge EFT effects, nuisances pulled far):
  smaller scan boxes from the new contours should remove most of them.
- Superseded by the fix: all v3 2D scans and plots (`full/` in mll and triple_diff; the old EOS 2D folders were deleted
  2026-11-20, in the EOS recycle bin), the "badly behaved mll scans" and FASTEXIT diagnostics below, and the profiled test2
  (`profiled_fits.md`). The 1D results (v3_pdfcut, no w11) are NOT affected.
- Also learnt: `plot_2d_scans.py` upsampled the grid with a cubic spline, drawing fake islands around isolated outliers
  (one bad point at 2dNLL 13727 among ~900 -> spline undershoot to -820). Use `runPlots.py 2` (mkEFTScan); the 2D chain
  now does (2026-11-19).
- Production on the fixed cards: launched 2026-11-19 17:48, `run_2d_perpair.sh` on `{mll,triple_diff}_w11export`
  (tmux `w11_mll`, `w11_td`; `chain.log` in each card folder). Done 2026-11-19 20:19: 351/351 pairs in both, no open
  boxes. Raw FASTEXIT lines in pass2 job logs: mll 2047 (229 logs), triple_diff 387 (63 logs); location vs contours
  still to check. Overlays: `plot_2d_scans.py` (no upsampling) in `propcorr_v3_2d/datacards/inc_mm/overlay_w11/`.
- Debug material: LLR `~/work/scratch_fastexit/` (EOS `2D_mll/debug/` deleted with the old 2D folders, 2026-11-20).

## Fits and results

### 1D (`propcorr_v3_pdfcut/datacards/inc_mm/<var>/`)
- `run_1d.sh` (copy: `analysis/combine_tools/run_1d_v3.sh`): createCombineJson, `createWS.py 1`, initial + scans, `runPlots_compare.py --compare-stat`, summaries (`--logscale`).
- All 4 variables x 27 operators converge. Stat+syst: 100 points. Stat-only: 1000 points (2026-11-17; 100 was too coarse for the narrow stat-only curves, old files in `stat_100pts/` and `scan_old_stat100pts/`).
- Plots in `<var>/scan/`, nuisance plots in `<var>/nuis_combine/` (drawNuisances.py) and `<var>/nuis_eft/` (drawNuisancesEFT.py).
- v3 vs v1_minnlo: high-mass operators tighter, several bulk operators 5-10% looser, entirely from the (less constrained) PDF nuisance (`minnlo_pdf_outliers.md`, PDF constraint check).

### 2D (`propcorr_v3_2d/datacards/inc_mm/{mll,triple_diff}/full/`) — SUPERSEDED (missing w11 step, see above); kept for the record
- All 351 pairs, both variables, per-pair scan boxes: `run_2d_perpair.sh` (coarse 50x50 -> boxes from the 95% extent -> fine 100x100). `full/final/` merged scans, `full/plots/` plots (`plot_2d_scans.py`, mkEFTScan colours).
- Triple diff closes every mll near-flat direction (all 300 distinct pairs closed and compact); exceptions: the exact identities cHQ1/cHQ3, cQl1/cQl3. Remaining strong (closed) correlations among cHDD, cHWB, cHl1, cHl3 (axis ratio ~7).
- Badly behaved mll scans (27 pairs, hand-checked): pairs of Z-coupling operators (cHj1, cHj3, cHu, cHd, cHQ1, cHbq, cHe, cHl1, cll1) and b-dipoles. Their mll effect is 92-99.7% a Z-peak normalisation change (`check_norm_degeneracy.py`) -> degenerate with each other and with normalisation nuisances -> near-flat, syst-dominated, fits land in different minima. 4x more points (40k) does not help; strict minimizer + per-point minimum helps a little (`merge_min_scans.py`). Not caused by propcorr (v9 templates even more collinear).
- Negative yields (checked 2026-11-19, pass2 condor logs): combine prints "Number of events is negative" / "underflow to -nan" (always mll bin 0) in 190/351 mll pairs (25 of the 27 bad ones; worst cHd_cHbq, cbBRe_cHbq, cbWRe_cHbq, cbBRe_cHd), 53/351 in triple_diff with far fewer counts. The NOMINAL prediction (sm+lin+quad+mixed, all bins) never goes below 0.63x(sm+bkg) anywhere in the scan boxes, for any pair. At the warning points the nuisances have run away to |theta| ~ 1e2-1e7 (median QCDScale 9e3, alphaS 6e3, PDF 5e3 sigma): the minimizer diverges along the near-flat direction and only then hits negative yields. So: not too few points, not negative templates, but runaway fits. Possible fix to test: bound the nuisances (e.g. +-5 sigma via --setParameterRanges) on a few bad pairs.
- Nuisance-block test (2026-11-19, Giacomo's method), mll cHd_cHbq, full pass2 10k grid, exact pass2 command, only
  --freezeParameters changed (`~/work/scratch_fastexit/run.sh`). FASTEXIT counts: all floating 284; no nuisances 0;
  PDF frozen 274; PSWeight frozen 267; alphaS frozen 309; theory frozen (exp floating) 198; exp frozen (theory floating) 20;
  only lumi / only mu SFs / only PU+prefire / only rochester+tt_ptrw floating: 0 each; lumi+muSF 5; muSF+PU+pref 0;
  lumi+muSF+PU+pref 47; all exp except lumi 131. My conclusion at the time ("no single culprit, collective
  degeneracy") was WRONG: the test never froze tt_ptrw alone with everything else floating (see next point).
- Background-removal test (2026-11-19, Giacomo), same pair/grid, one frozen rateParam per background (0 = removed
  from model and Asimov; `~/work/scratch_fastexit/nobkg_run.sh`): all bkg 276; no bkg 1; DYtt only 2; TT only 152;
  all but TT 1; TT only + QCDScale frozen 227; + PDF frozen 234; + tt_ptrw frozen 11.
  ORIGINAL setup with only tt_ptrw frozen: 284 -> 1 FASTEXIT. => culprit = tt_ptrw (TT only), in combination with the
  other nuisances (alone it gives 0). Template looks healthy: symmetric, smooth, +2.1% (bin 0) to +11% (last bin).
  Bin 0 composition: sm 86%, DYtt 8.1%, TT 4.3%, rest <1%. Why tt_ptrw drives the runaway: not understood yet.
- With hindsight (same day): freezing tt_ptrw (or PSWeight+tt_ptrw) removed most logged errors but NOT the stripes;
  -v 3 showed that the logged errors mostly recover (status 0) while the stripe points are silent, continuous ridges of
  the NLL. Root cause = the missing w11 step above, not tt_ptrw/PSWeight.

### 2D w11export quality check (2026-11-20/21, `check_2d_scans.py`, output `full/check_2d.{txt,csv}`)
- "inside 95%" = neighbour median 2dNLL < 5.99. mll: 12997 missing points (2324 inside, 236 pairs), 915 error status,
  20 outliers inside. triple_diff: 2192 missing (64 inside), 223 error, 99 outliers inside (20 pairs, ~5 units above
  neighbours: cHl3 x cld/ced, cHQ1 x clu/cje/ceu) = the islands in the overlays.
- mll missing points (2026-11-21, `~/work/scratch_missing/`): failed minimizations, silently not written by combine
  (RandStartPt commits a point only if minimize() is ok; at -v 0 nothing is logged). Not jobs/chunks: spread uniformly
  over the 200-point chunks; 811/2324 isolated, 805 with one missing neighbour. They cluster at specific values of one
  operator, whatever the partner: cHl1 grid row 62 (k=0.048) in 9 pairs (also as y axis), cHj1 62, cHe 42/45, cHl3 39.
  Reproduced on cHl1_cld 6225 (k_cHl1=0.048, k_cld=-0.039, 2dNLL ~3.7 from neighbours): Migrad strategy 0 status 3,
  FVAL 0, EDM 4.6e5 after 56 calls, "No improvement in line search"; fallback fails the same way.
  Strategy 1 (main + fallback): still fails. Same point alone with --setParameters k at the point: converges
  (EDM 5e-8). Whole chunk 6200-6399 with start k_cHl1=0.048,k_cld=0: 200/200 points, BUT the 2dNLL surface changes
  by O(1) w.r.t. production (not a constant offset; also rows that fitted fine before, e.g. row 63 j=40 1.57 -> ~0.5,
  and points next to the minimum move by +-0.17) => the fit result depends on the starting POI values
  (= the state in which the snapshot is taken, `snap` is taken after nll0 = nll.getVal() in MultiDimFit::doGrid).
  Production (start k=0) is self-consistent with strategy 1 to ~0.1. Not understood; to settle before re-running points.

### 2D w11export refit of bad points (2026-11-21) — ON HOLD, ask Giacomo whether to pursue
- Single-point refits of the bad points inside 95% (missing / error / outlier from `check_2d_scans.py --list-out`), same
  command as production (same snapshot, start k=0) except the minimizer; `refit_2d_points.py prepare|merge`. Merge keeps,
  per grid node, the lower 2dNLL of production and refit; `final/` repointed to `full/<tag>/merged/` (previous targets in
  `full/<tag>/final_before.txt`, production pass2/pass3 files untouched). Then runPlots, overlays, EOS (same 3 folders,
  overwritten). Chains: `inc_mm/refit_chain.sh`, `refit2_chain.sh` (logs next to them).
- Settings tested on single points (`~/work/scratch_refit_test/`, round2/):
  - production (strategy 0, fallback Migrad0:0.2) reproduces the failures.
  - F = strategy 2 + fallbacks Migrad0:0.2 -> Simplex0:0.1 -> Migrad1:0.2: fixes most outliers and missing points.
    No effect: MINIMIZER_analytic off, --cminPreScan.
  - G = F + --cminDefaultMinimizerTolerance=0.01, fallbacks tight first: Migrad1:0.01 -> Simplex0:0.01 -> Migrad2:0.05
    -> Migrad0:0.2. Reason: in F the loose Migrad0:0.2 fallback returns a "valid" false minimum, so the better ones never
    run. G recovered all 12 test points that F lost (values = neighbour median within ~0.05).
- Round 1 (F, `full/refit/`): mll 2410 points (249 pairs), triple_diff 168 (29 pairs); 0 failed jobs.
  mll: 2305 recovered, 60 lowered, 26 not lower. triple_diff: 64 recovered, 99 lowered, 5 not lower.
- Round 2 (G, `full/refit2/`, cQl1_cQl3 excluded: flat by cQl1 = cQl3): mll 36 points, triple_diff 23.
  mll: 30 recovered, 6 lowered. triple_diff: 7 recovered, 14 lowered, 2 not lower.
- After round 2 (inside 95%; "error status" counts in `check_2d_refit2.txt` refer to the pre-refit pass2 logs, ignore):

  | | mll before -> after | triple_diff before -> after |
  |---|---|---|
  | missing | 2324 -> 0 | 64 -> 2 |
  | outliers > 1 | 20 -> 12 (8 in cQl1_cQl3) | 99 -> 23 (7 pairs) |

- Left, not understood: (a) isolated triple_diff spikes ~+5 (cHQ1_clu i,j=77,40: 9.52 vs neighbours 4.5; cHl1_ced
  21,57: 10.05 vs 5.4) identical under 8 settings (F, G, tolerance, Combined, preScan, no analytic, no stepSize,
  no FITTER_BOUND); minimizer never converges (EDM 0.16 -> 0.009), only the loose last fallback returns a value.
  Warm start from a neighbour is awkward: with -t -1 a different snapshot changes the Asimov. (b) mll cbWRe_clj1: two
  interleaved branches (~3.5 and ~5.5-7) over a band of points; the outlier check catches only part of it.
  (c) relation to the start-value dependence above (missing points section): the refits use the production start, so
  they are consistent with production, but the dependence itself is open.
- Rejected submission: the first prepare wrote a scan.sub still queueing pass2's joblist (path spelled /grid_mnt vs
  ~/work); 1404 jobs removed while idle, none ran; only condor event lines appended to `pass2/condor_scan_bundled/logs/*.log`.
  Fixed (scan.sub written explicitly, asserts no pass2 path).

### Profiled fits
- In progress on triple_diff, see `profiled_fits.md`. test2 used the old v3_2d card (w11 without variations):
  redo on `triple_diff_w11export`.

## Template checks
- Positivity (`check_positivity_morphing.py`): every full template (sm, w1, wm1, w11) >= 0 in every bin, all variations.
- quad = (w1+wm1)/2 - sm < 0 for cHj3 and cll1 in the Z-peak mll bin (-3e-4 to -4e-4 of SM), also at LHE level and only in propcorr: generator-level, from `NPall<=2` allowing two EFT insertions in the amplitude (quad = |A1|^2 + 2 Re(A0* A2)). Plot `checks/quad_terms_mll.png` (`plot_quad_terms.py`). Question to Giacomo/Fabian open (`todo.md`).
- 2D quadratic form PSD per bin (`check_quad_psd.py`): violations only from the cHj3/cll1 Z-peak bin.
- Negative background bins: DYtt ~ -0.1 events in 3 triple_diff bins, TT ~ -0.03 in 2 (open question, `todo.md`).

## EOS (https://aldufour.web.cern.ch/RECO3D_propcorr_v3/)
- top: `eft_summary_multivar.*`
- `<var>/`: `eft_summary_two_panel.*`, `scan/` (1D), `dist/`, `nuis_combine/`, `nuis_eft/`, `prefit/` (datacard prefit, Asimov, no EFT, +1sigma shift per nuisance; `plot_prefit_datacard.py`, 2026-11-19)
- `2D_mll/`, `2D_triple_diff/` (runPlots.py, 351 each), `2D_mll_vs_triple_diff/` (overlay_<pair>, 351): fixed `*_w11export`
  cards (2026-11-20; uploaded as `*_w11`, renamed after deleting the old folders, which are in the EOS recycle bin);
  overwritten 2026-11-21 with the refit grids (rounds 1+2, see "refit of bad points")
- `checks/`: `quad_terms_mll.png`, `mll_2d_10k_vs_40k_vs_strict.png`
- Push: `eos_tools/push_v3_to_eos.sh` (1D + summaries, `--overwrite`), 2D via `run_2d_perpair.sh`; all through `pushtoEOS_parallel.py`.

## Environment gotcha
`dy_combine_morphing` (`analysis/combine_tools/env_llr_morphing.sh`) used to put stale May copies of runScans.py /
createWS.py / createCombineJson.py / createJson.py / runPlots.py (CMSSW_spritz `tools/`) before the repo in PATH.
~/.bashrc functions masked this for typed commands, but nice/xargs/condor/subprocess got the old copies.
Fixed 2026-11-17: repo `combine_tools` first in PATH.
