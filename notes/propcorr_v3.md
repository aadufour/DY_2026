# propcorr v3 — reference version (2026-11-12 to 11-17)

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
| Cross terms | `w11_*` for all 351 pairs, NOMINAL ONLY (no systematic variations on the mixed term) |

## Productions (LLR `spritz_fabian/configs/`)

| dir | what | how |
|---|---|---|
| `propcorr_v3_pdfcut/` | 1D analysis: all samples, sm/w1/wm1 with all variations | `config.py` = propcorr_v1 + runner_pdfcut; 3000 jobs long queue; merge, `pp_minnlo.py`, cards, `run_1d.sh`, `summaries.sh`, `nuis_plots.sh` |
| `propcorr_v3_eft2d/` | signal-only run with the 351 cross terms | `analysis/spritz/config_propcorr_v3_eft2d.py` + `runner_2d.py` (batched EFT fill); `combine_2Dscans.md` "2D production v3" |
| `propcorr_v3_2d/` | 2D datacards = v3_pdfcut + w11 from the eft2d run | `analysis/spritz/merge_eft2d.py --bkg propcorr_v3_pdfcut --sig propcorr_v3_eft2d --out propcorr_v3_2d`, MiNNLO postproc, spritz-cards-eft |

Datacards: `<dir>/datacards/inc_mm/<var>/{datacard.txt, shapes.root, jsonComb.json, metadata.json}`.
v3_pdfcut and v3_2d shapes are identical except for the w11 templates.

## Fits and results

### 1D (`propcorr_v3_pdfcut/datacards/inc_mm/<var>/`)
- `run_1d.sh` (copy: `analysis/combine_tools/run_1d_v3.sh`): createCombineJson, `createWS.py 1`, initial + scans, `runPlots_compare.py --compare-stat`, summaries (`--logscale`).
- All 4 variables x 27 operators converge. Stat+syst: 100 points. Stat-only: 1000 points (2026-11-17; 100 was too coarse for the narrow stat-only curves, old files in `stat_100pts/` and `scan_old_stat100pts/`).
- Plots in `<var>/scan/`, nuisance plots in `<var>/nuis_combine/` (drawNuisances.py) and `<var>/nuis_eft/` (drawNuisancesEFT.py).
- v3 vs v1_minnlo: high-mass operators tighter, several bulk operators 5-10% looser, entirely from the (less constrained) PDF nuisance (`minnlo_pdf_outliers.md`, PDF constraint check).

### 2D (`propcorr_v3_2d/datacards/inc_mm/{mll,triple_diff}/full/`)
- All 351 pairs, both variables, per-pair scan boxes: `run_2d_perpair.sh` (coarse 50x50 -> boxes from the 95% extent -> fine 100x100). `full/final/` merged scans, `full/plots/` plots (`plot_2d_scans.py`, mkEFTScan colours).
- Triple diff closes every mll near-flat direction (all 300 distinct pairs closed and compact); exceptions: the exact identities cHQ1/cHQ3, cQl1/cQl3. Remaining strong (closed) correlations among cHDD, cHWB, cHl1, cHl3 (axis ratio ~7).
- Badly behaved mll scans (27 pairs, hand-checked): pairs of Z-coupling operators (cHj1, cHj3, cHu, cHd, cHQ1, cHbq, cHe, cHl1, cll1) and b-dipoles. Their mll effect is 92-99.7% a Z-peak normalisation change (`check_norm_degeneracy.py`) -> degenerate with each other and with normalisation nuisances -> near-flat, syst-dominated, fits land in different minima. 4x more points (40k) does not help; strict minimizer + per-point minimum helps a little (`merge_min_scans.py`). Not caused by propcorr (v9 templates even more collinear).

### Profiled fits
- In progress on triple_diff, see `profiled_fits.md`.

## Template checks
- Positivity (`check_positivity_morphing.py`): every full template (sm, w1, wm1, w11) >= 0 in every bin, all variations.
- quad = (w1+wm1)/2 - sm < 0 for cHj3 and cll1 in the Z-peak mll bin (-3e-4 to -4e-4 of SM), also at LHE level and only in propcorr: generator-level, from `NPall<=2` allowing two EFT insertions in the amplitude (quad = |A1|^2 + 2 Re(A0* A2)). Plot `checks/quad_terms_mll.png` (`plot_quad_terms.py`). Question to Giacomo/Fabian open (`todo.md`).
- 2D quadratic form PSD per bin (`check_quad_psd.py`): violations only from the cHj3/cll1 Z-peak bin.
- Negative background bins: DYtt ~ -0.1 events in 3 triple_diff bins, TT ~ -0.03 in 2 (open question, `todo.md`).

## EOS (https://aldufour.web.cern.ch/RECO3D_propcorr_v3/)
- top: `eft_summary_multivar.*`
- `<var>/`: `eft_summary_two_panel.*`, `scan/` (1D), `dist/`, `nuis_combine/`, `nuis_eft/`
- `2D_mll/`, `2D_triple_diff/` (scan_<pair>, 351 each), `2D_mll_vs_triple_diff/` (overlay_<pair>, 351)
- `checks/`: `quad_terms_mll.png`, `mll_2d_10k_vs_40k_vs_strict.png`
- Push: `eos_tools/push_v3_to_eos.sh` (1D + summaries, `--overwrite`), 2D via `run_2d_perpair.sh`; all through `pushtoEOS_parallel.py`.

## Environment gotcha
`dy_combine_morphing` (`analysis/combine_tools/env_llr_morphing.sh`) used to put stale May copies of runScans.py /
createWS.py / createCombineJson.py / createJson.py / runPlots.py (CMSSW_spritz `tools/`) before the repo in PATH.
~/.bashrc functions masked this for typed commands, but nice/xargs/condor/subprocess got the old copies.
Fixed 2026-11-17: repo `combine_tools` first in PATH.
