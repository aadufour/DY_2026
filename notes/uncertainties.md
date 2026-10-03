# Systematic uncertainties — how they are built and applied

Scope: lumi, QCD scale, PDF (+ alpha_s), muon efficiency SFs (reco, ID+iso, trigger).
Reference setup: `eft_bkg_fullsyst_v9` (2018, `inc_mm`, mll). `propcorr_v1` uses the
same config except the EFT sample names (see "Config check" at the bottom).

Code references:
- nuisance definitions: `analysis/spritz/config_v9.py` (`nuisances = {...}`)
- event weights: `analysis/spritz/runner.py` (weight product ~l.353)
- Up/Down construction + k-factor: `analysis/spritz/post_process.py`
- spritz modules (on LLR): `spritz_fabian/src/spritz/modules/{theory_unc,lepton_sf,trigger_sf}.py`

## Common pipeline

1. Per event: nominal weight + alternative weights (one per variation).
2. Per bin: spritz fills one histogram per variation; postproc reduces them to
   `histo_<proc>_<syst>Up/Down` according to the nuisance `kind`:
   - `weight`: take the Up/Down histogram directly
   - `envelope`: per-bin max / min over the variations
   - `square`: nom ± sqrt(sum_i (X_i - X_nom)^2)
   - `stdev`: nom ± std over the variations
3. k-factor (MiNNLO / LO sm, computed from nominal only) multiplies nominal AND all
   Up/Down of sm, w1_*, wm1_* → relative variations are unchanged by the k-factor.
4. Combine: each `shape` nuisance = one theta ~ N(0,1), shared by all processes it is
   listed for. Split into a normalisation part (asymmetric lnN, kappa = N_Up/N_nom,
   N_Down/N_nom) and a shape part (vertical morphing, quadratic for |theta|<1, linear
   beyond).

Correction vs uncertainty: a correction changes the nominal prediction (e.g. muon SF,
k-factor); its uncertainty is the systematic that enters the fit as a nuisance.

## Luminosity — `lumi`

- `lnN 1.0084` on all MC processes (backgrounds, sm, all EFT templates).
- Pure normalisation, one Gaussian-constrained (log-normal) nuisance.
- Origin: `data/common/lumi.json` in spritz, `Full2018v9.rel_unc = 1.0084` (same file
  that gives `tot` = 59561.26 pb^-1). Per-year table: 2016 HIPM/noHIPM 1.0123,
  2017 1.0082, 2018 1.0084. Giacomo's per-year DY configs read it as
  `lumi_unc = lumis[year]["rel_unc"]`; our config hardcodes the same number
  (originally copied from Giacomo's `dy-genlevel`). Verified 2026-10-03.
- TODO: confirm with Giacomo / current LUM POG Run 2 table that these are the
  recommended values (older recommendation was 2.5% for 2018).
- TODO: read `rel_unc` from the json instead of hardcoding (needed for Run 2).
- Plot label says 59.74 fb^-1, normalisation uses 59.56 fb^-1 → make consistent.
- Run 2: a single `lumi` name across years = fully correlated; the official scheme
  has correlated + uncorrelated components.

## QCD scale — `QCDScale`

- Source: `LHEScaleWeight` in each sample's NanoAOD; `theory_unc` registers
  `QCDScale_i = weight × LHEScaleWeight[i]`.
- `kind: envelope` → Up/Down = per-bin max/min over the selected weights.
  - backgrounds (TT, WW, Single_Top): 7 weights (drop the two anti-correlated
    muR/muF combinations)
  - EFT samples (8 weights stored): indices [0,1,3,5,7]
- At LO with QCD=0 there is no alpha_s in the matrix element → muR has no effect.
  Whatever the weight ordering, the chosen indices include muF = 0.5 and muF = 2, so
  for sm/EFT the envelope is effectively a 2-point muF variation.
- Applied to: sm + all EFT, DYtt, TT, WW, Single_Top. Not GGToLL, WZ, ZZ (no LHE
  scale weights available).
- Reviewer question to prepare: SM is normalised to NNLO (MiNNLO k-factor) but its
  scale uncertainty is LO muF only (no muR, likely larger than NNLO scale unc.).

## PDF — `PDFweight`

- PDF set: NNPDF3.1 NNLO, LHAPDF 325300 (`NNPDF31_nnlo_as_0118_mc_hessian_pdfas`),
  `ErrorType: symmhessian+as`, 103 members:
  - mem 0: central value
  - mem 1–100: symmetric Hessian eigenvectors
  - mem 101 / 102: central value with alpha_s(mZ) = 0.116 / 0.120
- Verified chain (2026-10-03):
  1. gridpack (MG5) writes per-event weights for all 103 members into the LHE
  2. EFT NanoAOD branch `LHEPdfWeight`: "LHA IDs 325300 - 325402", n = 103
  3. `theory_unc`: `PDFWeight_i = weight × LHEPdfWeight[:, i]` (index to index)
  4. config: members 0–100 → `PDFweight`, 101/102 → `alphaS`
  → PDF and alpha_s uncertainties on sm/EFT come entirely from our own SMEFTsim
  samples (nothing from MiNNLO / other samples).
- `kind: square`: sigma = sqrt(sum_{i=1}^{100} (X_i - X_0)^2), Up/Down = nom ± sigma
  (symmetric). This is the correct prescription for a symmetric Hessian set.
  - NOT "replicas split into up/down by sign": the sign of a Hessian eigenvector is
    arbitrary, so sorting deviations by sign is meaningless (and underestimates).
  - Single nuisance → all bins move coherently. Rigorous alternative: 100 nuisances,
    one per eigenvector (keeps the true bin-to-bin correlations).
- Applied to: sm + all EFT, DYtt, TT, WW. Not Single_Top (unreliable LHEPdfWeight),
  GGToLL, WZ, ZZ.

## alpha_s — `alphaS`

- Members 101 / 102 (alpha_s = 0.116 / 0.120, i.e. ±0.002), separate nuisance from PDF.
- `kind: envelope` → Up/Down = per-bin max/min of the two.
  - Cleaner: Up = mem 102, Down = mem 101 (alpha_s has a physical sign; an envelope
    could flip the meaning of theta=+1 between bins). In practice identical here:
    mem 102 is above mem 101 in every bin.
  - PDF4LHC recommends ±0.0015 (scale deviations by 0.75); cosmetic at ~1%.
- Size on sm: ~+1% / -1.2%→0% across mll, mild shape.
- Applied to: sm + all EFT, DYtt.

## Why the theory ratios are identical across sm / w1 / wm1

Each EFT template weight = w_evt × r_k (EFT reweighting factor). The PDF/scale weight
ratio for an event depends only on x1, x2, muF (same kinematics for all templates), not
on the Wilson coefficients → Up/nom is the same for every template, to all digits.
All templates share one nuisance per systematic (100% correlated by construction).

## Muon efficiency scale factors — `mu_reco`, `mu_idiso`, `mu_trig`

### What they are
- Muon POG (Physics Object Group) tag-and-probe on Z→mumu, run on 2018 data and DY MC,
  per muon in (eta, pT) bins: SF = eps_data / eps_MC.
- Factorised chain, each step conditional on the previous:
  eps_tot = eps_reco × eps_ID|reco × eps_iso|ID × eps_trig|ID,iso → SFs multiply.
- Year-specific (alignment, chamber conditions, trigger menu: 2018 = IsoMu24). Each
  year's `cfg.json` points to its own SF file (`cfg["leptonSF"]`, correctionlib).
- The correctionlib map returns `nominal`, `stat`, `syst` per (eta, pT) bin; the
  uncertainty is entirely from the POG.

### How they are applied (spritz_fabian version — the one used for v9)
- Nominal: every MC event (all backgrounds, sm, every EFT template; not data):
  `weight = genWeight × puWeight × topPtWeight × RecoSF × TightSF × prefireWeight × TriggerSFweight_2l`
  (runner.py ~l.357), then EFT template k: `weight × r_k` (runner.py ~l.385).
- SFs depend only on the reconstructed muons → identical for all templates of an
  event → factor out of LIN/QUAD consistently.
- No envelope: `kind: weight`, the POG gives a single ±1 sigma per muon.

| nuisance | per-muon SF | per-muon sigma | Up / Down |
|---|---|---|---|
| mu_reco | `Muon_RecoSF` | sqrt(stat^2 + syst^2) | SF ± sigma, both muons |
| mu_idiso | `Muon_IdSF_<WP>` × `Muon_IsoSF_<WP>` | sqrt(SF_id^2 sig_iso^2 + SF_iso^2 sig_id^2) | SF ± sigma, both muons |
| mu_trig | `Muon_TriggerSF_tightId` (trig-matched muons) | sqrt(stat^2 + syst^2) | event SF ± event sigma |

- eta clamped to ±2.4, pT clamped to ≥15 GeV (≥26 for trigger), no upper pT clamp →
  muons above the last pT bin of the map get the last bin's SF and sigma.
- Up/Down: SF of every muon in every event shifted together → one nuisance, fully
  correlated across muons, bins and processes (conservative for the stat part).
- Example (ID+iso plateau): SF = 0.990, sigma = 0.0012 → 0.990^2 = 0.9801 vs
  0.9912^2 = 0.9825 → +0.24%.

### Size on sm (Up/nom - 1)
- mu_reco: flat 0.13%
- mu_idiso: 1.0% (50–60), 0.5% (60–80), ~0.25% plateau above
- mu_trig: 0.01–0.06% (artificially small, see below)

### Open points
- Trigger SF combination for 2 matched muons uses `1 - (1-sf1)(1-sf2)`: SFs plugged in
  where efficiencies belong. Correct:
  SF_evt = [1-(1-eps1_data)(1-eps2_data)] / [1-(1-eps1_MC)(1-eps2_MC)].
  Error propagated as sqrt(((1-sf2) sig1)^2 + ((1-sf1) sig2)^2) → suppressed ~50×.
  Numerically small (correct ~0.08% vs code ~0.01% for dsf = 0.5%), but formally
  wrong. Same formula in Giacomo's version.
- High pT: Z tag-and-probe has few muons above ~200 GeV; mll goes to 3 TeV. Check
  `muWP` (high-pT ID with momentum-binned SFs?) and the extrapolation uncertainty.
- Same-data use: SFs are measured on 2018 Z→mumu, i.e. largely our Z peak. Fine for
  the nominal (efficiency is a detector property, Z yield cancels), but the fit can
  constrain mu_idiso/mu_reco/mu_trig from our 38M-event Z peak → double counting.
  Check post-fit pulls/constraints.
- Full Run 2: per-year SFs; decide correlation across years (stat uncorrelated,
  syst correlated).

## MC statistics (for reference)

`autoMCStats 10 0 1`: Barlow-Beeston-lite, backgrounds only. sm/EFT templates have
`noStat=True` (bin variances zeroed) and `includeSignal=0`: they are reweightings of
the same events, so independent per-template stat nuisances would be wrong.

## Config check: v9 vs propcorr_v1

- Local `config_v9.py` vs `config_propcorr_v1.py`: identical except EFT sample names
  (`DYSMEFTsim_LO_mll_*` → `DYSMEFTsim_LO_propcorr_mll_*`). Same nuisances, binning,
  backgrounds, runner path.
- Local git: runner.py unchanged since 2026-07-21 (before both productions);
  post_process.py / make_cards.py changed 2026-09-07 only to add cross-term (`w11_*`)
  handling.
- propcorr NanoAOD `LHEPdfWeight`: "LHA IDs 325300 - 325402", n = 103 (verified
  2026-10-03) → same PDF/alpha_s chain as v9.
- LLR copies of runner.py, post_process.py, make_cards.py and both config.py files
  have the same md5 as the local repo (verified 2026-10-03).
- v9 `datacards/inc_mm/mll/datacard.txt` vs propcorr_v1
  `datacards_single/inc_mm/mll/datacard.txt`: same processes, same nuisance names, same
  nuisance lines (type + per-process values), same `autoMCStats 10 0 1`.
→ Everything in this note applies to propcorr_v1 as well.
