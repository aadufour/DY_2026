# Systematic uncertainties — how they are built and applied

Scope: lumi, QCD scale, PDF (+ alpha_s), parton shower, muon efficiency SFs (reco, ID+iso, trigger),
muon momentum (Rochester), pileup, L1 prefiring, top pT reweighting, MC statistics.
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

## Background samples (2018 UL)

| Background | Sample(s) | Generator / order | sigma [pb] |
|---|---|---|---|
| DYtt | DYJetsToTauTau_M-50_AtLeastOneEorMuDecay | POWHEG MiNNLO + Pythia8 + Photos (NNLO QCD) | 1164 |
| TT | TTTo2L2Nu | POWHEG NLO + Pythia8 | 89.3 |
| Single_Top | s-ch (aMC@NLO 4f, leptonic), t-ch top/antitop (POWHEG 5f), tW top/antitop (POWHEG, no fully had.) | NLO + Pythia8 | 3.4 / 134 / 21.7 ... |
| WW | WWTo2L2Nu | POWHEG NLO + Pythia8 | 11.1 |
| WZ | WZ_TuneCP5_13TeV-pythia8 (inclusive) | Pythia8 LO | 27.6 |
| ZZ | ZZ_TuneCP5_13TeV-pythia8 (inclusive) | Pythia8 LO | 12.2 |
| GGToLL | GGToMuMu_Pt-5, 5 mll bins × {El-El, Inel-El, Inel-Inel} | CepGen + LPair (+ Pythia6 for dissociation) | e.g. 0.25 / 0.37 / 0.56 (50–200) |

- DY mumu MiNNLO (DYJetsToMuMu_*, 11 mass bins) is not in the card: replaced by sm + EFT,
  used only for the k-factor.
- All backgrounds normalised from MC: sigma × 59.56 fb^-1. No control regions or
  rateParams.
- Source: config_v9.py (datasets/samples), samples_fabian.json.

## Nuisance × process matrix (datacard, identical for v9 and propcorr_v1)

| Nuisance | GGToLL | Single_Top | TT | WW | WZ | ZZ | DYtt | sm + EFT |
|---|---|---|---|---|---|---|---|---|
| QCDScale | – | ✓ | ✓ | ✓ | – | – | ✓ | ✓ |
| PDFweight | – | – | ✓ | ✓ | – | – | ✓ | ✓ |
| alphaS | – | – | – | – | – | – | ✓ | ✓ |
| PSWeight | – | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| mu_reco, mu_idiso, mu_trig | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| rochester_stat, rochester_syst | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| PU, prefireWeight | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| tt_ptrw | – | – | ✓ | – | – | – | – | – |
| lumi | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| MC stat | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – |

Rule: experimental nuisances on every process; theory variations wherever the sample
provides the generator weights. Gaps:

| Gap | Reason |
|---|---|
| GGToLL: no QCDScale/PDF/alpha_s/PS | CepGen/LPair: photon fluxes (form factors / structure functions), no PDFs, no Pythia8 shower weights |
| WZ, ZZ: no QCDScale/PDF/alpha_s | Pythia8-only LO samples, no LHE weights |
| Single_Top: no PDF/alpha_s | unreliable LHEPdfWeight in its NanoAOD (config comment) |
| TT, WW, Single_Top: no alpha_s | CHOICE: alphaS configured only for DY-type samples + EFT; TT/WW have members 101/102 → one-line change to add |
| tt_ptrw: TT only | correction defined for ttbar |
| MC stat: not on sm/EFT | same events reweighted (covariance treatment in progress) |

Reviewer questions:
- No cross-section lnN on any background (WZ, ZZ, GGToLL have only lumi + experimental
  normalisation freedom; TT/WW only scale/PDF). Same in dy-eft-2018 → framework
  convention, but expect the question (few %–10% for VV, more for gamma gamma → mu mu).
- Nonprompt (fake) muons missing: no W+jets/QCD samples, MC leptons gen-matched to
  prompt (`skip_genmatching: False`), fake template not in the card yet (Fabian's fake
  templates, 1 Oct; same-sign region `inc_mm_ss` exists in the config).
- Inclusive LO WZ/ZZ samples.

## Comparison with Giacomo's `dy-eft-2018` config (2026-10-04)

Path: `/grid_mnt/data__data.polcms/cms/adufour/spritz_giacomo/configs/dy-eft-2018/config.py`

Same in both (→ framework convention, nothing to defend):

| Nuisance | Ours | dy-eft-2018 |
|---|---|---|
| lumi | lnN 1.0084, all MC | lnN lumi_unc (= 1.0084), all MC |
| MC stat | autoMCStats 10 0 1 | identical |
| PU, prefireWeight | all MC | all MC |
| mu_reco, mu_trig | all MC | all MC |
| rochester_syst | sets 2–4, quadrature | sets 2–4, quadrature |
| tt_ptrw | TT only | TT only |
| background sigma lnN | none | none |

Different:

| Item | Ours | dy-eft-2018 | Comment |
|---|---|---|---|
| Muon ID/iso | merged mu_idiso | separate mu_id, mu_iso | framework version (Fabian's vs Giacomo's lepton_sf) |
| rochester_stat | yes (100 replicas) | no | his module needs do_stat=True; ours more complete |
| Theory | QCDScale, PDF, alpha_s, PS (where weights exist) | none of these | instead NLO_EW and N3LO_QCD correction uncertainties, DY only (DYmm_NNLO, DYtt for EW) |
| b-tag, jet PU-ID | none | btagSF_sf, btagSF_eff, puidSF | his selection uses jets (b-veto against ttbar) |
| Fakes | none | data-driven: fakes_param (fit) + fakes_model | the gap in ours |

Takeaways for SMP-V:
1. Fakes: most visible gap → "in progress".
2. b-veto: TT reaches ~25% of the total at 200–500 GeV in ours; a b-veto would cut TT
   and its uncertainties (adds b-tag SF nuisances). Discuss with Giacomo.
3. Higher-order DY: reference corrects DY to N3LO QCD + NLO EW with uncertainties on
   the corrections; ours stops at the NNLO QCD MiNNLO k-factor (meeting 28 Sept item).

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

## Parton shower — `PSWeight`

Slide one-liner: "PSWeight: Pythia8 ISR/FSR scale variations (×2, ×½), envelope of the
4 weights, all MC except γγ → μμ."

### What it is
- MG5 sample is pure LO (q qbar → mu mu, no extra partons). Pythia8 (CP5) adds:
  - ISR: radiation from the incoming partons → recoil → gives the dimuon pT_ll
  - FSR: QCD radiation from outgoing coloured partons
  - hadronisation, underlying event
- QED radiation off the muons is done by Photos++, not by the Pythia shower.
- The shower evaluates alpha_s at a scale ~ emission pT; the choice is arbitrary (like
  muR at fixed order) → vary it.

### No (QCD) FSR in LO DY
- The hard final state is two colourless muons → no QCD FSR from the hard process.
- Pythia's FSR still acts on partons created by ISR (an ISR gluon showers further into
  a jet) → changes ISR-jet structure, barely touches the muons (at most tiny isolation
  effects).
- QED FSR off the muons (radiative tail below the Z peak) is real but done by Photos;
  PS weights vary only the QCD alpha_s scale → not covered by PSWeight.
- → In practice PSWeight is an ISR variation for our signal; FSR weights should give
  ratios ≈ 1. Also the answer to "why merge ISR and FSR": the FSR part is ~empty.

### How the variations are built
- Pythia computes per-event weights on the fly (same shower history, alpha_s scale in
  ISR or FSR × 2 or × ½). NanoAOD branch `PSWeight`, 4 entries:
  [0] ISR×2 FSR×1, [1] ISR×1 FSR×2, [2] ISR×½ FSR×1, [3] ISR×1 FSR×½.
- `theory_unc`: `PSWeight_i = weight × PSWeight[:, i]` (only if 4 entries).
- `kind: envelope` → Up = per-bin max, Down = per-bin min over the 4 → ONE nuisance for
  ISR + FSR (index order irrelevant).

### Applied to
- All processes except GGToLL (no PS weights): DYtt, TT, WW, WZ, ZZ, Single_Top, sm,
  every EFT template. One nuisance, fully correlated (all Pythia8 CP5).
- EFT: PS weight belongs to the event's shower history, not to the WCs → identical
  for all templates of an event → same ratios across sm/w1/wm1.

### Size on sm (Up/nom - 1)
| mll [GeV] | 50–60 | 60–80 | 80–100 | 100–120 | 120–180 | 180–270 | 270–700 | 700–1000 | 1000–3000 |
|---|---|---|---|---|---|---|---|---|---|
| PSWeight | 3.75% | 0.16% | 0.19% | 0.56% | 0.27–0.32% | 0.68–0.83% | 0.13–0.16% | 0.80% | 1.6–1.7% |

- Shower is ~unitary → effects only through kinematics/acceptance. At LO+PS all of
  pT_ll comes from ISR → harder/softer recoil → muon pT shifts → threshold crossings.
  mll barely moves; cos theta* (Collins–Soper) built to be insensitive to pT_ll.
- 50–60 GeV (3.75%): muons at ~25–30 GeV, right at the 26 GeV trigger threshold.
  Largest single experimental-type effect.
- Middle range 0.1–0.8%, irregular → partly statistical.
- Above 1 TeV ~1.6%: probably partly statistical (event-by-event weight fluctuations,
  reduced effective SM statistics at high mass).
- Up/nom > 1 always: envelope max.

### Reviewer questions
- ISR and FSR merged in one envelope; common practice = two nuisances (Up = ×2,
  Down = ×½). Envelope loses the physical direction (theta=+1 can be ISR-up in one bin,
  ISR-down in another). Answer: FSR is ~empty for DY, so merged ≈ ISR alone.
- LO+PS signal: no real emission in the ME → the shower does all of pT_ll → larger PS
  uncertainty than an NLO-matched sample (MiNNLO). Motivation for generating
  `p p > l+ l- j` (meeting 23 Sept).
- One nuisance shared across generators (MG5 LO, MiNNLO, POWHEG): justified by the
  common shower (Pythia8 CP5).

## Why the theory ratios are (practically) identical across sm / w1 / wm1

- Per event the weight factorises: w_{k,i} = w_gen × r_k × p_i, with r_k one of the
  406 EFT weights (LHEReweightingWeight) and p_i one of the 103 PDF weights
  (LHEPdfWeight). Only 406 + 103 numbers are stored; spritz forms the products on the
  fly (theory_unc varies `weight`, the template fill multiplies by r_k, runner.py
  ~l.385) → no need for 406 × 103 weights.
- Exact at LO: dsigma = f(x1,muF) f(x2,muF) × |M(c)|^2. p_i depends only on x1, x2,
  muF (not on the WCs); r_k only on the hard kinematics and couplings (not on the PDF).
  Same for muF scale weights. Would break only if the EFT changed alpha_s-dependent
  pieces of the ME (muR at higher order); not the case with QCD=0.
- Per bin: Up_k/nom_k = sum_e(w r_k p_i) / sum_e(w r_k) = average of p_i weighted by
  template k's events → identical across templates only to the extent p_i does not
  correlate with r_k inside a bin (e.g. u ubar vs d dbar: different PDF ratios AND
  different EFT sensitivity). Measured: identical to 4 decimals → practically
  identical, not exactly.
- Each template gets its own Up/Down (filled with its own weights); nothing is
  "computed on SM and applied to EFT".
- All templates share one nuisance per systematic (100% correlated by construction).

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

## Muon momentum (Rochester) — `rochester_stat`, `rochester_syst`

Slide one-liner: "Rochester corrections (Muon POG): muon momentum scale (data + MC) and
resolution (MC) corrections. Uncertainties from 100 statistical replicas (std. dev.)
and 3 systematic variations (quadrature); applied as pT variations → bin migration."

### The problem
- Muon pT from track curvature in the 3.8 T field. Biases:
  - tracker misalignment: additive curvature bias, opposite sign for mu+/mu-,
    delta(pT)/pT ∝ pT → grows with momentum
  - B-field and material/energy-loss mismodelling: mostly multiplicative scale
  - all depend on charge, eta, phi
- MC resolution slightly better than data → simulated Z peak too narrow.

### The correction
- Muon POG standard, text file per year (RoccoR-type), derived from 2018 Z→mumu.
- Scale (data + MC): curvature corrected multiplicatively + additively in
  (charge, eta, phi) bins so that <1/pT> and the Z mass match expectations
  (data: `kScaleDT`).
- Resolution (MC only): `kSpreadMC` if gen-matched (stretches reco - gen), otherwise
  `kSmearMC` (random Gaussian, depends on nTrackerLayers; fixed seed 0).
- NOT a weight: replaces the muon pT (runner.py ~l.182-184) before selection, trigger
  matching, SFs and mll. Applied to data and all MC (backgrounds, sm, EFT). For EFT:
  the event moves to its (new) bin carrying all its template weights → consistent.
- Code: `spritz/modules/rochester.py` (Fabian's `varyRochester`/`correctRochester`
  always computes the 100 stat members; Giacomo's needs `do_stat=True`).

### The uncertainties
- Rochester file: sets (s) with members (m).
  - set 1: 100 members, correction re-derived with statistical fluctuations of the
    Z sample
  - sets 2, 3, 4: alternative choices in the derivation (systematics)
  - nominal in the code: set 5; variations = SF_i / SF_0 × SF_5 (shift relative to
    set 0, applied on top of set 5)
- Each variation re-runs selection + filling with the varied pT:

| nuisance | input | kind | Up / Down | analogue |
|---|---|---|---|---|
| rochester_stat | 100 stat members (set 1) | stdev | nom ± std over 100 | MC-replica PDFs |
| rochester_syst | sets 2, 3, 4 | square | nom ± sqrt(sum (X_i - X_nom)^2) | Hessian PDFs |

- Up/nom ≥ 1 by construction (nom + |sigma|).
- pT changes → events migrate across mll bin edges and in/out of the selection
  (pT thresholds, trigger matching). Largest where the spectrum is steep near an edge.
- Applied to all processes (datacard: 1.0 in every column), one nuisance each.

### Size on sm (Up/nom - 1)
| mll [GeV] | 50–60 | 60–80 | 80–100 | 100–120 | 120–140 | 140–700 | 700–1000 | 1000–1500 | 1500–3000 |
|---|---|---|---|---|---|---|---|---|---|
| stat | 0.19% | 0.04% | 0.00% | 0.06% | 0.05% | 0.02–0.05% | 0.13% | 0.29% | 0.16% |
| syst | 0.19% | 0.04% | 0.01% | 0.25% | 0.36% | 0.01–0.09% | 0.19% | 0.37% | 0.09% |

- 80–100 ≈ 0: a few-hundred-MeV shift keeps the Z peak inside the bin.
- 100–140 (syst up to 0.36%): steep falling edge above the peak; upper edge (100) is
  closer to mZ than the lower one (80): 8.8 vs 11.2 GeV → more migration across 100.
  Scale effect → syst > stat.
- 50–60: ~25 GeV muons at the pT thresholds → in/out of the selection.
- Above 700 GeV: misalignment bias grows with pT; corrections derived from ~30–60 GeV
  Z muons → least constrained.
- 1500–3000: wide bin, migrations stay inside.

### Open points
- Why is the nominal set 5 and not set 0? Check the Rochester file README.
- High pT: Rochester recommended up to ~200 GeV (to my knowledge); above that the
  Muon POG high-pT scale treatment (Generalized Endpoint, TuneP). mll reaches 3 TeV →
  same question as for the muon SFs.

## Pileup — `PU`

### What it is
- ~30–60 simultaneous pp interactions per bunch crossing. For muons the main effect is
  on isolation (pileup energy in the cone → fewer muons pass the relative-iso cut).
- MC is simulated with a fixed pileup scenario (UL18 mixing / premix library), which
  does not exactly match the 2018 data distribution.

### The correction
- Per event: w_PU(n) = P_data(n) / P_MC(n), with n = `Pileup.nTrueInt` (true mean
  number of interactions the MC event was simulated with).
- P_data: from the luminosity measurement (per-lumisection recorded lumi × sigma_inel,
  summed over the 2018 golden JSON). Central LUM POG product, standard for all 2018 UL
  analyses — nothing is extracted from our events (no circularity, unlike muon SFs).
- P_MC: the official UL18 mixing profile, built into the correction file. Valid for
  our private SMEFTsim samples because they were premixed with the standard UL18
  scenario (DIGIPremix step).
- Code: `spritz/modules/puweight.py` (identical in Fabian's and Giacomo's versions),
  correctionlib `cfg["puWeights"]` / `cfg["puWeightsKey"]`,
  `evaluate(nTrueInt, "nominal" | "up" | "down")`. Year-specific.

### The uncertainty
- Only unknown in P_data: sigma_inel = 69.2 mb ± 4.6% (CMS standard).
- Up/Down = weights recomputed with sigma_inel × 1.046 / × 0.954 (data profile shifts
  by ~±1.5 interactions on average). Provided directly by the file → `kind: weight`,
  no envelope.

### How it is applied
- Nominal: `puWeight` is in the event weight of every MC event (runner.py ~l.359):
  all 7 backgrounds (GGToLL, Single_Top, TT, WW, WZ, ZZ, DYtt), sm, every EFT template.
  Depends only on the event's nTrueInt → identical for all templates of an event →
  factors out of LIN/QUAD.
- Up/Down: `puWeight` swapped for `puWeight_PU_up/down`, all histograms refilled. One
  nuisance `PU` shared and fully correlated across all processes (config:
  `bkg_samples + DYll + eft_samples`; datacard: 1.0 in every column).

### Size on sm (Up/nom - 1)
- −0.1% to −0.5% in most bins, −1.1% at 50–60 GeV, +0.96% at 1500–3000 GeV.
- Sign: weights average to ~1 over the MC sample, so the total is unchanged; only the
  selection efficiency changes. Up = more pileup → lower iso efficiency → fewer events.
- Largest at 50–60 GeV: isolation is relative, so the same pileup energy matters more
  for ~25 GeV muons.
- Jumpy pattern and the sign flip in the last bin are statistical: Up weights are
  largest in the high-pileup tail where MC is sparse → noisy Up template.

### Open points
- Reviewer question: noisy PU template (smoothing? negligible impact?), check its
  post-fit pull/constraint.
- Not verified directly: which file/key is used (cfg.json not found in the config dir;
  the runner reads it from the job working directory). Low priority.
- Optional: compare `Pileup.nTrueInt` of a SMEFTsim file with an official UL18 sample
  to confirm the same pileup scenario.

## L1 prefiring — `prefireWeight`

### What it is
- L1 assigns each trigger primitive to a bunch crossing. Occasionally a primitive is
  assigned to BX-1 → L1 fires on the (empty) previous crossing → trigger rules forbid
  an L1 accept right after → the real event in BX0 is lost ("prefiring").
- Data only: the simulation does not model the timing problem → MC slightly
  overestimates the yield.
- Sources: ECAL endcaps (2.0 < |eta| < 3.0, radiation-induced timing drift; 2016–2017,
  fixed for 2018) and the muon system (wrong-BX muon primitives; all years). For our
  2018 muon analysis only the muon component matters.

### The correction
- w_prefire = prod_i (1 - p_i): probability that no object in the event prefires.
  p_i from prefiring maps in (eta, pT) (and phi for muons), measured centrally on
  "unprefirable" events (crossings where a BX-1 trigger was impossible by the trigger
  rules). Not measured on our events.
- Computed at NanoAOD production and stored in `L1PreFiringWeight` (`Nom`, `Up`, `Dn`).
  Our private SMEFTsim NanoAODs (standard NanoAODv9 config) have it like the official
  backgrounds.

### The uncertainty
- Up/Dn from the uncertainties of the prefiring maps, stored directly in NanoAOD
  → `kind: weight`, no envelope.

### How it is applied
- runner.py ~l.199-212: `prefireWeight = L1PreFiringWeight.Nom` in the event weight of
  every MC event (all 7 backgrounds, sm, every EFT template). Object-based → identical
  for all templates of an event → factors out of LIN/QUAD.
- Up/Down: `L1PreFiringWeight.Up/.Dn`, histograms refilled. One nuisance
  `prefireWeight`, shared and fully correlated across all processes (datacard: 1.0 in
  every column).
- If a sample had no `L1PreFiringWeight` branch the runner would set the weight to 1
  and register no variation; all our samples have it.

### Size on sm (Up/nom - 1)
- Flat −0.15% in every mll bin.
- Negative: Up = higher prefiring probability → lower weight → fewer events.
- Flat: only the muon component is active in 2018, roughly pT-independent, and the
  muon eta distribution barely changes with mll.

## Top pT reweighting — `tt_ptrw`

Slide one-liner: "TT reweight, applied only to TT bkg: corrects the excess of high-pT
tops in the NLO MC; uncertainty = full size of the correction."

### The problem
- TT sample: `TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8` (official UL18), POWHEG NLO QCD
  + Pythia8 (samples_fabian.json ~l.580). Not MiNNLO: MiNNLO (NNLO QCD) is only the DY
  sample used for the k-factor.
- Top pT spectrum is softer in data than in the NLO MC (too many high-pT tops,
  ~10% at a few hundred GeV). NNLO QCD (+ NLO EW) describes data → missing higher
  orders in the NLO MC.
- Relevance: harder tops → harder leptons → higher mll. TT is ~0.6% of the total
  overall but up to ~25% around 200–500 GeV.

### The correction
- TOP PAG recommendation, function of the generator-level (last-copy) top pT
  (`spritz/modules/tt_reweight.py`, functionally identical in Fabian's and Giacomo's
  versions):
  w(pT) = 0.103 exp(-0.0118 pT) - 0.000134 pT + 0.973
  (w = 1.08, 1.02, 0.99, 0.96, 0.92 at pT = 0, 50, 100, 200, 400 GeV) → raises low-pT,
  lowers high-pT tops, roughly normalisation-preserving.
- Event weight: sqrt(w(pT_top) × w(pT_antitop)); 1 if the event lacks a top or antitop.
- Central recommendation (top measurements + theory), not fitted to our events.

### The uncertainty
- Uncertainty = full size of the correction: Up = w + |w-1|, Down = w - |w-1|.
  For w < 1: Up = 1 (no reweighting), Down = 2w - 1 (correction applied twice).
  Built per event → `kind: weight`, no envelope.

### How it is applied
- Only TT: `top_pt_rwgt: True` only on the `TTTo2L2Nu` dataset (config_v9.py ~l.101).
  All other samples (incl. Single_Top, sm, EFT) get topPtWeight = 1 (runner.py ~l.177).
- Nominal: `topPtWeight` in the event weight (runner.py ~l.360).
- Datacard: `tt_ptrw` has 1.0 only in the TT column. Only process-specific
  experimental-type nuisance.

### Size (propcorr_v1 shapes)
| mll [GeV] | 180–220 | 220–270 | 270–350 | 350–500 | 500–700 |
|---|---|---|---|---|---|
| TT Up/nom | +3.3% | +3.8% | +4.4% | +5.3% | +6.4% |
| TT / total | 21% | 24% | 25% | 23% | 17% |
| effect on total | 0.7% | 0.9% | 1.1% | 1.2% | 1.1% |

- TT Up/nom rises from +2.1% (50–60) to +10.9% (1500–3000); TT/total peaks at ~25%
  (270–350) and falls below 5% above 1 TeV.
- ~1% on the total prediction in 200–700 GeV, comparable to PSWeight there.
- Up/nom > 1 everywhere: most TT events have w < 1, so Up (= no reweighting) gives
  more TT.

## MC statistics — `autoMCStats`

Slide / talk line: "Then we have the Monte Carlo statistical uncertainty, currently
applied only to the backgrounds: the SM and EFT templates are reweightings of the same
events, so their statistical fluctuations are fully correlated and can't be treated as
independent bins. A proper covariance-matrix treatment is in the next steps."
→ add "full covariance-matrix treatment of SM/EFT template MC stat" to next steps.

### The problem
- Templates are built from a finite number of weighted MC events → each bin's
  prediction has sigma_b = sqrt(sum_i w_i^2).
- Different from data statistics (Poisson term of the likelihood, the "Stat only"
  curve): MC stat is a systematic on the prediction itself.

### Combine treatment (Barlow-Beeston lite)
- No Up/Down histograms: Combine reads the bin errors in shapes.root.
- Card: `inc_mm_mll autoMCStats 10 0 1` = threshold 10, include-signal 0, hist-mode 1.
- Per bin: total included yield nu_b, error sigma_b, n_eff = (nu_b / sigma_b)^2.
  - n_eff > 10: one Gaussian nuisance per bin (`prop_bin<b>`) on the total MC yield
    ("lite": one per bin for all processes, minimised analytically → fast)
  - n_eff <= 10: one Poisson nuisance per process in that bin (`prop_bin<b>_<proc>`)

### Applied to: backgrounds only
- sm + 54 EFT templates excluded twice: include-signal = 0 (signal by process ID <= 0)
  and `noStat=True` in the config (bin variances zeroed in shapes.root;
  config_v9.py ~l.182-188).
- Why: sm, w1, wm1 are the same events reweighted → fully correlated fluctuations.
  Independent per-template treatment would add them in quadrature and break the
  LIN/QUAD algebra. Correct treatment = full covariance matrix between templates
  (not implemented yet; see meetings.md, 1 Oct).

### Consequences / open points
- The SM/EFT MC stat uncertainty is currently NOT in the fit (DY ~99% of the yield).
- 3D fit: SMEFTsim MC stat dominates in large |cos theta*| / |yll| bins (Fabian,
  1 Oct) → covariance treatment or more statistics needed.
- Effective SM statistics: samples generated at the base point (all ops = 1) and
  reweighted to SM → at high mll (EFT-dominated base) SM weights are small and spread,
  n_eff can be much lower than the generated 3.5M per mll bin.
- TODO: quantify for the 1D mll fit (sm MC stat vs expected data stat, n_eff per bin)
  from histos.root (variances still present there):
  `python3 -c "import uproot, numpy as np; h=uproot.open('<propcorr_v1>/histos.root')['inc_mm/mll/histo_sm']; v=h.values(); e=np.sqrt(h.variances()); print(100*e/v, 100/np.sqrt(v), (v/e)**2)"`

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
