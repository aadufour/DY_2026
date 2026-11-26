## Questions for Giacomo (from 2026-11-25, after the SMP-26-014 preapproval)
- [ ] propcorr: quad = (w1+wm1)/2 - sm < 0 for cHj3 and cll1 in the Z-peak mll bin (q/sm ~ -4e-4, ~0.5-0.8% of lin).
      Identical in propcorr_v1, v1_minnlo, v3; absent in eft_bkg_fullsyst_v9 (no propagator correction).
      All full templates (sm, w1, wm1, w11) >= 0. Expected from the propagator/width expansion at O(c^2), or a gridpack issue?
      (check: analysis/combine_tools/check_positivity_morphing.py)
      Also at LHE level (post-fix propcorr cache): mll [80-100] quad/SM = -3.9e-4 (cHj3), -3.0e-4 (cll1), ~22 sigma MC stat;
      all 20 Z-peak triple_diff bins negative -> generator-level, not reco/morphing (2026-11-17).
- [ ] Negative background bins (DYtt ~ -0.1 ev in 3 triple_diff bins, TT ~ -0.03 in 2): clip to 0 / small positive and document, or leave as is?
- [x] Remaining bad 2D grid points: Giacomo 2026-11-25 "per la tesi direi che va bene così" -> leave as is.
      (Refit rounds 1+2, propcorr_v3.md "refit of bad points"; left: isolated stuck spikes, e.g. cHl1_ced, strict
      minimizers never converge there, only the loose last fallback returns a value.)
      Still open (minor): adopt setting G (strategy 2, tol 0.01, tight fallbacks first) as runScans.py default?
- [ ] Lin-only / quad-only fits (relevance of missing dim-8, validity of the linear approx, sm+lin vs full):
      formulas (Giacomo's whiteboard + N-operator generalisation) in notes/lin_quad_only.md. No combine model yet;
      Giacomo will look at it soon (2026-11-25) -> pick up his version if it comes, else implement.
- [ ] EFT-template MC noise (binning_study.md, N80 nuisance study): most of the N80 1D gain is an artefact of noisy
      EFT linear templates (reweighting from all-ops=1, n_eff 0.2-48% per slice), not covered by autoMCStats.
      Options: smooth/fit the EFT templates, min EFT n_eff per bin, or EFT template MC stat in the likelihood.

## Analysis, open
- [ ] EFT-template MC noise: check how much the v3 triple_diff (280 bins) moves (smoothed EFT templates).
- [ ] sm template weight outliers: LHE 200-400 slice has one event with SM/gen weight x782 (n_eff 1%).
      Check the reco sm template (n_eff per bin, spikes) in v3 and N80.
- [ ] Profiled fits on the full set: combine cQl1+cQl3 -> cQl13, cHQ1+cHQ3 -> cHQ13 (identical templates: use the
      cQl1/cHQ1 templates, drop cQl3/cHQ3 and their w11; check w11_cHQ1_cHQ3 = sm + 2*lin + 4*quad).
      set16 on triple_diff_w11export done 2026-11-23 (profiled_fits.md, EOS profiled_triple_diff/): ratios 1.1-2.1.
      3 exactly flat linear directions (cHl1, cll1, cHe, cHj3, cHl3, cHDD, cHWB): all 25 cannot be profiled together.
- [ ] FASTEXIT: after the w11 fix only at 2dNLL >> 100 (huge EFT effects, nuisances pulled far), see propcorr_v3.md.
      If needed, Giacomo's "usual method": run with NO nuisances, then add them back in blocks. Also costhetastar / rapll_abs?
- [ ] Summary plot: fix for operators with two minima (yll, costheta) (Govoni)
- [ ] MC stat uncertainty in the fits: already implemented, just check it is on (13 Oct)

## Done (recent)
- [x] w11 cross terms for 2D (propcorr_v3_2d, mll + triple_diff) + exportNuisances step (`*_w11export`), 2D scans,
      overlays, EOS (2026-11-19/20); bad-point refits (2026-11-21).
- [x] Prefit plots for mll + triple_diff: RECO3D_propcorr_v3/<var>/prefit/ (2026-11-19)
- [x] Same-binning N80 chain run + compared with v3 (2026-11-21, binning_study.md, EOS RECO3D_same_binning_N80/)
- [x] Profiled test2 redone on triple_diff_w11export (unchanged) + set16 (2026-11-23)
- [x] AN propcorr section (§7.3) + one-jet gridpack paragraph in §7.3.1 (2026-11-24)

# THESIS PLAN (2026-11-03)
Analysis work is over (talk before pre-approval done). Goal: final presentation + thesis (read by Govoni).
Only the points discussed with the profs, then writing. No precision uncertainties.
Own samples, 2018 only, no full Run 2.

Decided with Giacomo (2026-11-05):
- own samples, 2018 only: OK
- NO higher-order corrections (N3LO QCD / NLO EW) and NO fakes: would need understanding + explaining in the thesis
- NO correlated MC stat on the signal (same reason); keep noStat on the EFT templates, state it as a limitation
- theory uncertainties: PORT THE MiNNLO WAY (per-variation k = MiNNLO_V / LO_sm_V, as Fabian's renormalize_hists):
  SM part at NNLO, only EFT/SM at LO. Postproc-only change -> post_process_minnlo.py

Analysis: see "Analysis, open" at the top.

Plots (only what is visually wrong; spritz-eft-plot list below)
- [ ] axis labels: mll and 3D Events/GeV, rapll/costheta Events (Roberto)
- [ ] stat uncertainty band (another colour, separate from syst) (Roberto)

Writing only
- [x] propcorr section ("never done before"): AN §7.3
- [ ] binning variation: mention as possibility, don't pursue
- [ ] cHQ1/cHQ3 exact degeneracy (b-quark-only process): one line

Optional (discussed, decide)
- [ ] K-factor SMEFTsim vs SMEFTatNLO for one operator (Govoni: "interessante, mai fatto") - heaviest
- [ ] injection test in combine (--setParameters c=x) - easy

=====================================================================
PARKED: everything below is NOT for the thesis (kept for reference)
=====================================================================

- [x] combine fits still says 138 fb^-1 -> to change (just visual!)


----

- [x] (primo) introduce syst one  by one: we can work out why the lines wiggle. prima tolgo teoriche (pdf, qcd scale), poi exp. controllo che la parametrizzazione sia giusta: /grid_mnt/data__data.polcms/cms/adufour/CMSSW_spritz/CMSSW_14_1_0_pre4/src/tools/checks
probabile da aggiornare per morphing.
checkPositibity assume shapes che si chiamano "sm", "quad_" mentre tu hai "sm", "w1_" "wm1" etc...

...non mi torna come mai debba essere necessario modificare queste cose... comunque ho capito che il problema è QCD
- [ ] tolgo QCD scale dai fit.
- [ ] alphaS entra?

----

- [ ] check SM of this plot against by SMEFTsim SM (reweighting weight)(l'ha fatto fiacomo nella '/Users/albertodufour/Desktop/3DY_Fabian_Samples_LO - Presentazioni Google.pdf')


---
- [ ] ```spritz-eft-plot```:
    - [x] bottom panel for 3D (why is there data points fixed at 0 in the unblind region, screws up everything. This could be solved just by changing the range but I first want to check wether its just a plotting bug or actual selection bug);
    - [ ] rapll/costheta are NOT events/GeV (just Events)
    - [ ] finally move the legend for the bottom panel
    - [ ] adjust the title for tripl diff (too high up)
    - [ ] put a better legend for triple diff (maybe common, with the title)
    - [ ] triple diff: move the legend with the range for costhetastar and yll omn top right (now its on top of the peak)



---
roberto
- [ ] ceu vs clu (fattore 3?)
- [ ] make better plots for the ratios of propcorr/baseline
(stat band and Events/GeV labels: in "Plots" above)

-----

- [ ] FABIAN UPDATED THE REPO?? (anche no)


-----
Meeting notes moved to notes/meetings.md (2026-10-03)

-----
Uncertainties, open points for SMP-V (2026-10-03, details in notes/uncertainties.md)
- [ ] lumi: confirm 0.84% for 2018 (spritz data/common/lumi.json rel_unc) with Giacomo / current LUM POG table (older recommendation 2.5%)
- [ ] lumi: read rel_unc from lumi.json instead of hardcoding 1.0084 (config_v9.py, config_propcorr_v1.py)
- [ ] lumi label: plots say 59.74 fb^-1, normalisation uses 59.56 fb^-1 -> make consistent on slides/plots
- [x] COMBINE slide: PDF text is wrong ("103 replicas divided up/down") -> NNPDF3.1 symmetric Hessian, 100 eigenvectors summed in quadrature, alphaS separate (101/102)
- [ ] mu_trig: 2-muon formula 1-(1-sf1)(1-sf2) uses SFs instead of efficiencies, uncertainty suppressed ~50x (same in Giacomo's version)
- [ ] high-pT muons: check muWP and SF validity above ~200 GeV (mll up to 3 TeV, SFs clamped to last pT bin)
- [ ] fit: check post-fit pulls/constraints of mu_reco/mu_idiso/mu_trig (Z peak can constrain SFs measured on the same Z events)
- [ ] alphaS: Up = mem 102, Down = mem 101 instead of envelope (same result now, cleaner)
- [ ] Rochester: why is the nominal set 5 instead of set 0?
- [ ] Run 2: lumi correlation scheme across years, per-year muon SF correlations

-----
Final version of the analysis (after SMP-V, 2026-10-04, details in notes/uncertainties.md)
- [ ] fakes: add data-driven nonprompt-muon estimate (Fabian's fake templates, same-sign region inc_mm_ss) with fakes_param + fakes_model nuisances, as in spritz_giacomo/configs/dy-eft-2018
- [ ] b-veto against ttbar (TT ~25% of total at 200-500 GeV) + btagSF_sf / btagSF_eff / puidSF nuisances -> discuss with Giacomo
- [ ] higher-order DY: N3LO QCD + NLO EW corrections with their uncertainties (as in dy-eft-2018) instead of only the NNLO QCD MiNNLO k-factor
- [ ] background cross-section lnN uncertainties (WZ, ZZ, GGToLL have none; TT/WW only scale/PDF)
- [ ] alphaS on TT (members 101/102 available, ±3%; WW has none: replica set 320900)
- [ ] WW PDF: stored set 320900 is MC replicas → use std. dev. (kind stdev, members 1-100), not quadrature (~10x overestimate now)
- [ ] ST s-channel: LHE scale/PDF weights off by factor 2 → rescale x2 or exclude from QCDScale (ask Giacomo)
- [ ] Single_Top PDF/alphaS: tW weights broken → flat lnN or borrow ttbar relative variation
- [ ] background xsec: WZ/ZZ normalised at LO (27.59 / 12.17 pb vs NLO ~47 / ~16.5) → fix or add lnN; verify WW (11.09) and tW (21.7) against XSDB
- [ ] (maybe) switch VV to the exclusive NLO samples as in spritz_giacomo/configs/dy-2018 (Fabian already updated them): WWTo2L2Nu, WZTo3LNu, WZTo2Q2L, ZZTo4L, ZZTo2L2Nu, ZZTo2Q2L (POWHEG / aMC@NLO) → NLO shape + norm, QCD scale/PDF/alphaS available (alphaS not for WW, ZZTo2L2Nu). Would also solve the WZ/ZZ-at-LO item. Giacomo (2026-10-05): current inclusive Pythia setup is consistent, no contradiction
- [ ] MC stat of sm/EFT templates: full covariance-matrix treatment; first quantify n_eff per bin for the 1D mll fit (histos.root)
- [ ] PSWeight: split ISR and FSR into separate nuisances?
