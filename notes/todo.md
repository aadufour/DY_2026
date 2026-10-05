- [x] combine fits still says 138 fb^-1 -> to change (just visual!)

---
- [ ] WRITE ON THE NOTE?? For sure this part on propcorr ("never done by anyone before")

----

- [x] (primo) introduce syst one  by one: we can work out why the lines wiggle. prima tolgo teoriche (pdf, qcd scale), poi exp. controllo che la parametrizzazione sia giusta: /grid_mnt/data__data.polcms/cms/adufour/CMSSW_spritz/CMSSW_14_1_0_pre4/src/tools/checks
probabile da aggiornare per morphing.
checkPositibity assume shapes che si chiamano "sm", "quad_" mentre tu hai "sm", "w1_" "wm1" etc...

...non mi torna come mai debba essere necessario modificare queste cose... comunque ho capito che il problema è QCD
- [ ] tolgo QCD scale dai fit.
- [ ] alphaS entra?

----------
- [ ]

------------




- [ ] Introduce mixed terms fo 2D fits
    - cannot MERGE them!
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

2. [ ] cHQ1 and cHQ3 are the same!


-------


3. [x] Plots with SAME BINNING FOR EACH VARIABLE! 

---


4. [ ] Scan with ONLY LINEAR/QUADRATIC (maybe very close to c=0) -> check giacomos model on mattermost

---
roberto
- [ ] ceu vs clu (fattore 3?)
- [ ] add stat uncertainty band (another color, unrelated to the systematic)
- [ ] mll e 3D sono events/GeV!!
- [ ] make better plots for the ratios of propcorr/baseline

-----

- [ ] FABIAN UPDATED THE REPO?? (anche no)

-----
- [ ] need to check for validity of the linear approx? e.g. checking sm+lin vs full. chiedo a giacomo

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
