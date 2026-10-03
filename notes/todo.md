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
- [ ] COMBINE slide: PDF text is wrong ("103 replicas divided up/down") -> NNPDF3.1 symmetric Hessian, 100 eigenvectors summed in quadrature, alphaS separate (101/102)
- [ ] mu_trig: 2-muon formula 1-(1-sf1)(1-sf2) uses SFs instead of efficiencies, uncertainty suppressed ~50x (same in Giacomo's version)
- [ ] high-pT muons: check muWP and SF validity above ~200 GeV (mll up to 3 TeV, SFs clamped to last pT bin)
- [ ] fit: check post-fit pulls/constraints of mu_reco/mu_idiso/mu_trig (Z peak can constrain SFs measured on the same Z events)
- [ ] alphaS: Up = mem 102, Down = mem 101 instead of envelope (same result now, cleaner)
- [ ] Rochester: why is the nominal set 5 instead of set 0?
- [ ] Run 2: lumi correlation scheme across years, per-year muon SF correlations
