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







-------------
Govoni 23 Sept
- perche queste variabili? non si potrebbe prima selezionare op -> scelgo variabili migliori per sens? [possibile sviluppo]

- REWEIGHTING: se w1 w wm1 non fossero reweighted non potrei fare algebra -> lin quad
- TopU3l: è la terza generazione quark singled out o solo top?
- MINNLO è NNTLO QCD e NLO QCD? da capire
- il k factor da minnlo non dovrebbe propagare incertezza QCD?
- cos'è tt_ptrw?

- variare binning: citare come possibilità ma non approfondire
- studio iniezione segnale BSM: so fittare? injection study
- bias study: cosa succede se si stima male il fondo (piu tecnico)

- fix summary plot per quando ho due minimi (e.g. yll e ctheta buggati)

- profilazione: lascio variare gli altri coefficienti durante un fit (o una parte)