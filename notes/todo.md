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
Govoni 23 Sept +. risposte giacomo 28 sept (dopo il !)
- perche queste variabili? non si potrebbe prima selezionare op -> scelgo variabili migliori per sens? [possibile sviluppo]
! viene da paper teorico. è così semplice che xsec full può essere completamente descritta da queste tre. c'è ridondanza.

- REWEIGHTING: se w1 w wm1 non fossero reweighted non potrei fare algebra -> lin quad
- TopU3l: è la terza generazione quark singled out o solo top?
- MINNLO è NNLO QCD e NLO EWK? da capire
! NLO EWK is not embedded in the matrix element but comes from pythia (ISR, FSR): missing virtual contributions
- il k factor da minnlo non dovrebbe propagare incertezza QCD?
! either take the minnlo or LO, free to choose. MINNLO better (smaller). take the minnlo and add on top a little uncertainty (lognormal). It would be nice to add the generation generate p p > l l j (adding real emission part of the diadgrams) -> prediction gets closer to LO. would need to generate new gridpacks.
- cos'è tt_ptrw?

- variare binning: citare come possibilità ma non approfondire
! si possono fare i due estremi: 1. stesso humero assoluto di bin tra single var e 3D (mll può essere più sensibile ma quando guardo in 2D PUò ESSERE MEGLIO, meno correlazione. quando faccio i fit profiled (e poi globali. quella è lòa fogura di merito principale). si rivsolvono e.g. flat direction) o 2. la variabili


- studio iniezione segnale BSM: so fittare? injection study
! abbastanza facile (comando di combine mettendo WC a zero, --setparameter a x invece di 0)
- bias study: cosa succede se si stima male il fondo (piu tecnico)



- fix summary plot per quando ho due minimi (e.g. yll e ctheta buggati)

- profilazione: lascio variare gli altri coefficienti durante un fit (o una parte).
! si può fare, tutto già implementato. Posso scegliere sottogruppo (servono i fit 2D).
devo escludere operatori flat finche il set di operatori è closed e faccio un profiled.

PRIMA COSA PROFILED. Imparo le direzioni che non posso costringere con questo approccio, qual è il set massimo di operatori che posso metter insieme (PCA?? unico modo per non buttare operatori). Per esempio sarebbe bello dare un'idea con 2018 delle flat 

Devo lavorare su incertezza. I k factor sono uguali tra Sm e EFT? Posso prenddere un operatore che c'è sia in SMEFTsim e SMEFTatNLO (sempre MG5 ma NLO). voglio controllare che tra i due il k factor rimanga uguale lin(op) (????????????->mattermost 15:18). se sono uguali posso applicare anche a EFT. ho il confronto solo a NLO pero lo estrapolo anche A NNLO! giustifico il fatto di usare k factor basato su SM per ordini successivi. Mai stato fatto, interessante.

- Fare full run2 è un esercizio di stile. Non ci guadagno granché



------------
28 sept giacomo fabian
- covariance matrix?
- last slide of fabian (bottom formula)
- we will maybe add a lognormal to take into account mismodeling uncertainty. or smt else?


Fabian timeline:
- defending beginning of feb.
- finish thesis end of nov.


- missing all reviews (they all require impact plots)
    ->ptv (o btv?), jme and muon
- then preapproval
    -> would be nice to have nlo EW + NNLO QCD SMEFTsim, smeftatnlo (not corr matrix). 1D results and maybe 1 2D examples
    -> maybe SMP general oct 16 (or 1 nov?)

- maybe SMPV on oct6? check on thursday