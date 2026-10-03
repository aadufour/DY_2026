# Meeting notes

## 2026-09-23 — Govoni (risposte Giacomo 28 sept, dopo il !)

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

## 2026-09-28 — Giacomo, Fabian

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

## 2026-10-01 — Triple diff meeting (Giacomo, Fabian, Raffaella)

- fabian question: is there a better way to deal with uncertainties besides the up/down envelope? For PDF is not necessary as they are all statistical replicas so we can sum the squares. For QCD (muf, mur) we will be asked to split mur and muf so ok (idk giacomo said this, did not really understand). If we do not have up/down but just 1 variation we can symmetrize it!

fabian presentation: /Users/albertodufour/Desktop/DrellYan/update_fabian_1oct.pdf

If we had hessian eigenvariations we should keep all the 100? we can then symmetrize it and have 100 shapes up/down. We can divide bin by bin but it becomes a huge datacard

- latest fabian push: fake templates. we have to add smt to the config, today slide 4 ()
- we can see that for now the MC stats is dominating for smeftsim, especially in large ctheta and yll bins. should be fixed with the new generation but maybe we will have to generate more for those bins.
- slide 9: what are these super constraints? to implement full covariance matrix (not yet implemented) we will have to set bin errors to zero (automcstats) [also on SM]; otherwise it will assume they are all independent and they will be added in quadrature and we will have giant errors (algebra doesnt close...)


- to check the wiggly lines: we can remove qcd and see what happens.
tuesday 6
