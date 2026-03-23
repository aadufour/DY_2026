## COMBINE
20/03/2026
cfr presentazione Giacomo

Likelihood = Prod(canali, spazio fasi) * Prod(bin) {Likelihood(mu, theta | data) * prior(theta)}

spazi fasi devono esser blocchi statisticamente indipendenti (e.g. mll<50, mll>50 o anno 2016,2017,... o canali che non si parlano tipo jet=0 e jet=1,...)



theta nuisance params: teoriche (QCD scale [ossia scale di rinormalizzaione e fattorizzazione], valore alpha strong, PDF,...) o sperimentali(trigger, muoni, ricostruzione, jet, JES e jet energy resolution -> incertezze sulla met, pileup reweighting, jet pilup ID)
Inizialmente io introdurrò teoriche (non ho ancora MC)

datacard hardcoded: devo definire un certo numero di variabili per combine (quanti nuisance params, quanti "bin" ossia phase space,...)

(....)

In fondo rate parameter: nuisance ma con flat prior (lognormal e.g. ce l'ha gauss): stiam data driven della xs di un fondo

In fondo autoMCstat: incertezza statsitica del MC (==barlow biston).


obiettivo: fare data card per DY EFT.
Il problema dei fit EFT è che le interferenze possono essere negativi (LIN). Combine crasha se trovo valori neg (PDF non possono essere neg)
trick riscrivendo in termini di roba ppos def (pag.11 slides, bottom right) -> scrivo datacard cosi


devo linkargli i root file con le distribuzioni (inizio con mll, poi 3d unrolled)


codice che crea datacard data osservabile e operatore.
datacard diverse per operatore
devo costruire i ROOT file con le distribuzioni (i.e. con gli istogrammi delle mie variabili tirate fuori da madgraph)








(cd /Users/albertodufour/combine) in realtà non serve??
text2workspace.py /Users/albertodufour/code/DY2026/datacard_test/datacard.txt -P HiggsAnalysis.CombinedLimit.AnomalousCouplingEFTNegative:analiticAnomalousCouplingEFTNegative -o model_ced.root --PO addToCompleteOperators=ced --PO eftOperators=ced --X-pack-asympows --optimize-simpdf-constraints=cms --X-optimizeMHDependency=fixed --X-allow-no-signal --X-allow-no-background

text2workspace.py trasforms my datacard.txt (text) into a roofit workspace (not human readable:likelihoods,...)
inputs the datacard as argv[1]
-P is the physics model to be used
-o is output with the physics model


1. datacard
2. workspace : p (argv[1] è il humero di operatori, se 2 fa tutte le combo di 2 op.)
3. fit:
    - initial fit
    - likelihood scan

    runScans.py <n_operatori> <action (= initial or scan)>
    initial fa il best fit (minimizza chi2 una volta e trova minimo globale di tutti i paramteri definiti nella datacard: PoI (params of interest) come k_ced e nuisances tipo lumi)
    per il likelihood scan fisso a il parameter of interest



!!TODO: DEBUG
runScans.py 1 scan --doSplitPoints=10
DEBUGGED: qualcosa a che fare con multiprocessing su macos, Claude fixed it
-> runScan_fixed.py


poi da capire runPlots.py

TODO:
    likelihood scan per ogni operatore (singolarmente)
        -> faccio una datacard con
        process sm quad_ced am_lin_quad_ced quad_cHDD lin_quad_cHDD ....
        poi seleziono con createWS.py --PO eftoperator....


    capire a quali siamo più sensibili (cfr. tesi Bulla p.67)

di fatto questo è un sensitivity scan, riproduco quell'articolo che mi ha mandato Giacomo...












(aggiunti a .zshrc: sono eseguibili, dentro hanno una ref a che pyhton usa [sarà un problema])

1. createJSon.py --datacard <path/to/datacard.txt>
    -> metadata.json
2. createCombineJSon.py --datacard <path/to/datacards.txt>
3. createWS.py [1,2,3] (num operatori): wrapper for text2workspace.py



#23/03/2026
debuggato tutto. Recap:

1. build_datacard.py
    -> fa histograms.root e costruisce datacard.txt
2. createJSon.py --datacard <path/to/datacard.txt>
    -> metadata.json !attenzione: da qui dentro si possono cambiare i limiti per il likelihood scan
3. createCombineJSon.py --datacard <path/to/datacards.txt>
    -> jsonComb.json
4. createWS.py <n_op>
5. runScans_fixed.py <n_op> initial
6. runScans_fixed.py <n_op> scan --doSplitPoints=10
7. runPlots.py <n_op>