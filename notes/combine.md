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



da aggiungere pesi misti
shape a datacard con sm_lin_quad_mixed_op1_op2


provo a generare gli eventi sulle macchine llr: vorrei systematics (sbatti lhapdf)





## apptainer
spritz: framework di analisi
src/runners: codice principale
modules -> leptoni,...

flessibile: si configura co configs file

sample definiti in data
ogni anno c'è una production chain diversa
samples.json: dati e MCs (nostro anno 2018)
non avrò "fake" bkg (data driven), solo bkg MC

!! cfr cmsweb das per vedere i datacenter con i file...


5 stadi per girarlo

1. config.py
2. ricavo file e replicas: spritz-fileset
3. fare i "chunks" (si raggruppano i file in base al numero di eventi): spritz-chunks
4. fare i jobs: spritz-batch -dr (dr è dry run, non sottomette a condor)
! per girare in locale devo girare runner.py per ogni batch, lungo -> ./run_local.sh (da modificare i path!!)
    Apptainer> for i in 0 1 2 3 4 5 6 7; do ./run_local.sh $i; done
5. collezionare output: spritz-merge
6. postprocessing: spritz-postproc
7/8: plots e datacard: spritz-plots e spritz-datacard


grid proxy cerificate

## apptainer
apptainer shell -B /etc/grid-security/certificates:/etc/grid-security/certificates -B /cvmfs <PATH/A/spritz-env.sif>
apptainer shell -B /etc/grid-security/certificates:/etc/grid-security/certificates -B /cvmfs spritz-env.sif


/home/llr/cms/adufour/spritz/configs/vbfz-2018/config.py

voms-proxy-init --rfc --voms cms -valid 192:00


echo $X509_USER_PROXY


Apptainer> ls condor
job_0  job_1  job_2  job_3  job_4  job_5  job_6  job_7	runner.py  run.sh  submit.jdl
Apptainer> for i in 0 1 2 3 4 5 6 7; do ./run_local.sh $i; done





# generating events to work on systematics on llr (14/04/2026)

DY SMEFT LHE → Cache → Datacard workflow
─────────────────────────────────────────
Repo: /grid_mnt/data__data.polcms/cms/adufour/DY_2026
MG5:  /grid_mnt/data__data.polcms/cms/adufour/MG5/mg5amcnlo
venv: /grid_mnt/data__data.polcms/cms/adufour/dy_venv
cache output: MG5/.../CACHE/lhe_cache.pkl

Activate env:
  source /grid_mnt/data__data.polcms/cms/adufour/dy_venv/bin/activate
  export SETUPTOOLS_USE_DISTUTILS=stdlib
  export PATH=/grid_mnt/data__data.polcms/cms/adufour/dy_venv/bin:$PATH

Event generation (7 bins from DYSMEFTMll50_120 to DYSMEFTMll1000_3000):
  - reweight_card.dat must be in Cards/ before generate_events -f
  - reweight card: analysis/gridpack/gridpack_misc/cards/reweight_card.dat
  - 406 reweight points: 1 SM + 27×(+1) + 27×(-1) + 351 pairs
  - PDF: lhapdf / 303600 (NNPDF31_nnlo_as_0118)
  - numpy==1.26.4 required for f2py (numpy 2.x breaks it)

Analysis:
  1. python3 analysis/combine/build_cache.py     → lhe_cache.pkl
  2. python3 analysis/combine/build_datacard.py --op cHDD --lumi 59740
     → histograms.root + datacard.txt
  3. createWS.py / runScans_fixed.py / runPlots.py (see notes/combine.md)




# new setup
using combine through CMSSW
cd /grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/combine
source env_llr.sh



cd /grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/combine

# steps 1-2: analysis python
source env_llr.sh analysis
python3 build_cache.py
python3 build_datacard.py --op cHDD

# steps 3-7: combine
source env_llr.sh
python3 createJson.py --datacard datacard.txt
python3 createCombineJson.py --datacard datacard.txt
python3 createWS.py 1
python3 runScans.py 1 initial
python3 runScans.py 1 scan



#------------
# new workflow
# Steps 1-2: analysis
dy_analysis
build_cache.py
build_datacard.py --op cHDD

# Steps 3-7: combine
dy_combine
createJson.py --datacard datacard.txt
createCombineJson.py --datacard datacard.txt
createWS.py 1
runScans.py 1 initial
runScans.py 1 scan
