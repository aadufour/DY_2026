## running MadGraph5
Installed locally via git clone.

To run MG5
cd /Users/albertodufour/MG5/mg5amcnlo
./bin/mg5_aMC

    
    generate p p > e+ e-
    output filename
    launch

#changing UFOs: SMEFTsim
Multiple sub-UFOs
create a link to a specific one
ln -s SMEFTsim/UFO_models/SMEFTsim_topU3l_MwScheme_UFO SMEFTsim_topU3l_MwScheme_UFO

    generate p p > e+ e- NP<=1     #NP<=1 means at most one EFT operator insertion per diagram: SM + linear interference terms and quadratic!!. 
    output DY_SMEFT
    launch

-----------------------
###(12/03/2026)
## CHOICE OF PARAMETERS
NP<=1 means we can have New Physics (includes EFT operators. NP=(=)0 is only SM)
NP^2=2 only selects EFT squared amplitude: ONLY diagrams with at least one NP operator (no SM). If NP^2<2(facciamo finta che sia la sintassi) we also consider interferences.
Suppose the SM contributes with 8 diagrams and EFTs. with 40, SM suqred amplitudes will have 8, NP^2=2 will have 40, the interference has 48 (must consider all diagrams) (interference scales like C, squared like C^2).


## RESTRICTION CARD TABLE
Only modify NON CP violating operators:
    ###################################
    ## INFORMATION FOR SMEFT
    ###################################
    Block SMEFT

DO NOT TOUCH (for now):
    ###################################
    ## INFORMATION FOR SMEFTCPV
    ###################################
    Block SMEFTcpv

Change one operator at a time from 0.00000... to 9.99999...e-01 (set to one: this is because 1.00000... gives problems with madgraph)
If the operator does not contribute to any graph it does not plot any.
The idea is to prepare restriction cards for each operator individually and make a table (just yes or no).
You can also run commands with madgraph as such: ./bin/mg5_aMC check_cjj11.txt
where check_cjj.txt is just
    import model SMEFTsim_topU3l_MwScheme_UFO-all_contributing_massless
    generate p p > mu+ mu- SMHLOOP=0 NP=1 (NP^2==2)
    output test_cjj

!! syntax:
    import model SMEFTsim_topU3l_MwScheme_UFO-cjj_massless
    the dash - and then the name of the restiction card






-----------------------------
13/03
REWEIGHTING
Idea: just generate events one time and assign a different weight to each event, depending on which oeprator is turned on (the new weight is basically |M_old|^2/|M_new|^2, essentially it tells you the pprbability of that event happening with the modified restrictions, e.g. turning off al operators [SM]: if some event could nevwer happen in standard model the weight will be zero). Saves a lot of time!
1. generate p p > mu+ mu- with ALL operators set to 0.99999 and output to myOUTPUTDIR
2. make Cards/reweight_card.dat and introduce a launch for all different reweights you will need to calculate
    this card will show something like (!! I(w/ Claude) made a script that produces it):
    """
    #first rewight: everything to zero [SM]
    launch --rwgt_name=SM
        set SMEFT 5 0 
        set SMEFT 9 0 
        set SMEFT 26 0 
        set SMEFT 28 0 
        set SMEFT 30 0 
        set SMEFT 103 0 
        set SMEFT 104 0
        set SMEFT 105 0 
        set SMEFT 107 0 
        set SMEFT 108 0 
        set SMEFT 109 0 
        set SMEFT 115 0 
        set SMEFT 117 0 
        set SMEFT 119 0
        set SMEFT 121 0 


        (...)
    """
and do this for each SINGLE operator (CRUCIAL: do not put all operators, we do not care about it, not useful)
3. cd path/to/myOUTPUTDIR
    ./bin/generate_events (make sure that the reweight card is called reweight_card.dat)

3. the output will be a single LHE file with all of the reweight info (generating sample only once!, then rescaling according to the cross sections, which is a much simple integral, compared to generating everything again)

One interesting thing might be to do a diff of the distributions of the weights.

To deal with TWO operators we have to calculate the weights for each turned on individually and both on at the same time, along with SM.
To deal with 3 or more its still only needed to reweight combinations of two 8this is because the cross section, and therefore the weights, are the matrix element squared, so we only have double products.

!!! It is crucial to have all operators set to 1 (0.99999....) in the restriction card, so that the output knows that it has to take into account every operator. these will be later turned on and off based on the reweight card.
TIP: on the reweight card you can write 1 and -1 without the damn 9s...




!!! I now work with two conda venvs

dy_analysis for pylhe and mlphep
MG5_2_9_18 for MG5


!!!!!!!!!!! watch out: am now working with a very specific version of madgraph and python: NEVER use the folder MG5_old OR MG5_env 



#------------------------------
16/03/2026
## TRIPLE DIFFERENTIAL
UNROLLING HIST 3D
3D (or even 2d) histograms cannot be fitted (technical limit): we can unroll them onto a 1D histogram so they can be fitted [what kind of fit??]
The next step wold be to plot the chi2 (with obs=sm, coupling C=0) (triple diff) vs C: it is a kind od sensitivity test

Ratio plot



run_card.dat
genero eventi per diversi intervalli mll
copio DY_all -> DY_all_min_max [GeV] cambiando run_card con taglio in mll
nell'analisi merge lhe (PESATI!!!)
scelgo [50, 120, 200, 400, 600, 800, 1000, 3000]

For tmrw:
- run /Users/albertodufour/code/DY2026/triple_diff_test/run_generate_events.sh to generate events in all DY_all_min_max
- manually export lhe files (?) in myLHE
- run /Users/albertodufour/code/DY2026/triple_diff_test/rw_triple_diff.py so it merges and plots the unrolled histo

10k events are ok for each bin (so standard)
there should be mìfuller tails and stuff (hopefully no empty bins anywhere)


Ptl1/2, etal1/2, phi1/2, ptz, dr(l1,l2), Deltaphi, DeltaEta

fare più stats su 800-1000 e 1000-3000: 100k eventi.

!!!!!!!!!!!!
Qusto andrà sicuramente in tesi!!





#
stampo anche plot triple diff con solo lin e uno con solo quad (senza log, si vede meglio la shape e le incertezze) (sempre unrolled)


## fixing stuff
import model SMEFTsim_topU3l_MwScheme_UFO-<operator>_massless (syntax UFO-[restrict]_operator_massless[.dat])
define p = p b b~          (needed to cfr with giacomos model: bs are massive so they are not in p, for some reason)
generate p p > l+ l- QCD=0