1. [x] Propapagtor correction! especially important when you are on shell!
    - [x] build gridpack with propagator correction: I just updated restrict all massless so 1 1.000000e+00 # linearPropCorrections. actually we need to put it back in the tar.gz
    - [x] validate the propagator correction: does it work with 10k events?
    - [x] make nanoaods
    - [x] analysis and fits
    - [x] comparisons

--
- [x] combine fits still says 138 fb^-1 -> to change (just visual!)

---
- [ ] WRITE ON THE NOTE?? For sure this part on propcorr ("never done by anyone before")

----

- [ ] introduce syst one  by one: we can work out why the lines wiggle


- [ ] check SM of this plot against by SMEFTsim SM (reweighting weight)


---
- [ ] ```spritz-eft-plot```:
    - [x] bottom panel for 3D (why is there data points fixed at 0 in the unblind region, screws up everything. This could be solved just by changing the range but I first want to check wether its just a plotting bug or actual selection bug);
    - [ ] rapll/costheta are NOT events/GeV (just Events)
    - [ ] finally move the legend for the bottom panel
    - [ ] adjust the title for tripl diff (too high up)
    - [ ] put a better legend for triple diff (maybe common, with the title)
    - [ ] triple diff: move the legend with the range for costhetastar and yll omn top right (now its on top of the peak)



---

2. [ ] cHQ1 and cHQ3 are the same???


-------


3. [x] Plots with SAME BINNING FOR EACH VARIABLE! 

---


4. [ ] Scan with ONLY LINEAR/QUADRATIC (maybe very close to c=0) -> check giacomos model on mattermost
---

5. [ ] 2D operator scan
 - [x] mll
 - [ ] other vars



---
roberto
- [ ] ceu vs clu (fattore 3?)
- [ ] add stat uncertainty band (another color, unrelated to the systematic)
- [ ] mll e 3D sono events/GeV!!
- [ ] make better plots for the ratios of propcorr/baseline

-----

- [ ] FABIAN UPDATED THE REPO??







































-------
tmp: 2D fits

operators appearing in a "degenerate/flat fit" pair from build_correlation_matrix.py - box/metadata.json range likely still wrong for these
- [x] cbWRe
- [x] cbBRe
- [ ] clj1: NO!
- [ ] cHj1: NO!
- [ ] ced: NO!
- [ ] cje: NO!
- [ ] clu: NO!
- [ ] cld: NO!
- [ ] cHQ1: NO!
- [ ] cHQ3: NO!
- [ ] cHu: NO!
- [ ] cHd: NO!
- [x] cHbq
- [ ] clj3: NO!
- [ ] cHe: NO!
- [x] cQl1
- [x] cQl3
- [x] cbe
- [x] cQe
- [x] cbl
- [ ] ceu: NO!

(didn't work for fits)

IS UNCORRELATED:
- [x] cbWRe_cbBRe: PDF AND PNG DO NOT MATCH!!
- [x] cbWRe_clj1
- [x] cbBRe_clj1
- [x] cHj1_clj1
- [ ] cHj1_ced
- [x] cHj1_cje
- [x] cHj1_clu
- [x] cHj1_cld
- [ ] cHQ1_clj1
- [ ] cHQ1_cld
- [ ] cHQ3_clj1
- [ ] cHQ3_cld
- [x] cHu_clj1
- [ ] cHu_cje
- [x] cHu_clu
- [ ] cHu_cld
- [x] cHd_cHbq
- [x] cHd_clj1
- [x] cHd_cje
- [x] cHd_clu
- [x] cHd_cld
- [x] cHbq_clj1
- [x] cHbq_clj3
- [x] cHbq_ced
- [ ] cHe_clj1
- [ ] cHe_cje
- [ ] cHe_clu
- [ ] cHe_cld
- [x] clj1_cQl1
- [x] clj1_cQl3
- [ ] clj1_cbe
- [x] clj1_cQe
- [ ] clj1_cbl
- [ ] clj3_cQl1
- [ ] clj3_cQl3
- [ ] clj3_cQe
- [ ] clj3_cbl
- [ ] cQl1_ceu
- [ ] cQl1_ced
- [x] cQl1_cje
- [x] cQl1_clu
- [ ] cQl1_cld
- [ ] cQl3_ceu
- [ ] cQl3_ced
- [x] cQl3_cje
- [x] cQl3_clu
- [ ] cQl3_cld
- [ ] ceu_cQe
- [ ] ced_cbe
- [ ] ced_cQe
- [ ] ced_cbl
- [ ] cbe_cje
- [x] cbe_clu
- [ ] cbe_cld
- [x] cje_cQe
- [ ] cje_cbl
- [ ] cQe_clu
- [ ] cQe_cld
- [ ] clu_cbl
- [ ] cld_cbl