1. [x] Propapagtor correction! especially important when you are on shell!
    - [x] build gridpack with propagator correction: I just updated restrict all massless so 1 1.000000e+00 # linearPropCorrections. actually we need to put it back in the tar.gz
    - [x] validate the propagator correction: does it work with 10k events?
    - [x] make nanoaods
    - [x] analysis and fits
    - [x] comparisons

--
- [ ] combine fits still says 138 fb^-1 -> to change (just visual!)

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