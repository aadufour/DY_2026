# Closing Strategy — DY SMEFT Run2

*Last updated: 2026-09-29. Based on meetings with Govoni (23 Sept), Giacomo (28 Sept), Fabian (28 Sept).*

---

## Timeline

| Milestone | Target |
|-----------|--------|
| Fabian thesis done | end of Nov 2026 |
| Fabian defense | beginning of Feb 2027 |
| SMP general meeting | ~Oct 16 or Nov 1 (TBC) |
| SMPV | ~Oct 6 (check Thursday) |
| Pre-approval | after all reviews pass |

---

## 0. 1-jet study (already in progress)

- [x] Generate all 28 gridpacks (4 sets × 7 bins) with MLM matching (`ickkw=1`, `xqcut=15`)
- [ ] Integrate 1j samples into spritz pipeline (spritz-merge, spritz-postproc, spritz-postproc-eft, spritz-cards-eft)
- [ ] Compare propcorr / nonpropcorr / 0j+1j predictions for mll
- [ ] Motivation: real emission from `p p > l l j` brings LO prediction closer to NLO (as suggested by Giacomo in 23 Sept meeting)

---

## 1. Uncertainties — highest priority

### 1a. K-factor for EFT (key cross-check)
- [ ] Check whether SM K-factor (MINNLO/LO ratio) is the same for EFT operators in SMEFTsim vs SMEFTatNLO (NLO MG5).
  - Pick one operator present in both SMEFTsim and SMEFTatNLO.
  - Compare LO and NLO K-factors for the linear EFT term.
  - If equal → justify applying SM K-factor also to EFT (never done, scientifically interesting).
  - This would validate using NNLO QCD K-factor from MINNLO for EFT terms.
  - Giacomo: mattermost thread 15:18 on 28 Sept.

### 1b. Lognormal mismodeling uncertainty
- [ ] Govoni/Giacomo discussed adding a lognormal to account for mismodeling.
  - Decide whether to add it and for which templates.
  - Implement in spritz-cards-eft.

### 1c. QCD scale uncertainties in fits
- [ ] Remove QCD scale systematics from fits (they are causing wiggles).
  - `alphaS` — does it enter? Check.
- [ ] Theoretical systematics (PDF, QCD scale): introduce one by one, isolate the source of wiggles.
  - Check `checkPositivity` script: `/grid_mnt/.../CMSSW_spritz/.../tools/checks` — needs update for morphing (currently assumes `sm`, `quad_` naming, not `sm`, `w1_`, `wm1`).

---

## 2. Fit improvements

### 2a. Profiled fits — first priority among fits
- [ ] Profiled fits: let other WCs float freely during 1D fits.
  - Already implemented in combine, just needs to be run.
  - Identify the maximum set of operators that can be simultaneously constrained (no flat directions in the subset).
  - Use PCA to avoid discarding operators.
  - Exclude flat operators for now; include once the set is closed.
  - Give an example with 2018 data.

### 2b. Mixed terms for 2D fits
- [ ] Introduce mixed (interference) terms for 2D fits.
  - **Cannot merge them** with current approach — needs dedicated handling.

### 2c. Linear/quadratic-only scans
- [ ] Scan with ONLY linear terms (near c=0) and ONLY quadratic terms.
  - Check Giacomo's model on Mattermost.
  - Useful for understanding where the linear approximation breaks down.

### 2d. Linear approximation validity check
- [ ] Check validity of linear approx: compare `sm+lin` vs `full` prediction.
  - Ask Giacomo for the reference / expected scale.

---

## 3. Plots and cosmetics

### 3a. spritz-eft-plot fixes
- [ ] `rapll` and `costheta` axes: label as "Events" not "Events/GeV"
- [ ] `mll` and 3D: label as "Events/GeV" ✓ (already correct)
- [ ] Bottom panel: fix data points stuck at 0 in unblind region (plotting bug or selection bug?)
- [ ] Bottom panel: move legend
- [ ] Triple diff: fix title position (too high)
- [ ] Triple diff: better legend for range labels (costheta*, yll — move to top right, away from peak)
- [ ] Fix summary plot for operators with two minima (yll and costheta currently bugged)

### 3b. propcorr/baseline ratio plots
- [ ] Make better ratio plots: propcorr / baseline for each bin.

### 3c. Stat uncertainty band
- [ ] Add a separate stat uncertainty band (different color from systematics).

### 3d. Combine luminosity label
- [ ] Update combine output to show correct lumi (currently shows 138 fb⁻¹ — visual only, not a physics bug).

### 3e. cHQ1 / cHQ3 degeneracy note
- [ ] Document that cHQ1 and cHQ3 are exactly degenerate (b-quark-only process → physics, not a bug).
  - Already confirmed, just needs a line in the AN.

---

## 4. Reviews (blocking pre-approval)

All reviews require **impact plots**. These need to be produced before the review submissions.

- [ ] PTV (or BTV?) review
- [ ] JME review
- [ ] Muon review

---

## 5. Pre-approval target

For pre-approval, Giacomo/Fabian want:
- [ ] NLO EW + NNLO QCD SMEFTsim
- [ ] SMEFTatNLO (no correlation matrix)
- [ ] 1D results for all operators
- [ ] At least 1–2 2D examples (profiled)

---

## 6. Analysis note (AN)

- [ ] Write section on propcorr ("never done before by anyone in CMS DY") — Govoni specifically flagged this.
- [ ] Include the linear approximation validity check.
- [ ] Include the K-factor EFT validation study (if completed).

---

## 7. "Nice to have" (explicitly de-prioritized)

- Full Run2 combination: "exercise in stile, not much to gain" (Giacomo).
- Injection study / bias study: feasible but not blocking.
- Binning variation: mention as a possibility in the AN, do not pursue.
- FABIAN REPO update: also no (anche no).
- ceu vs clu factor of 3 check (Roberto): low priority.
