# Propcorr section for AN-2024/110: open questions and doubts

Drafts: `propcorr_section_local.tex` / `results_propcorr_local.tex` (with `% src:`
comments, refs in words; compile with `main.tex` + `propcorr.bib`) and the `_AN.tex`
copies for Overleaf (no comments, real `\ref`/`\cite`). Edit the local ones and
redo the `_AN` copies from them. AN reference version: v5,
13 Oct 2026.

## To ask Giacomo

Blocking for the section:

1. ~~**Edit access** to the Overleaf project (currently viewer only).~~ -> answered 4 Nov, see Settled.
2. ~~Placement~~ -> answered 4 Nov, see Settled.
3. ~~Nominal templates~~ -> answered 5 Nov: 9.3 is with propcorr, see Settled.
4. **Missing SMEFTsim LO documentation.** For SMEFTatNLO the AN has m_ll-binned gridpacks,
   NanoAOD production and reco validation (§7.4-7.6). For SMEFTsim LO there is only the LHE
   validation (§7.2, with two test gridpacks), even though SMEFTsim LO results appear in §9
   (Figs. 110/111). Should I write those parts too? Material exists: 7 m_ll bins, reweight
   card with 406 points, CRAB NanoAOD chain (500k/bin, patched `GenWeightsTableProducer`),
   spritz v9 reco setup.

Can wait:

- ~~Giacomo's large production~~ -> answered 5 Nov, see Settled. Still open, low priority:
  does it use our propcorr gridpacks (one sentence in §7.3 samples paragraph)?
5. ~~Reco-level propcorr results~~ -> yes, new §9.4 (4 Nov). Still to ask: which
   variable(s) in 9.4 (triple_diff only, or all four as makeSummaryPropcorr produces), and
   are the 9.3 figures from our `config_v9` setup? (`config_propcorr_v1` = v9 with only the
   EFT signal swapped, so 9.4 compares like with like only if 9.3 is v9.)
6. **SMEFTatNLO and propcorr.** If the main result is NLO-based, the NLO templates
   (probably, not verified) don't include propagator corrections. Comparison, systematic,
   something else?
7. **First (invalid) production** with `linearPropCorrections = 0`: mention it, or only
   describe the correct production? Depends on whether those results were ever shown.
   Paragraph is commented out in the draft.
8. **Plots.** Giacomo asked for "physics at linear and quadratic level", so the draft now
   has all 15 operators, LIN and QUAD (two figures, like Figs. 57/58), plus the cQe control.
   Cross-section table replaced by one sentence (5 Nov); the numbers, for reference:
   50-120 1188.67 -> 1088.54 pb (0.916); 120-200 13.187 -> 13.138 (0.996);
   200-400 4.287 -> 4.291 (1.001); 400-600 2.076 -> 2.076 (1.000);
   600-800 1.565 -> 1.569 (1.003); 800-1000 1.240 -> 1.239 (0.999);
   1000-3000 4.197 -> 4.204 (1.002). Propcorr from LHE <init> (XERRUP < 0.3%), baseline
   from samples.json (its integration error not recorded).
9. ~~Credit for the patched model~~ -> dropped 5 Nov (not needed in the AN).

## Physics / technical doubts

10. ~~Truncation of the expansion~~ -> settled by the LHE plots (5 Nov): the quadratic
    templates change with propcorr and become negative at the pole for cll1, cHl3, cHj3,
    so QUAD contains the dim-6-vertex x Z' interference (|Z'|^2 removed by NPprop^2<=2).
    Now stated in the draft, with the chain of sources in the src comment. Still worth
    asking Fabian whether this truncation is intended. A negative QUAD also means
    SM + c LIN + c^2 QUAD can go negative at large |c| at the pole: possibly related to
    the negative-yield fit failures seen in wide scans (not checked).
11. **Validity of the linearisation.** (The guide's "accurate up to O(10%)" statement is
    for the W, Fig. 2 right of the guide; removed from the draft on 5 Nov as not about the Z.) Do we need an SM+lin vs full check for the operators
    with the largest delta Gamma_Z, or is citing the SMEFTsim linearised method enough?
    (Already an open item in `notes/todo.md` / `notes/closing_strategy.md`.)
12. **b-quark mass.** The AN's §7 uses the `_b_massless` model; our production uses
    `restrict_all_massless.dat`. Its `MB` value on LLR is NOT verified: the repo copy
    `analysis/restriction_card/restrict_all_massless.dat` is the stale card (it still has
    `linearPropCorrections = 0`), and `MB = 4.18` was read from Fabian's card, not ours. Check on
    LLR: `grep -A3 "Block MASS" <model>/restrict_all_massless.dat` and the `5 ` line.
    Hint from the LHE plots (5 Nov): the baseline (no propcorr) linear templates of
    cbWRe/cbBRe are small but non-zero (about -0.09 / -0.05 pb/bin at the pole), which
    suggests m_b != 0 in our production, unlike the AN's `_b_massless` model. The dipoles O_bW, O_bB enter
    delta Gamma_Z only through a term proportional to m_b. Which treatment is the reference?
    Also §7.2 explains the vanishing O_bW/O_bB effect with flavour-changing currents, while
    our explanation is the m_b chirality flip: the wording should be made consistent.
13. **O_HQ^(1) / O_HQ^(3) degeneracy.** Identical templates event by event (only the
    b quark is accessible, so only C_HQ^(1) + C_HQ^(3) enters, in both the vertex and
    delta Gamma_Z). Exact physics, not a bug. How to treat it in the fit: keep both, merge
    into one parameter, or fix one?
14. **cHj3 outlier.** In the LHE comparison, cHj3 shows a much larger propcorr/baseline
    difference than the other operators (max event weight ~24x the baseline one). Not
    understood yet.

## Self-checks before submitting

15. Real AN labels for Table 20, Eq. 58, §5.1.2, §7.2, §7.5 (marked `TODO label`), and the
    citation key of the SMEFTsim 3.0 guide.
16. Baseline cross sections in the xsec table come from the earlier `samples.json`: confirm
    they are from the same restriction card and m_ll bins (same proc cards apart from the
    `NP<=1` line).
17. **Operator count: O_He missing from Table 20.** AN §7 says "only 26" operators affect
    neutral DY at LO, and Table 20 lists 26; our configs use 27. The one missing from
    Table 20 is O_He (`cHe`, Z coupling to right-handed leptons): it is also not the
    O_HQ^(1)/O_HQ^(3) degeneracy, both of which are listed. Almost certainly an omission in
    the AN: O_He modifies the Z->ll vertex at tree level, the AN's own Table 21 maps
    `cHe -> cpe = cpmu = cpta`, and the NLO text counts 27. It is also one of the 15
    propcorr operators. Point out to Giacomo (add O_He to Table 20, 26 -> 27).
    Also ask: the AN (l. 909-910) says the gridpacks include "all operators listed in
    Table 20 and Table 21" -- did their SMEFTsim LO gridpack actually include O_He?
    (Verified visually on AN v5 PDF pages 84-85.)
18. ~~Selection~~ (selection sentence cut from the draft, 5 Nov; the plots still use the cuts below). The validation plots use pT > 25/10 GeV, |eta| < 2.4; §7.2 uses
    pT > 10 GeV, |eta| < 3.0. Align, or state the difference.
19. ~~Check the draft's statements against the plots~~ -> done 5 Nov; descriptions in §7.3
    written from the regenerated plots (values in the src comments).
    Tooling note: `run_compare_all_ops.sh` defaults to the high-stats baseline cache, which
    with `--no-normalize` comes out ~10x too high (compare_ops.py divides by a fixed
    N_GEN = 100k, the 50-120 bin has a 1M top-up). Use the default 100k cache, or fix the
    normalisation, before using that runner again.
20. ~~Impact subsection~~ -> moved to §9.4 (`results_propcorr_local.tex`), figures still to produce.
21. **"Only the Z propagator is corrected"** in neutral DY. Sourced in principle (SMEFTsim has
    dummies only for Z, W, h, t, §6.1; no photon dummy), but worth confirming on the actual
    diagrams. On LLR, for one propcorr process directory:
    `grep -l -E "\bz1\b" SubProcesses/P*/matrix*.ps` and the same for `w1`, `h1`, `t1`
    (expected: only z1 matches).
22. **Production restriction card** in the draft listing (`linearPropCorrections = 1`) is from
    the LLR output pasted during the bug fix. Update the repo copy so it matches, or the next
    reader of the repo will find the wrong value.
23. ~~Stale numbers in AN §7.2 text~~ -> wrong claim, corrected 5 Nov: "Equation 55" is a
    `\ref{eq:eft_parametrization_1op}`, the one-operator parametrization N(k) = SM + k Lin +
    k^2 Quad in terms of w(0), w(+-1). Our sections now cite that one (no new label needed);
    (58) only restates the algebra. "Figure 60" for Figs. 57/58 not checked in the source.

## Settled (already in the draft)

- **9.3 and 9.4 (Giacomo, 5 Nov):** the 9.3 expected limits use the propcorr samples (his
  production) with all uncertainties. 9.4 only has to show the effect of the propagator
  corrections, not a full study: our v9 vs propcorr_v1 fits are fine even without MC stat
  on the signal. Draft uses the two plots extracted from the SMP-V Keynote
  (`figures/propcorr/`). Ratio plot regenerated 5 Nov (variable shown inside the axes, lumi
  label); same values as 9 Sep, so both plots are v9 vs propcorr_v1/datacards_single.
  Still to do: explain cHd widening and the dipoles
  moving (cbWRe 0.70, cbBRe 1.28) or leave them unexplained.
- **Length (5 Nov):** §7.3 cut to two subsubsections (implementation and generation;
  linear and quadratic templates), cross-section table replaced by one sentence, W
  linearisation sentence dropped. cQe control figure kept separate (user's choice).
- 138 fb^-1 label on the 9.3 plots vs "2018 only" in the text: not worth raising.

- **Placement (Giacomo, 4 Nov):** new subsection right after §7.2 showing the parametrization
  with propagator corrections at linear and quadratic level (`propcorr_section_local.tex`, now
  titled "EFT parametrization with propagator corrections at LO LHE level"), and a new §9.4
  "Expected limits with propagator corrections" (`results_propcorr_local.tex`, same style as 9.3).

- **Only the width is corrected, not the mass.** We use the {m_W, m_Z, G_F} input scheme
  (`MwScheme` UFO), so m_Z and m_W are inputs and get no EFT shift (delta m_Z = 0 by
  construction). SMEFTsim also switches off its m_W shift parameter whenever
  `linearPropCorrections` is non-zero. At LO in neutral-current DY only the Z propagator
  appears, so the whole correction is delta Gamma_Z. (Draft §7.3, third paragraph.)
- Corrections to the original brief: the generation uses `NPall<=2 NPprop^2<=2`, not
  `NP<=1`; the per-operator LIN/QUAD figures are Figs. 57/58 (not 47/48); the LIN/QUAD
  algebra is Eq. 58 (restated); the definition our text cites is Eq. 55.
- Opening sentence reworded so it states only SMEFTsim's default (vertex insertions
  only), not what the AN's templates contain.
- **Sourcing pass (3 Nov).** Every paragraph of the draft has a `% src:` comment. Removed my
  own derivations (the -2 delta Gamma/Gamma peak formula and the interference-shape formula);
  replaced with guide Eq. (6.6) at delta m = 0, Eq. (6.10) and the on-shell statement of §6.
  Operator list from `dWZ` in the patched `parameters.py` plus `dGf`/`dMZ2` (= guide Table 4,
  topU3l row), cross-checked with the EFT2Obs auto-detect script. O(10%) linearisation
  statement is the guide's, made for the W: the draft says so.
