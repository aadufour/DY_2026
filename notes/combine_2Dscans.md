# 2D EFT Scans — Missing Cross-Term Investigation & Fix

*(as of 2026-08-18, propcorr_v1)*

---

## The problem

Every "double" (2-operator) EFT fit run so far (`createWS.py 2`, `runScans.py 2 initial/scan`, all
351 pairs) was built from datacards containing only per-operator templates
(`sm`, `w1_<op>`, `wm1_<op>`). The **cross/interference term between the two
floating operators was never included** — every 2D fit has implicitly been
assuming the mixed term is zero, silently, with no error anywhere in the
chain. This note documents how that was found, confirmed, and fixed.

See also: `notes/combine.md` (1D/2D scan tooling), `notes/spritz.md`
(spritz pipeline), `~/code/TEX/TESI_DY/Chapters/Chapter3.tex` (thesis
derivation this whole investigation is cross-checked against).

---

## The physics

The cross section is a degree-2 polynomial in the Wilson coefficients
(thesis Eq. 3.x, Chapter3.tex:138-143):

```
σ(c) = σ_SM + Σ_α c_α σ_α^int + Σ_{α≤β} c_α c_β σ_αβ^quad
```

The sum over `α≤β` includes the 351 off-diagonal cross terms (α≠β), not
just the 27 diagonal quadratic terms. But the formula actually used to build
the Combine datacard (thesis Eq. 3.x, Chapter3.tex:246-250) only has a
single operator index:

```
N_k(c_α) = N_k^SM + c_α N_k^{int,α} + c_α² N_k^{quad,α}
```

— confirmed by the very next line: *"Fits are performed with one Wilson
coefficient floating at a time, with all others fixed to their SM value of
zero."* The 406-configuration generation methodology (SM + 27 linear + 27
quadratic + 351 cross, Chapter3.tex:213-224) correctly derives how to
extract every one of these components from the reweighted MC sample — the
cross terms were designed for and are generated, they just were never
plumbed through to Combine.

---

## Diagnosis chain (confirmed, file by file)

1. **`analysis/spritz/make_cards.py`** (`is_signal` check, originally
   line ~62) — only recognized `sm`, `w1_*`, `wm1_*` as signal processes. No
   cross-term naming pattern anywhere.

2. **`analysis/spritz/post_process.py`** (`apply_kfactor`, originally line
   ~182) — same gap: `eft_prefixes` hardcoded to `sm`/`w1_*`/`wm1_*` only.
   This file is otherwise generic (just histograms whatever `config.py`'s
   `samples` dict lists), so the real question was one level up.

3. **`analysis/spritz/config_propcorr_v1.py`** (`subsamples_eft`, lines
   74-81) — `rwgt = "events.LHEReweightingWeight"`, only extracts indices
   `0` (SM), `1..27` (w1), `28..54` (wm1). Comment even said "27 operators,
   matching LHE weight indices 1..27" with no mention of anything beyond 54.

4. **The raw MG5 reweight card** —
   `.../genproductions/.../DYSMEFTMll50_120_propcorr_reweight_card.dat`.
   Confirmed via `grep -c "^launch"` = **406** total configurations. The
   351 pair blocks (e.g. `launch --rwgt_name=cHWB_cbBRe` setting
   `SMEFT 9 1.0` and `SMEFT 25 1.0` simultaneously) are genuinely present —
   **the raw per-event weights already exist**, they're just never read.
   Index mapping, verified two independent ways (position in
   `itertools.combinations(OPERATORS, 2)`, and launch-block line-number
   arithmetic): `index = 55 + itertools.combinations(OPERATORS, 2).index((op_i, op_j))`
   for `op_i` before `op_j` in the canonical `OPERATORS` list order.
   Spot-checked: `cHDD_cHWB` → 55, `cHDD_cbWRe` → 56, `cHWB_cbBRe` → 82
   (both methods agree exactly).

5. **The physics model** —
   `HiggsAnalysis.AnalyticAnomalousCoupling.AnomalousCouplingMorphing_comb`
   (CMSSW, `CMSSW_14_1_0_pre4/src/HiggsAnalysis/AnalyticAnomalousCoupling/python/AnomalousCouplingMorphing_comb.py`).
   **Already fully supports the mixed term** — no changes needed here at all:
   - `doParametersOfInterest()` unconditionally builds a
     `func_sm_linear_quadratic_mixed_{bin}_{opi}_{opj}` RooFormulaVar
     (`"SM + Lin_i + Lin_j + Quad_i + Quad_j + 2*M_ij"`) whenever more than
     one operator is active for a bin — this has been sitting in every 2D
     workspace already built, just wired to nothing.
   - `getYieldScale_()` looks for a process literally named
     `w11_<op_i>_<op_j>` (either operator order accepted) to attach to that
     function. Since no such process has ever existed in any datacard, the
     lookup simply never matches — not "scaled to 0", just never a summand
     in the total yield at all. Functionally identical outcome, confirming
     the diagnosis precisely.
   - The search is scoped to the *currently active* operator subset only
     (`self.numOperators = len(active_ops)`), so a 2D build only ever needs
     the one relevant `w11_<opi>_<opj>` process, never all 351 at once.

---

## The fix

**`createWS.py` and the physics model need zero changes.** Everything is
upstream, in the three files that produce the datacard's raw process list:

1. **`analysis/spritz/config_propcorr_v1_crossterms.py`** (NEW file — a
   full copy of `config_propcorr_v1.py`, not an in-place edit, since this
   needs its own gridmount config dir / reprocessing run kept separate from
   the original production for comparison). Adds:
   - 351 more `subsamples_eft` entries: `w11_<opi>_<opj>` pulling
     `LHEReweightingWeight[:, 55 + combinations-index]`
   - 351 more `samples` dict entries (same naming pattern as `w1_`/`wm1_`)
   - 351 more `colors` entries (cosmetic)
   - `eft_samples` extended to include all `w11_*` names, which
     automatically wires them into every existing systematic-variation
     block (QCDscale, PDF, alphaS, PSWeight, rochester stat/syst, mu
     reco/idiso/trig, PU, prefireWeight) via the dict-comprehensions
     already present — no separate nuisance bookkeeping needed.

2. **`analysis/spritz/post_process.py`** (edited in place — generic,
   config-agnostic, purely additive change, safe for any config without
   `w11_*` histograms too) — `apply_kfactor`'s `eft_prefixes` now also
   includes all 351 `w11_<opi>_<opj>` names, so the k-factor rescaling
   (MiNNLO/SMEFTsim-LO-SM normalisation matching) applies to the cross-term
   histograms consistently with `sm`/`w1_*`/`wm1_*`.

3. **`analysis/spritz/make_cards.py`** (edited in place, same
   additive/safe reasoning) — `is_signal` now also recognizes
   `sample_name.startswith("w11_")`, so cross-term histograms get written
   into `shapes.root` as `histo_w11_<opi>_<opj>` and added to the
   datacard's process list with a signal index — exactly the process name
   `getYieldScale_` is already looking for.

All three files syntax-checked (`python3 -m py_compile`) after editing.
Index-mapping logic independently verified against the actual reweight card
(see diagnosis chain above).

---

## Deployment (LLR) — not yet run

Reprocessing needed (raw NanoAOD already has all 406
`LHEReweightingWeight` columns — **no new gridpack/LHE/NanoAOD
production required**, this is a rehistogramming pass over already-produced
files):

```bash
spritz-shell
mkdir -p /grid_mnt/data__data.polcms/cms/adufour/spritz_fabian/configs/propcorr_v1_crossterms
cp /grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/spritz/config_propcorr_v1_crossterms.py \
   /grid_mnt/data__data.polcms/cms/adufour/spritz_fabian/configs/propcorr_v1_crossterms/config.py
cd /grid_mnt/data__data.polcms/cms/adufour/spritz_fabian/configs/propcorr_v1_crossterms
mkdir -p data
```

**Open question before proceeding** — does `propcorr_v1` need a
`fileset.json` per-subsample patch (like the older `v7` pipeline in
`notes/spritz.md`), or does the internal `"subsamples": subsamples_eft`
mechanism in the dataset registration mean the *existing*
`propcorr_v1/data/fileset.json` (same physical files) can just be reused
as-is? Check:
```bash
grep -c "_w1_\|_wm1_" /grid_mnt/data__data.polcms/cms/adufour/spritz_fabian/configs/propcorr_v1/data/fileset.json
```
- `0` → reuse existing fileset.json (`cp .../propcorr_v1/data/fileset.json data/`), skip `spritz-fileset`
- nonzero → needs the same kind of patch as the documented `v7` pattern, extended for the 351 `w11_` pairs too

Then:
```bash
spritz-chunks
spritz-batch-llr
exit
cd condor && condor_submit submit.jdl
# ... wait for jobs ...
spritz-shell
cd /grid_mnt/data__data.polcms/cms/adufour/spritz_fabian/configs/propcorr_v1_crossterms
spritz-merge
spritz-postproc-eft
spritz-cards-eft
# then the usual combine workflow (createWS.py 2, runScans.py 2 initial/scan, ...)
# now automatically picks up and correctly uses the w11_ cross terms - no changes needed there.
```

Given this reprocesses the full dataset across all 7 mll bins, worth
checking the first job's output (or a small `njobs` test) before trusting
the whole batch — matches the same caution that caught the negative-yield
truncation issue earlier in the 2D-scan boundary work (see
`notes/combine.md`, "2D (Double) EFT Scans" chapter).
