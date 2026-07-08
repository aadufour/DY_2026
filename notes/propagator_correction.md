# Propagator correction — implementation notes (Jul 2026)

Goal: add propagator-correction contributions to the DY SMEFT gridpack
("especially important when you are on-shell" — see todo.md). Fabian gave
me a patched UFO model to base this on.

---

## 1. What Fabian actually gave me

Folder: `SMEFTsim_topU3l_MwScheme_propagatorhack_UFO/` (kept in this repo as
a reference — **not** something that needs to be installed/used directly,
just the source of 4 files).

It is **not** a full new model. Only 4 files differ from the stock
`SMEFTsim_topU3l_MwScheme_UFO` model:
- `coupling_orders.py` — adds a new coupling order `NPprop`
  (`expansion_order = 99`) and `NPall` (`expansion_order = 99`, didn't
  exist before)
- `couplings.py`
- `restrict_massless.dat`
- `restrict_SMlimit_massless.dat`

Everything else (`propagators.py`, `parameters.py`, `particles.py`,
`vertices.py`, and crucially all ~130 per-operator `restrict_<op>_massless.dat`
reweight cards) is untouched and must be preserved.

Fabian's instructions:
1. Replace those 4 files inside the existing `SMEFTsim_topU3l_MwScheme_UFO`
   model (wherever it's actually used for generation).
2. In the proc card, change the process definition from
   `NP<=1` → `NPall<=2 NPprop^2<=2`.

---

## 2. Where the model actually needed patching

Two separate places, both required for the production pipeline (the local
MG5 install at `/grid_mnt/.../MG5/mg5amcnlo/models/...` is a third copy,
used only for quick interactive checks — **not required for production**,
and as of writing is still unpatched since it's not on the critical path).

### a) The genproductions gridpack-build tarball
```
/grid_mnt/data__data.polcms/cms/adufour/genproductions/bin/MadGraph5_aMCatNLO/cards/SMEFTsim_topU3l_MwScheme_UFO.tar.gz
```
This is what `gridpack_generation.sh` unpacks into a fresh MG5 install on
each condor worker node (see line ~260: `tar -axvf .../SMEFTsim_topU3l_MwScheme_UFO.tar.gz --strip-components=5 -C $WORKDIR/$MGBASEDIRORIG/models/`).

**Gotcha:** the tarball's internal path
(`Users/albertodufour/MG5_2_9_18/mg5amcnlo/models/SMEFTsim_topU3l_MwScheme_UFO/...`)
has exactly 5 leading directory segments before the model folder, matching
the hardcoded `--strip-components=5`. This looks like garbage (it's a leaked
Mac path from whoever originally packed it) but is load-bearing — do **not**
flatten it when repacking, or extraction breaks. `gridpack_generation.sh` is
shared across other processes in that `cards/` dir too, so don't edit the
strip-components number there either; just keep repacking with the same
5-segment structure.

Patch procedure used (done in place, original saved as `.tar.gz.bak`):
```bash
cd .../MadGraph5_aMCatNLO/cards
tar xzf SMEFTsim_topU3l_MwScheme_UFO.tar.gz
MODEL=Users/albertodufour/MG5_2_9_18/mg5amcnlo/models/SMEFTsim_topU3l_MwScheme_UFO
cp <repo>/SMEFTsim_topU3l_MwScheme_propagatorhack_UFO/{coupling_orders.py,couplings.py,restrict_massless.dat,restrict_SMlimit_massless.dat} $MODEL/
tar czf SMEFTsim_topU3l_MwScheme_UFO.tar.gz Users
rm -rf Users
```
Confirmed applied via `diff` against the coworker's reference files (shows
`NPall`/`NPprop` present) and against the `.tar.gz.bak` original (shows the
diff only where expected).

### b) Proc cards (all 7 mll bins)
`import model SMEFTsim_topU3l_MwScheme_UFO-all_massless` stays the same —
`all_massless` is my own reweight-base restriction card (all ops ≈ 1), not
one of Fabian's 4 changed files, and doesn't need to change since operator
indices/names are unaffected by the patch.

Changed line 39/40 in every `*_proc_card.dat`:
```
generate p p > mu+ mu- QCD=0 SMHLOOP=0 NP<=1
```
→
```
generate p p > mu+ mu- QCD=0 SMHLOOP=0 NPall<=2 NPprop^2<=2
```
via `sed -i 's/NP<=1$/NPall<=2 NPprop^2<=2/'` across all 7
`DYSMEFTMll*_proc_card.dat` files.

---

## 3. Card/output organization (to avoid clobbering the baseline)

Renamed everything with a `_propcorr` suffix so it can never collide with
the existing non-propagator-corrected production:

- Cards: `cards/DY_SMEFT_Gridpacks/DYSMEFTMll<bin>_propcorr/` (moved out of
  a scratch `propagator_corr/` subfolder up to sit alongside the baseline
  bin folders — `gridpack_job.sh`'s path convention expects
  `cards/DY_SMEFT_Gridpacks/$1` directly, no subfolder). `output` line
  inside each proc card updated to match (`output DYSMEFTMll<bin>_propcorr`).
- Built gridpack tarballs: `DYSMEFTMll<bin>_propcorr_slc7_amd64_gcc700_CMSSW_10_6_19_tarball.tar.xz`
  in `.../MadGraph5_aMCatNLO/`.
- Gridpack-build submit: custom `gridpack_job_propcorr.sh` +
  `gridpacks_propcorr.submit` (8 cores, `T3Queue=long`, matching the
  original working recipe — found and rejected two stale variants first:
  `gridpacks_1core.submit`/`gridpack_job.sh` used only 1 core and pointed
  at a `SYST/` cards subfolder that doesn't exist for this bin set; the
  `old/gridpack_job_bkp.sh` also had the same stale `SYST/` path).
- LHE event-generation output + logs (kept **out of the git repo**):
  `/grid_mnt/data__data.polcms/cms/adufour/LHE/propcorr/DYSMEFTMll<bin>_propcorr/{lhe_output/,logs/}`
  — sibling to the existing `LHE/SYST_slc7/` and `LHE/old/`.
- Cache: `/grid_mnt/data__data.polcms/cms/adufour/LHE/propcorr/CACHE/lhe_cache_propcorr.pkl`

---

## 4. Gridpack build (condor)

7 jobs, one per bin, via `gridpacks_propcorr.submit` (8 cores, 20G,
`T3Queue=long`). All 7 finished with exit code 0, tarballs ~47-48MB each.

**Confirmed the patch actually took effect** (not just "ran without
crashing"): the build log shows real diagram generation under the new
syntax —
```
INFO: Generating Feynman diagrams for Process: b b~ > mu+ mu- NPall<=2 QCD=0 SMHLOOP=0 NPprop^2<=2 @1
...
Generated helas calls for 4 subprocesses (129 diagrams) in 0.365 s
```
160 total diagrams (vs whatever the plain `NP<=1` count was — not logged
for direct comparison, would be worth doing via the diagram-isolation check
below if I want a cleaner number).

## 5. Event generation (condor, parallelized)

Got a script pair from a colleague (Kirill, via Biriukov) that uses
HTCondor's `transfer_input_files` to ship each job its own private copy of
the gridpack tarball into an ephemeral per-job scratch directory — this is
what actually allows true parallelism (MadEvent can't have two jobs sharing
one `process/madevent/` dir, since it uses a `RunWeb` lock file; the first,
naive attempt at "more jobs" would have needed manual per-job directory
bookkeeping, this sidesteps that entirely).

Adapted as `launch_lhe_propcorr.sh` (generic, takes `$TARBALL`/`$SEEDID`/`$SAVEFILE`
as env vars) + one `.sub` file per bin
(`launch_lhe_<bin>_propcorr.sub`), each with a `seq` range of 10 unique
seeds → `queue SampleId from seq <base+1> <base+10> |`. 10 jobs × 10000
events = 100k events/bin, 70 jobs total, fully parallel.

Result: 70/70 files produced, `output_<seed>.lhe`, 653MB each (matches
1000-event test scaled ×10: 66MB×10≈660MB). One job took ~2h44m wall time
but only used ~1 of the requested 8 cores effectively (`Cpus: Usage 0.99/8`)
— worth investigating if this needs to be faster in future, not chased
down this time.

## 6. Cache building

`build_cache_new.py`/`build_cache_parallel.py` expect one merged LHE file
per bin. Rather than physically merging the 10 files/bin (~45GB total,
wasteful I/O), wrote `analysis/combine_tools/build_cache_propcorr.py` —
adapted so each bin's worker loops over its 10 files internally and
concatenates in memory (parses weight IDs once from the first file per bin,
since headers are identical across a bin's files).

Run: `python3 build_cache_propcorr.py --nodoubles --workers 16 --eta-max 2.4 --pt-lead 25 --pt-sub 10`
(note: `--workers 16` doesn't actually help — task granularity is one
task per bin, so max useful parallelism = 7 regardless of `--workers`).

Result: 660,151 events kept (of 700,000 read) after fiducial cuts. Cut
efficiency is strongly mass-dependent in a physically sensible way — the
pT cut (leptons >25/10 GeV) removes ~10% of events in the 50-120 GeV bin
(soft leptons near threshold) but is essentially irrelevant above 400 GeV
(plenty of energy to go around); η cut is roughly flat (~2-5%) across all
bins, consistent with pure geometric acceptance. PDF member count (103) and
scale-variation key count (2, MUF-only envelope) match the documented
baseline convention exactly (see lhe.md) — confirms the patch didn't
disturb the standard systematics machinery.

`build_datacard_new.py` needs **no code changes** — `--cache`/`--output`/
`--datacard` are already CLI args, just point `--cache` at
`lhe_cache_propcorr.pkl`. `N_GEN_PER_SAMPLE = 100_000` hardcoded constant
is still correct (same per-bin event count as baseline). Not yet run for
propcorr as of writing.

---

## 7. Physics validation attempted so far

Compared normalized mll shape (0.5 GeV bins, 80-101.5 GeV window) between
`lhe_cache_propcorr.pkl` and the baseline `lhe_cache_parallel.pkl`:

- First tried `w_SM` (all Wilson coefficients = 0) — result was
  statistics-dominated noise, and in hindsight probably the wrong quantity:
  if the propagator correction is an EFT-induced effect (tied to a
  coupling order that only activates with nonzero Wilson coefficients), it
  may vanish identically at the pure SM point.
- Tried `cll1` operator-on weight — got identical results to `w_SM` in
  both samples. **Not a bug**: `cll1` is a purely-leptonic 4-fermion
  contact operator with no tree-level diagram for quark-initiated DY, so
  reweighting to `cll1=1` trivially reproduces the SM prediction. Wrong
  operator to probe with.
- Tried `cHDD` (Higgs-doublet kinetic term — directly touches Z/W
  propagator normalization, more relevant candidate) with a proper pull/χ²
  test (Poisson errors from both samples): **χ²/dof = 1.44 over 43 bins**
  — a mild (~2σ) excess above pure noise, not a strong/conclusive
  detection either way. Inconclusive at current statistics (100k
  events/bin, ~10k events per fine bin at the peak).

**Strongest confirmation so far is not statistical but structural**: the
`NPprop` coupling order exists in the patched model, the proc card
correctly shows `NPall<=2 NPprop^2<=2` (not the old `NP<=1`), and the
process genuinely generated 160 real diagrams under the new syntax (not
zero, not silently falling back to old behavior).

---

## 8. Open items / possible next steps

- Diagram-level isolation check (proposed, not yet run since it turned out
  to need the local MG5 install patched too, which isn't otherwise
  needed): compare diagram counts for `NPprop^2==0` vs `NPprop^2==2` in an
  interactive `mg5_aMC` session — would give an unambiguous count of
  exactly how many diagrams involve the propagator correction, cleaner
  than the statistical shape test.
- Run `build_datacard_new.py --cache .../lhe_cache_propcorr.pkl` to get
  actual histograms/datacard for the propcorr sample.
- If a real shape effect needs to be resolved cleanly, current statistics
  (100k events/bin) are probably not enough — the χ² test above is right
  at the edge of sensitivity for a percent-level effect.
- The 2h44m-per-job / ~1-core-effective-usage event generation timing is
  unexplained and unoptimized — fine for one production run, might matter
  if this needs to be repeated often.
