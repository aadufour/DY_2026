# 1-Jet DY SMEFTsim Gridpacks

## Motivation

Added `p p > mu+ mu- j` to the existing 0-jet DY SMEFTsim process to include 1-jet final states. This brings the LO simulation closer to NLO accuracy (LO multileg / LO+1j approach). The operator scan is not redone with these samples — they are produced for methodology studies.

## Process cards

Two sets of gridpacks were produced:

### Propcorr (`linearPropCorrections = 1`)

```
generate p p > mu+ mu- QCD=0 SMHLOOP=0 NPall<=2 NPprop^2<=2
add process p p > mu+ mu- j SMHLOOP=0 NPall<=2 NPprop^2<=2
```

### Non-propcorr (`linearPropCorrections = 0`)

```
generate p p > mu+ mu- QCD=0 SMHLOOP=0 NP<=1
add process p p > mu+ mu- j SMHLOOP=0 NP<=1
```

Notes:
- `QCD=0` kept for the 0-jet process (pure EW)
- `QCD=0` dropped for the +jet process (jet requires QCD vertex)
- `j` defined as `g u c d s u~ c~ d~ s~` (no b quarks in jet)
- `ptj > 10 GeV` added to run card
- mll bins: 50-120, 120-200, 200-400, 400-600, 600-800, 800-1000, 1000-3000 GeV

## MLM matching

MLM jet matching is required to avoid double-counting between the 0-jet and 1-jet processes: without it, Pythia can shower the 0-jet process and produce a jet, overlapping with the explicit 1-jet matrix element.

Parameters added to all 1j run cards (both propcorr and nonpropcorr):

```
  1     = ickkw            ! 0 no matching, 1 MLM
  1.0   = alpsfact         ! scale factor for QCD emission vx
  False = chcluster        ! cluster only according to channel diag
  5     = asrwgtflavor     ! highest quark flavor for a_s reweight
  False = auto_ptj_mjj     ! Automatic setting of ptj and mjj if xqcut >0
  15.0  = xqcut            ! minimum kt jet measure between partons
```

`xqcut = 15 GeV` (generation ptj cut is 10 GeV; xqcut is somewhat arbitrary in the 10–15 GeV range). The corresponding Pythia GEN config should use `JetMatching:qCut = 25` (convention: qCut = xqcut + 10).

**Note**: the 0-jet gridpacks do NOT have MLM matching (`ickkw=0`). Matching is only applied when running the combined 0j+1j sample in Pythia.

## Model: linearPropCorrections

Controlled via `restrict_all_massless.dat` inside the SMEFTsim model tarball:
- Propcorr: `linearPropCorrections = 1`
- Non-propcorr: `linearPropCorrections = 0`

On the cluster, two tarballs are kept:
- `/grid_mnt/.../models/SMEFTsim_topU3l_MwScheme_UFO.tar.gz` — currently non-propcorr (0)
- `/grid_mnt/.../models/SMEFTsim_topU3l_MwScheme_UFO.tar.gz.propcorr` — propcorr (1)

**Important**: swap the tarball before submitting gridpack jobs depending on which version you need.

## Cards location (on cluster)

```
/grid_mnt/data__data.polcms/cms/adufour/genproductions/bin/MadGraph5_aMCatNLO/cards/DY_SMEFT_Gridpacks/
    DYSMEFTMll{bin}_1j_propcorr/   ← propcorr cards
    DYSMEFTMll{bin}_1j/            ← non-propcorr cards
```

## Gridpack tarballs (on cluster)

```
/grid_mnt/data__data.polcms/cms/adufour/genproductions/bin/MadGraph5_aMCatNLO/
    DYSMEFTMll{bin}_1j_propcorr_slc7_amd64_gcc700_CMSSW_10_6_19_tarball.tar.xz
    DYSMEFTMll{bin}_1j_slc7_amd64_gcc700_CMSSW_10_6_19_tarball.tar.xz   ← pending
```

Propcorr tarballs: all 7 bins complete (produced 2026-09-29, with MLM matching). ✓
Non-propcorr tarballs: all 7 bins complete (produced 2026-09-29, with MLM matching). ✓

## Submit scripts

```
gridpacks_1j_propcorr.submit   ← propcorr 1j
gridpacks_1j.submit            ← non-propcorr 1j
```

Both use `gridpack_job_propcorr.sh` as executable, `T3Queue = long`, 8 CPUs, 20 GB RAM.

**Submission order**: nonpropcorr first (model tarball is currently nonpropcorr). After nonpropcorr jobs finish, swap the model tarball to propcorr, then submit propcorr jobs:

```bash
# swap to propcorr
cd /grid_mnt/.../models/
mv SMEFTsim_topU3l_MwScheme_UFO.tar.gz SMEFTsim_topU3l_MwScheme_UFO.tar.gz.nonpropcorr
mv SMEFTsim_topU3l_MwScheme_UFO.tar.gz.propcorr SMEFTsim_topU3l_MwScheme_UFO.tar.gz
condor_submit gridpacks_1j_propcorr.submit
```

## Known issues

### Non-propcorr first submission failed (batch 815280)

The old non-propcorr cards (from `old/DYSMEFTMll*/`) contained `set max_t_for_channel 99` in the proc card. This command does not exist in MG5 v2.6.5 used on the cluster (`InvalidCmd` error), causing MG5 to abort before generating the output directory. Fix: removed the line from all 7 non-propcorr 1j proc cards with `sed -i '/set max_t_for_channel 99/d'`, cleaned up failed work directories, and resubmitted.

### `set zerowidth_tchannel True` and `set loop_color_flows False`

The old non-propcorr cards also contain these lines (absent from propcorr cards). These were not removed — they appear to be accepted by MG5 v2.6.5 (no error in propcorr logs where they are absent, but they were present in the original non-propcorr cards that previously worked for the 0-jet samples).

### Non-propcorr 1j run card issues (fixed 2026-09-29)

The old nonpropcorr 1j run cards (from `.bak`) had two problems:
1. `sde_strategy = 1` — not supported by MG5 v2.6.5, causes Fortran compilation failure (`Symbol 'sde_strategy' has no IMPLICIT type` in `setrun.f`). MG5 injects this into `run_card.inc` during the pilot run for 1j processes.
2. Wrong PDF: `nn23lo1`/230000 (NNPDF2.3 LO) instead of `lhapdf`/325300 (NNPDF3.1 NLO).

Fix: replaced all 7 nonpropcorr 1j run cards with copies of the propcorr 1j run cards (identical physics settings, just different proc card and model tarball distinguish the two).

### Non-propcorr model tarball had wrong internal path (fixed 2026-09-29)

The non-propcorr `SMEFTsim_topU3l_MwScheme_UFO.tar.gz` was originally packed from a Mac local path, giving internal structure `Users/albertodufour/MG5_2_9_18/mg5amcnlo/models/SMEFTsim_topU3l_MwScheme_UFO/...` instead of the correct `SMEFTsim_topU3l_MwScheme_UFO/...`. MG5 could not find `restrict_all_massless.dat`. Fixed by extracting, moving the directory to the correct relative location, and repacking. The propcorr tarball (`.propcorr` suffix) always had the correct structure.
