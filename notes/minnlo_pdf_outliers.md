# MiNNLO DY: unphysical PDF weights at high mass (2026-11-08 to 11-11)

## Problem
With the MiNNLO theory variations (`post_process_minnlo.py`: per-variation k = DYll_i / sm_i), the PDF
nuisance on `sm` and the EFT templates blows up at high mll:

| mll bin | MG5 LO cards (v1) | MiNNLO cards |
|---|---|---|
| < 1 TeV | <= 3% | same |
| 1-1.5 TeV | ~3% | ~6% |
| 1.5-3 TeV | ~4% | +30% / -28% |

triple_diff high-mass bins: 100-350%.

## How the PDF nuisance is built
Each event fills H_nom with w and H_i with w * r_i, r_i = LHEPdfWeight[i]/LHEPdfWeight[0], i = 1..100,
same binning. Per bin: sigma = sqrt(sum_i (H_i - H_nom)^2), Up/Down = H_nom +- sigma (`kind: square`).
One event with r_i ~ 1e5 dominates H_i in its bin, hence sigma. The MiNNLO k-factor copies DYll's relative
variation per member onto `sm` and the EFT templates.

## Cause
Single events in `DYJetsToMuMu_M-2000toInf_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos`
(RunIISummer20UL18NanoAODv9, one file `.../2550000/5424027A-33C9-BC44-9BBE-20FF3E40ED4B.root`, 300k events):

| run:lumi:event | member | ratio | LHE mll | reco mll* | share of the sample's 1.5-3 TeV yield moved |
|---|---|---|---|---|---|
| 1:130:86524 | 78 | +1.05e5 | 2581 | 2564 | -80% |
| 1:141:94231 | 52 | -9.9e4 | 3078 | 2997 | +48% |
| 1:219:146025 | 73 | -1.6e4 | 2029 | 2681 | -32% |

(*two leading muons, no analysis selection; sign flips from negative genWeight.)
In DYll (all mass samples) the last bin moves by -18.8% (78), +11.7% (52), -7.5% (73); median member 0.4%.

Not a binning or generation-range problem: no overflow folding in spritz (events with reco mll > 3 TeV are
dropped), and these events sit at 2-3 TeV, x ~ 0.15-0.25, where the PDFs are constrained.

Fabian (`dy-2018`) uses the same sample with the same `square` treatment but his binning stops at 500 GeV;
Giacomo's `dy-eft-2018` comments out the MiNNLO bins from 700to800 up. Neither looked at this mass range.

## Veto test (`configs/propcorr_v1_veto/`)
`runner_pdfveto.py` = runner.py + veto of the 3 events. Only `DYmm_M-2000toInf` reprocessed (3 chunks),
swapped into a copy of v1's merged pkl, then `post_process_minnlo.py` (4 workers).
Last mll bin PDF quadrature on DYll: **24.3% -> 6.8%**; bins below 1.5 TeV unchanged.
Remaining: member 58 gives -4.3% in 1-1.5 TeV (from `M-1000to1500`).

## Full scan: it is a tail, not 3 events
`scan_pdf_weight_tail.py`, first file (up to 300k events) per sample, max_i |r_i - 1| per event:

| sample | > 10x | > 1e3 | > 1e4 | max shift of one member's sample yield from events > 10x |
|---|---|---|---|---|
| 10to50 | 0.00% | 0 | 0 | 0.1% |
| 100to200 | 0.00% | 1 | 0 | 0.8% |
| 200to400 | 0.00% | 0 | 0 | 0.4% |
| 400to500 | 0.01% | 0 | 0 | 0.2% |
| 500to700 | 0.01% | 3 | 2 | 5.4% (m51) |
| 700to800 | 0.02% | 0 | 0 | 0.2% |
| 800to1000 | 0.03% | 3 | 0 | 1.7% |
| 1000to1500 | 0.08% | 6 | 0 | 3.1% (m58) |
| 1500to2000 | 0.31% | 1 | 0 | 17.7% (m83), only 34k events read |
| 2000toInf | 1.25% | 82 | 9 | 68.3% (m78) |

50to100 not scanned (different dataset name).

9 events > 1e4 in 2000toInf, traced through our runner's selection (copy of runner.py keeping only them):

| run:lumi:event | member | ratio | selection | analysis mll |
|---|---|---|---|---|
| 1:141:94231 | 52 | -9.9e4 | passes (vetoed) | ~2998 |
| 1:130:86524 | 78 | +1.0e5 | passes (vetoed) | ~2604 |
| 1:219:146025 | 73 | -1.6e4 | passes (vetoed) | 1939 or 2066 |
| 1:443:296350 | 81 | -1.4e4 | passes | ~1730 |
| 1:110:73328 | 3 | +1.1e4 | passes | ~1949 |
| 1:336:224715 | 76 | -1.7e4 | passes | 1939 or 2066 |
| 1:300:200146 | 74 | +1.2e4 | passes | ~3600 (outside bins) |
| 1:16:10077 | 10 | -5.1e4 | fails lepton cut | - |
| 1:35:23216 | 92 | +4.3e4 | fails l2Tight | - |

(trace printed event ids sorted but mll in file order: 146025/224715 pairing ambiguous, both in 1.5-3 TeV.)
So 3 more ~1e4 events are in the 1.5-3 TeV bin after the veto; each moves its member by only ~1-2%, hidden
among normal members, part of the remaining 6.8%. 296350's quick reco mll (5.5 TeV) differs from the analysis
one (~1.7 TeV) because the analysis picks a different muon pair.

## Plots (`plots_minnlo_compare/`, LLR `propcorr_v1_minnlo/pdf_check/`, `propcorr_v1_veto/pdf_check/`)
- `nuis_PDFweight_sm_mll_mg5.png`, `nuis_PDFweight_sm_mll_minnlo.png`: the nuisance combine sees (drawNuisances.py --label)
- `pdf_members_DYll_mll.png`: 100 members, before veto (52, 58, 78 highlighted)
- `pdf_members_DYll_mll_veto.png`: after veto; highlighting is wrong (30, 35 are normal members), to redo with explicit members
- `pdf_event_weights.png`: per event, ratio vs member index for the 3 vetoed events vs the central 99% band

## Scripts
- `analysis/spritz/scan_pdf_weight_tail.py`: standalone scan (uproot/awkward/numpy), thresholds, worst events with LHE/reco mll and bin
- `analysis/spritz/plot_pdf_event_weights.py`: per-event member ratios
- `analysis/combine_tools/plot_pdf_members.py`: per-member binned variation
- `analysis/combine_tools/drawNuisances.py`: `--label` added
- `analysis/spritz/runner_pdfveto.py`: runner + 3-event veto (test only)

## Open
- scan 50to100 and full 1500to2000
- general rule instead of a hand-picked veto, e.g. drop events with any |r| > 100, all DY mass samples (condor)
- redo the after-veto members plot with explicit members
- rebuild mll card + fits from the veto/cut histos (node was loaded by another user)

## Update 2026-11-12
- Combine nuisance on `sm` (MiNNLO card, mll last bin): +-29.0% (MG5 card: +-3.4%). The 24.3% above is DYll in
  histos.root; why `sm` gets 29% instead of 24% through the k-factor is not checked.
- Full scan of all DY samples: `propcorr_v1_veto/pdf_check/scan_all.txt` (50to100 = `DYJetsToMuMu`: clean).
  1500to2000 has 2 events at ~1e5 (155212 m93, 101904 m55) + 146952 (8.9e3): all three FAIL the selection
  (lepton cut, l2Tight, earlier), so irrelevant for us.
- Member plot with only 52/73/78 labelled by event: `pdf_members_DYll_mll_events.png`.
- Message sent to Giacomo and Fabian (2026-11-12) with the 4 plots; waiting for their answer.
- Next time: scan only events passing the analysis selection (runner copy dumping per-event max ratio in
  inc_mm, on condor) and/or an approximate `--select` in the standalone script.
