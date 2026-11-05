# LLR workspace map (surveyed 2026-11-05)

Access: `ssh llrcms` (proxy jump via llrgate01, port 22) -> llruicms01.in2p3.fr.

## Rules
- `$HOME = /home/llr/cms/adufour`: configs only (.bashrc, .ssh, .condor, .t3). Never work here.
- Work area: `~/work` -> `/grid_mnt/data__data.polcms/cms/adufour` (`cdwork`).
- Scratch: a subdir of `~/work`, removed afterwards. Never /tmp.
- Disk `/grid_mnt/data__data.polcms` is shared and 95% full (13 TB free of 236 TB).

## Shell environment (~/.bashrc)
- `spritz-shell`: apptainer `~/work/spritz-env.sif` (Python 3.10), PYTHONPATH -> `spritz_fabian/src`
- `spritz-shell-giacomo`: same sif, `SPRITZ_GIACOMO=1` -> `spritz_giacomo/src` + conda env `~/.conda/envs/spritz`
- `dy_analysis`, `dy_combine`, `dy_combine_morphing`: source `DY_2026/analysis/combine_tools/env_llr*.sh`
- PATH includes `DY_2026/analysis/spritz` and `DY_2026/analysis/combine_tools` (repo scripts are the commands)
- `tlog <name>`: tmux session logging to `./tmux-<name>.log` (that's where the tmux-*.log files come from)
- `t3stat`: condor queue status

## Work area (~/work), sizes
| dir | size | what |
|---|---|---|
| `spritz_fabian/` | 360G | old spritz env used for all the thesis results; outputs in `configs/` |
| `LHE/` | 137G | LHE validation: `baseline`, `propcorr`, `SYST_slc7`, `compare_ops_plots` (old/ + bugged deleted 2026-11-05) |
| `spritz_giacomo/` | ~5G | spritz on branch pr6 (PR #6 route, parked) |
| `gridpacks/` | 41G | `SYST_slc7/DYSMEFTMll*` (7 mll bins, current) + `old/` (incl. `propcorr_BUGGED`) |
| `spritz/` | 23G | original spritz checkout |
| `genproductions/` | 19G | gridpack production |
| `MG5/` | 11G | MadGraph |
| `CMSSW_14_1_0_pre4/`, `CMSSW_spritz/` | 0.7G each | combine (HiggsAnalysis) + `tools/` (checks, plotters, systematic_tools) |
| `combine/` | 0.6G | old cards (`spritz_cards/v6`), `validation/` |
| `DY_2026/` | 0.4G | this repo (git, main) |
| `3DYGeneration/` | 41M | crab/condor submission of NanoAOD production per year/order |
| `analysis_venv/`, `dy_venv/` | 0.3G | python venvs |
| `spritz-env.sif` | 1.5G | spritz container (copy also in $HOME) |

## Analysis outputs (`spritz_fabian/configs/`)
| config | contents |
|---|---|
| `propcorr_v1/` | **thesis baseline.** `condor/results_merged_new.pkl` (160 MB, merged, input for postproc), `histos.root`, `datacards_single/` + `datacards_double/` (inc_mm/{mll,rapll_abs,costhetastar,triple_diff}: datacard + 1D fits + summary plots), `no_PDF/`, `no_QCD/`, `no_theory/` (systematics studies), `plots*`, `summary_plots/` |
| `propcorr_v2/` | v1 + 351 cross terms; `condor/` deleted 2026-10-08, only a small `histos.root` + config left |
| `eft_bkg_fullsyst_v9/`, `_v5/` | non-propcorr predecessors, each with full `condor/` (3000 jobs) |

`spritz_giacomo/configs/`: `propcorr_new/` (merged pkl backup 1.2 GB; `merge_condor/` deleted 2026-11-05), `propcorr_pr6/` (config only), plus Giacomo's/Fabian's reference configs (`dy-eft-2018`, `dy-2018`, ...).

## Samples
Own propcorr NanoAODs on EOS GRIF:
`root://eos.grif.fr:1094//eos/grif/cms/llr/store/user/aldufour/3DY_SMEFTsim_LO_propcorr/...` (7 mll bins).
Reading them needs a valid grid proxy (`voms-proxy-init`); `~/work/proxy.pem` is from 2026-10-07.

## Open
- MiNNLO postproc (`post_process_minnlo.py`) can run directly on `propcorr_v1/condor/results_merged_new.pkl`.
- Checking `nLHEScaleWeight` on the own samples needs a fresh proxy.
- Disk cleanup candidates (ask before deleting): the `condor/job_*` dirs of v5/v9.
