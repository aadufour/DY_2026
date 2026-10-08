# spritz_giacomo setup notes

*Started 2026-09-29. Setting up Giacomo's `giacomo-pr-cleanup` branch on the LLR cluster.*

---

## What Giacomo's branch adds (fstaeg/spritz PR#4)

- `spritz-merge --condor`: merge step now runs on condor instead of locally
- **Megahisto runner** (`runner_3DY_eft_full_morphing_megahisto.py`): splits 406 EFT templates + 82,621 covariance pairs into batches sharing one `hist.Hist`. ~14× faster than original eval-based runner (30s vs 428s per chunk). Validated bit-for-bit identical to original.
- Corrupt file handling in fileset
- SMEFTatNLO config
- `renorm_samples`: bin-by-bin normalization of SMEFTsim LO templates to MiNNLO NNLO SM prediction
- Pickle histograms

Fabian also pushed to this branch: SM+EFT normalization (`renorm_samples`), pickle histograms,
fakes in postproc (`fakes_dict`), and (5-6 Oct, `623dade`..`95fca23`): EFT points and covariances in
separate histograms (memory), operator index lookup tool (`get_rw_idx`), updated `dy-eft-2018` config,
Giacomo's 1j LO EFT samples (0+1 jet merged).

Updating: `git status -sb && git fetch --all && git log --oneline HEAD..origin/giacomo-pr-cleanup`, then
`git diff --stat HEAD origin/giacomo-pr-cleanup` (check `spritz.yaml`: if changed, the conda env may need updates).
`pip install -e` means the conda env picks up pulled code immediately.

---

## spritz_giacomo directory

Created by rsync from `spritz_fabian` (code only, no job outputs, 36MB):

```bash
rsync -av --exclude='configs/' /grid_mnt/.../spritz_fabian/ /grid_mnt/.../spritz_giacomo/
cd /grid_mnt/.../spritz_giacomo
git stash
git checkout -b giacomo origin/giacomo-pr-cleanup
```

Location: `/grid_mnt/data__data.polcms/cms/adufour/spritz_giacomo`

---

## Config in use: `configs/propcorr_new` (2026-10-07)

Source: `analysis/spritz/config_propcorr_2018_megahisto_full.py` in DY_2026, copied to
`spritz_giacomo/configs/propcorr_new/config.py`. spritz_giacomo at `95fca23` (giacomo-pr-cleanup).

**Base = Fabian's `dy-eft-2018`** (backgrounds, nuisances, corrections, renorm, fakes) +
**EFT = Giacomo's megahisto setup** (27 operators, 406 points: sm, w1, wm1, w11).

- EFT samples: `DYMuMu_LO_EFT_SMEFTsim_propcorr_mll<bin>_Photos_startingOne`, 7 mll bins, GRIF EOS
- `LHEReweightingWeight` indices hardcoded (sm=0, w1=1-27, wm1=28-54, w11=55-405, `itertools.combinations` order);
  cross-checked = Fabian's `get_rw_idx` lookup tool for all 406 points in all 7 bins (0 mismatches)
- **No covariance pairs** (MC stat via autoMCStats, `noStat` on EFT templates, as in propcorr v1/v2)
- Naming as in Fabian's dy-eft-2018: datasets `DYmm_mll<bin>`, `DYmm_MiNNLO_M-<bin>`; samples `DYmm_<point>`, `DYmm_MiNNLO`
  (process names are applied only at postproc: can be renamed later without rerunning condor)
- **K-factor**: `renorm_samples` target `DYmm_MiNNLO`, reference `DYmm_sm` (bin-by-bin, in `spritz-postproc`)
- **`DYmm_MiNNLO` has `exclude_from_datacard`**: renormalized EFT templates replace it in the cards (no DY double counting)
- HO corrections on MiNNLO: N3LO QCD (3D, mll 30-3000) + NLO EW (mll 50-3000); nuisances also on EFT templates
- Region `bveto_mm` (b-veto) + `bveto_mm_ss` (fakes source); fakes from SS data
- **Blinding**: `data_mll_max: 500` on `bveto_mm` (data cut by runner, MC up to 3000); SS region `mll < 500` for all
- Theory systs on (`do_theory_variations: True`): QCDscale, PDF, alphaS, PS on bkgs + MiNNLO + EFT (as in v1/v2)
- triple_diff: mll=[50,70,80,90,100,110,120,150,200,300,500,700,1000,1500,3000], costheta 4 bins, rapll 4 bins
- njobs=6100

### Runner: `analysis/spritz/runner_3DY_eft_blind.py`
= Fabian's `runner_3DY_eft.py` (post `623dade`: EFT points and covariances in separate batch histograms)
+ 3 lines: data blinding via region key `data_mll_max` (same mechanism as our old `runner.py`).
Config points to it directly in DY_2026 (no copy into spritz_giacomo needed).
If Fabian updates `runner_3DY_eft.py`: re-copy as `runner_3DY_eft_fabian.py` and re-apply the 3 lines.

### Theory uncertainties: how they are combined (vs propcorr v2)
- Runner: each dataset gets theory variations from its own weights (EFT = MG5 LO weights, MiNNLO = its own).
- `post_process.renormalize_hists`: for each MiNNLO variation V, `k_V = MiNNLO_V / LO_sm_V`, `EFT_V = k_V * EFT_LO_V`.
  So `DYmm_sm` variation = MiNNLO (NNLO) variation exactly; EFT points = MiNNLO_V x (EFT/SM ratio at LO).
- **In propcorr v2 the full prediction had MG5-LO theory unc.; now only the EFT/SM ratio does.**
- Caveat: if a variation exists for MiNNLO but not for `DYmm_sm`, `k_var` is silently reused from the previous
  variation. Not an issue now (all MiNNLO nuisances are also on the EFT samples).

---

## Python environment

The old `spritz-env.sif` has Python 3.10. Giacomo's branch requires Python 3.12+ (f-string nested quote syntax in `batch.py`). Giacomo uses a mamba environment directly.

**Solution**: create a conda env from `spritz.yaml` in the repo, inside the apptainer (conda 24.9 is available in the sif). Two packages had to be removed from the yaml before creation:
- `correctionlib==0.1.dev155+...` — local dev version, not on PyPI
- `rucio-clients==35.5.0` — build fails (`pkg_resources` missing); installed separately after

```bash
# inside apptainer
sed -i '/correctionlib/d; /rucio/d' /tmp/spritz_giacomo.yaml
conda env create -f /tmp/spritz_giacomo.yaml
/home/llr/cms/adufour/.conda/envs/spritz/bin/pip install -e /grid_mnt/.../spritz_giacomo
/home/llr/cms/adufour/.conda/envs/spritz/bin/pip install rucio-clients
```

---

## bashrc changes

All tracked in `~/work/DY_2026/.bashrc`. Key sections:

### 1. Conditional PYTHONPATH / SPRITZ_PATH / PATH

```bash
if [ -z "$SPRITZ_GIACOMO" ]; then
    export PYTHONPATH=/grid_mnt/.../spritz_fabian/src:$PYTHONPATH
    export SPRITZ_PATH=/grid_mnt/.../spritz_fabian
else
    export PYTHONPATH=/grid_mnt/.../spritz_giacomo/src:$PYTHONPATH
    export SPRITZ_PATH=/grid_mnt/.../spritz_giacomo
    export PATH=/home/llr/cms/adufour/.conda/envs/spritz/bin:$PATH  # py3.12 first
fi
```

### 2. spritz-shell-giacomo alias

Sets `SPRITZ_GIACOMO=1` so the block above activates the giacomo paths, and binds the giacomo `samples.json`:

```bash
alias spritz-shell-giacomo='SPRITZ_GIACOMO=1 apptainer exec \
  -B /etc/grid-security/certificates:/etc/grid-security/certificates \
  -B /cvmfs -B /grid_mnt \
  -B /grid_mnt/.../spritz_giacomo/data/Full2018v9/samples/samples.json:/opt/spritz/data/Full2018v9/samples/samples.json \
  /grid_mnt/.../spritz-env.sif bash --rcfile ~/.bashrc'
```

### 3. `~/.local/bin` only added for fabian shell

```bash
if [ -z "$SPRITZ_GIACOMO" ]; then
    export PATH=$HOME/.local/bin:$PATH
fi
```

---

## Condor on LLR T3

`spritz-batch` must be run via the wrapper **`spritz-batch-llr`** (in `analysis/spritz/`, on PATH), inside `spritz-shell-giacomo`.
`condor_submit: command not found` at the end of spritz-batch is expected (no condor in the container); then submit
outside apptainer: `cd condor && condor_submit submit.jdl`.

The wrapper patches `submit.jdl` and `run.sh`:
- T3 queue (`T3Queue = short`, `WNTag = el9`, `include : /opt/exp_soft/cms/t3/t3queue |`)
- `request_memory=8192` (default 2048 is too low: jobs reach 4-7.7 GB)
- `run.sh`: apptainer + proxy copied from `/grid_mnt/.../proxy.pem`; with `SPRITZ_GIACOMO` set, uses the
  conda py3.12 (`/home/llr/cms/adufour/.conda/envs/spritz/bin/python`), needed for the runner's py3.12 syntax

**How spritz stores outputs**: `transfer_output_remaps = "results.pkl = $(Folder)/chunks_job.pkl"` — each job
overwrites its input `chunks_job.pkl` with the results (no separate results.pkl). Each chunk has a `result` field;
the runner skips chunks already done, so **resubmitting `submit.jdl` only redoes failed/missing chunks**.
`chunks_job.pkl` is zlib-compressed: read with `spritz.framework.framework.read_chunks`, not plain pickle.
Runner catches errors per chunk: **exit code 0 does not mean success** — check `err.txt` for `Error for chunk`.
`out.txt`/`err.txt` only arrive when the job ends (no live progress).

### propcorr_new run, 2026-10-07 (cluster 816402)

- Fileset: 6100 jobs x 5 chunks, 125k-401k events/job, exactly 1000 EFT events/job
- **1st attempt (cluster 816324) aborted**: proxy had expired -> every remote read failed with
  `XRootD error: [FATAL] Auth failed: No protocols left to try`. Fix: renew proxy, resubmit.
  Leftover `err.txt` from the removed cluster still show Auth failed (all written at 14:48, many nodes) — ignore.
- Memory: default `request_memory=2048` -> jobs held (`over cgroup memory limit`, usage 4.0-5.7 GB).
  Fixed live with `condor_qedit 816402 RequestMemory 6144 && condor_release 816402`, then 8192 (usage up to 7.7 GB):
  `condor_qedit -constraint 'ClusterId==816402 && (JobStatus==5 || JobStatus==1)' RequestMemory 8192 && condor_release 816402`
- Wall time: full jobs take ~1.5-1.7 h (max seen ~100 min); the `short` queue limit is **2 h (job killed)**.
  Jobs that only redid failed chunks of the 1st attempt took ~18 min (misleading early estimate).
  Runner writes results only at the end -> a job killed at 2 h loses all its chunks.
- Queues: only `short` and `long`. `long` is heavily throttled (~9 of our jobs ran at once vs ~600 on short):
  don't move many jobs there; use it only for the few that get killed on short.
  Queue is a job attribute: `condor_qedit -constraint '...' T3Queue '"long"'`.
- Output: **503 GB** in `condor/` (24,416 files). `/grid_mnt/data__data.polcms` was 98% full (5.2 TB free).
  Result size per job follows the EFT mass bin of its EFT chunk (compression of empty bins):
  mll50_120 / 120_200 up to ~113 MB, 200_400 ~86-93 MB, 400_600 ~72-84 MB.
- Done 2026-10-08 morning (~18 h after submission), plus retries:
  - 260 jobs missing chunks (all missing their EFT chunk): GRIF `TimeoutError: Operation expired` (transient,
    GRIF overloaded by many jobs reading the same EFT files), leftover Auth failed, few GRIF `unauthorized identity`.
    GRIF read test with proxy OK -> retried with `submit_retry.jdl` (only those folders), all recovered.
  - 7 jobs held for memory before the fix -> released with 10240 MB.
- **Jobs 6059-6099 (41) contain no EFT chunk** (6059 EFT chunks -> one each in job_0..job_6058): their result
  file is only ~2 MB but they are complete. Don't mistake them for failures.
- Final check: 0 untouched inputs (`find job_* -name chunks_job.pkl -size -500k`), no complete-size job with
  `Error for chunk`.
- Proxy valid until 2026-10-15.

### Checking which DONE jobs are actually incomplete

Condor exit 0 is not enough. Size of `chunks_job.pkl` (input ~10 kB, jobs w/o EFT ~2 MB, with EFT 70-113 MB):
```bash
condor_q $USER -af TransferInput | grep -o 'job_[0-9]*' | sort -u > queued.txt
find job_* -name chunks_job.pkl -printf "%s %h\n" > sizes.txt
awk '$1<1e7 {print $2}' sizes.txt | grep -vxF -f queued.txt > incomplete.txt   # then inspect, see below
```
Inspect contents (inside `spritz-shell-giacomo`; only small files, fast): chunks with `result == {}` are missing.
```python
from spritz.framework.framework import read_chunks
for c in read_chunks("job_N/chunks_job.pkl"):
    print(c["data"]["dataset"], len(c["result"]), (c.get("error") or "")[:80])
```
Errors in `err.txt` from unreadable files are `Exception: ('Error, could not read any of the filenames ...` with the
XRootD errors inside one string: grep for `Operation expired|Auth failed|unauthorized identity|Unable to open file`.
Retry only some folders:
```bash
sed "s/^queue 1 Folder in .*/queue 1 Folder in $(paste -sd, incomplete.txt)/" submit.jdl > submit_retry.jdl && condor_submit submit_retry.jdl
```

### Merge on condor (`spritz-merge --condor`), 2026-10-08

Run **from the config dir** (`configs/propcorr_new`, not `condor/`: it appends `condor/` itself), inside
`spritz-shell-giacomo`:
```bash
spritz-merge --condor -dr      # prepares merge_condor/ (100 groups of ~61 jobs, --njobs to change), no submit
```
The generated files are for **CERN lxplus** and must be patched for LLR before submitting:
- `run.sh`: Fabian's afs proxy path + bare `python` -> replaced by apptainer + conda py3.12 + spritz_giacomo
  (no proxy needed: the worker reads `job_*/chunks_job.pkl` by absolute path on `/grid_mnt`, no input transfer)
- `special_start.sh`: pointed to `spritz_fabian/src` -> spritz_giacomo
- `submit.jdl`: remove `MY.SingularityImage`, `use_x509userproxy`, `+JobFlavour`; add T3 queue lines;
  `request_memory=65536` (default 2048; groups need 15-50 GB, see below)
```bash
cd merge_condor && G=/grid_mnt/data__data.polcms/cms/adufour/spritz_giacomo && cat > run.sh <<EOF
#!/bin/bash
time apptainer exec -B /cvmfs -B /grid_mnt \\
    --env SPRITZ_PATH=$G --env PYTHONPATH=$G/src:$G \\
    /grid_mnt/data__data.polcms/cms/adufour/spritz-env.sif \\
    /home/llr/cms/adufour/.conda/envs/spritz/bin/python merge_worker.py .
EOF
chmod +x run.sh && echo "export PYTHONPATH=$G/src:\$PYTHONPATH" > special_start.sh && sed -i -e '/MY.SingularityImage/d' -e '/use_x509userproxy/d' -e '/+JobFlavour/d' -e 's/^request_memory *=.*/request_memory=65536/' -e '/^log/a T3Queue = short\nWNTag   = el9\ninclude : /opt/exp_soft/cms/t3/t3queue |' submit.jdl
condor_submit submit.jdl     # outside apptainer
```
Each merge job overwrites `merge_condor/job_N/chunks_job.pkl` with `{"real_results": merged}`. When all are done:
`spritz-merge --merge` (local) -> `results_merged_new.pkl`. **Watch memory**: the local final merge holds all
7 EFT mass bins x 406 points x ~146 variations at once (possibly 10+ GB RAM on the login node); check the size
of the merged group files first. TODO: turn the patches into a `spritz-merge-llr` wrapper.

Status: merge submitted 2026-10-08 (cluster 816414, 100 jobs).

**Merge memory is large** (usage 15-45+ GB per job, 12 GB request was not enough):
- `add_dict` is correct (lists of batch histograms are summed element-wise, any key), memory is real.
- Groups are **interleaved** (merge job k gets `job_k, job_k+100, job_k+200, ...`), so every group holds all
  7 EFT mass bins -> more merge jobs (`--njobs`) would NOT reduce memory.
- Per EFT mass bin in RAM: triple_diff 576 bins (with flow) x 407 points x 3 region slots x 147 syst x 16 B
  ~1.65 GB + 1D variables ~0.6 GB = ~2.3 GB; x 7 bins ~16 GB, plus ~2-3 copies during `add_dict`.
- Fix: `condor_qedit -constraint 'ClusterId==816414 && JobStatus==5' RequestMemory 65536 && condor_release 816414`
  (48 GB still too little for one group). Worker nodes have 96-128 GB, so 64 GB slots exist (fewer per node).
  Use 65536 directly next time.
- Login node (llruicms01) has 187 GB (163 available): the final local `spritz-merge --merge` (expected
  ~50 GB+) should fit; run it at a quiet moment with `/usr/bin/time -v` to record the peak (needed for postproc).
- Next runs: booking EFT histograms only for `triple_diff` (`eft_variables`) would cut ~25% of EFT memory.

### Useful monitoring commands

```bash
condor_q -totals
condor_q 816402 -hold -af HoldReason | cut -c1-140 | sort | uniq -c
condor_q 816402 -run -af '(time()-EnteredCurrentStatus)/60' | sort -n | tail -3       # longest running (min)
# wall time of finished jobs (min), from the `time` line at the end of err.txt
find job_* -name err.txt -newermt "2026-10-07 14:50" -exec grep -h "^real" {} + | awk '{split($2,a,"m"); print a[1]+0}' | sort -n | tail
# failed chunks per dataset
grep -h "Error for chunk" job_*/err.txt | grep -o "'dataset': '[^']*'" | sort | uniq -c
/opt/exp_soft/cms/t3/t3stat -q          # jobs per user/queue/status
```
Avoid `condor_history -constraint ...` without `-limit`: it scans the whole schedd history and hangs for minutes.

---

## Next steps

1. ~~Batch jobs~~ done 2026-10-08, all 6100 complete
2. Merge: wait for the 100 `merge_condor` jobs (check holds/memory, `err.txt`), then `spritz-merge --merge`
3. Write `spritz-merge-llr` wrapper with the patches above
4. `spritz-postproc`: check histogram names to confirm renorm matches variations at nuisance level (QCDscaleUp etc.)
5. `spritz-cards-eft`; ask Giacomo whether it works without `covariance.root`
6. Process names are `DYmm_<point>`: our combine tools (`AnomalousCouplingMorphing_comb`) look up `w11_<i>_<j>`,
   so either rename samples to bare names (postproc only) or adapt the tools
7. Future runs: book EFT histograms only for `triple_diff` (`eft_variables` option) to cut output size;
   needs a check that postproc/renorm tolerate missing EFT 1D histograms
8. Later: covariance pairs (Giacomo's MC-stat covariance for morphing); renorm now handles them in Fabian's config
