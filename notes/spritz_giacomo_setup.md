# spritz_giacomo setup notes

*Started 2026-09-29. Setting up Giacomo's `giacomo-pr-cleanup` branch on the LLR cluster.*

---

## What Giacomo's branch adds (fstaeg/spritz PR#4)

- `spritz-merge --condor`: merge step now runs on condor instead of locally
- **Megahisto runner** (`runner_3DY_eft_full_morphing_megahisto.py`): splits 406 EFT templates + 82,621 covariance pairs into batches sharing one `hist.Hist`. ~14× faster than original eval-based runner (30s vs 428s per chunk). Validated bit-for-bit identical to original.
- Corrupt file handling in fileset
- SMEFTatNLO config
- `renorm_samples`: bin-by-bin normalization of SMEFTsim LO templates to MiNNLO NNLO SM prediction (`target: DYmm_NNLO`, `reference: DYmm_LO_sm`) — currently only in `dy-eft-2018`, not yet in megahisto config
- Pickle histograms

Fabian also pushed to this branch: SM+EFT normalization (`renorm_samples`), pickle histograms.

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

## Config in use

`configs/dy-eft-full-smeftsim-propcorr-2018-megahisto`

- 2018 only, SMEFTsim propcorr NanoAODs from Giacomo (`root://eos.grif.fr:1094//eos/grif/cms/llr//store/user/gboldrin/3DY_SMEFTsim_LO/`)
- 27 operators, 406 EFT points, 82,621 covariance pairs
- Triple-diff binning: mll=[50,120,200,400,600,800,1000,3000], costheta=[-1,0,1], rapll=[0,1.2,2.4] (coarse, to be refined)
- njobs=6100, megahisto runner
- Only lumi + stat nuisances. **No `renorm_samples` yet** (to be added later from `dy-eft-2018`)

---

## bashrc changes

Two changes made to `~/.bashrc` (tracked in `~/work/DY_2026/.bashrc`):

### 1. spritz-shell-giacomo alias

Uses `SPRITZ_GIACOMO=1` env var to switch paths inside `--rcfile ~/.bashrc`:

```bash
alias spritz-shell-giacomo='SPRITZ_GIACOMO=1 apptainer exec \
  -B /etc/grid-security/certificates:/etc/grid-security/certificates \
  -B /cvmfs -B /grid_mnt \
  -B /grid_mnt/.../spritz_giacomo/data/Full2018v9/samples/samples.json:/opt/spritz/data/Full2018v9/samples/samples.json \
  /grid_mnt/.../spritz-env.sif bash --rcfile ~/.bashrc'
```

### 2. Conditional PYTHONPATH / SPRITZ_PATH (lines 48-54)

```bash
if [ -z "$SPRITZ_GIACOMO" ]; then
    export PYTHONPATH=/grid_mnt/.../spritz_fabian/src:$PYTHONPATH
    export SPRITZ_PATH=/grid_mnt/.../spritz_fabian
else
    export PYTHONPATH=/grid_mnt/.../spritz_giacomo/src:$PYTHONPATH
    export SPRITZ_PATH=/grid_mnt/.../spritz_giacomo
fi
```

### 3. Added `~/.local/bin` to PATH

```bash
export PATH=$HOME/.local/bin:$PATH
```

---

## spritz installation inside apptainer

The `.sif` has Python 3.10. Giacomo's branch requires >=3.10 but uses f-string nested quote syntax that only works in 3.12+. Need to use the sif's own pip:

```bash
/usr/local/envs/spritz/bin/pip install -e /grid_mnt/.../spritz_giacomo --user
```

Scripts installed to `~/.local/bin/` (now on PATH).

---

## Open issue: Python 3.10 f-string syntax error

`batch.py` line 168 uses:
```python
command += f"cp {get_fw_path()}/data/{an_dict["year"]}/cfg.json {job_dir}/; "
```

This nested double-quote syntax only works in Python 3.12+. The sif has Python 3.10.15. Giacomo likely uses a newer sif. Waiting for his reply on Mattermost.

**Fix**: change `an_dict["year"]` → `an_dict['year']` in that f-string (trivial, one character).

---

## Next steps

1. Resolve Python version issue with Giacomo
2. Run `spritz-fileset configs/dy-eft-full-smeftsim-propcorr-2018-megahisto`
3. Run `spritz-batch configs/dy-eft-full-smeftsim-propcorr-2018-megahisto`
4. After jobs finish: `spritz-merge`, `spritz-postproc`, `spritz-postproc-eft`, `spritz-cards-eft`
5. Add `renorm_samples` from `dy-eft-2018` to megahisto config before postproc
6. Consider finer triple-diff binning (costheta 5 bins, rapll 4 bins) — discuss with Giacomo
