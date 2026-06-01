# LHE Pipeline — DY SMEFT Analysis Notes

---

## Overview

The LHE pipeline provides a parton-level (detector-free) validation of the EFT
constraints. It reads MadGraph LHE files directly, builds a cache of event
weights, then fills histograms for combine.

```
LHE files  →  build_cache_new.py  →  lhe_cache_new.pkl
                                           │
                                    build_datacard_new.py
                                           │
                              histograms.root + datacard.txt
                                           │
                          createCombineJson.py --binname quad_
                          createWS_lhe.py 1
                                           │
                              runScans → scan_*.root
```

---

## 1. LHE File Structure

A Les Houches Event (LHE) file is a plain-text XML-like file produced by
MadGraph. It has two main sections: the **header** and the **event block**.

### 1.1 Header

The header contains a `<initrwgt>` block that lists all weight variants the
generator computed for each event. Each variant is tagged with a unique numeric
`id` and attributes describing the physics hypothesis:

```xml
<header>
  <initrwgt>
    <weightgroup name="sm_weights">
      <weight id="1"> MUR=1 MUF=1 PDF=325300 </weight>   <!-- nominal SM -->
    </weightgroup>

    <weightgroup name="scale_variation" combine="envelope">
      <weight id="2"> MUR=1   MUF=0.5 PDF=325300 </weight>
      <weight id="3"> MUR=1   MUF=2   PDF=325300 </weight>
      <weight id="4"> MUR=0.5 MUF=1   PDF=325300 </weight>
      <weight id="5"> MUR=2   MUF=1   PDF=325300 </weight>
      <weight id="6"> MUR=0.5 MUF=0.5 PDF=325300 </weight>
      <weight id="7"> MUR=2   MUF=2   PDF=325300 </weight>
      <!-- (MUR=2,MUF=0.5) and (MUR=0.5,MUF=2) excluded: anti-correlated -->
    </weightgroup>

    <weightgroup name="PDF_325300" combine="hessian">
      <weight id="10"> MUR=1 MUF=1 PDF=325300 </weight>  <!-- member 0 = central -->
      <weight id="11"> MUR=1 MUF=1 PDF=325301 </weight>  <!-- member 1 -->
      ...
      <weight id="112"> MUR=1 MUF=1 PDF=325402 </weight> <!-- member 102 -->
    </weightgroup>

    <weightgroup name="SMEFT_reweighting">
      <weight id="200"> SM </weight>
      <weight id="201"> cHDD </weight>       <!-- cHDD=+1, all others=0 -->
      <weight id="202"> minuscHDD </weight>  <!-- cHDD=-1, all others=0 -->
      <weight id="203"> cHWB </weight>
      <weight id="204"> minuscHWB </weight>
      ...
      <weight id="500"> cHDD_cHWB </weight>  <!-- cHDD=+1, cHWB=+1 (cross term) -->
      ...
    </weightgroup>
  </initrwgt>
</header>
```

Key attributes on each `<weight>` tag:

| Attribute | Meaning |
|-----------|---------|
| `MUR` | μR/μ₀ — renormalization scale ratio |
| `MUF` | μF/μ₀ — factorization scale ratio |
| `PDF` | LHAPDF set ID |
| `DYN_SCALE` | Dynamic scale flag — **always skipped** |
| `ALPSFACT` | αs emission variation — **always skipped** |

`parse_weight_ids()` in `build_cache_new.py` reads this header by regex and
returns the relevant weight IDs grouped by type.

### 1.2 Events

```xml
<event>
  <eventinfo weight="+1.234e-03" ... />
  <!-- status=1: final state particles -->
  <particle id="-13" status="1" px="12.3" py="-5.6" pz="44.1" e="46.2" />  <!-- μ⁺ -->
  <particle id="13"  status="1" px="-12.3" py="5.6" pz="55.0" e="56.3" />  <!-- μ⁻ -->
  <weights>
    <wgt id="1">   1.234e-03 </wgt>  <!-- SM nominal -->
    <wgt id="2">   1.180e-03 </wgt>  <!-- MUF=0.5 -->
    <wgt id="3">   1.290e-03 </wgt>  <!-- MUF=2   -->
    ...
    <wgt id="200"> 1.234e-03 </wgt>  <!-- SMEFT SM -->
    <wgt id="201"> 1.356e-03 </wgt>  <!-- cHDD c=+1 -->
    <wgt id="202"> 1.112e-03 </wgt>  <!-- cHDD c=-1 -->
    ...
  </weights>
</event>
```

The **event weight** (`eventinfo.weight`) is the nominal generation weight in
pb. For unweighted (flat-weight) generation:

```
|eventinfo.weight| = σ_total / N_generated   (pb)
→  sum over all events = σ_total
```

The `<wgt>` values inside `<weights>` are reweighted versions of the same
phase-space point for alternative physics hypotheses. They are in the same
units (pb) and carry the same phase-space Jacobian — only the matrix element
squared changes.

### 1.3 SMEFT weight naming convention in our gridpack

| Key in `event.weights` | Meaning |
|------------------------|---------|
| `'SM'` | All Wilson coefficients = 0 |
| `'cHDD'` | cHDD = +1, all others = 0 |
| `'minuscHDD'` | cHDD = −1, all others = 0 |
| `'cHDD_cHWB'` | cHDD = +1 and cHWB = +1 simultaneously |

The order of operators in cross-weight keys (`op1_op2` vs `op2_op1`) is not
guaranteed — `build_cache_new.py` tries both and uses whichever is found.

---

## 2. Cache Building (`build_cache_new.py`)

### 2.1 Why a cache?

Reading the LHE files event-by-event with `pylhe` is slow (~10–30 min for all
7 mll bins combined). The cache pre-extracts all kinematics and weights into
compact numpy arrays that load in a few seconds for repeated use.

### 2.2 Kinematic variables computed

For each event passing the mll cut (50–3000 GeV):

**Invariant mass** m_ll:
```python
p = p1 + p2   # 4-vector sum (px, py, pz, E)
mll = sqrt(E² - px² - py² - pz²)
```

**Rapidity** |y_Z|:
```python
y = |0.5 × ln((E + pz) / (E - pz))|
```

**Collins-Soper angle** cos θ*:
The lepton angle in the Z rest frame relative to the beam bisector. Sensitive
to the forward-backward asymmetry A_FB:
```python
beta  = pz_Z / E_Z          # Z boost along beam
gamma = E_Z / mll
pz1_boosted = gamma × (pz_l1 - beta × E_l1)   # boost lepton to Z rest frame
costheta = pz1_boosted / |p_l1_boosted|
```

### 2.3 Cache structure

```python
cache = {
    'mll':           float64[N],          # invariant mass (GeV)
    'rap':           float64[N],          # |y_Z|
    'cstar':         float64[N],          # cos θ*
    'w_SM':          float64[N],          # SMEFT SM weight (pb)
    'xwgt':          float64[N],          # raw eventinfo.weight
    'w_p1':          {op: float64[N]},    # weight at c=+1 per operator
    'w_m1':          {op: float64[N]},    # weight at c=−1 per operator
    'w_pp':          {(op1,op2): float64[N]},  # cross weights (pairs)
    'w_scale':       {id: float64[N]},    # MUF scale variation weights
    'w_pdf_central': float64[N],          # standalone PDF-central weight
    'pdf_325300':    float64[N×103],      # all 103 PDF members (5-flav)
}
```

### 2.4 Scale variations — MUF only, why not MUR

At **Leading Order (LO)** the matrix element for qq̄→Z has exactly one power
of αs (from the Z production vertex... actually DY is αs⁰ at Born level,
EW process). More precisely:

- **μR** controls the running of αs. At LO there are no QCD loops, so αs
  appears at most as an overall constant. Varying μR has no logarithmic
  sensitivity — the dependence cancels at NLO and is numerically negligible
  at strict LO.
- **μF** controls where to split collinear radiation between the PDF and the
  hard scattering. Even at LO, changing μF changes which PDFs are evaluated and
  at what scale, producing a genuine uncertainty through DGLAP evolution.

**Consequence:** we keep only MUF-only variations (`MUR=1, MUF=0.5` and
`MUR=1, MUF=2`) — a **2-point envelope** instead of the 6- or 8-point envelope
used at NLO.

### 2.5 PDF set — 5-flavour only, no 4-flavour

The gridpack was generated in the **5-flavour scheme**: the b quark is treated
as massless and included as a parton inside the proton. This corresponds to
LHAPDF set **325300** (NNPDF31_nnlo_as_0118_mc_hessian_pdfas).

The 4-flavour scheme (set 325500) treats the b as massive and excludes
b-initiated diagrams. Using a 4-flavour PDF with a 5-flavour matrix element
is scheme-inconsistent. It is therefore dropped entirely from `build_cache_new.py`.

### 2.6 Checkpoint / resume

Processing all 7 LHE files takes ~30 min. After each file completes, a
checkpoint is saved. Re-running the script automatically resumes:

```bash
python3 build_cache_new.py             # full run or resume
python3 build_cache_new.py --nodoubles # skip operator-pair cross weights (faster)
python3 build_cache_new.py --nevents 5000  # test run with subset of events
```

---

## 3. Datacard and Histogram Building (`build_datacard_new.py`)

### 3.1 EFT decomposition

The EFT cross section for a single operator with Wilson coefficient c:

```
σ(c) = σ_SM  +  c · σ_lin  +  c² · σ_quad
```

The three components are extracted analytically from the three raw templates:

```python
w_lin  = 0.5 × (w_p1 − w_m1)          # linear / interference term
w_quad = 0.5 × (w_p1 + w_m1) − w_SM   # quadratic (pure EFT²) term
w_slq  = w_SM + C×w_lin + C²×w_quad   # SM+lin+quad at reference C=1
```

**Why w_lin can be negative:** the interference between the EFT amplitude and
the SM amplitude can be destructive, giving σ_lin < 0 for some operators and
some phase-space regions. Combine cannot handle negative histogram bins.
`AnomalousCouplingEFTNegative_comb` solves this by storing `quad_{op}` and
`sm_lin_quad_{op}` (always ≥ 0 for reasonable C) and reconstructing the linear
term internally.

For two operators simultaneously (cross-term for 2D fits):
```python
w_inter = w_pp[op1,op2] − w_p1[op1] − w_p1[op2] + w_SM
w_mixed = C1 × C2 × w_inter
```

### 3.2 Normalization fix

The cache stores un-normalized weights. Each weight ≈ σ_total per event (not
σ/N_gen), so `sum(w_SM) ≈ N_gen × σ` instead of `σ`.

The fix is applied immediately after loading:

```python
N_gen = len(w_SM)       # total number of generated events in cache
w_SM  = w_SM / N_gen    # sum(w_SM) now ≈ σ_total (pb)
# same rescaling applied to w_p1, w_m1, w_pp, w_scale, w_pdf_central, pdf_arr
```

After the fix:
```
sum(w_SM) × LUMI ≈ σ_total × LUMI ≈ N_expected_events  ✓
```

Without this fix `sum(w_SM) × LUMI ≈ N_gen × N_expected ≈ 10¹²`, which causes
MINUIT to lose numerical precision and fail to converge.

### 3.3 Histogram filling

```python
h.fill(mll_arr, weight=weights * LUMI)
```

Each bin accumulates:
```
bin_i = Σ_{events in bin i} (w_event / N_gen) × LUMI  ≈  σ_bin × LUMI  =  N_expected_in_bin
```

This is a **yield histogram** (expected event counts), which is what combine
expects.

### 3.4 Binning

The mll binning matches the RECO analysis (`config_dy_smeftsim_v7/v8`) exactly:

```
50–75 GeV:    5 GeV steps   (pre-Z, 5 bins)
76–105 GeV:   2 GeV steps   (Z peak, 15 bins — resolves the lineshape)
106–119 GeV:  5 GeV steps   (post-Z shoulder, 3 bins)
120–3000 GeV: variable      (120, 150, 200, 250, 300, 400, 600, 800, 1000, 1500, 3000)
```

**Total: 34 bins.** Matching the RECO binning makes the LHE and RECO likelihood
scan constraints directly comparable bin-by-bin.

Note: finer binning = more bins with few events at LHE level (~100k events
total spread over 34 bins). The 2 GeV Z-peak bins will each have O(1k) events,
which is enough for smooth templates but autoMCStats will add uncertainties
there.

### 3.5 Systematic shape nuisances

**MUF scale (2-point envelope):**
```python
scale_hists = [make_hist(w_nom × (w_scale[k] / w_pdf_central)) for k in scale_keys]
h_up.value[b] = max(scale_hists[:][b])    # bin-wise maximum
h_dn.value[b] = min(scale_hists[:][b])    # bin-wise minimum
```
The ratio `w_scale[k] / w_pdf_central` reweights from the nominal to the
alternative scale point. `w_pdf_central` (not `w_SM`) is used as denominator
because the scale weight group in MadGraph is defined relative to the PDF
central member weight, not the SMEFT-reweighted SM weight.

**PDF (asymmetric quadrature sum):**
For `mc_hessian` sets, each non-central member is a ±1σ eigenvector. The PDF
uncertainty in each histogram bin b is:
```python
dev = replica_values[:, b] - nominal[b]    # deviation of each replica
σ_up[b] = sqrt(Σ dev²  for dev > 0)       # upward quadrature
σ_dn[b] = sqrt(Σ dev²  for dev < 0)       # downward quadrature
```
This is the **quadrature sum** (not the RMS), because eigenvectors are
orthogonal and their contributions add in quadrature.

### 3.6 Datacard process indices

`AnomalousCouplingEFTNegative_comb` uses a specific sign convention:

```
sm                     → index −1    (SM reference)
quad_{op1}             → index −2
sm_lin_quad_{op1}      → index −3
quad_{op2}             → index −4
sm_lin_quad_{op2}      → index −5
...
sm_lin_quad_mixed_...  → large negative indices
```

All indices are negative — combine treats all processes as signal-type EFT
templates. This is different from the RECO morphing pipeline where `sm` has
index +1 (background).

### 3.7 Usage

```bash
dy_analysis
cd /grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/combine/5_flav_bins

python3 /grid_mnt/.../combine_tools/build_datacard_new.py \
    --all_op \
    --output histograms.root \
    --datacard datacard.txt
# --lumi defaults to 59740 pb^-1 (2018 full)
```

Then follow the LHE combine workflow in `notes/combine.md`.

---

## 4. Files Reference

| File | Location | Purpose |
|------|----------|---------|
| `build_cache_new.py` | `analysis/combine_tools/` | Build `lhe_cache_new.pkl` |
| `build_datacard_new.py` | `analysis/combine_tools/` | Build histograms + datacard |
| `build_cache_syst.py` | `analysis/combine_tools/` | Old cache (6-point scale, 4+5 flav PDF) |
| `build_datacard_reco_bins.py` | `analysis/combine_tools/` | Old datacard (needs `--pdf-flavour`) |
| `lhe_cache_new.pkl` | `LHE/SYST_slc7/CACHE/` | Active cache (MUF-only, 5-flav PDF) |
| `lhe_cache_syst.pkl` | `LHE/SYST_slc7/CACHE/` | Old cache (kept for reference) |
| LHE files | `LHE/SYST_slc7/DYSMEFTMll{lo}_{hi}/` | Raw MadGraph output, 7 mll bins |
