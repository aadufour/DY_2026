"""
build_cache_pdf.py

Read all DY_pdf_LO_HI LHE files, extract kinematics + systematics weights,
save to a pickle cache. Supports checkpoint/resume: after each file the
accumulated arrays are saved so the script can be interrupted and restarted.

Cache structure:
    mll        : (N,)        dilepton invariant mass [GeV]
    rap        : (N,)        |dilepton rapidity|
    cstar      : (N,)        cos(theta*) in dilepton rest frame
    w_central  : (N,)        nominal weight  (MUR=1, MUF=1, PDF member 0, ID 45)
    w_pdf      : (100, N)    PDF replica weights (IDs 46-145, members 1-100)
    w_scale    : (8, N)      scale variation weights (IDs 1,6,11,16,25,30,35,40)
    scale_ids  : list[str]   ordered list of scale weight IDs (matches w_scale rows)

Usage (on llrcms):
    python3 build_cache_pdf.py
"""

import os
import pickle
import warnings

import numpy as np
import pylhe

# ── Config ────────────────────────────────────────────────────────────────────

MLL_EDGES = [50, 120, 200, 400, 600, 800, 1000, 3000]

LHE_FILES = [
    f"/home/llr/cms/adufour/MG5/mg5amcnlo/DY_pdf_{lo}_{hi}/Events/run_01/unweighted_events.lhe.gz"
    for lo, hi in zip(MLL_EDGES[:-1], MLL_EDGES[1:])
]

CACHE_FILE      = "/home/llr/cms/adufour/DY_2026/analysis_pdf/cache_pdf.pkl"
CHECKPOINT_FILE = "/home/llr/cms/adufour/DY_2026/analysis_pdf/cache_pdf_checkpoint.pkl"

CENTRAL_ID = '45'
PDF_IDS    = [str(i) for i in range(46, 146)]    # replicas 1-100

# All 8 non-nominal static scale variations (all 9 MUR/MUF combinations
# from {0.5,1,2}^2, excluding the nominal MUR=1 MUF=1 which is CENTRAL_ID)
SCALE_IDS = ['1', '6', '11', '16', '25', '30', '35', '40']

# ── Kinematic helpers ─────────────────────────────────────────────────────────

def mll(p1, p2):
    p = np.array(p1) + np.array(p2)
    return np.sqrt(max(p[3]**2 - sum(p[i]**2 for i in range(3)), 0.0))

def rap(p1, p2):
    p = np.array(p1) + np.array(p2)
    E, pz = p[3], p[2]
    return abs(0.5 * np.log((E + pz) / (E - pz)))

def cstar(p1, p2):
    p1 = np.array(p1); p2 = np.array(p2)
    p = p1 + p2
    E, pz = p[3], p[2]
    mass = mll(p1, p2)
    beta = pz / E
    gamma = E / mass
    pz1_b = gamma * (p1[2] - beta * p1[3])
    p1mag = np.sqrt(p1[0]**2 + p1[1]**2 + pz1_b**2)
    return pz1_b / p1mag

# ── Load checkpoint if it exists ──────────────────────────────────────────────

if os.path.exists(CHECKPOINT_FILE):
    print(f"Resuming from checkpoint: {CHECKPOINT_FILE}")
    with open(CHECKPOINT_FILE, 'rb') as f:
        ckpt = pickle.load(f)
    acc_mll      = ckpt['mll']
    acc_rap      = ckpt['rap']
    acc_cstar    = ckpt['cstar']
    acc_central  = ckpt['w_central']
    acc_pdf      = ckpt['w_pdf']      # (100, N_so_far)
    acc_scale    = ckpt['w_scale']    # (8,   N_so_far)
    done_files   = ckpt['completed_files']
    print(f"  {len(done_files)} file(s) done, {acc_mll.shape[0]:,} events loaded\n")
else:
    acc_mll     = np.empty(0, dtype=np.float64)
    acc_rap     = np.empty(0, dtype=np.float64)
    acc_cstar   = np.empty(0, dtype=np.float64)
    acc_central = np.empty(0, dtype=np.float64)
    acc_pdf     = np.empty((100, 0), dtype=np.float64)
    acc_scale   = np.empty((8,   0), dtype=np.float64)
    done_files  = []

# ── Process each LHE file ─────────────────────────────────────────────────────

for lhe_file in LHE_FILES:
    if lhe_file in done_files:
        print(f"Skipping (done): {os.path.basename(os.path.dirname(os.path.dirname(lhe_file)))}")
        continue

    print(f"\nReading {lhe_file}")
    if not os.path.exists(lhe_file):
        print("  WARNING: file not found, skipping.")
        continue

    buf_mll     = []
    buf_rap     = []
    buf_cstar   = []
    buf_central = []
    buf_pdf     = [[] for _ in PDF_IDS]    # 100 lists
    buf_scale   = [[] for _ in SCALE_IDS] # 8 lists

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for i, ev in enumerate(pylhe.read_lhe_with_attributes(lhe_file)):
            if (i + 1) % 5000 == 0:
                print(f"  {i + 1} events processed")

            leptons = [p for p in ev.particles
                       if int(p.status) == 1 and abs(int(p.id)) in {11, 13}]
            if len(leptons) < 2:
                continue

            lm = next((p for p in leptons if int(p.id) > 0), leptons[0])
            lp = next((p for p in leptons if int(p.id) < 0), leptons[1])
            v_lm = [lm.px, lm.py, lm.pz, lm.e]
            v_lp = [lp.px, lp.py, lp.pz, lp.e]

            m = mll(v_lm, v_lp)
            if not (MLL_EDGES[0] <= m <= MLL_EDGES[-1]):
                continue

            buf_mll.append(m)
            buf_rap.append(rap(v_lm, v_lp))
            buf_cstar.append(cstar(v_lm, v_lp))
            buf_central.append(ev.weights[CENTRAL_ID])
            for j, k in enumerate(PDF_IDS):
                buf_pdf[j].append(ev.weights[k])
            for j, k in enumerate(SCALE_IDS):
                buf_scale[j].append(ev.weights[k])

    n = len(buf_mll)
    print(f"  {n:,} events kept")

    # Convert to numpy and concatenate into accumulators
    acc_mll     = np.concatenate([acc_mll,     np.array(buf_mll,     dtype=np.float64)])
    acc_rap     = np.concatenate([acc_rap,     np.array(buf_rap,     dtype=np.float64)])
    acc_cstar   = np.concatenate([acc_cstar,   np.array(buf_cstar,   dtype=np.float64)])
    acc_central = np.concatenate([acc_central, np.array(buf_central, dtype=np.float64)])
    acc_pdf     = np.concatenate([acc_pdf,     np.array(buf_pdf,     dtype=np.float64)], axis=1)
    acc_scale   = np.concatenate([acc_scale,   np.array(buf_scale,   dtype=np.float64)], axis=1)

    # Save checkpoint
    done_files.append(lhe_file)
    ckpt = {
        'mll':             acc_mll,
        'rap':             acc_rap,
        'cstar':           acc_cstar,
        'w_central':       acc_central,
        'w_pdf':           acc_pdf,
        'w_scale':         acc_scale,
        'completed_files': done_files,
    }
    with open(CHECKPOINT_FILE, 'wb') as f:
        pickle.dump(ckpt, f)
    mb = os.path.getsize(CHECKPOINT_FILE) / 1e6
    print(f"  Checkpoint saved ({mb:.1f} MB, {len(done_files)}/{len(LHE_FILES)} files done)")

# ── Save final cache ──────────────────────────────────────────────────────────

print(f"\nTotal events: {acc_mll.shape[0]:,}")
print(f"w_pdf shape:  {acc_pdf.shape}   (100 replicas × N events)")
print(f"w_scale shape:{acc_scale.shape} (8 variations × N events)")

cache = {
    'mll':       acc_mll,
    'rap':       acc_rap,
    'cstar':     acc_cstar,
    'w_central': acc_central,
    'w_pdf':     acc_pdf,      # (100, N)
    'w_scale':   acc_scale,    # (8, N)
    'scale_ids': SCALE_IDS,    # ordered list so you know which row is which
}

os.makedirs(os.path.dirname(CACHE_FILE), exist_ok=True)
with open(CACHE_FILE, 'wb') as f:
    pickle.dump(cache, f)

mb = os.path.getsize(CACHE_FILE) / 1e6
print(f"Cache saved to {CACHE_FILE}  ({mb:.1f} MB)")

if os.path.exists(CHECKPOINT_FILE):
    os.remove(CHECKPOINT_FILE)
    print("Checkpoint removed.")
