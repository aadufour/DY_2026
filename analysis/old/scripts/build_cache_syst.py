"""
build_cache_syst.py

Like build_cache.py but for the SYST_slc7 LHE files, which additionally
contain PDF systematic weights and QCD scale weights.

New cache keys vs. build_cache.py:
  pdf_325300   ndarray [n_events, 103]  NNPDF31_nnlo_as_0118_mc_hessian_pdfas
  pdf_325500   ndarray [n_events, 101]  NNPDF31_nnlo_as_0118_nf_4_mc_hessian
  scale_weights ndarray [n_events, 9]  standard MUR/MUF ∈ {0.5,1,2}^2 (no DYN_SCALE)

Weight key map (from LHE header):
  scale_weights columns → keys 1001,1006,1011,1016,1021,1026,1031,1036,1041
    (MUR,MUF): (1,1),(2,1),(0.5,1),(1,2),(2,2),(0.5,2),(1,0.5),(2,0.5),(0.5,0.5)
  pdf_325300 columns   → keys 1048–1150  (MemberID 0–102)
  pdf_325500 columns   → keys 1151–1251  (MemberID 0–100)
"""

import os
import pickle
import warnings
from itertools import combinations
import numpy as np
import pylhe

# ── Config ────────────────────────────────────────────────────────────────────
MLL_BIN_EDGES = [50, 120, 200, 400, 600, 800, 1000, 3000]

LHE_FILES = [
    f"/grid_mnt/data__data.polcms/cms/adufour/LHE/SYST_slc7/DYSMEFTMll{lo}_{hi}/unweighted_events.lhe"
    for lo, hi in zip(MLL_BIN_EDGES[:-1], MLL_BIN_EDGES[1:])
]

MLL_LO = 50.0
MLL_HI = 3000.0

OPERATORS = [
    'cHDD' , 'cHWB', 'cbWRe',
    'cbBRe', 'cHj1', 'cHQ1' ,
    'cHj3' , 'cHQ3', 'cHu'  ,
    'cHd'  , 'cHbq', 'cHl1' ,
    'cHl3' , 'cHe' , 'cll1' ,
    'clj1' , 'clj3', 'cQl1' ,
    'cQl3' , 'ceu' , 'ced'  ,
    'cbe'  , 'cje' , 'cQe'  ,
    'clu'  , 'cld' , 'cbl'
]

OP_PAIRS = list(combinations(OPERATORS, 2))

# 9 standard scale combinations: (MUR, MUF) without DYN_SCALE
SCALE_KEYS = ['1001', '1006', '1011', '1016', '1021', '1026', '1031', '1036', '1041']
SCALE_LABELS = [
    'mur1_muf1', 'mur2_muf1', 'mur05_muf1',
    'mur1_muf2', 'mur2_muf2', 'mur05_muf2',
    'mur1_muf05', 'mur2_muf05', 'mur05_muf05',
]

PDF_325300_KEYS = [str(i) for i in range(1048, 1151)]  # 103 members (MemberID 0–102)
PDF_325500_KEYS = [str(i) for i in range(1151, 1252)]  # 101 members (MemberID 0–100)

CACHE_FILE      = "/grid_mnt/data__data.polcms/cms/adufour/LHE/SYST_slc7/CACHE/lhe_cache_syst.pkl"
CHECKPOINT_FILE = "/grid_mnt/data__data.polcms/cms/adufour/LHE/SYST_slc7/CACHE/lhe_cache_syst_checkpoint.pkl"

# ── Kinematic functions ────────────────────────────────────────────────────────
def mll(p1, p2):
    p = np.array(p1) + np.array(p2)
    return np.sqrt(max(p[3]**2 - sum(p[i]**2 for i in range(3)), 0.0))

def rap(p1, p2):
    p  = np.array(p1) + np.array(p2)
    E, pz = p[3], p[2]
    return np.abs(0.5 * np.log((E + pz) / (E - pz)))

def cstar(p1, p2):
    p1 = np.array(p1); p2 = np.array(p2)
    p  = p1 + p2
    E, pz = p[3], p[2]
    mass  = mll(p1, p2)
    beta  = pz / E
    gamma = E / mass
    pz1_b = gamma * (p1[2] - beta * p1[3])
    p1mag = np.sqrt(p1[0]**2 + p1[1]**2 + pz1_b**2)
    return pz1_b / p1mag

# ── Load checkpoint if it exists ──────────────────────────────────────────────
if os.path.exists(CHECKPOINT_FILE):
    print(f"Resuming from checkpoint: {CHECKPOINT_FILE}")
    with open(CHECKPOINT_FILE, "rb") as f:
        ckpt = pickle.load(f)
    acc_mll          = ckpt["mll"]
    acc_rap          = ckpt["rap"]
    acc_cstar        = ckpt["cstar"]
    acc_w_SM         = ckpt["w_SM"]
    acc_w_p1         = ckpt["w_p1"]
    acc_w_m1         = ckpt["w_m1"]
    acc_w_pp         = ckpt["w_pp"]
    acc_pdf_325300   = ckpt["pdf_325300"]
    acc_pdf_325500   = ckpt["pdf_325500"]
    acc_scale        = ckpt["scale_weights"]
    completed_files  = ckpt["completed_files"]
    print(f"  {len(completed_files)} file(s) already done: {[os.path.basename(f) for f in completed_files]}")
    print(f"  {len(acc_mll):,} events loaded from checkpoint\n")
else:
    acc_mll          = np.empty(0, dtype=np.float64)
    acc_rap          = np.empty(0, dtype=np.float64)
    acc_cstar        = np.empty(0, dtype=np.float64)
    acc_w_SM         = np.empty(0, dtype=np.float64)
    acc_w_p1         = {op: np.empty(0, dtype=np.float64) for op in OPERATORS}
    acc_w_m1         = {op: np.empty(0, dtype=np.float64) for op in OPERATORS}
    acc_w_pp         = {pair: np.empty(0, dtype=np.float64) for pair in OP_PAIRS}
    acc_pdf_325300   = np.empty((0, 103), dtype=np.float64)
    acc_pdf_325500   = np.empty((0, 101), dtype=np.float64)
    acc_scale        = np.empty((0, 9),   dtype=np.float64)
    completed_files  = []

os.makedirs(os.path.dirname(CACHE_FILE), exist_ok=True)

# ── Collect ────────────────────────────────────────────────────────────────────
for lhe_file in LHE_FILES:
    if lhe_file in completed_files:
        print(f"Skipping (already done): {os.path.basename(lhe_file)}")
        continue

    print(f"\nReading {lhe_file}")
    if not os.path.exists(lhe_file):
        print("  WARNING: not found, skipping.")
        continue

    buf_mll        = []
    buf_rap        = []
    buf_cstar      = []
    buf_w_SM       = []
    buf_w_p1       = {op: [] for op in OPERATORS}
    buf_w_m1       = {op: [] for op in OPERATORS}
    buf_w_pp       = {pair: [] for pair in OP_PAIRS}
    buf_pdf_325300 = []   # list of 103-element lists
    buf_pdf_325500 = []   # list of 101-element lists
    buf_scale      = []   # list of 9-element lists

    pp_keys = {}

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        events = pylhe.read_lhe_with_attributes(lhe_file)

        for i, event in enumerate(events):
            if (i + 1) % 5000 == 0:
                print(f"  {i + 1} events processed")

            leptons = [
                p for p in event.particles
                if int(p.status) == 1 and abs(int(p.id)) in {11, 13}
            ]
            if len(leptons) < 2:
                continue

            lm = next((p for p in leptons if int(p.id) > 0), leptons[0])
            lp = next((p for p in leptons if int(p.id) < 0), leptons[1])
            v_lm = [lm.px, lm.py, lm.pz, lm.e]
            v_lp = [lp.px, lp.py, lp.pz, lp.e]

            m  = mll(v_lm, v_lp)
            y  = rap(v_lm, v_lp)
            cs = cstar(v_lm, v_lp)

            if not (MLL_LO <= m <= MLL_HI):
                continue

            if not pp_keys:
                for op1, op2 in OP_PAIRS:
                    pp_keys[(op1, op2)] = (
                        f'{op1}_{op2}' if f'{op1}_{op2}' in event.weights
                        else f'{op2}_{op1}'
                    )

            buf_mll.append(m)
            buf_rap.append(y)
            buf_cstar.append(cs)
            buf_w_SM.append(event.weights['SM'])
            for op in OPERATORS:
                buf_w_p1[op].append(event.weights[op])
                buf_w_m1[op].append(event.weights[f'minus{op}'])
            for pair in OP_PAIRS:
                buf_w_pp[pair].append(event.weights[pp_keys[pair]])

            buf_pdf_325300.append([event.weights[k] for k in PDF_325300_KEYS])
            buf_pdf_325500.append([event.weights[k] for k in PDF_325500_KEYS])
            buf_scale.append([event.weights[k] for k in SCALE_KEYS])

    n_kept = len(buf_mll)
    print(f"  {n_kept} events kept from this file")

    acc_mll   = np.concatenate([acc_mll,   np.array(buf_mll,   dtype=np.float64)])
    acc_rap   = np.concatenate([acc_rap,   np.array(buf_rap,   dtype=np.float64)])
    acc_cstar = np.concatenate([acc_cstar, np.array(buf_cstar, dtype=np.float64)])
    acc_w_SM  = np.concatenate([acc_w_SM,  np.array(buf_w_SM,  dtype=np.float64)])
    for op in OPERATORS:
        acc_w_p1[op] = np.concatenate([acc_w_p1[op], np.array(buf_w_p1[op], dtype=np.float64)])
        acc_w_m1[op] = np.concatenate([acc_w_m1[op], np.array(buf_w_m1[op], dtype=np.float64)])
    for pair in OP_PAIRS:
        acc_w_pp[pair] = np.concatenate([acc_w_pp[pair], np.array(buf_w_pp[pair], dtype=np.float64)])
    acc_pdf_325300 = np.concatenate([acc_pdf_325300, np.array(buf_pdf_325300, dtype=np.float64).reshape(-1, 103)], axis=0)
    acc_pdf_325500 = np.concatenate([acc_pdf_325500, np.array(buf_pdf_325500, dtype=np.float64).reshape(-1, 101)], axis=0)
    acc_scale      = np.concatenate([acc_scale,      np.array(buf_scale,      dtype=np.float64).reshape(-1, 9)  ], axis=0)

    completed_files.append(lhe_file)
    ckpt = {
        "mll":             acc_mll,
        "rap":             acc_rap,
        "cstar":           acc_cstar,
        "w_SM":            acc_w_SM,
        "w_p1":            acc_w_p1,
        "w_m1":            acc_w_m1,
        "w_pp":            acc_w_pp,
        "pdf_325300":      acc_pdf_325300,
        "pdf_325500":      acc_pdf_325500,
        "scale_weights":   acc_scale,
        "completed_files": completed_files,
    }
    with open(CHECKPOINT_FILE, "wb") as f:
        pickle.dump(ckpt, f)
    size_mb = os.path.getsize(CHECKPOINT_FILE) / 1e6
    print(f"  Checkpoint saved ({size_mb:.1f} MB, {len(completed_files)}/{len(LHE_FILES)} files done)")

print(f"\nTotal: {len(acc_mll):,} events")

# ── Save final cache ───────────────────────────────────────────────────────────
cache = {
    'mll':           acc_mll,
    'rap':           acc_rap,
    'cstar':         acc_cstar,
    'w_SM':          acc_w_SM,
    'w_p1':          acc_w_p1,
    'w_m1':          acc_w_m1,
    'w_pp':          acc_w_pp,
    'pdf_325300':    acc_pdf_325300,   # [n_events, 103]
    'pdf_325500':    acc_pdf_325500,   # [n_events, 101]
    'scale_weights': acc_scale,        # [n_events, 9]  — see SCALE_LABELS for column order
    'scale_labels':  SCALE_LABELS,
}

with open(CACHE_FILE, 'wb') as f:
    pickle.dump(cache, f)

size_mb = os.path.getsize(CACHE_FILE) / 1e6
print(f"Cache saved to {CACHE_FILE}  ({size_mb:.1f} MB)")

if os.path.exists(CHECKPOINT_FILE):
    os.remove(CHECKPOINT_FILE)
    print("Checkpoint removed.")
