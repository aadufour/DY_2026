"""
rw_build_cache.py

Used for rw_triple_diff.py when merging all DY_all_min_max LHE files.

Read all LHE files once, extract kinematics and all operator weights,
save to a pickle. Subsequent runs of rw_triple_diff.py load the pickle
instead of re-parsing the LHE for faster loading.

Supports checkpointing: progress is saved after each LHE file, so the
script can be interrupted and restarted without re-processing completed files.

Memory strategy: Python lists are used only within each file (small, temporary).
After each file they are converted to numpy arrays and concatenated into the
accumulated arrays. This keeps long-lived storage compact (~10x less memory
than keeping everything as Python lists).

Usage:
    python3 rw_build_cache.py
"""

import os
import pickle
import warnings
from itertools import combinations
import numpy as np
import pylhe

# ── Config (must match rw_triple_diff.py) ─────────────────────────────────────
MLL_BIN_EDGES = [50, 120, 200, 400, 600, 800, 1000, 3000]
# LHE_FILES = [
#     f"/Users/albertodufour/MG5_2_9_18/mg5amcnlo/DY_all_{lo}_{hi}/myLHE/unweighted_events.lhe"
#     for lo, hi in zip(MLL_BIN_EDGES[:-1], MLL_BIN_EDGES[1:])
# ]

# updated for MG5 llr paths
LHE_FILES = [
    f"/grid_mnt/data__data.polcms/cms/adufour/MG5/mg5amcnlo/DYSMEFTMll{lo}_{hi}/Events/run_01/unweighted_events.lhe.gz"
    for lo, hi in zip(MLL_BIN_EDGES[:-1], MLL_BIN_EDGES[1:])
]






MLL_LO = 50.0
MLL_HI = 3000.0

# updated with the "b quark" operators
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

# CACHE_FILE      = "/Users/albertodufour/code/DY2026/analysis/lhe_cache.pkl"
# CHECKPOINT_FILE = "/Users/albertodufour/code/DY2026/analysis/lhe_cache_checkpoint.pkl"

# updated for MG5 llr paths
CACHE_FILE      = "/grid_mnt/data__data.polcms/cms/adufour/MG5/mg5amcnlo/CACHE/lhe_cache.pkl"
CHECKPOINT_FILE = "/grid_mnt/data__data.polcms/cms/adufour/MG5/mg5amcnlo/CACHE/lhe_cache_chekpoint.pkl"

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
# Accumulated arrays are numpy (compact). Python lists only live within one file.
if os.path.exists(CHECKPOINT_FILE):
    print(f"Resuming from checkpoint: {CHECKPOINT_FILE}")
    with open(CHECKPOINT_FILE, "rb") as f:
        ckpt = pickle.load(f)
    acc_mll         = ckpt["mll"]
    acc_rap         = ckpt["rap"]
    acc_cstar       = ckpt["cstar"]
    acc_w_SM        = ckpt["w_SM"]
    acc_w_p1        = ckpt["w_p1"]
    acc_w_m1        = ckpt["w_m1"]
    acc_w_pp        = ckpt["w_pp"]
    completed_files = ckpt["completed_files"]
    print(f"  {len(completed_files)} file(s) already done: {[os.path.basename(f) for f in completed_files]}")
    print(f"  {len(acc_mll):,} events loaded from checkpoint\n")
else:
    acc_mll         = np.empty(0, dtype=np.float64)
    acc_rap         = np.empty(0, dtype=np.float64)
    acc_cstar       = np.empty(0, dtype=np.float64)
    acc_w_SM        = np.empty(0, dtype=np.float64)
    acc_w_p1        = {op: np.empty(0, dtype=np.float64) for op in OPERATORS}
    acc_w_m1        = {op: np.empty(0, dtype=np.float64) for op in OPERATORS}
    acc_w_pp        = {pair: np.empty(0, dtype=np.float64) for pair in OP_PAIRS}
    completed_files = []

# ── Collect ────────────────────────────────────────────────────────────────────
for lhe_file in LHE_FILES:
    if lhe_file in completed_files:
        print(f"Skipping (already done): {os.path.basename(lhe_file)}")
        continue

    print(f"\nReading {lhe_file}")
    if not os.path.exists(lhe_file):
        print("  WARNING: not found, skipping.")
        continue

    # Temporary Python lists for this file only
    buf_mll   = []
    buf_rap   = []
    buf_cstar = []
    buf_w_SM  = []
    buf_w_p1  = {op: [] for op in OPERATORS}
    buf_w_m1  = {op: [] for op in OPERATORS}
    buf_w_pp  = {pair: [] for pair in OP_PAIRS}

    # Pre-resolve cross-weight keys once from the first event
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

            # Resolve cross-weight key names once
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

    n_kept = len(buf_mll)
    print(f"  {n_kept} events kept from this file")

    # Convert buffers to numpy and concatenate into accumulators
    acc_mll   = np.concatenate([acc_mll,   np.array(buf_mll,   dtype=np.float64)])
    acc_rap   = np.concatenate([acc_rap,   np.array(buf_rap,   dtype=np.float64)])
    acc_cstar = np.concatenate([acc_cstar, np.array(buf_cstar, dtype=np.float64)])
    acc_w_SM  = np.concatenate([acc_w_SM,  np.array(buf_w_SM,  dtype=np.float64)])
    for op in OPERATORS:
        acc_w_p1[op] = np.concatenate([acc_w_p1[op], np.array(buf_w_p1[op], dtype=np.float64)])
        acc_w_m1[op] = np.concatenate([acc_w_m1[op], np.array(buf_w_m1[op], dtype=np.float64)])
    for pair in OP_PAIRS:
        acc_w_pp[pair] = np.concatenate([acc_w_pp[pair], np.array(buf_w_pp[pair], dtype=np.float64)])

    # ── Save checkpoint (numpy arrays — compact) ──────────────────────────────
    completed_files.append(lhe_file)
    ckpt = {
        "mll":             acc_mll,
        "rap":             acc_rap,
        "cstar":           acc_cstar,
        "w_SM":            acc_w_SM,
        "w_p1":            acc_w_p1,
        "w_m1":            acc_w_m1,
        "w_pp":            acc_w_pp,
        "completed_files": completed_files,
    }
    with open(CHECKPOINT_FILE, "wb") as f:
        pickle.dump(ckpt, f)
    size_mb = os.path.getsize(CHECKPOINT_FILE) / 1e6
    print(f"  Checkpoint saved ({size_mb:.1f} MB, {len(completed_files)}/{len(LHE_FILES)} files done)")

print(f"\nTotal: {len(acc_mll):,} events")

# ── Save final cache ───────────────────────────────────────────────────────────
cache = {
    'mll':   acc_mll,
    'rap':   acc_rap,
    'cstar': acc_cstar,
    'w_SM':  acc_w_SM,
    'w_p1':  acc_w_p1,
    'w_m1':  acc_w_m1,
    'w_pp':  acc_w_pp,
}

with open(CACHE_FILE, 'wb') as f:
    pickle.dump(cache, f)

size_mb = os.path.getsize(CACHE_FILE) / 1e6
print(f"Cache saved to {CACHE_FILE}  ({size_mb:.1f} MB)")

# ── Remove checkpoint once final cache is written ─────────────────────────────
if os.path.exists(CHECKPOINT_FILE):
    os.remove(CHECKPOINT_FILE)
    print("Checkpoint removed.")
