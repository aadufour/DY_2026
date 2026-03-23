"""
rw_build_cache.py

Used for rw_triple_diff.py when merging all DY_all_min_max LHE files.

Read all LHE files once, extract kinematics and all operator weights,
save to a pickle. Subsequent runs of rw_triple_diff.py load the pickle
instead of re-parsing the LHE for faster loading.

Usage:
    python3 rw_build_cache.py
"""

import os
import pickle
import warnings
import numpy as np
import pylhe

# ── Config (must match rw_triple_diff.py) ─────────────────────────────────────
MLL_BIN_EDGES = [50, 120, 200, 400, 600, 800, 1000, 3000]
LHE_FILES = [
    f"/Users/albertodufour/MG5_2_9_18/mg5amcnlo/DY_all_{lo}_{hi}/myLHE/unweighted_events.lhe"
    for lo, hi in zip(MLL_BIN_EDGES[:-1], MLL_BIN_EDGES[1:])
]

MLL_LO = 50.0
MLL_HI = 3000.0

OPERATORS = [
    'cHDD', 'cHWB', 'cHj1',
    'cHj3', 'cHu',  'cHd',
    'cHl1', 'cHl3', 'cHe',
    'cll1', 'clj1', 'clj3',
    'ceu',  'ced',  'cje',
    'clu',  'cld',
]

CACHE_FILE = "/Users/albertodufour/code/DY2026/triple_diff_test/lhe_cache.pkl"

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

# ── Collect ────────────────────────────────────────────────────────────────────
mll_vals   = []
rap_vals   = []
cstar_vals = []
w_SM       = []
w_p1       = {op: [] for op in OPERATORS}
w_m1       = {op: [] for op in OPERATORS}

for lhe_file in LHE_FILES:
    print(f"\nReading {lhe_file}")
    if not os.path.exists(lhe_file):
        print("  WARNING: not found, skipping.")
        continue

    n_before = len(mll_vals)
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

            mll_vals.append(m)
            rap_vals.append(y)
            cstar_vals.append(cs)
            w_SM.append(event.weights['SM'])
            for op in OPERATORS:
                w_p1[op].append(event.weights[op])
                w_m1[op].append(event.weights[f'minus{op}'])

    print(f"  {len(mll_vals) - n_before} events kept from this file")

print(f"\nTotal: {len(mll_vals)} events")

# ── Save ───────────────────────────────────────────────────────────────────────
cache = {
    'mll':   np.array(mll_vals),
    'rap':   np.array(rap_vals),
    'cstar': np.array(cstar_vals),
    'w_SM':  np.array(w_SM),
    'w_p1':  {op: np.array(v) for op, v in w_p1.items()},
    'w_m1':  {op: np.array(v) for op, v in w_m1.items()},
}

with open(CACHE_FILE, 'wb') as f:
    pickle.dump(cache, f)

size_mb = os.path.getsize(CACHE_FILE) / 1e6
print(f"Cache saved to {CACHE_FILE}  ({size_mb:.1f} MB)")
