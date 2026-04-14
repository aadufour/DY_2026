"""
build_cache.py

Read all DYSMEFTMll LHE files once, extract kinematics, all SMEFT operator
weights, and QCD scale / PDF systematic weights. Save everything to a pickle.

Subsequent runs of build_datacard.py load the pickle instead of re-parsing
the LHE files.

Supports checkpointing: progress is saved after each LHE file so the script
can be interrupted and restarted without re-processing completed files.

Systematic weights extracted
─────────────────────────────
  Scale variations : all MUR/MUF combinations at the central PDF (303600).
                     Key pattern: MUR=X_MUF=Y_PDF=303600
  PDF replicas     : central scale (MUR=1, MUF=1) at each PDF replica member.
                     Key pattern: MUR=1_MUF=1_PDF=<replica_id>

These key patterns match standard MadGraph 2.9.x systematics output.
On the first event the script prints all detected keys so you can verify.

Usage:
    python3 build_cache.py
"""

import os
import pickle
import re
import warnings
from itertools import combinations

import numpy as np
import pylhe

# ── Config ─────────────────────────────────────────────────────────────────────

MLL_BIN_EDGES = [50, 120, 200, 400, 600, 800, 1000, 3000]

LHE_FILES = [
    f"/grid_mnt/data__data.polcms/cms/adufour/MG5/mg5amcnlo/DYSMEFTMll{lo}_{hi}/Events/run_01/unweighted_events.lhe.gz"
    for lo, hi in zip(MLL_BIN_EDGES[:-1], MLL_BIN_EDGES[1:])
]

MLL_LO = 50.0
MLL_HI = 3000.0

OPERATORS = [
    'cHDD', 'cHWB', 'cbWRe', 'cbBRe', 'cHj1', 'cHQ1',
    'cHj3', 'cHQ3', 'cHu',  'cHd',   'cHbq', 'cHl1',
    'cHl3', 'cHe',  'cll1', 'clj1',  'clj3', 'cQl1',
    'cQl3', 'ceu',  'ced',  'cbe',   'cje',  'cQe',
    'clu',  'cld',  'cbl',
]

OP_PAIRS = list(combinations(OPERATORS, 2))

CACHE_FILE      = "/grid_mnt/data__data.polcms/cms/adufour/MG5/mg5amcnlo/CACHE/lhe_cache.pkl"
CHECKPOINT_FILE = "/grid_mnt/data__data.polcms/cms/adufour/MG5/mg5amcnlo/CACHE/lhe_cache_checkpoint.pkl"

# ── Systematic key patterns ────────────────────────────────────────────────────
# MadGraph 2.9.x systematics module writes weight IDs like:
#   MUR=0.5_MUF=0.5_PDF=303600   ← scale variation, central PDF
#   MUR=1.0_MUF=1.0_PDF=303601   ← central scale, PDF replica 1
CENTRAL_PDF = 303600

# Scale: any MUR/MUF combo at central PDF
_RE_SCALE = re.compile(
    r'^MUR=[\d.]+_MUF=[\d.]+_PDF=%d$' % CENTRAL_PDF, re.IGNORECASE)

# PDF replicas: central scale, any PDF id != central
_RE_PDF = re.compile(
    r'^MUR=1[._]?0*_MUF=1[._]?0*_PDF=(\d+)$', re.IGNORECASE)

# ── Kinematic functions ────────────────────────────────────────────────────────

def mll(p1, p2):
    p = np.array(p1) + np.array(p2)
    return np.sqrt(max(p[3]**2 - sum(p[i]**2 for i in range(3)), 0.0))

def rap(p1, p2):
    p = np.array(p1) + np.array(p2)
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

# ── Load checkpoint ────────────────────────────────────────────────────────────

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
    acc_w_scale     = ckpt.get("w_scale", {})
    acc_w_pdf       = ckpt.get("w_pdf",   {})
    completed_files = ckpt["completed_files"]
    print(f"  {len(completed_files)} file(s) done, {len(acc_mll):,} events loaded\n")
else:
    acc_mll         = np.empty(0, dtype=np.float64)
    acc_rap         = np.empty(0, dtype=np.float64)
    acc_cstar       = np.empty(0, dtype=np.float64)
    acc_w_SM        = np.empty(0, dtype=np.float64)
    acc_w_p1        = {op:   np.empty(0, dtype=np.float64) for op   in OPERATORS}
    acc_w_m1        = {op:   np.empty(0, dtype=np.float64) for op   in OPERATORS}
    acc_w_pp        = {pair: np.empty(0, dtype=np.float64) for pair in OP_PAIRS}
    acc_w_scale     = {}   # {key: array}, populated on first event
    acc_w_pdf       = {}   # {key: array}, populated on first event
    completed_files = []

# ── Main loop ──────────────────────────────────────────────────────────────────

for lhe_file in LHE_FILES:
    if lhe_file in completed_files:
        print(f"Skipping (already done): {os.path.basename(lhe_file)}")
        continue

    print(f"\nReading {lhe_file}")
    if not os.path.exists(lhe_file):
        print("  WARNING: not found, skipping.")
        continue

    # ── Per-file buffers (Python lists — cheap to append) ──────────────────
    buf_mll   = []
    buf_rap   = []
    buf_cstar = []
    buf_w_SM  = []
    buf_w_p1  = {op:   [] for op   in OPERATORS}
    buf_w_m1  = {op:   [] for op   in OPERATORS}
    buf_w_pp  = {pair: [] for pair in OP_PAIRS}
    buf_w_scale = {}   # {key: list}, filled once SCALE_KEYS are known
    buf_w_pdf   = {}   # {key: list}, filled once PDF_KEYS are known

    # Sentinel flags
    pp_keys    = {}    # cross-weight key resolution (resolved once)
    SCALE_KEYS = list(acc_w_scale.keys()) if acc_w_scale else None
    PDF_KEYS   = list(acc_w_pdf.keys())   if acc_w_pdf   else None
    syst_init  = (SCALE_KEYS is not None)  # already resolved from checkpoint

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        events = pylhe.read_lhe_with_attributes(lhe_file)

        for i, event in enumerate(events):
            if (i + 1) % 5000 == 0:
                print(f"  {i + 1} events processed")

            # ── Select final-state charged leptons ──────────────────────────
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

            wkeys = event.weights

            # ── Resolve SMEFT cross-weight keys (once per file) ─────────────
            if not pp_keys:
                for op1, op2 in OP_PAIRS:
                    pp_keys[(op1, op2)] = (
                        f'{op1}_{op2}' if f'{op1}_{op2}' in wkeys
                        else f'{op2}_{op1}'
                    )

            # ── Detect and initialise systematic keys (once globally) ───────
            if not syst_init:
                all_keys   = list(wkeys.keys())
                SCALE_KEYS = [k for k in all_keys if _RE_SCALE.match(k)]
                PDF_KEYS   = [
                    k for k in all_keys
                    if _RE_PDF.match(k)
                    and int(_RE_PDF.match(k).group(1)) != CENTRAL_PDF
                ]
                print(f"  Detected {len(SCALE_KEYS)} scale variation keys")
                print(f"  Detected {len(PDF_KEYS)} PDF replica keys")
                if SCALE_KEYS:
                    print(f"    Scale keys : {SCALE_KEYS}")
                if PDF_KEYS:
                    print(f"    PDF keys (first 3): {PDF_KEYS[:3]} ...")
                if not SCALE_KEYS:
                    print("  WARNING: no scale keys matched — check weight IDs in LHE")
                if not PDF_KEYS:
                    print("  WARNING: no PDF replica keys matched — check weight IDs in LHE")
                buf_w_scale = {k: [] for k in SCALE_KEYS}
                buf_w_pdf   = {k: [] for k in PDF_KEYS}
                syst_init   = True

            # ── Fill buffers ────────────────────────────────────────────────
            buf_mll.append(m)
            buf_rap.append(y)
            buf_cstar.append(cs)
            buf_w_SM.append(wkeys['SM'])

            for op in OPERATORS:
                buf_w_p1[op].append(wkeys[op])
                buf_w_m1[op].append(wkeys[f'minus{op}'])

            for pair in OP_PAIRS:
                buf_w_pp[pair].append(wkeys[pp_keys[pair]])

            for k in SCALE_KEYS:
                buf_w_scale[k].append(wkeys.get(k, wkeys['SM']))

            for k in PDF_KEYS:
                buf_w_pdf[k].append(wkeys.get(k, wkeys['SM']))

    n_kept = len(buf_mll)
    print(f"  {n_kept} events kept from this file")

    # ── Convert and concatenate into accumulators ───────────────────────────
    acc_mll   = np.concatenate([acc_mll,   np.array(buf_mll,   dtype=np.float64)])
    acc_rap   = np.concatenate([acc_rap,   np.array(buf_rap,   dtype=np.float64)])
    acc_cstar = np.concatenate([acc_cstar, np.array(buf_cstar, dtype=np.float64)])
    acc_w_SM  = np.concatenate([acc_w_SM,  np.array(buf_w_SM,  dtype=np.float64)])

    for op in OPERATORS:
        acc_w_p1[op] = np.concatenate([acc_w_p1[op], np.array(buf_w_p1[op], dtype=np.float64)])
        acc_w_m1[op] = np.concatenate([acc_w_m1[op], np.array(buf_w_m1[op], dtype=np.float64)])

    for pair in OP_PAIRS:
        acc_w_pp[pair] = np.concatenate([acc_w_pp[pair], np.array(buf_w_pp[pair], dtype=np.float64)])

    for k in SCALE_KEYS:
        prev = acc_w_scale.get(k, np.empty(0, dtype=np.float64))
        acc_w_scale[k] = np.concatenate([prev, np.array(buf_w_scale[k], dtype=np.float64)])

    for k in PDF_KEYS:
        prev = acc_w_pdf.get(k, np.empty(0, dtype=np.float64))
        acc_w_pdf[k] = np.concatenate([prev, np.array(buf_w_pdf[k], dtype=np.float64)])

    # ── Save checkpoint ─────────────────────────────────────────────────────
    completed_files.append(lhe_file)
    ckpt = {
        "mll":             acc_mll,
        "rap":             acc_rap,
        "cstar":           acc_cstar,
        "w_SM":            acc_w_SM,
        "w_p1":            acc_w_p1,
        "w_m1":            acc_w_m1,
        "w_pp":            acc_w_pp,
        "w_scale":         acc_w_scale,
        "w_pdf":           acc_w_pdf,
        "completed_files": completed_files,
    }
    with open(CHECKPOINT_FILE, "wb") as f:
        pickle.dump(ckpt, f)
    size_mb = os.path.getsize(CHECKPOINT_FILE) / 1e6
    print(f"  Checkpoint saved ({size_mb:.1f} MB, {len(completed_files)}/{len(LHE_FILES)} files done)")

print(f"\nTotal: {len(acc_mll):,} events")
print(f"  Scale variations : {len(acc_w_scale)} keys")
print(f"  PDF replicas     : {len(acc_w_pdf)} keys")

# ── Save final cache ───────────────────────────────────────────────────────────
cache = {
    'mll':     acc_mll,
    'rap':     acc_rap,
    'cstar':   acc_cstar,
    'w_SM':    acc_w_SM,
    'w_p1':    acc_w_p1,
    'w_m1':    acc_w_m1,
    'w_pp':    acc_w_pp,
    'w_scale': acc_w_scale,
    'w_pdf':   acc_w_pdf,
}

with open(CACHE_FILE, 'wb') as f:
    pickle.dump(cache, f)

size_mb = os.path.getsize(CACHE_FILE) / 1e6
print(f"Cache saved → {CACHE_FILE}  ({size_mb:.1f} MB)")

if os.path.exists(CHECKPOINT_FILE):
    os.remove(CHECKPOINT_FILE)
    print("Checkpoint removed.")
