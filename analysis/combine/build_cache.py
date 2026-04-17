"""
build_cache.py

Read all DYSMEFTMll LHE files once, extract kinematics, all SMEFT operator
weights, and QCD scale / PDF systematic weights. Save everything to a pickle.

Systematic weight IDs are resolved by parsing the LHE header once:
  Scale : no DYN_SCALE, PDF=303600, MUR/MUF in {0.5,1,2}x{0.5,1,2},
          excluding the two anti-correlated extremes (0.5,2) and (2,0.5)
          and the central point (1,1) — standard 6-point envelope.
  PDF   : MUR=1, MUF=1, PDF != 303600, no DYN_SCALE — replicas 1..100.
  Central (nominal) : MUR=1, MUF=1, PDF=303600, no DYN_SCALE — used as
          the base weight for ratio computation in build_datacard.py.

Usage:
    python3 build_cache.py
    python3 build_cache.py --nodoubles   # skip operator-pair weights (reweight card without pairs)
"""

import argparse
import gzip
import os
import pickle
import re
import warnings
from itertools import combinations

import numpy as np
import pylhe

# ---- Argument parsing --------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument('--nodoubles', action='store_true',
                    help='Skip operator-pair (quadratic cross) weights — use when reweight card has no pairs')
args = parser.parse_args()

SKIP_PAIRS = args.nodoubles

# ---- Config ------------------------------------------------------------------------------------------------------------------------------------------

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

OP_PAIRS = [] if SKIP_PAIRS else list(combinations(OPERATORS, 2))

CACHE_FILE      = "/grid_mnt/data__data.polcms/cms/adufour/MG5/mg5amcnlo/CACHE/lhe_cache.pkl"
CHECKPOINT_FILE = "/grid_mnt/data__data.polcms/cms/adufour/MG5/mg5amcnlo/CACHE/lhe_cache_checkpoint.pkl"

CENTRAL_PDF = 303600

# ---- Parse LHE header to resolve weight IDs -----------------------

def parse_weight_ids(lhe_file):
    """
    Read the LHE header and return:
      scale_ids   : list of string IDs for the 6-point scale envelope
      pdf_ids     : list of string IDs for PDF replicas
      central_id  : string ID of MUR=1, MUF=1, PDF=303600 (no DYN_SCALE)
    """
    open_fn = gzip.open if lhe_file.endswith('.gz') else open
    header  = []
    with open_fn(lhe_file, 'rt') as f:
        for line in f:
            header.append(line)
            if '</header>' in line or '<event' in line:
                break
    header_text = ''.join(header)

    # Parse all <weight id="X" ...> tags
    pattern = re.compile(
        r'<weight\s+id=["\'](\d+)["\'][^>]*>', re.IGNORECASE)
    attr_re = {
        'MUR':       re.compile(r'MUR=["\']?([\d.]+)',       re.IGNORECASE),
        'MUF':       re.compile(r'MUF=["\']?([\d.]+)',       re.IGNORECASE),
        'PDF':       re.compile(r'PDF=["\']?(\d+)',          re.IGNORECASE),
        'DYN_SCALE': re.compile(r'DYN_SCALE=["\']?([\d.]+)', re.IGNORECASE),
    }

    weights = {}
    for line in header_text.splitlines():
        m = pattern.search(line)
        if not m:
            continue
        wid = m.group(1)
        info = {}
        for attr, rx in attr_re.items():
            am = rx.search(line)
            info[attr] = am.group(1) if am else None
        weights[wid] = info

    scale_ids  = []
    pdf_ids    = []
    central_id = None

    for wid, info in weights.items():
        mur = float(info['MUR']) if info['MUR'] else None
        muf = float(info['MUF']) if info['MUF'] else None
        pdf = int(info['PDF'])   if info['PDF'] else None
        dyn = info['DYN_SCALE']

        if dyn is not None:
            continue   # skip dynamic scale choices

        # Central point
        if mur == 1.0 and muf == 1.0 and pdf == CENTRAL_PDF:
            central_id = wid
            continue

        # Scale variations: central PDF, exclude anti-correlated extremes
        if pdf == CENTRAL_PDF and mur is not None and muf is not None:
            if not (mur == 0.5 and muf == 2.0) and not (mur == 2.0 and muf == 0.5):
                scale_ids.append(wid)

        # PDF replicas: central scale, non-central PDF member
        if mur == 1.0 and muf == 1.0 and pdf is not None and pdf != CENTRAL_PDF:
            pdf_ids.append(wid)

    return scale_ids, pdf_ids, central_id

# ---- Kinematic functions ------------------------------------

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

# ---- Load checkpoint -----------------------------------------

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
    acc_w_scale       = ckpt.get("w_scale",       {})
    acc_w_pdf         = ckpt.get("w_pdf",         {})
    acc_w_pdf_central = ckpt.get("w_pdf_central", np.empty(0, dtype=np.float64))
    completed_files   = ckpt["completed_files"]
    print(f"  {len(completed_files)} file(s) done, {len(acc_mll):,} events loaded\n")
else:
    acc_mll           = np.empty(0, dtype=np.float64)
    acc_rap           = np.empty(0, dtype=np.float64)
    acc_cstar         = np.empty(0, dtype=np.float64)
    acc_w_SM          = np.empty(0, dtype=np.float64)
    acc_w_p1          = {op:   np.empty(0, dtype=np.float64) for op   in OPERATORS}
    acc_w_m1          = {op:   np.empty(0, dtype=np.float64) for op   in OPERATORS}
    acc_w_pp          = {pair: np.empty(0, dtype=np.float64) for pair in OP_PAIRS}
    acc_w_scale       = {}
    acc_w_pdf         = {}
    acc_w_pdf_central = np.empty(0, dtype=np.float64)
    completed_files   = []

# ---- Main loop ---------------------------------------------------

for lhe_file in LHE_FILES:
    if lhe_file in completed_files:
        print(f"Skipping (already done): {os.path.basename(lhe_file)}")
        continue

    print(f"\nReading {lhe_file}")
    if not os.path.exists(lhe_file):
        print("  WARNING: not found, skipping.")
        continue

    # Resolve systematic weight IDs from header
    scale_ids, pdf_ids, central_id = parse_weight_ids(lhe_file)
    print(f"  Scale IDs ({len(scale_ids)}): {scale_ids}")
    print(f"  PDF IDs   ({len(pdf_ids)}):  {pdf_ids[:5]} ...")
    print(f"  Central ID: {central_id}")

    # Per-file buffers
    buf_mll     = []
    buf_rap     = []
    buf_cstar   = []
    buf_w_SM    = []
    buf_w_p1    = {op:   [] for op   in OPERATORS}
    buf_w_m1    = {op:   [] for op   in OPERATORS}
    buf_w_pp    = {pair: [] for pair in OP_PAIRS}
    buf_w_scale       = {k: [] for k in scale_ids}
    buf_w_pdf         = {k: [] for k in pdf_ids}
    buf_w_pdf_central = []

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

            wkeys = event.weights

            if not pp_keys and not SKIP_PAIRS:
                for op1, op2 in OP_PAIRS:
                    pp_keys[(op1, op2)] = (
                        f'{op1}_{op2}' if f'{op1}_{op2}' in wkeys
                        else f'{op2}_{op1}'
                    )

            buf_mll.append(m)
            buf_rap.append(y)
            buf_cstar.append(cs)
            buf_w_SM.append(wkeys['SM'])

            for op in OPERATORS:
                buf_w_p1[op].append(wkeys[op])
                buf_w_m1[op].append(wkeys[f'minus{op}'])

            if not SKIP_PAIRS:
                for pair in OP_PAIRS:
                    buf_w_pp[pair].append(wkeys[pp_keys[pair]])

            for k in scale_ids:
                buf_w_scale[k].append(wkeys.get(k, wkeys['SM']))

            if central_id is not None:
                buf_w_pdf_central.append(wkeys.get(central_id, wkeys['SM']))

            for k in pdf_ids:
                buf_w_pdf[k].append(wkeys.get(k, wkeys['SM']))

    n_kept = len(buf_mll)
    print(f"  {n_kept} events kept")

    # Convert and accumulate
    acc_mll   = np.concatenate([acc_mll,   np.array(buf_mll,   dtype=np.float64)])
    acc_rap   = np.concatenate([acc_rap,   np.array(buf_rap,   dtype=np.float64)])
    acc_cstar = np.concatenate([acc_cstar, np.array(buf_cstar, dtype=np.float64)])
    acc_w_SM  = np.concatenate([acc_w_SM,  np.array(buf_w_SM,  dtype=np.float64)])

    for op in OPERATORS:
        acc_w_p1[op] = np.concatenate([acc_w_p1[op], np.array(buf_w_p1[op], dtype=np.float64)])
        acc_w_m1[op] = np.concatenate([acc_w_m1[op], np.array(buf_w_m1[op], dtype=np.float64)])

    for pair in OP_PAIRS:
        acc_w_pp[pair] = np.concatenate([acc_w_pp[pair], np.array(buf_w_pp[pair], dtype=np.float64)])

    for k in scale_ids:
        prev = acc_w_scale.get(k, np.empty(0, dtype=np.float64))
        acc_w_scale[k] = np.concatenate([prev, np.array(buf_w_scale[k], dtype=np.float64)])

    for k in pdf_ids:
        prev = acc_w_pdf.get(k, np.empty(0, dtype=np.float64))
        acc_w_pdf[k] = np.concatenate([prev, np.array(buf_w_pdf[k], dtype=np.float64)])

    acc_w_pdf_central = np.concatenate([acc_w_pdf_central, np.array(buf_w_pdf_central, dtype=np.float64)])

    # Save checkpoint
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
        "w_pdf_central":   acc_w_pdf_central,
        "completed_files": completed_files,
    }
    with open(CHECKPOINT_FILE, "wb") as f:
        pickle.dump(ckpt, f)
    size_mb = os.path.getsize(CHECKPOINT_FILE) / 1e6
    print(f"  Checkpoint saved ({size_mb:.1f} MB, {len(completed_files)}/{len(LHE_FILES)} files done)")

print(f"\nTotal: {len(acc_mll):,} events")
print(f"  Scale variations : {len(acc_w_scale)} keys")
print(f"  PDF replicas     : {len(acc_w_pdf)} keys")

# Save final cache
cache = {
    'mll':     acc_mll,
    'rap':     acc_rap,
    'cstar':   acc_cstar,
    'w_SM':    acc_w_SM,
    'w_p1':    acc_w_p1,
    'w_m1':    acc_w_m1,
    'w_pp':    acc_w_pp,
    'w_scale':       acc_w_scale,
    'w_pdf':         acc_w_pdf,
    'w_pdf_central': acc_w_pdf_central,
}

with open(CACHE_FILE, 'wb') as f:
    pickle.dump(cache, f)

size_mb = os.path.getsize(CACHE_FILE) / 1e6
print(f"Cache saved -> {CACHE_FILE}  ({size_mb:.1f} MB)")

if os.path.exists(CHECKPOINT_FILE):
    os.remove(CHECKPOINT_FILE)
    print("Checkpoint removed.")
