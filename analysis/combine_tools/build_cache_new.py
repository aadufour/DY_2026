#!/usr/bin/env python3
"""
build_cache_new.py
==================
Simplified version of build_cache_syst.py for LO generation with b quark in proton.

Changes vs build_cache_syst.py:
  - Scale variations: MUF-only envelope (MUR fixed to 1).
    At LO there are no virtual/real corrections so MUR does not appear
    in the matrix element — varying it has no effect.
  - PDF: only NNPDF31_nnlo_as_0118_mc_hessian_pdfas (325300, 5-flavour).
    The 4-flavour set (325500) is dropped because the events were generated
    with b quark in the proton (5-flavour scheme).

Usage:
    python3 build_cache_new.py
    python3 build_cache_new.py --nodoubles
    python3 build_cache_new.py --nevents 20000
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
                    help='Skip operator-pair (quadratic cross) weights')
parser.add_argument('--nevents', type=int, default=None,
                    help='Maximum events to read per LHE file (default: all)')
parser.add_argument('--eta-max', type=float, default=None,
                    help='Maximum |eta| for both leptons (default: no cut). '
                         'CMS muon acceptance: 2.4')
parser.add_argument('--pt-lead', type=float, default=None,
                    help='Minimum pT of leading lepton in GeV (default: no cut). '
                         'Typical CMS value: 25 GeV')
parser.add_argument('--pt-sub', type=float, default=None,
                    help='Minimum pT of subleading lepton in GeV (default: no cut). '
                         'Typical CMS value: 10 GeV')
args = parser.parse_args()

SKIP_PAIRS = args.nodoubles
MAX_EVENTS = args.nevents
ETA_MAX    = args.eta_max
PT_LEAD    = args.pt_lead
PT_SUB     = args.pt_sub

do_fiducial = any(x is not None for x in [ETA_MAX, PT_LEAD, PT_SUB])
if do_fiducial:
    print("Fiducial cuts:")
    if ETA_MAX  is not None: print(f"  |eta| < {ETA_MAX}")
    if PT_LEAD  is not None: print(f"  pT(lead) > {PT_LEAD} GeV")
    if PT_SUB   is not None: print(f"  pT(sub)  > {PT_SUB} GeV")
else:
    print("No fiducial cuts applied (parton-level inclusive)")

# ---- Config ------------------------------------------------------------------

MLL_BIN_EDGES = [50, 120, 200, 400, 600, 800, 1000, 3000]

LHE_FILES = [
    f"/grid_mnt/data__data.polcms/cms/adufour/LHE/SYST_slc7/DYSMEFTMll{lo}_{hi}/unweighted_events.lhe"
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

CACHE_FILE      = "/grid_mnt/data__data.polcms/cms/adufour/LHE/SYST_slc7/CACHE/lhe_cache_new.pkl"
CHECKPOINT_FILE = "/grid_mnt/data__data.polcms/cms/adufour/LHE/SYST_slc7/CACHE/lhe_cache_new_checkpoint.pkl"

CENTRAL_PDF = 325300   # 5-flavour central PDF

# ---- Parse LHE header --------------------------------------------------------

def parse_weight_ids(lhe_file):
    """
    Return:
      scale_ids       list of IDs for MUF-only envelope (MUR=1, MUF!=1,
                      PDF=CENTRAL_PDF, no DYN_SCALE)
      pdf_325300_ids  ordered list of 103 IDs (MemberID 0-102, MUR=MUF=1)
      central_id      ID of MUR=1, MUF=1, PDF=CENTRAL_PDF standalone nominal
    """
    open_fn = gzip.open if lhe_file.endswith('.gz') else open
    header = []
    with open_fn(lhe_file, 'rt') as f:
        for line in f:
            header.append(line)
            if '</header>' in line or '<event' in line:
                break
    header_text = ''.join(header)

    pattern = re.compile(r'<weight\s+id=["\'](\d+)["\'][^>]*>', re.IGNORECASE)
    attr_re = {
        'MUR':       re.compile(r'MUR=["\']?([\d.]+)',       re.IGNORECASE),
        'MUF':       re.compile(r'MUF=["\']?([\d.]+)',       re.IGNORECASE),
        'PDF':       re.compile(r'PDF=["\']?(\d+)',          re.IGNORECASE),
        'DYN_SCALE': re.compile(r'DYN_SCALE=["\']?([\d.]+)', re.IGNORECASE),
        'ALPSFACT':  re.compile(r'ALPSFACT=["\']?([\d.]+)',  re.IGNORECASE),
    }

    weights = {}
    order   = []
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
        order.append(wid)

    scale_ids      = []
    pdf_325300_ids = []
    central_id     = None

    for wid in order:
        info = weights[wid]
        mur = float(info['MUR']) if info['MUR'] else None
        muf = float(info['MUF']) if info['MUF'] else None
        pdf = int(info['PDF'])   if info['PDF'] else None
        dyn = info['DYN_SCALE']

        if dyn is not None or info['ALPSFACT'] is not None:
            continue

        # Standalone nominal
        if mur == 1.0 and muf == 1.0 and pdf == CENTRAL_PDF and central_id is None:
            central_id = wid
            continue

        # MUF-only scale envelope: MUR fixed to 1, MUF varies, central PDF
        # At LO only MUF matters — MUR has no effect on the matrix element.
        if pdf == CENTRAL_PDF and mur == 1.0 and muf is not None and muf != 1.0:
            scale_ids.append(wid)

        # PDF 325300 members (5-flavour, 103 members)
        if mur == 1.0 and muf == 1.0 and pdf is not None:
            if CENTRAL_PDF <= pdf <= CENTRAL_PDF + 102:
                pdf_325300_ids.append(wid)

    return scale_ids, pdf_325300_ids, central_id

# ---- Kinematic functions -----------------------------------------------------

def mll(p1, p2):
    p = np.array(p1) + np.array(p2)
    return np.sqrt(max(p[3]**2 - sum(p[i]**2 for i in range(3)), 0.0))

def rap(p1, p2):
    p = np.array(p1) + np.array(p2)
    E, pz = p[3], p[2]
    return np.abs(0.5 * np.log((E + pz) / (E - pz)))

def eta(p):
    """Pseudorapidity of a single particle 4-vector [px, py, pz, E]."""
    px, py, pz, e = p
    pmag = np.sqrt(px**2 + py**2 + pz**2)
    if pmag == abs(pz):   # collinear with beam — infinite eta
        return np.inf
    return 0.5 * np.log((pmag + pz) / (pmag - pz))

def pt(p):
    """Transverse momentum of a single particle 4-vector [px, py, pz, E]."""
    px, py, pz, e = p
    return np.sqrt(px**2 + py**2)

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

# ---- Load checkpoint ---------------------------------------------------------

if os.path.exists(CHECKPOINT_FILE):
    print(f"Resuming from checkpoint: {CHECKPOINT_FILE}")
    with open(CHECKPOINT_FILE, "rb") as f:
        ckpt = pickle.load(f)
    acc_mll           = ckpt["mll"]
    acc_rap           = ckpt["rap"]
    acc_cstar         = ckpt["cstar"]
    acc_w_SM          = ckpt["w_SM"]
    acc_xwgt          = ckpt.get("xwgt", np.empty(0, dtype=np.float64))
    acc_w_p1          = ckpt["w_p1"]
    acc_w_m1          = ckpt["w_m1"]
    acc_w_pp          = ckpt["w_pp"]
    acc_w_scale       = ckpt.get("w_scale",       {})
    acc_w_pdf_central = ckpt.get("w_pdf_central", np.empty(0, dtype=np.float64))
    acc_pdf_325300    = ckpt.get("pdf_325300",    np.empty((0, 103), dtype=np.float64))
    completed_files   = ckpt["completed_files"]
    print(f"  {len(completed_files)} file(s) done, {len(acc_mll):,} events loaded\n")
else:
    acc_mll           = np.empty(0, dtype=np.float64)
    acc_rap           = np.empty(0, dtype=np.float64)
    acc_cstar         = np.empty(0, dtype=np.float64)
    acc_w_SM          = np.empty(0, dtype=np.float64)
    acc_xwgt          = np.empty(0, dtype=np.float64)
    acc_w_p1          = {op:   np.empty(0, dtype=np.float64) for op   in OPERATORS}
    acc_w_m1          = {op:   np.empty(0, dtype=np.float64) for op   in OPERATORS}
    acc_w_pp          = {pair: np.empty(0, dtype=np.float64) for pair in OP_PAIRS}
    acc_w_scale       = {}
    acc_w_pdf_central = np.empty(0, dtype=np.float64)
    acc_pdf_325300    = np.empty((0, 103), dtype=np.float64)
    completed_files   = []

os.makedirs(os.path.dirname(CACHE_FILE), exist_ok=True)

# ---- Main loop ---------------------------------------------------------------

for lhe_file in LHE_FILES:
    if lhe_file in completed_files:
        print(f"Skipping (already done): {os.path.basename(lhe_file)}")
        continue

    print(f"\nReading {lhe_file}")
    if not os.path.exists(lhe_file):
        print("  WARNING: not found, skipping.")
        continue

    scale_ids, pdf_325300_ids, central_id = parse_weight_ids(lhe_file)
    n_pdf300 = len(pdf_325300_ids)
    print(f"  MUF-only scale IDs ({len(scale_ids)}): {scale_ids}")
    print(f"  PDF 325300 ({n_pdf300} members): IDs {pdf_325300_ids[0]}–{pdf_325300_ids[-1]}")
    print(f"  Central ID: {central_id}")

    buf_mll           = []
    buf_rap           = []
    buf_cstar         = []
    buf_w_SM          = []
    buf_xwgt          = []
    buf_w_p1          = {op:   [] for op   in OPERATORS}
    buf_w_m1          = {op:   [] for op   in OPERATORS}
    buf_w_pp          = {pair: [] for pair in OP_PAIRS}
    buf_w_scale       = {k:    [] for k    in scale_ids}
    buf_w_pdf_central = []
    buf_pdf_325300    = []

    pp_keys = {}

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        events = pylhe.read_lhe_with_attributes(lhe_file)

        for i, event in enumerate(events):
            if MAX_EVENTS is not None and i >= MAX_EVENTS:
                break
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

            # ---- Fiducial cuts (optional) ------------------------------------
            if do_fiducial:
                pt_lm  = pt(v_lm);   pt_lp  = pt(v_lp)
                eta_lm = eta(v_lm);  eta_lp = eta(v_lp)
                pt_lead_val = max(pt_lm, pt_lp)
                pt_sub_val  = min(pt_lm, pt_lp)
                if ETA_MAX is not None:
                    if abs(eta_lm) > ETA_MAX or abs(eta_lp) > ETA_MAX:
                        continue
                if PT_LEAD is not None:
                    if pt_lead_val < PT_LEAD:
                        continue
                if PT_SUB is not None:
                    if pt_sub_val < PT_SUB:
                        continue
            # ------------------------------------------------------------------

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
            buf_xwgt.append(event.eventinfo.weight)

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

            buf_pdf_325300.append([wkeys.get(k, wkeys['SM']) for k in pdf_325300_ids])

    n_kept = len(buf_mll)
    print(f"  {n_kept} events kept")

    acc_mll   = np.concatenate([acc_mll,   np.array(buf_mll,   dtype=np.float64)])
    acc_rap   = np.concatenate([acc_rap,   np.array(buf_rap,   dtype=np.float64)])
    acc_cstar = np.concatenate([acc_cstar, np.array(buf_cstar, dtype=np.float64)])
    acc_w_SM  = np.concatenate([acc_w_SM,  np.array(buf_w_SM,  dtype=np.float64)])
    acc_xwgt  = np.concatenate([acc_xwgt,  np.array(buf_xwgt,  dtype=np.float64)])

    for op in OPERATORS:
        acc_w_p1[op] = np.concatenate([acc_w_p1[op], np.array(buf_w_p1[op], dtype=np.float64)])
        acc_w_m1[op] = np.concatenate([acc_w_m1[op], np.array(buf_w_m1[op], dtype=np.float64)])

    for pair in OP_PAIRS:
        acc_w_pp[pair] = np.concatenate([acc_w_pp[pair], np.array(buf_w_pp[pair], dtype=np.float64)])

    for k in scale_ids:
        prev = acc_w_scale.get(k, np.empty(0, dtype=np.float64))
        acc_w_scale[k] = np.concatenate([prev, np.array(buf_w_scale[k], dtype=np.float64)])

    acc_w_pdf_central = np.concatenate([acc_w_pdf_central, np.array(buf_w_pdf_central, dtype=np.float64)])

    arr_325300 = np.array(buf_pdf_325300, dtype=np.float64).reshape(-1, n_pdf300)
    acc_pdf_325300 = np.concatenate([acc_pdf_325300, arr_325300], axis=0)

    completed_files.append(lhe_file)
    ckpt = {
        "mll":             acc_mll,
        "rap":             acc_rap,
        "cstar":           acc_cstar,
        "w_SM":            acc_w_SM,
        "xwgt":            acc_xwgt,
        "w_p1":            acc_w_p1,
        "w_m1":            acc_w_m1,
        "w_pp":            acc_w_pp,
        "w_scale":         acc_w_scale,
        "w_pdf_central":   acc_w_pdf_central,
        "pdf_325300":      acc_pdf_325300,
        "completed_files": completed_files,
    }
    with open(CHECKPOINT_FILE, "wb") as f:
        pickle.dump(ckpt, f)
    size_mb = os.path.getsize(CHECKPOINT_FILE) / 1e6
    print(f"  Checkpoint saved ({size_mb:.1f} MB, {len(completed_files)}/{len(LHE_FILES)} files done)")

print(f"\nTotal: {len(acc_mll):,} events")
print(f"  MUF scale variations : {len(acc_w_scale)} keys")
print(f"  pdf_325300 shape     : {acc_pdf_325300.shape}")

# ---- Save final cache --------------------------------------------------------

cache = {
    'mll':           acc_mll,
    'rap':           acc_rap,
    'cstar':         acc_cstar,
    'w_SM':          acc_w_SM,
    'xwgt':          acc_xwgt,
    'w_p1':          acc_w_p1,
    'w_m1':          acc_w_m1,
    'w_pp':          acc_w_pp,
    'w_scale':       acc_w_scale,
    'w_pdf_central': acc_w_pdf_central,
    'pdf_325300':    acc_pdf_325300,   # [n_events, 103]
    # Selection applied when building this cache
    'cuts': {
        'eta_max':  ETA_MAX,
        'pt_lead':  PT_LEAD,
        'pt_sub':   PT_SUB,
        'mll_lo':   MLL_LO,
        'mll_hi':   MLL_HI,
    },
}

with open(CACHE_FILE, 'wb') as f:
    pickle.dump(cache, f)

size_mb = os.path.getsize(CACHE_FILE) / 1e6
print(f"Cache saved -> {CACHE_FILE}  ({size_mb:.1f} MB)")

if os.path.exists(CHECKPOINT_FILE):
    os.remove(CHECKPOINT_FILE)
    print("Checkpoint removed.")
