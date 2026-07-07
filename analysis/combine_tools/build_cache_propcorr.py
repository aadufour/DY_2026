#!/usr/bin/env python3
"""
build_cache_propcorr.py
========================
Variant of build_cache_parallel.py for the propagator-corrected production,
where each mll bin has 10 separate LHE files (from parallel condor jobs)
instead of one merged file.

Each worker reads all of its bin's files in one pass and concatenates the
events internally, rather than requiring a pre-merged unweighted_events.lhe.
Header-derived weight IDs (scale/PDF) are parsed once from the first file
per bin, since all files of a bin share the same run_card/reweight_card.

Usage:
    python3 build_cache_propcorr.py
    python3 build_cache_propcorr.py --eta-max 2.4 --pt-lead 25 --pt-sub 10
    python3 build_cache_propcorr.py --nodoubles --nevents 5000
    python3 build_cache_propcorr.py --workers 4   # limit parallelism
"""

import argparse
import gzip
import glob
import os
import pickle
import re
import warnings
from itertools import combinations
from multiprocessing import Pool

import numpy as np
import pylhe

# ---- Argument parsing --------------------------------------------------------

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--nodoubles', action='store_true',
                    help='Skip operator-pair (quadratic cross) weights')
parser.add_argument('--nevents', type=int, default=None,
                    help='Maximum events to read per bin (default: all)')
parser.add_argument('--eta-max', type=float, default=None,
                    help='Maximum |eta| for both leptons (CMS muon acceptance: 2.4)')
parser.add_argument('--pt-lead', type=float, default=None,
                    help='Minimum pT of leading lepton in GeV (typical CMS: 25)')
parser.add_argument('--pt-sub', type=float, default=None,
                    help='Minimum pT of subleading lepton in GeV (typical CMS: 10)')
parser.add_argument('--workers', type=int, default=None,
                    help='Number of parallel workers (default: one per bin)')
args = parser.parse_args()

SKIP_PAIRS = args.nodoubles
MAX_EVENTS = args.nevents
ETA_MAX    = args.eta_max
PT_LEAD    = args.pt_lead
PT_SUB     = args.pt_sub
do_fiducial = any(x is not None for x in [ETA_MAX, PT_LEAD, PT_SUB])

# ---- Config ------------------------------------------------------------------

MLL_BIN_EDGES = [50, 120, 200, 400, 600, 800, 1000, 3000]

LHE_BASE = "/grid_mnt/data__data.polcms/cms/adufour/LHE/propcorr"

# one (tag, [file list]) pair per bin
LHE_BINS = [
    (
        f"DYSMEFTMll{lo}_{hi}_propcorr",
        sorted(glob.glob(f"{LHE_BASE}/DYSMEFTMll{lo}_{hi}_propcorr/lhe_output/output_*.lhe"))
    )
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

CACHE_DIR       = "/grid_mnt/data__data.polcms/cms/adufour/LHE/propcorr/CACHE"
CACHE_FILE      = os.path.join(CACHE_DIR, "lhe_cache_propcorr.pkl")
PER_FILE_DIR    = os.path.join(CACHE_DIR, "per_file")

CENTRAL_PDF = 325300

# ---- Header parser -----------------------------------------------------------

def parse_weight_ids(lhe_file):
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
        'MUR':       re.compile(r'MUR=["\']?([\d.]+)',        re.IGNORECASE),
        'MUF':       re.compile(r'MUF=["\']?([\d.]+)',        re.IGNORECASE),
        'PDF':       re.compile(r'PDF=["\']?(\d+)',           re.IGNORECASE),
        'DYN_SCALE': re.compile(r'DYN_SCALE=["\']?([\d.]+)', re.IGNORECASE),
        'ALPSFACT':  re.compile(r'ALPSFACT=["\']?([\d.]+)',   re.IGNORECASE),
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
        if info['DYN_SCALE'] is not None or info['ALPSFACT'] is not None:
            continue
        if mur == 1.0 and muf == 1.0 and pdf == CENTRAL_PDF and central_id is None:
            central_id = wid
            continue
        # MUF-only scale envelope (MUR fixed to 1 — LO, MUR has no effect)
        if pdf == CENTRAL_PDF and mur == 1.0 and muf is not None and muf != 1.0:
            scale_ids.append(wid)
        # 5-flavour PDF members
        if mur == 1.0 and muf == 1.0 and pdf is not None:
            if CENTRAL_PDF <= pdf <= CENTRAL_PDF + 102:
                pdf_325300_ids.append(wid)

    return scale_ids, pdf_325300_ids, central_id

# ---- Kinematic helpers -------------------------------------------------------

def _mll(p1, p2):
    p = p1 + p2
    return np.sqrt(max(p[3]**2 - p[0]**2 - p[1]**2 - p[2]**2, 0.0))

def _rap(p1, p2):
    p = p1 + p2
    E, pz = p[3], p[2]
    return abs(0.5 * np.log((E + pz) / (E - pz)))

def _cstar(p1, p2):
    p     = p1 + p2
    E, pz = p[3], p[2]
    mass  = _mll(p1, p2)
    beta  = pz / E
    gamma = E / mass
    pz1_b = gamma * (p1[2] - beta * p1[3])
    p1mag = np.sqrt(p1[0]**2 + p1[1]**2 + pz1_b**2)
    return pz1_b / p1mag

def _eta(p):
    px, py, pz, e = p
    pmag = np.sqrt(px**2 + py**2 + pz**2)
    if pmag == abs(pz):
        return np.inf
    return 0.5 * np.log((pmag + pz) / (pmag - pz))

def _pt(p):
    return np.sqrt(p[0]**2 + p[1]**2)

# ---- Per-bin worker -----------------------------------------------------------

def process_bin(bin_spec):
    """
    Read all LHE files belonging to one mll bin, apply cuts, accumulate arrays.
    Returns a dict of numpy arrays or loads from per-bin checkpoint if present.
    """
    tag, lhe_files = bin_spec
    ckpt = os.path.join(PER_FILE_DIR, f"{tag}.pkl")

    if os.path.exists(ckpt):
        print(f"[{tag}] Loading from checkpoint")
        with open(ckpt, "rb") as f:
            return pickle.load(f)

    print(f"[{tag}] Starting ({len(lhe_files)} files)")
    if not lhe_files:
        print(f"[{tag}] WARNING: no files found, skipping")
        return None

    # headers are identical across all files of a bin (same run_card/reweight_card)
    scale_ids, pdf_325300_ids, central_id = parse_weight_ids(lhe_files[0])
    n_pdf = len(pdf_325300_ids)
    print(f"[{tag}] scale IDs={scale_ids}  PDF members={n_pdf}  central={central_id}")

    buf_mll           = []
    buf_rap           = []
    buf_cstar         = []
    buf_w_SM          = []
    buf_xwgt          = []
    buf_w_p1          = {op: [] for op in OPERATORS}
    buf_w_m1          = {op: [] for op in OPERATORS}
    buf_w_pp          = {pair: [] for pair in OP_PAIRS}
    buf_w_scale       = {k: [] for k in scale_ids}
    buf_w_pdf_central = []
    buf_pdf_325300    = []

    pp_keys   = {}
    n_read    = 0
    n_kept    = 0
    n_cut_eta = 0
    n_cut_pt  = 0

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)

        for lhe_file in lhe_files:
            if MAX_EVENTS is not None and n_read >= MAX_EVENTS:
                break

            events = pylhe.read_lhe_with_attributes(lhe_file)

            for event in events:
                if MAX_EVENTS is not None and n_read >= MAX_EVENTS:
                    break
                n_read += 1
                if n_read % 10000 == 0:
                    print(f"[{tag}] {n_read} events read, {n_kept} kept")

                leptons = [
                    p for p in event.particles
                    if int(p.status) == 1 and abs(int(p.id)) in {11, 13}
                ]
                if len(leptons) < 2:
                    continue

                lm = next((p for p in leptons if int(p.id) > 0), leptons[0])
                lp = next((p for p in leptons if int(p.id) < 0), leptons[1])
                v_lm = np.array([lm.px, lm.py, lm.pz, lm.e])
                v_lp = np.array([lp.px, lp.py, lp.pz, lp.e])

                m  = _mll(v_lm, v_lp)
                if not (MLL_LO <= m <= MLL_HI):
                    continue

                # Fiducial cuts
                if do_fiducial:
                    if ETA_MAX is not None:
                        if abs(_eta(v_lm)) > ETA_MAX or abs(_eta(v_lp)) > ETA_MAX:
                            n_cut_eta += 1
                            continue
                    if PT_LEAD is not None or PT_SUB is not None:
                        pt_lm = _pt(v_lm); pt_lp = _pt(v_lp)
                        pt_l  = max(pt_lm, pt_lp); pt_s = min(pt_lm, pt_lp)
                        if PT_LEAD is not None and pt_l < PT_LEAD:
                            n_cut_pt += 1
                            continue
                        if PT_SUB is not None and pt_s < PT_SUB:
                            n_cut_pt += 1
                            continue

                y  = _rap(v_lm, v_lp)
                cs = _cstar(v_lm, v_lp)

                wkeys = event.weights

                if not pp_keys and not SKIP_PAIRS:
                    for op1, op2 in OP_PAIRS:
                        pp_keys[(op1, op2)] = (
                            f'{op1}_{op2}' if f'{op1}_{op2}' in wkeys else f'{op2}_{op1}'
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
                n_kept += 1

    print(f"[{tag}] Done: {n_read} read, {n_kept} kept"
          + (f", {n_cut_eta} cut by eta, {n_cut_pt} cut by pT" if do_fiducial else ""))

    result = {
        'mll':           np.array(buf_mll,   dtype=np.float64),
        'rap':           np.array(buf_rap,   dtype=np.float64),
        'cstar':         np.array(buf_cstar, dtype=np.float64),
        'w_SM':          np.array(buf_w_SM,  dtype=np.float64),
        'xwgt':          np.array(buf_xwgt,  dtype=np.float64),
        'w_p1':          {op: np.array(buf_w_p1[op], dtype=np.float64) for op in OPERATORS},
        'w_m1':          {op: np.array(buf_w_m1[op], dtype=np.float64) for op in OPERATORS},
        'w_pp':          {pair: np.array(buf_w_pp[pair], dtype=np.float64) for pair in OP_PAIRS},
        'w_scale':       {k: np.array(buf_w_scale[k], dtype=np.float64) for k in scale_ids},
        'w_pdf_central': np.array(buf_w_pdf_central, dtype=np.float64),
        'pdf_325300':    np.array(buf_pdf_325300, dtype=np.float64).reshape(-1, n_pdf),
    }

    os.makedirs(PER_FILE_DIR, exist_ok=True)
    with open(ckpt, "wb") as f:
        pickle.dump(result, f)
    size_mb = os.path.getsize(ckpt) / 1e6
    print(f"[{tag}] Checkpoint saved ({size_mb:.1f} MB)")

    return result

# ---- Main --------------------------------------------------------------------

if __name__ == "__main__":
    os.makedirs(CACHE_DIR, exist_ok=True)

    print(f"Processing {len(LHE_BINS)} mll bins in parallel "
          f"({sum(len(f) for _, f in LHE_BINS)} LHE files total)")
    for tag, files in LHE_BINS:
        print(f"  [{tag}] {len(files)} files found")
    if do_fiducial:
        print(f"Fiducial cuts: eta_max={ETA_MAX}, pt_lead={PT_LEAD}, pt_sub={PT_SUB}")
    else:
        print("No fiducial cuts (parton-level inclusive)")
    print()

    n_workers = args.workers or len(LHE_BINS)
    with Pool(processes=n_workers) as pool:
        results = pool.map(process_bin, LHE_BINS)

    # Filter failed bins
    results = [r for r in results if r is not None]
    if not results:
        raise RuntimeError("No results — check LHE_BASE / bin directories.")

    print(f"\nMerging {len(results)} bin results ...")

    def concat(arrays):
        return np.concatenate(arrays)

    cache = {
        'mll':           concat([r['mll']   for r in results]),
        'rap':           concat([r['rap']   for r in results]),
        'cstar':         concat([r['cstar'] for r in results]),
        'w_SM':          concat([r['w_SM']  for r in results]),
        'xwgt':          concat([r['xwgt']  for r in results]),
        'w_p1':          {op: concat([r['w_p1'][op] for r in results]) for op in OPERATORS},
        'w_m1':          {op: concat([r['w_m1'][op] for r in results]) for op in OPERATORS},
        'w_pp':          {pair: concat([r['w_pp'][pair] for r in results]) for pair in OP_PAIRS},
        'w_scale':       {},
        'w_pdf_central': concat([r['w_pdf_central'] for r in results]),
        'pdf_325300':    np.concatenate([r['pdf_325300'] for r in results], axis=0),
        'cuts': {
            'eta_max': ETA_MAX,
            'pt_lead': PT_LEAD,
            'pt_sub':  PT_SUB,
            'mll_lo':  MLL_LO,
            'mll_hi':  MLL_HI,
        },
    }

    # Merge scale weight dicts (keys may differ per file if header varies)
    all_scale_keys = set()
    for r in results:
        all_scale_keys.update(r['w_scale'].keys())
    for k in all_scale_keys:
        arrays = [r['w_scale'][k] for r in results if k in r['w_scale']]
        cache['w_scale'][k] = concat(arrays)

    n_total = len(cache['mll'])
    print(f"Total events in cache : {n_total:,}")
    print(f"pdf_325300 shape      : {cache['pdf_325300'].shape}")
    print(f"Scale variation keys  : {len(cache['w_scale'])}")

    with open(CACHE_FILE, 'wb') as f:
        pickle.dump(cache, f)
    size_mb = os.path.getsize(CACHE_FILE) / 1e6
    print(f"\nCache saved -> {CACHE_FILE}  ({size_mb:.1f} MB)")
