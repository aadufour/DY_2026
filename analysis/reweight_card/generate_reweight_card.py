#!/usr/bin/env python3
"""
Generate MG5 reweight_card.dat for SMEFT DY analysis.

Base run: all active operators set to 0.9999...
This card produces:
  1   SM          (all C = 0)
  N   single +1   (one C = +1, rest = 0)
  N   single -1   (one C = -1, rest = 0)
  N*(N-1)/2 pairs (two C = +1, rest = 0)

The active operator list is read dynamically from operator_scan.html
(green rows = contributing operators).

Writes to Cards/reweight_card.dat for every mass-window directory.
"""

import os
import re
from itertools import combinations

# ── HTML source of contributing operators ──────────────────────────────
HTML_PATH = '/Users/albertodufour/code/DY2026/analysis/restriction_card/operator_scan.html'

# ── Full SMEFT block index mapping (SMEFTatNLO ordering) ───────────────
NAME_TO_INDEX = {
    'cG':         1,  'cW':         2,  'cH':         3,  'cHbox':      4,
    'cHDD':       5,  'cHG':        6,  'cHW':        7,  'cHB':        8,
    'cHWB':       9,  'cuHRe':     10,  'ctHRe':     11,  'cdHRe':     12,
    'cbHRe':     13,  'cuGRe':     14,  'ctGRe':     15,  'cuWRe':     16,
    'ctWRe':     17,  'cuBRe':     18,  'ctBRe':     19,  'cdGRe':     20,
    'cbGRe':     21,  'cdWRe':     22,  'cbWRe':     23,  'cdBRe':     24,
    'cbBRe':     25,  'cHj1':      26,  'cHQ1':      27,  'cHj3':      28,
    'cHQ3':      29,  'cHu':       30,  'cHt':       31,  'cHd':       32,
    'cHbq':      33,  'cHudRe':    34,  'cHtbRe':    35,  'cjj11':     36,
    'cjj18':     37,  'cjj31':     38,  'cjj38':     39,  'cQj11':     40,
    'cQj18':     41,  'cQj31':     42,  'cQj38':     43,  'cQQ1':      44,
    'cQQ8':      45,  'cuu1':      46,  'cuu8':      47,  'ctt':       48,
    'ctu1':      49,  'ctu8':      50,  'cdd1':      51,  'cdd8':      52,
    'cbb':       53,  'cbd1':      54,  'cbd8':      55,  'cud1':      56,
    'ctb1':      57,  'ctd1':      58,  'cbu1':      59,  'cud8':      60,
    'ctb8':      61,  'ctd8':      62,  'cbu8':      63,  'cutbd1Re':  64,
    'cutbd8Re':  65,  'cju1':      66,  'cQu1':      67,  'cju8':      68,
    'cQu8':      69,  'ctj1':      70,  'ctj8':      71,  'cQt1':      72,
    'cQt8':      73,  'cjd1':      74,  'cjd8':      75,  'cQd1':      76,
    'cQd8':      77,  'cbj1':      78,  'cbj8':      79,  'cQb1':      80,
    'cQb8':      81,  'cjQtu1Re':  82,  'cjQtu8Re':  83,  'cjQbd1Re':  84,
    'cjQbd8Re':  85,  'cjujd1Re':  86,  'cjujd8Re':  87,  'cjujd11Re': 88,
    'cjujd81Re': 89,  'cQtjd1Re':  90,  'cQtjd8Re':  91,  'cjuQb1Re':  92,
    'cjuQb8Re':  93,  'cQujb1Re':  94,  'cQujb8Re':  95,  'cjtQd1Re':  96,
    'cjtQd8Re':  97,  'cQtQb1Re':  98,  'cQtQb8Re':  99,  'ceHRe':    100,
    'ceWRe':    101,  'ceBRe':    102,  'cHl1':     103,  'cHl3':     104,
    'cHe':      105,  'cll':      106,  'cll1':     107,  'clj1':     108,
    'clj3':     109,  'cQl1':     110,  'cQl3':     111,  'cee':      112,
    'ceu':      113,  'cte':      114,  'ced':      115,  'cbe':      116,
    'cje':      117,  'cQe':      118,  'clu':      119,  'ctl':      120,
    'cld':      121,  'cbl':      122,  'cle':      123,  'cledjRe':  124,
    'clebQRe':  125,  'cleju1Re': 126,  'cleQt1Re': 127,  'cleju3Re': 128,
    'cleQt3Re': 129,
}


def read_contributing_operators(html_path):
    """
    Parse operator_scan.html and return an ordered dict {index: name}
    for every green (contributing) row.
    """
    with open(html_path) as f:
        content = f.read()

    # Each green row looks like:
    #   <tr><td ...>NAME</td><td ... background:#d4edda ...>
    pattern = r'<tr><td[^>]*>(\w+)</td><td[^>]*background:#d4edda[^>]*>'
    names = re.findall(pattern, content)

    operators = {}
    missing = []
    for name in names:
        if name in NAME_TO_INDEX:
            operators[NAME_TO_INDEX[name]] = name
        else:
            missing.append(name)

    if missing:
        print(f"WARNING: no SMEFT index known for: {missing}")
        print("  Add them to NAME_TO_INDEX in this script.")

    # Sort by block index so the card is deterministic
    return dict(sorted(operators.items()))


# ── Load operators from HTML ────────────────────────────────────────────
OPERATORS = read_contributing_operators(HTML_PATH)
print(f"Loaded {len(OPERATORS)} contributing operators from {HTML_PATH}")

# ── Mass window directories to write cards for ─────────────────────────
MG5_BASE  = '/Users/albertodufour/MG5_2_9_18/mg5amcnlo'
MLL_EDGES = [50, 120, 200, 400, 600, 800, 1000, 3000]
# OUTPUTS   = [
#     os.path.join(MG5_BASE, f'DY_all_{lo}_{hi}', 'Cards', 'reweight_card.dat')
#     for lo, hi in zip(MLL_EDGES[:-1], MLL_EDGES[1:])
# ]
OUTPUTS = ['/Users/albertodufour/code/DY2026/analysis/reweight_card/reweight_card.dat']


def format_block(rwgt_name, on_indices, value=1.0):
    lines = [f'launch --rwgt_name={rwgt_name}']
    for idx in OPERATORS:
        val = value if idx in on_indices else 0
        lines.append(f'    set SMEFT {idx} {val}')
    return '\n'.join(lines)


n       = len(OPERATORS)
n_pairs = n * (n - 1) // 2
HEADER  = f"""\
#*************************************************************************
#                          Reweight Module                               *
#                Matrix-Element reweighting at LO/NLO                    *
#  Mattelaer Olivier                                    arxiv:1607.00763 *
#*************************************************************************
#
# Base run: all active operators set to 0.9999... (effectively C=1)
#
# This card generates {1 + 2*n + n_pairs} reweight points:
#   1      SM           (all C = 0)
#   {n:<6d} single +1    (one C = +1, rest = 0)
#   {n:<6d} single -1    (one C = -1, rest = 0)
#   {n_pairs:<6d} pairs        (two C = +1, rest = 0)
#
# Operators read from: {HTML_PATH}
#
# CP-odd operators (SMEFTcpv block) are untouched throughout.
#
#*************************************************************************

change mode NLO

"""

# ── Build blocks (same for every mass window) ──────────────────────────
blocks = []
op_items = list(OPERATORS.items())

blocks.append(format_block('SM', set()))

for idx, name in op_items:
    blocks.append(format_block(name,           {idx}, value= 1.0))
for idx, name in op_items:
    blocks.append(format_block(f'minus{name}', {idx}, value=-1.0))

for (idx1, name1), (idx2, name2) in combinations(op_items, 2):
    blocks.append(format_block(f'{name1}_{name2}', {idx1, idx2}))

content = HEADER + '\n\n'.join(blocks) + '\n'

n_single = len(OPERATORS)
n_pairs  = len(blocks) - 1 - 2 * n_single
print(f"1 SM  +  {n_single} (+1)  +  {n_single} (-1)  +  {n_pairs} pairs  =  {len(blocks)} blocks")

# ── Write to each directory that exists ───────────────────────────────
for path in OUTPUTS:
    cards_dir = os.path.dirname(path)
    if not os.path.isdir(cards_dir):
        print(f"  SKIP (Cards/ dir not found): {path}")
        continue
    with open(path, 'w') as f:
        f.write(content)
    print(f"  Written: {path}")
