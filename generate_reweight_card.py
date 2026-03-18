#!/usr/bin/env python3
"""
Generate MG5 reweight_card.dat for SMEFT DY analysis.

Base run: all active operators set to 0.9999...
This card produces:
  1   SM          (all C = 0)
  17  single +1   (one C = +1, rest = 0)
  17  single -1   (one C = -1, rest = 0)
  136 pairs       (two C = +1, rest = 0)
  ─────────────────────────────────────────────
  171 launch blocks total

Writes to Cards/reweight_card.dat for every mass-window directory.
"""

import os
from itertools import combinations

# Active operators in the base run (SMEFT block index : name)
OPERATORS = {
    5:   'cHDD',
    9:   'cHWB',
    26:  'cHj1',
    28:  'cHj3',
    30:  'cHu',
    32:  'cHd',
    103: 'cHl1',
    104: 'cHl3',
    105: 'cHe',
    107: 'cll1',
    108: 'clj1',
    109: 'clj3',
    113: 'ceu',
    115: 'ced',
    117: 'cje',
    119: 'clu',
    121: 'cld',
}

# Mass window directories to write cards for
MG5_BASE  = '/Users/albertodufour/MG5_2_9_18/mg5amcnlo'
MLL_EDGES = [50, 120, 200, 400, 600, 800, 1000, 3000]
OUTPUTS   = [
    os.path.join(MG5_BASE, f'DY_all_{lo}_{hi}', 'Cards', 'reweight_card.dat')
    for lo, hi in zip(MLL_EDGES[:-1], MLL_EDGES[1:])
]


def format_block(rwgt_name, on_indices, value=1.0):
    lines = [f'launch --rwgt_name={rwgt_name}']
    for idx in OPERATORS:
        val = value if idx in on_indices else 0
        lines.append(f'    set SMEFT {idx} {val}')
    return '\n'.join(lines)


HEADER = """\
#*************************************************************************
#                          Reweight Module                               *
#                Matrix-Element reweighting at LO/NLO                    *
#  Mattelaer Olivier                                    arxiv:1607.00763 *
#*************************************************************************
#
# Base run: all active operators set to 0.9999... (effectively C=1)
#
# This card generates 171 reweight points:
#   1   SM           (all C = 0)
#   17  single +1    (one C = +1, rest = 0)
#   17  single -1    (one C = -1, rest = 0)
#   136 pairs        (two C = +1, rest = 0)
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
    blocks.append(format_block(name,          {idx}, value= 1.0))
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
