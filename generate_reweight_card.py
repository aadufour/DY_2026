#!/usr/bin/env python3
"""
Generate MG5 reweight_card.dat for SMEFT DY analysis.

Base run: all active operators set to 0.9999...
This card produces:
  1   SM weight          (all C = 0)
  17  single-operator    (one C = 1, rest = 0)
  136 two-operator pairs (two C = 1, rest = 0)
  ─────────────────────────────────────────────
  154 launch blocks total
"""

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

OUTPUT = '/Users/albertodufour/MG5_2_9_18/mg5amcnlo/DY_all/Cards/reweight_card.dat'


def format_block(rwgt_name, on_indices, value=1.0):
    """Format a single launch block, setting on_indices to value and the rest to 0."""
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
# This card generates 154 reweight points:
#   1   SM          (all C = 0)
#   17  single-op   (one C = 1, rest = 0)
#   136 pairs       (two C = 1, rest = 0)
#
# CP-odd operators (SMEFTcpv block) are untouched throughout.
#
#*************************************************************************

change mode NLO

"""

blocks = []

# ── 1. SM ──────────────────────────────────────────────────────────────
blocks.append(format_block('SM', set()))

# ── 2. Individual operators ────────────────────────────────────────────
for idx, name in OPERATORS.items():
    blocks.append(format_block(name, {idx}))

# ── 3. All two-operator pairs ──────────────────────────────────────────
op_items = list(OPERATORS.items())
for (idx1, name1), (idx2, name2) in combinations(op_items, 2):
    blocks.append(format_block(f'{name1}_{name2}', {idx1, idx2}))

# ── Write ──────────────────────────────────────────────────────────────
with open(OUTPUT, 'w') as f:
    f.write(HEADER)
    f.write('\n\n'.join(blocks))
    f.write('\n')

n_single = len(OPERATORS)
n_pairs  = len(blocks) - 1 - n_single
print(f"Written to {OUTPUT}")
print(f"  1 SM  +  {n_single} single-op  +  {n_pairs} pairs  =  {len(blocks)} blocks total")
