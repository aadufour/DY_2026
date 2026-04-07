#!/usr/bin/env python3
"""
write_reweight_card.py

Generate and write the MG5 reweight_card.dat to one or more explicit output paths.
Accepts --output arguments so it works on any machine regardless of local MG5
directory layout.

Operators included:
  17 confirmed (green in operator scan)
  +9 b quark ops: cbBRe, cHQ1, cHQ3, cHbq, cQl1, cQl3, cbe, cQe, cbl
  = 26 total

  1 SM  +  26 (+1)  +  26 (-1)  +  325 pairs  =  378 blocks

Usage (llrcms, single DY_all):
    python3 write_reweight_card.py \
        --output /home/llr/cms/adufour/MG5/mg5amcnlo/DY_all/Cards/reweight_card.dat

Usage (dry run, print to stdout):
    python3 write_reweight_card.py --stdout
"""

import argparse
import os
from itertools import combinations

# SMEFT block index : operator name
OPERATORS = {
    # --- confirmed green ---
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
# b quark ops
    25:  'cbBRe',
    27:  'cHQ1',
    29:  'cHQ3',
    33:  'cHbq',
    110: 'cQl1',
    111: 'cQl3',
    116: 'cbe',
    118: 'cQe',
    122: 'cbl',
}

HEADER = """\
#*************************************************************************
#                          Reweight Module                               *
#                Matrix-Element reweighting at LO/NLO                    *
#  Mattelaer Olivier                                    arxiv:1607.00763 *
#*************************************************************************
#
# Base run: all active operators set to 0.9999... (effectively C=1)
#
# Operators: 17 confirmed (green) + 9 uncertain (yellow) = 26 total
#
# This card generates 378 reweight points:
#   1   SM           (all C = 0)
#   26  single +1    (one C = +1, rest = 0)
#   26  single -1    (one C = -1, rest = 0)
#   325 pairs        (two C = +1, rest = 0)
#
# CP-odd operators (SMEFTcpv block) are untouched throughout.
#
#*************************************************************************

change mode NLO

"""


def format_block(rwgt_name, on_indices, value=1.0):
    lines = [f'launch --rwgt_name={rwgt_name}']
    for idx in OPERATORS:
        val = value if idx in on_indices else 0
        lines.append(f'    set SMEFT {idx} {val}')
    return '\n'.join(lines)


def build_content():
    op_items = list(OPERATORS.items())
    blocks = []
    blocks.append(format_block('SM', set()))
    for idx, name in op_items:
        blocks.append(format_block(name, {idx}, value=1.0))
    for idx, name in op_items:
        blocks.append(format_block(f'minus{name}', {idx}, value=-1.0))
    for (idx1, name1), (idx2, name2) in combinations(op_items, 2):
        blocks.append(format_block(f'{name1}_{name2}', {idx1, idx2}))
    n = len(blocks)
    n_ops = len(OPERATORS)
    print(f"1 SM  +  {n_ops} (+1)  +  {n_ops} (-1)  +  {n - 1 - 2*n_ops} pairs  =  {n} blocks")
    return HEADER + '\n\n'.join(blocks) + '\n'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Write MG5 reweight_card.dat.")
    parser.add_argument('--output', metavar='PATH', action='append', default=[],
                        help="Output path (can be repeated for multiple directories)")
    parser.add_argument('--stdout', action='store_true',
                        help="Print card to stdout instead of writing files")
    args = parser.parse_args()

    content = build_content()

    if args.stdout:
        print(content)
    elif not args.output:
        parser.error("Provide at least one --output path, or use --stdout")
    else:
        for path in args.output:
            cards_dir = os.path.dirname(path)
            if not os.path.isdir(cards_dir):
                print(f"  SKIP (Cards/ dir not found): {path}")
                continue
            with open(path, 'w') as f:
                f.write(content)
            print(f"  Written: {path}")
