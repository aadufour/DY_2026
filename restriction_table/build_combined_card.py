#!/usr/bin/env python3
"""
Build a single restriction card with all contributing operators turned on.

Usage:
    python build_combined_card.py --operators cHj1 cHQ1 cHe cHl1 cHl3 ...
    python build_combined_card.py --from-html operator_scan.html
    python build_combined_card.py --from-txt contributing.txt
"""

import argparse
import os
import re

params = [
    [' 1', 'cG'],
    [' 2', 'cW'],
    [' 3', 'cH'],
    [' 4', 'cHbox'],
    [' 5', 'cHDD'],
    [' 6', 'cHG'],
    [' 7', 'cHW'],
    [' 8', 'cHB'],
    [' 9', 'cHWB'],
    [' 10', 'cuHRe'],
    [' 11', 'ctHRe'],
    [' 12', 'cdHRe'],
    [' 13', 'cbHRe'],
    [' 14', 'cuGRe'],
    [' 15', 'ctGRe'],
    [' 16', 'cuWRe'],
    [' 17', 'ctWRe'],
    [' 18', 'cuBRe'],
    [' 19', 'ctBRe'],
    [' 20', 'cdGRe'],
    [' 21', 'cbGRe'],
    [' 22', 'cdWRe'],
    [' 23', 'cbWRe'],
    [' 24', 'cdBRe'],
    [' 25', 'cbBRe'],
    [' 26', 'cHj1'],
    [' 27', 'cHQ1'],
    [' 28', 'cHj3'],
    [' 29', 'cHQ3'],
    [' 30', 'cHu'],
    [' 31', 'cHt'],
    [' 32', 'cHd'],
    [' 33', 'cHbq'],
    [' 34', 'cHudRe'],
    [' 35', 'cHtbRe'],
    [' 36', 'cjj11'],
    [' 37', 'cjj18'],
    [' 38', 'cjj31'],
    [' 39', 'cjj38'],
    [' 40', 'cQj11'],
    [' 41', 'cQj18'],
    [' 42', 'cQj31'],
    [' 43', 'cQj38'],
    [' 44', 'cQQ1'],
    [' 45', 'cQQ8'],
    [' 46', 'cuu1'],
    [' 47', 'cuu8'],
    [' 48', 'ctt'],
    [' 49', 'ctu1'],
    [' 50', 'ctu8'],
    [' 51', 'cdd1'],
    [' 52', 'cdd8'],
    [' 53', 'cbb'],
    [' 54', 'cbd1'],
    [' 55', 'cbd8'],
    [' 56', 'cud1'],
    [' 57', 'ctb1'],
    [' 58', 'ctd1'],
    [' 59', 'cbu1'],
    [' 60', 'cud8'],
    [' 61', 'ctb8'],
    [' 62', 'ctd8'],
    [' 63', 'cbu8'],
    [' 64', 'cutbd1Re'],
    [' 65', 'cutbd8Re'],
    [' 66', 'cju1'],
    [' 67', 'cQu1'],
    [' 68', 'cju8'],
    [' 69', 'cQu8'],
    [' 70', 'ctj1'],
    [' 71', 'ctj8'],
    [' 72', 'cQt1'],
    [' 73', 'cQt8'],
    [' 74', 'cjd1'],
    [' 75', 'cjd8'],
    [' 76', 'cQd1'],
    [' 77', 'cQd8'],
    [' 78', 'cbj1'],
    [' 79', 'cbj8'],
    [' 80', 'cQb1'],
    [' 81', 'cQb8'],
    [' 82', 'cjQtu1Re'],
    [' 83', 'cjQtu8Re'],
    [' 84', 'cjQbd1Re'],
    [' 85', 'cjQbd8Re'],
    [' 86', 'cjujd1Re'],
    [' 87', 'cjujd8Re'],
    [' 88', 'cjujd11Re'],
    [' 89', 'cjujd81Re'],
    [' 90', 'cQtjd1Re'],
    [' 91', 'cQtjd8Re'],
    [' 92', 'cjuQb1Re'],
    [' 93', 'cjuQb8Re'],
    [' 94', 'cQujb1Re'],
    [' 95', 'cQujb8Re'],
    [' 96', 'cjtQd1Re'],
    [' 97', 'cjtQd8Re'],
    [' 98', 'cQtQb1Re'],
    [' 99', 'cQtQb8Re'],
    [' 100', 'ceHRe'],
    [' 101', 'ceWRe'],
    [' 102', 'ceBRe'],
    [' 103', 'cHl1'],
    [' 104', 'cHl3'],
    [' 105', 'cHe'],
    [' 106', 'cll'],
    [' 107', 'cll1'],
    [' 108', 'clj1'],
    [' 109', 'clj3'],
    [' 110', 'cQl1'],
    [' 111', 'cQl3'],
    [' 112', 'cee'],
    [' 113', 'ceu'],
    [' 114', 'cte'],
    [' 115', 'ced'],
    [' 116', 'cbe'],
    [' 117', 'cje'],
    [' 118', 'cQe'],
    [' 119', 'clu'],
    [' 120', 'ctl'],
    [' 121', 'cld'],
    [' 122', 'cbl'],
    [' 123', 'cle'],
    [' 124', 'cledjRe'],
    [' 125', 'clebQRe'],
    [' 126', 'cleju1Re'],
    [' 127', 'cleQt1Re'],
    [' 128', 'cleju3Re'],
    [' 129', 'cleQt3Re'],
]

HERE = os.path.dirname(os.path.abspath(__file__))


def load_contributing_from_html(html_path):
    """Extract contributing operator names from the scan HTML output."""
    with open(html_path) as f:
        content = f.read()
    # Rows with checkmark have background:#d4edda (green)
    names = re.findall(r'<td[^>]*>([^<]+)</td>\s*<td[^>]*background:#d4edda[^>]*>', content)
    return [n.strip() for n in names]


def load_contributing_from_txt(txt_path):
    """Read operator names one per line (blank lines and # comments ignored)."""
    names = []
    with open(txt_path) as f:
        for line in f:
            line = line.split('#')[0].strip()
            if line:
                names.append(line)
    return names


def build_card(contributing_set, out_path):
    before_path = os.path.join(HERE, 'restrict_before.txt')
    after_path  = os.path.join(HERE, 'restrict_after.txt')
    with open(before_path) as f:
        contents_before = f.read()
    with open(after_path) as f:
        contents_after = f.read()

    with open(out_path, 'w') as f:
        f.write(contents_before)
        for param in params:
            idx, name = param[0], param[1]
            if name in contributing_set:
                f.write(f'   {idx} 9.999999e-01 # {name}\n')
            else:
                f.write(f'   {idx} 0 # {name}\n')
        f.write(contents_after)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Build combined restriction card for contributing SMEFT operators.')
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument('--operators', nargs='+', metavar='OP',
                     help='Operator names to enable, e.g. cHj1 cHl1 cHe ...')
    src.add_argument('--from-html', metavar='FILE',
                     help='Read contributing operators from scan_operators HTML output')
    src.add_argument('--from-txt', metavar='FILE',
                     help='Read contributing operators from a plain text file (one per line)')
    parser.add_argument('--out', default=os.path.join(HERE, 'restrict_all_contributing_massless.dat'),
                        help='Output restriction card path')
    args = parser.parse_args()

    if args.operators:
        contributing = set(args.operators)
    elif args.from_html:
        contributing = set(load_contributing_from_html(args.from_html))
    else:
        contributing = set(load_contributing_from_txt(args.from_txt))

    # Validate names
    known = {p[1] for p in params}
    unknown = contributing - known
    if unknown:
        print(f'WARNING: unknown operator names (will be ignored): {sorted(unknown)}')
        contributing -= unknown

    build_card(contributing, args.out)
    print(f'Written: {args.out}')
    print(f'Operators enabled ({len(contributing)}): {", ".join(sorted(contributing))}')
