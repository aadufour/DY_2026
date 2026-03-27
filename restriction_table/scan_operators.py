#!/usr/bin/env python3
"""
Scan all SMEFT operators for contributions to p p > l+ l-
Runs MG5 for each operator, parses diagram count from stdout,
and produces an HTML summary table.

Usage (run from MadGraph root):
    python scan_operators.py --mg5 ./bin/mg5_aMC
"""

import os
import re
import subprocess
import argparse
import shutil
import webbrowser

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

MODEL = 'SMEFTsim_topU3l_MwScheme_UFO'
PROCESS = 'p p > l+ l- QCD=0 SMHLOOP=0 NP=1 NP^2==2'

MG5_INPUT_TEMPLATE = """\
import model {model}-{name}_massless
define p = p b b~
generate {process}
display processes
quit
"""

MG5_DIAGRAMS_TEMPLATE = """\
import model {model}-{name}_massless
define p = p b b~
generate {process}
output {outdir}
quit
"""

def count_diagrams(stdout):
    """Parse mg5 stdout for diagram count. Returns int (0 if none found)."""
    total = 0
    for m in re.finditer(r'(\d+)\s+diagram', stdout):
        total += int(m.group(1))
    return total

def install_card(mg5root, name):
    """Copy the restriction card into the MadGraph model directory."""
    card_src = os.path.join(os.path.dirname(__file__), 'cards', f'restrict_{name}_massless.dat')
    model_dir = os.path.join(mg5root, 'models', MODEL)
    card_dst = os.path.join(model_dir, f'restrict_{name}_massless.dat')
    if not os.path.isdir(model_dir):
        raise FileNotFoundError(f'Model directory not found: {model_dir}')
    if not os.path.isfile(card_src):
        raise FileNotFoundError(f'Restriction card not found: {card_src}\nRun build_restrict.py first.')
    shutil.copy2(card_src, card_dst)
    return card_dst

def run_mg5(mg5_path, name, mg5root, debug=False, verbose=False):
    import tempfile
    install_card(mg5root, name)
    mg5_input = MG5_INPUT_TEMPLATE.format(model=MODEL, name=name, process=PROCESS)
    if verbose:
        print(f'\n--- mg5 commands for {name} ---')
        print(mg5_input)
        print('--- end commands ---\n')
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tmp:
        tmp.write(mg5_input)
        tmp_path = tmp.name
    try:
        result = subprocess.run(
            [mg5_path, tmp_path],
            capture_output=True,
            text=True,
        )
        stdout = result.stdout + result.stderr
    finally:
        os.unlink(tmp_path)
    if debug:
        print('\n--- mg5 output ---')
        print(stdout)
        print('--- end ---\n')
    n_diag = count_diagrams(stdout)
    return n_diag, stdout

def run_mg5_diagrams(mg5_path, name, mg5root, outdir=None, verbose=False):
    """Run MG5 with `output` to generate Feynman diagram HTML for a single operator."""
    import tempfile
    install_card(mg5root, name)

    if outdir is None:
        outdir = os.path.abspath(f'diagrams_{name}')

    # MG5 will refuse to overwrite without interaction — remove first
    if os.path.exists(outdir):
        shutil.rmtree(outdir)

    mg5_input = MG5_DIAGRAMS_TEMPLATE.format(model=MODEL, name=name, process=PROCESS, outdir=outdir)
    if verbose:
        print(f'\n--- mg5 commands for {name} (diagram mode) ---')
        print(mg5_input)
        print('--- end commands ---\n')

    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tmp:
        tmp.write(mg5_input)
        tmp_path = tmp.name
    try:
        subprocess.run([mg5_path, tmp_path], capture_output=False, text=True)
    finally:
        os.unlink(tmp_path)

    index = os.path.join(outdir, 'HTML', 'index.html')
    if os.path.isfile(index):
        print(f'\nDiagrams written to: {index}')
        webbrowser.open(f'file://{os.path.abspath(index)}')
    else:
        print(f'\nWarning: expected diagram index not found at {index}')
        print(f'Output directory contents: {os.listdir(outdir) if os.path.isdir(outdir) else "(missing)"}')

def write_html(results, outfile):
    contributing = [(n, d) for n, d, _ in results if d > 0]
    not_contributing = [(n, d) for n, d, _ in results if d == 0]

    rows = ''
    for name, n_diag, _ in results:
        if n_diag > 0:
            mark = '&#10003;'  # checkmark
            style = 'background:#d4edda; color:#155724;'
        else:
            mark = '&#10007;'  # cross
            style = 'background:#f8d7da; color:#721c24;'
        rows += f'  <tr><td style="padding:4px 12px">{name}</td><td style="padding:4px 12px; text-align:center; {style}">{mark} {n_diag}</td></tr>\n'

    html = f"""\
<!DOCTYPE html>
<html>
<head><meta charset="utf-8"><title>DY SMEFT operator scan</title>
<style>
  body {{ font-family: monospace; }}
  table {{ border-collapse: collapse; }}
  td {{ border: 1px solid #ccc; }}
  h2 {{ margin-top: 1.5em; }}
</style>
</head>
<body>
<h1>pp &rarr; l<sup>+</sup>l<sup>-</sup> &mdash; SMEFT operator scan</h1>
<p>{len(contributing)} contributing / {len(not_contributing)} not contributing</p>
<table>
  <tr><th style="padding:4px 12px">Operator</th><th style="padding:4px 12px">Diagrams</th></tr>
{rows}</table>
</body>
</html>
"""
    with open(outfile, 'w') as f:
        f.write(html)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mg5', required=True, help='Path to mg5_aMC executable')
    parser.add_argument('--mg5root', default=None, help='MadGraph root directory (default: inferred from --mg5)')
    parser.add_argument('--out', default='operator_scan.html', help='Output HTML file')
    parser.add_argument('--debug', action='store_true', help='Print raw mg5 output for first operator then exit')
    parser.add_argument('--verbose', action='store_true', help='Print the exact mg5 commands sent for each operator')
    parser.add_argument('--operators', nargs='+', metavar='OP',
                        help='Only scan these operators (e.g. --operators cHl1 cHl3 cHe)')
    parser.add_argument('--show-diagrams', metavar='OP',
                        help='Generate and open Feynman diagrams for a single operator (skips the scan)')
    args = parser.parse_args()

    mg5root = args.mg5root or os.path.dirname(os.path.dirname(os.path.abspath(args.mg5)))
    print(f'MadGraph root: {mg5root}')
    print(f'Model directory: {os.path.join(mg5root, "models", MODEL)}')

    # Diagram-view mode: generate output for a single operator and open HTML
    if args.show_diagrams:
        name = args.show_diagrams
        known = {p[1] for p in params}
        if name not in known:
            raise SystemExit(f'Unknown operator: {name}')
        print(f'Generating diagrams for {name} ...')
        run_mg5_diagrams(args.mg5, name, mg5root, verbose=args.verbose)
        raise SystemExit(0)

    # Filter operators if --operators was given
    selected = params
    if args.operators:
        requested = set(args.operators)
        selected = [p for p in params if p[1] in requested]
        missing = requested - {p[1] for p in selected}
        if missing:
            print(f'Warning: unknown operators ignored: {", ".join(sorted(missing))}')
        if not selected:
            raise SystemExit('No matching operators found.')

    results = []
    for i, param in enumerate(selected):
        name = param[1]
        print(f'[{i+1:3d}/{len(selected)}] {name} ... ', end='', flush=True)
        debug_this = args.debug and i == 0
        n_diag, stdout = run_mg5(args.mg5, name, mg5root, debug=debug_this, verbose=args.verbose)
        if args.debug and i == 0:
            raise SystemExit('Debug mode: exiting after first operator.')
        status = f'{n_diag} diagram(s)' if n_diag > 0 else 'no diagrams'
        print(status)
        results.append((name, n_diag, stdout))

    write_html(results, args.out)
    print(f'\nDone. Results written to {args.out}')

    print('\nContributing operators:')
    for name, n_diag, _ in results:
        if n_diag > 0:
            print(f'  {name}: {n_diag} diagram(s)')
