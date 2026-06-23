#!/usr/bin/env python3
"""
print_syst_table.py

Print QCD scale and PDF systematic uncertainties as % of SM yield per mll bin.

Usage:
    python print_syst_table.py --input histograms.root [--channel triple_DY]
"""

import argparse
import numpy as np
import uproot

parser = argparse.ArgumentParser()
parser.add_argument("--input",   required=True)
parser.add_argument("--channel", default="triple_DY")
args = parser.parse_args()

with uproot.open(args.input) as f:
    def get(name):
        key = f"{args.channel}/{name}"
        return f[key].to_numpy()[0].copy() if key in f else None

    sm      = get("sm")
    qcd_up  = get("sm_qcd_scaleUp")
    qcd_dn  = get("sm_qcd_scaleDown")
    pdf_up  = get("sm_pdfUp")
    pdf_dn  = get("sm_pdfDown")
    _, edges = f[f"{args.channel}/sm"].to_numpy()

if sm is None:
    raise RuntimeError("sm histogram not found")

bin_labels = [f"[{int(lo)},{int(hi)}]" for lo, hi in zip(edges[:-1], edges[1:])]

# asymmetric up/down as % of SM
qcd_up_pct = np.where(sm > 0, (qcd_up - sm) / sm * 100, 0.) if qcd_up is not None else np.zeros_like(sm)
qcd_dn_pct = np.where(sm > 0, (sm - qcd_dn) / sm * 100, 0.) if qcd_dn is not None else np.zeros_like(sm)
pdf_up_pct = np.where(sm > 0, (pdf_up - sm) / sm * 100, 0.) if pdf_up is not None else np.zeros_like(sm)
pdf_dn_pct = np.where(sm > 0, (sm - pdf_dn) / sm * 100, 0.) if pdf_dn is not None else np.zeros_like(sm)

tot_up_pct = np.sqrt(qcd_up_pct**2 + pdf_up_pct**2)
tot_dn_pct = np.sqrt(qcd_dn_pct**2 + pdf_dn_pct**2)

print(f"\nSystematic uncertainties as % of SM yield  [{args.channel}]\n")
header = f"{'mll bin [GeV]':<18}  {'N_SM':>10}  {'QCD up%':>8}  {'QCD dn%':>8}  {'PDF up%':>8}  {'PDF dn%':>8}  {'Tot up%':>8}  {'Tot dn%':>8}"
print(header)
print("-" * len(header))
for i, label in enumerate(bin_labels):
    print(f"{label:<18}  {sm[i]:>10.1f}  {qcd_up_pct[i]:>7.2f}%  {qcd_dn_pct[i]:>7.2f}%  "
          f"{pdf_up_pct[i]:>7.2f}%  {pdf_dn_pct[i]:>7.2f}%  "
          f"{tot_up_pct[i]:>7.2f}%  {tot_dn_pct[i]:>7.2f}%")
print()
