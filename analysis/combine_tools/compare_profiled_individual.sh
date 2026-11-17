#!/bin/bash
# Plot profiled vs individual 1D scans (mkEFTScan) and print 68%/95% intervals.
# Usage (in the profiled dir): compare_profiled_individual.sh OPS INDIVIDUAL_DIR [LABEL]
#   INDIVIDUAL_DIR = dir with higgsCombine.<op>.individual.MultiDimFit.mH125.root (1D, same datacard)
OPS=$1; IND=$2; LABEL=${3:-"Profiled"}; N=${OPS//,/_}
source ~/.bashrc >/dev/null 2>&1; dy_combine_morphing >/dev/null 2>&1
mkdir -p plots
for op in ${OPS//,/ }; do
  ln -sf $IND/higgsCombine.$op.individual.MultiDimFit.mH125.root ind_$op.root
  mkEFTScan.py higgsCombine.${N}__${op}_profiled.individual.MultiDimFit.mH125.root -p k_$op -maxNLL 10 -lumi 59.74 -cms -preliminary \
    -o plots/scan_profiled_$op -ff png pdf -ml "$LABEL" --others "ind_$op.root:2:2:Individual" > /dev/null 2>&1
done
python3 - "$OPS" "$N" <<'PY'
import sys, ROOT, numpy as np
ROOT.gROOT.SetBatch(True)
ops, N = sys.argv[1].split(","), sys.argv[2]
def iv(fn, op):
    d = ROOT.RDataFrame("limit", fn).Filter("quantileExpected>-0.5").AsNumpy([f"k_{op}", "deltaNLL"])
    x, z = d[f"k_{op}"], 2 * d["deltaNLL"]; z = z - z.min()
    return [(x[z < l].min(), x[z < l].max()) for l in (1.0, 3.84)], len(x), x[np.argmin(z)]
f = lambda t: f"[{t[0]:+.4f},{t[1]:+.4f}]"
print(f"{'op':6s} {'individual 68%':>20s} {'profiled 68%':>20s} {'ratio':>6s} | {'individual 95%':>20s} {'profiled 95%':>20s} {'ratio':>6s} | pts  best-fit")
for op in ops:
    (i1, i2), _, _ = iv(f"ind_{op}.root", op)
    (p1, p2), n, bf = iv(f"higgsCombine.{N}__{op}_profiled.individual.MultiDimFit.mH125.root", op)
    r = lambda p, i: (p[1] - p[0]) / (i[1] - i[0])
    print(f"{op:6s} {f(i1):>20s} {f(p1):>20s} {r(p1,i1):6.2f} | {f(i2):>20s} {f(p2):>20s} {r(p2,i2):6.2f} | {n:4d} {bf:+.4f}")
PY
