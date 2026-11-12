#!/usr/bin/env python3
"""
plot_pdf_members.py

PDF nuisance plot per Hessian member: relative variation (%) of each PDFWeight_i (i=1..100)
w.r.t. nominal, per bin, for one process in a spritz histos.root. Outlier members are highlighted.

Usage (dy_combine_morphing env, from the config dir):
    python3 plot_pdf_members.py --file histos.root --var mll --proc DYll [--outdir .]
"""
import argparse
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep

plt.style.use(hep.style.CMS)

p = argparse.ArgumentParser()
p.add_argument("--file", default="histos.root")
p.add_argument("--region", default="inc_mm")
p.add_argument("--var", default="mll")
p.add_argument("--proc", default="DYll")
p.add_argument("--ntop", type=int, default=3, help="highlight the ntop members with the largest |dev|/median in any bin")
p.add_argument("--members", default=None, help="comma-separated members to highlight (overrides --ntop)")
p.add_argument("--labels", default=None, help="semicolon-separated legend text per --members entry, e.g. event ids")
p.add_argument("--tag", default="", help="suffix for the output file name")
p.add_argument("--outdir", default=".")
p.add_argument("--label", default=None, help="extra bold text on the plot")
a = p.parse_args()

f = uproot.open(a.file)
base = f"{a.region}/{a.var}/histo_{a.proc}"
h = f[base]
n = h.values()
edges = h.axis().edges()
M = np.array([f[f"{base}_PDFWeight_{i}"].values() for i in range(1, 101)])
safe = np.where(n != 0, n, 1.0)
r = 100 * (M / safe - 1)                      # (100 members, nbins)
quad = np.sqrt((r**2).sum(0))
med = np.median(np.abs(r), axis=0)
score = (np.abs(r) / np.maximum(med, 0.1)).max(axis=1)
if a.members:
    out = [int(x) for x in a.members.split(",")]
else:
    out = sorted(int(i) + 1 for i in np.argsort(-score)[:a.ntop])
labs = a.labels.split(";") if a.labels else [""] * len(out)

x = edges if a.var != "triple_diff" else np.arange(len(n) + 1)
fig, ax = plt.subplots(figsize=(12, 7))
first = True
for i in range(100):
    if i + 1 not in out:
        hep.histplot(r[i], x, ax=ax, color="lightgray", lw=0.8, label="other PDF members" if first else None)
        first = False
for k, m in enumerate(out):
    hep.histplot(r[m - 1], x, ax=ax, color=f"C{k}", lw=2, label=f"PDFWeight_{m}" + (f" ({labs[k]})" if labs[k] else ""))
hep.histplot(quad, x, ax=ax, color="black", ls="--", lw=2, label="quadrature sum (nuisance)")
hep.histplot(-quad, x, ax=ax, color="black", ls="--", lw=2)
ax.axhline(0, color="gray", lw=0.8)
ax.set_ylabel(f"PDF member / nominal - 1 [%]")
ax.set_xlabel({"mll": r"$m_{\ell\ell}$ [GeV]", "triple_diff": "bin index"}.get(a.var, a.var))
if a.var == "mll":
    ax.set_xscale("log")
ax.text(0.03, 0.95, f"{a.proc}, {a.region}, Hessian members 1-100", transform=ax.transAxes, fontsize=16, va="top")
if a.label:
    ax.text(0.03, 0.87, a.label, transform=ax.transAxes, fontsize=20, fontweight="bold", va="top")
ax.legend(fontsize=14, loc="lower left")
hep.cms.label("Preliminary", data=False, ax=ax, fontsize=20)
o = f"{a.outdir}/pdf_members_{a.proc}_{a.var}{a.tag}"
for ext in ("pdf", "png"):
    fig.savefig(f"{o}.{ext}")
print("highlighted:", out)
print("last bin dev% of highlighted:", [round(float(r[m - 1, -1]), 1) for m in out], " median |dev|:", round(float(med[-1]), 2))
print("saved", o)
