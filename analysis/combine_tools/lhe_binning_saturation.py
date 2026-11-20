#!/usr/bin/env python3
"""
How the information of each 1D variable saturates with the number of bins N (LHE level, same inputs and caveats as
lhe_binning_check.py): stat-only linear expected 1 sigma per operator (lin^2 corrected for the LHE MC variance)
vs N, relative to the v3 triple_diff (280 bins), plus the worst-bin LHE MC stat of the SM.

Binnings for a given N:
  mll          : equal numbers of LHE events per bin (fine where the samples are dense: Z peak, slice edges),
                 with a minimum width of --min-width GeV (resolution floor)
  costhetastar : uniform in [-1, 1]
  rapll_abs    : uniform in [0, 2.4]

Usage: lhe_binning_saturation.py [-o outdir] [--Ns 10,20,40,80,140,280]
"""
import argparse, os, pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep

p = argparse.ArgumentParser()
p.add_argument("--cache", default="/grid_mnt/data__data.polcms/cms/adufour/LHE/propcorr/CACHE/lhe_cache_propcorr_parallel.pkl")
p.add_argument("-o", "--outdir", default="lhe_saturation")
p.add_argument("--Ns", default="5,10,20,40,80,140,280")
p.add_argument("--min-width", type=float, default=1.0)
a = p.parse_args()
os.makedirs(a.outdir, exist_ok=True)
plt.style.use(hep.style.CMS)

NGEN, LUMI = 100_000, 59740.0
c = pickle.load(open(a.cache, "rb"))
mll, y, cs = c["mll"], np.abs(c["rap"]), c["cstar"]
wsm = c["w_SM"] / NGEN * LUMI
OPS = [op for op in c["w_p1"] if op not in ("cHQ3", "cQl3")]  # identical to cHQ1, cQl1
lin = {op: 0.5 * (c["w_p1"][op] - c["w_m1"][op]) / NGEN * LUMI for op in OPS}
del c


def sigma(i, n):
    ok = (i >= 0) & (i < n)
    sm = np.bincount(i[ok], weights=wsm[ok], minlength=n)
    mc = np.sqrt(np.bincount(i[ok], weights=wsm[ok] ** 2, minlength=n)) / np.where(sm > 0, sm, 1)
    out = {}
    for op in OPS:
        L = np.bincount(i[ok], weights=lin[op][ok], minlength=n)
        V = np.bincount(i[ok], weights=lin[op][ok] ** 2, minlength=n)
        F = np.sum(np.clip(L ** 2 - V, 0, None) / np.where(sm > 0, sm, np.inf))
        out[op] = 1 / np.sqrt(F) if F > 0 else np.inf
    return out, float(mc.max()), float(sm.min())


def mll_edges(N):
    q = np.quantile(mll, np.linspace(0, 1, N + 1))
    q[0], q[-1] = 50.0, 3000.0
    e = [q[0]]
    for x in q[1:-1]:
        if x - e[-1] >= a.min_width:
            e.append(x)
    if 3000.0 - e[-1] < a.min_width:
        e.pop()
    return np.array(e + [3000.0])


# reference: v3 triple_diff (unroll order irrelevant for the sums)
em = np.array([50, 60, 80, 100, 120, 140, 180, 220, 270, 350, 500, 700, 1000, 1500, 3000.])
ec = np.array([-1.0, -0.6, -0.2, 0.2, 0.6, 1.0]); ey = np.array([0.0, 0.48, 0.96, 1.44, 2.4])
im, ic, iy = np.digitize(mll, em) - 1, np.digitize(cs, ec) - 1, np.digitize(y, ey) - 1
ok = (im >= 0) & (im < 14) & (ic >= 0) & (ic < 5) & (iy >= 0) & (iy < 4)
ref, _, _ = sigma(np.where(ok, (iy * 5 + ic) * 14 + im, -1), 280)

Ns = [int(x) for x in a.Ns.split(",")]
res = {}
for var in ("mll", "costhetastar", "rapll_abs"):
    for N in Ns:
        if var == "mll":
            e = mll_edges(N); v = mll
        elif var == "costhetastar":
            e = np.linspace(-1, 1, N + 1); v = cs
        else:
            e = np.linspace(0, 2.4, N + 1); v = y
        s, mcmax, smmin = sigma(np.digitize(v, e) - 1, len(e) - 1)
        res[(var, N)] = dict(s=s, mc=mcmax, smmin=smmin, n=len(e) - 1,
                             wmin=float(np.diff(e).min()))
        print(f"{var:13s} N={N:4d} (actual {len(e)-1:4d}) MC stat max {100*mcmax:5.1f}%  min SM/bin {smmin:10.1f}"
              f"  min width {np.diff(e).min():.3g}")

# ---- plot: sigma(N) / sigma(v3 triple_diff) per operator, one panel per variable
fig, axs = plt.subplots(1, 3, figsize=(24, 8), sharey=True)
cmap = plt.get_cmap("tab20")
for ax, var in zip(axs, ("mll", "costhetastar", "rapll_abs")):
    for k, op in enumerate(OPS):
        r = [res[(var, N)]["s"][op] / ref[op] for N in Ns]
        ax.plot([res[(var, N)]["n"] for N in Ns], r, marker="o", ms=4, color=cmap(k % 20),
                ls="-" if k < 20 else "--", label=op)
    ax.axhline(1, color="k", lw=1)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel(f"number of {var} bins")
    ax.set_title(var, fontsize=22)
axs[0].set_ylabel(r"$\sigma_{1D}(N)\ /\ \sigma$(v3 triple_diff, 280 bins)")
axs[2].legend(ncol=2, fontsize=11, loc="upper right")
fig.savefig(f"{a.outdir}/saturation.png", dpi=90, bbox_inches="tight")

# ---- table: sigma(N)/sigma(N=max) to show where each variable saturates (median over operators) + MC stat
lines = ["| var | " + " | ".join(f"N={N}" for N in Ns) + " |", "|---|" + "---|" * len(Ns)]
for var in ("mll", "costhetastar", "rapll_abs"):
    med = [np.median([res[(var, N)]["s"][op] / res[(var, Ns[-1])]["s"][op] for op in OPS]) for N in Ns]
    lines.append(f"| {var}: median sigma(N)/sigma({Ns[-1]}) | " + " | ".join(f"{m:.2f}" for m in med) + " |")
    lines.append(f"| {var}: worst MC stat | " + " | ".join(f"{100*res[(var, N)]['mc']:.0f}%" for N in Ns) + " |")
    lines.append(f"| {var}: min SM events/bin | " + " | ".join(f"{res[(var, N)]['smmin']:.3g}" for N in Ns) + " |")
lines.append("")
lines.append("best 1D at N=max vs v3 triple_diff (280), sigma ratio:")
lines.append("| op | " + " | ".join(("mll", "costhetastar", "rapll_abs")) + " |"); lines.append("|---|---|---|---|")
for op in OPS:
    lines.append(f"| {op} | " + " | ".join(f"{res[(v, Ns[-1])]['s'][op] / ref[op]:.2f}" for v in ("mll", "costhetastar", "rapll_abs")) + " |")
open(f"{a.outdir}/saturation.md", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
