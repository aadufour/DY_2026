#!/usr/bin/env python3
"""
plot_pdf_event_weights.py

Per-event PDF weights: for selected events (run:lumi:event), plot LHEPdfWeight[i]/LHEPdfWeight[0]
for the 100 Hessian members, against the central 99% band of all events in the same reco mll range.

Usage (spritz-env.sif, valid proxy, X509_CERT_DIR set):
    python3 plot_pdf_event_weights.py --file root://.../file.root \
        --events 1:141:94231 1:219:146025 1:130:86524 --mllmin 1500 --mllmax 3000 --outdir .
"""
import argparse
import numpy as np
import awkward as ak
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep

plt.style.use(hep.style.CMS)

p = argparse.ArgumentParser()
p.add_argument("--file", required=True)
p.add_argument("--events", nargs="+", required=True, help="run:lumi:event")
p.add_argument("--mllmin", type=float, default=1500)
p.add_argument("--mllmax", type=float, default=3000)
p.add_argument("--label", default="DYJetsToMuMu_M-2000toInf (MiNNLO, UL18)")
p.add_argument("--outdir", default=".")
a = p.parse_args()

t = uproot.open(a.file, timeout=120)["Events"]
ev = t.arrays(["run", "luminosityBlock", "event", "LHEPdfWeight", "Muon_pt", "Muon_eta", "Muon_phi"])

# reco mll from the two leading muons (quick selection, not the analysis one)
mu = ak.argsort(ev.Muon_pt, ascending=False)[:, :2]
pt, eta, phi = ev.Muon_pt[mu], ev.Muon_eta[mu], ev.Muon_phi[mu]
px, py, pz = ak.sum(pt * np.cos(phi), 1), ak.sum(pt * np.sin(phi), 1), ak.sum(pt * np.sinh(eta), 1)
E = ak.sum(pt * np.cosh(eta), 1)
mll = np.where(ak.to_numpy(ak.num(ev.Muon_pt) >= 2), np.sqrt(np.maximum(ak.to_numpy(E**2 - px**2 - py**2 - pz**2), 0)), -1)

W = ak.to_numpy(ev.LHEPdfWeight)
R = W[:, 1:101] / W[:, [0]]                          # Hessian members 1..100
sel = (mll > a.mllmin) & (mll < a.mllmax)
lo, hi = np.quantile(R[sel], [0.005, 0.995], axis=0)
idx = np.arange(1, 101)

fig, axes = plt.subplots(len(a.events), 1, figsize=(12, 3.8 * len(a.events)), sharex=True)
axes = np.atleast_1d(axes)
for ax, tag in zip(axes, a.events):
    run, lumi, evt = map(int, tag.split(":"))
    j = np.where((ak.to_numpy(ev.run) == run) & (ak.to_numpy(ev.luminosityBlock) == lumi) & (ak.to_numpy(ev.event) == evt))[0]
    if len(j) == 0:
        print("event not found", tag)
        continue
    r = R[j[0]]
    ax.fill_between(idx, np.abs(lo), np.abs(hi), step="mid", color="lightgray",
                    label=f"central 99% of events, {a.mllmin:.0f} < reco $m_{{\\ell\\ell}}$ < {a.mllmax:.0f} GeV")
    ax.step(idx, np.abs(r), where="mid", color="black", lw=1.5, label=f"event {tag}")
    k = int(np.argmax(np.abs(r - 1)))
    ax.annotate(f"member {k + 1}: ratio {r[k]:+.3g}", (k + 1, abs(r[k])), xytext=(10, -5),
                textcoords="offset points", fontsize=14, color="tab:red", va="top")
    ax.set_yscale("log")
    ax.set_ylim(1e-2, 1e6)
    ax.set_ylabel("|PDF member / nominal|", fontsize=16)
    ax.legend(fontsize=12, loc="upper left")
    ax.text(0.98, 0.95, f"reco $m_{{\\ell\\ell}}$ = {mll[j[0]]:.0f} GeV", transform=ax.transAxes,
            ha="right", va="top", fontsize=14)
    print(tag, "member", k + 1, "ratio", r[k], "| other members max |r-1|:", np.sort(np.abs(r - 1))[-2])
axes[-1].set_xlabel("Hessian member index")
hep.cms.label("Preliminary", data=False, ax=axes[0], fontsize=18)
axes[0].set_title(a.label, fontsize=14, loc="center", pad=28)
o = f"{a.outdir}/pdf_event_weights"
for ext in ("pdf", "png"):
    fig.savefig(f"{o}.{ext}", bbox_inches="tight")
print("saved", o)
