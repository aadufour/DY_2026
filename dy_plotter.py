"""
dy_plotter.py

General-purpose 1D histogram plotter for Drell-Yan LHE files.

Define which variables to plot in the VARIABLES list at the bottom of this
file.  Each variable is a lambda over a DYEvent, so new observables are a
one-liner.

Usage
-----
  python dy_plotter.py events.lhe
  python dy_plotter.py events.lhe --weight SM       # pick a named weight
  python dy_plotter.py events.lhe --out plots/      # output directory
  python dy_plotter.py events.lhe --max-events 5000

Observables available on DYEvent
---------------------------------
  ev.l1 / ev.l2          : leading / sub-leading lepton (FourVector, by pT)
  ev.lm / ev.lp          : negative / positive lepton
  ev.Z                   : dilepton system (FourVector)

FourVector properties
---------------------
  .pt  .eta  .phi  .mass  .rap  .e  .px  .py  .pz

DYEvent computed properties
----------------------------
  ev.mll           : dilepton invariant mass
  ev.yll           : dilepton rapidity
  ev.cos_theta_star: cos(theta*) — angle of l- vs beam in dilepton rest frame
  ev.delta_r       : DeltaR(l1, l2)
  ev.delta_phi     : |Delta phi(l1, l2)|  in [0, pi]
  ev.delta_eta     : |Delta eta(l1, l2)|
"""

import argparse
import math
import os
import warnings

import boost_histogram as bh
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import pylhe


# ── Four-vector ───────────────────────────────────────────────────────────────

class FourVector:
    """Minimal 4-vector with kinematic properties. Convention: (px, py, pz, E)."""

    __slots__ = ("px", "py", "pz", "e")

    def __init__(self, px, py, pz, e):
        self.px = float(px)
        self.py = float(py)
        self.pz = float(pz)
        self.e  = float(e)

    # ── arithmetic ────────────────────────────────────────────────────────────

    def __add__(self, other):
        return FourVector(
            self.px + other.px,
            self.py + other.py,
            self.pz + other.pz,
            self.e  + other.e,
        )

    # ── transverse ────────────────────────────────────────────────────────────

    @property
    def pt(self):
        return math.sqrt(self.px**2 + self.py**2)

    @property
    def phi(self):
        return math.atan2(self.py, self.px)

    # ── longitudinal / angular ────────────────────────────────────────────────

    @property
    def eta(self):
        p3 = math.sqrt(self.px**2 + self.py**2 + self.pz**2)
        if p3 == abs(self.pz):          # on-axis protection
            return math.copysign(1e10, self.pz)
        return 0.5 * math.log((p3 + self.pz) / (p3 - self.pz))

    @property
    def rap(self):
        return 0.5 * math.log((self.e + self.pz) / (self.e - self.pz))

    # ── invariant mass ────────────────────────────────────────────────────────

    @property
    def mass(self):
        m2 = self.e**2 - self.px**2 - self.py**2 - self.pz**2
        return math.sqrt(max(m2, 0.0))


# ── DY Event ──────────────────────────────────────────────────────────────────

class DYEvent:
    """
    Wraps one pylhe event and exposes Drell-Yan observables.

    Identifies the two final-state leptons and the dilepton system Z.
    """

    def __init__(self, lm: FourVector, lp: FourVector, weight: float):
        """
        Parameters
        ----------
        lm : FourVector  : negatively charged lepton  
        lp : FourVector  : positively charged lepton
        weight : float
        """
        self._lm     = lm
        self._lp     = lp
        self._Z      = lm + lp
        self._weight = weight

        # cache l1/l2 assignment by pT
        if lm.pt >= lp.pt:
            self._l1, self._l2 = lm, lp
        else:
            self._l1, self._l2 = lp, lm

    # ── particles ─────────────────────────────────────────────────────────────

    @property
    def lm(self) -> FourVector:
        """Negatively charged lepton."""
        return self._lm

    @property
    def lp(self) -> FourVector:
        """Positively charged lepton."""
        return self._lp

    @property
    def l1(self) -> FourVector:
        """Leading-pT lepton."""
        return self._l1

    @property
    def l2(self) -> FourVector:
        """Sub-leading-pT lepton."""
        return self._l2

    @property
    def Z(self) -> FourVector:
        """Dilepton 4-vector."""
        return self._Z

    @property
    def weight(self) -> float:
        return self._weight

    # ── dilepton observables ──────────────────────────────────────────────────

    @property
    def mll(self) -> float:
        return self._Z.mass

    @property
    def yll(self) -> float:
        return self._Z.rap

    @property
    def cos_theta_star(self) -> float:
        """
        cos(theta*): angle of l- relative to z in the
        dilepton rest frame, via a longitudinal boost.
        """
        Z = self._Z
        lm = self._lm
        mass  = Z.mass
        if mass == 0.0:
            return 0.0
        beta  = Z.pz / Z.e
        gamma = Z.e  / mass
        # boost lm along z into the dilepton rest frame
        pz_boosted = gamma * (lm.pz - beta * lm.e)
        p_mag      = math.sqrt(lm.px**2 + lm.py**2 + pz_boosted**2)
        if p_mag == 0.0:
            return 0.0
        return pz_boosted / p_mag

    # ── pair topology ─────────────────────────────────────────────────────────

    @property
    def delta_phi(self) -> float:
        """|DeltaPhi(l1, l2), mapped to [0, pi]."""
        dphi = abs(self._l1.phi - self._l2.phi)
        if dphi > math.pi:
            dphi = 2 * math.pi - dphi
        return dphi

    @property
    def delta_eta(self) -> float:
        return abs(self._l1.eta - self._l2.eta)

    @property
    def delta_r(self) -> float:
        return math.sqrt(self.delta_eta**2 + self.delta_phi**2)


# ── Variable definition ───────────────────────────────────────────────────────

class Variable:
    """
    Describes one 1D histogram.

    Parameters
    ----------
    name   : short identifier used for the output filename
    fn     : callable  DYEvent -> float
    n_bins : number of bins
    x_range: (xmin, xmax) tuple
    xlabel : xlabel
    """
    def __init__(self, name, fn, n_bins, x_range, xlabel):
        self.name   = name
        self.fn     = fn
        self.n_bins = n_bins
        self.x_range = x_range
        self.xlabel = xlabel


# ── Variables to plot ─────────────────────────────────────────────────────────
# Add / remove / edit entries here.

VARIABLES = [
    # ── leading lepton ────────────────────────────────────────────────────────
    Variable("pt_l1",    lambda ev: ev.l1.pt,           50, (0,   300),  r"$p_T^{\ell_1}$ [GeV]"),
    Variable("eta_l1",   lambda ev: ev.l1.eta,          50, (-5,  5),    r"$\eta_{\ell_1}$"),
    Variable("phi_l1",   lambda ev: ev.l1.phi,          50, (-math.pi, math.pi), r"$\phi_{\ell_1}$"),

    # ── sub-leading lepton ────────────────────────────────────────────────────
    Variable("pt_l2",    lambda ev: ev.l2.pt,           50, (0,   300),  r"$p_T^{\ell_2}$ [GeV]"),
    Variable("eta_l2",   lambda ev: ev.l2.eta,          50, (-5,  5),    r"$\eta_{\ell_2}$"),
    Variable("phi_l2",   lambda ev: ev.l2.phi,          50, (-math.pi, math.pi), r"$\phi_{\ell_2}$"),

    # ── Z / dilepton ──────────────────────────────────────────────────────────
    Variable("pt_Z",     lambda ev: ev.Z.pt,            50, (0,   300),  r"$p_T^Z$ [GeV]"),
    Variable("mll",      lambda ev: ev.mll,             50, (50,  200),  r"$m_{\ell\ell}$ [GeV]"),
    Variable("yll",      lambda ev: ev.yll,             50, (-5,  5),    r"$y_{\ell\ell}$"),
    Variable("cos_star", lambda ev: ev.cos_theta_star,  50, (-1,  1),    r"$\cos\theta^*$"),

    # ── pair topology ─────────────────────────────────────────────────────────
    Variable("delta_r",  lambda ev: ev.delta_r,         50, (0,   6),    r"$\Delta R(\ell_1,\ell_2)$"),
    Variable("delta_phi",lambda ev: ev.delta_phi,       50, (0,   math.pi), r"$|\Delta\phi(\ell_1,\ell_2)|$"),
    Variable("delta_eta",lambda ev: ev.delta_eta,       50, (0,   10),   r"$|\Delta\eta(\ell_1,\ell_2)|$"),
]


# ── LHE parsing ───────────────────────────────────────────────────────────────

def parse_lhe(lhe_path: str, weight_key: str | None, max_events: int | None):
    """
    Read a LHE file and return a list of DYEvent objects.

    weight_key : if None, falls back to the event-level weight (eventinfo.weight).
    """
    events = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        lhe_events = pylhe.read_lhe_with_attributes(lhe_path)

        for i, event in enumerate(lhe_events):
            if max_events is not None and i >= max_events:
                break
            if (i + 1) % 10_000 == 0:
                print(f"  {i + 1} events read…")

            # identify final-state leptons
            leptons = [
                p for p in event.particles
                if int(p.status) == 1 and abs(int(p.id)) in {11, 13}
            ]
            if len(leptons) < 2:
                continue

            # separate by charge: positive PDG id → negative lepton (e-=11, μ-=13)
            lm_raw = next((p for p in leptons if int(p.id) > 0), leptons[0])
            lp_raw = next((p for p in leptons if int(p.id) < 0), leptons[1])

            lm = FourVector(lm_raw.px, lm_raw.py, lm_raw.pz, lm_raw.e)
            lp = FourVector(lp_raw.px, lp_raw.py, lp_raw.pz, lp_raw.e)

            # weight selection
            if weight_key is not None and event.weights:
                w = float(event.weights[weight_key])
            else:
                w = float(event.eventinfo.weight)

            events.append(DYEvent(lm, lp, w))

    print(f"  {len(events)} events parsed from {lhe_path}")
    return events


# ── Histogram filling and plotting ────────────────────────────────────────────

def fill_and_plot(dy_events, variables, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    hep.style.use("CMS")

    weights = np.array([ev.weight for ev in dy_events])

    for var in variables:
        # evaluate observable for all events
        values = np.array([var.fn(ev) for ev in dy_events])

        # fill histogram
        h = bh.Histogram(
            bh.axis.Regular(var.n_bins, var.x_range[0], var.x_range[1]),
            storage=bh.storage.Weight(),
        )
        h.fill(values, weight=weights)

        # plot
        fig, ax = plt.subplots(figsize=(7, 5))
        hep.histplot(h, ax=ax, histtype="step", color="steelblue", linewidth=1.5)
        ax.set_xlabel(var.xlabel, fontsize=13)
        ax.set_ylabel("Weighted events", fontsize=13)
        ax.semilogy()
        hep.cms.label(ax=ax, llabel="", data=True, lumi=None)

        out_path = os.path.join(out_dir, f"{var.name}.png")
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  saved {out_path}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Plot 1D DY observables from a LHE file."
    )
    parser.add_argument("lhe", help="Path to input LHE file")
    parser.add_argument(
        "--weight", default=None,
        help="Named weight key to use (e.g. 'SM'). Defaults to event-level weight.",
    )
    parser.add_argument(
        "--out", default="dy_plots",
        help="Output directory for PNG files (default: dy_plots/)",
    )
    parser.add_argument(
        "--max-events", type=int, default=None,
        help="Stop after this many events (useful for quick tests).",
    )
    args = parser.parse_args()

    print(f"Reading {args.lhe}")
    dy_events = parse_lhe(args.lhe, args.weight, args.max_events)

    if not dy_events:
        print("No events found — check lepton PDG IDs or file path.")
        return

    print(f"\nFilling and plotting {len(VARIABLES)} variables…")
    fill_and_plot(dy_events, VARIABLES, args.out)
    print("\nDone.")


if __name__ == "__main__":
    main()
