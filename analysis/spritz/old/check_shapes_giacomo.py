#!/usr/bin/env python3

import uproot
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep

# --------------------------------------------------
# Configuration
# --------------------------------------------------

root_file = "histos.root"

for op in ['cqlm2', 'cql32', 'cqe2', 'cll1221', 'cpdc', 'cpwb', 'cpl2', 'c3pl1', 'c3pl2', 'cpmu', 'cpqmi', 'cpq3i', 'cpq3', 'cpqm', 'cpu', 'cpd']:
    hist_names = {
        "SM": "inc_mm/mll/nominal/histo_DYMuMu_sm",
        f"{op}": f"inc_mm/mll/nominal/histo_DYMuMu_{op}",
        f"{op}_m1": f"inc_mm/mll/nominal/histo_DYMuMu_{op}_m1",
        # Replace/add more if needed
        # "other": "path/to/hist"
    }

    # --------------------------------------------------
    # Style
    # --------------------------------------------------

    hep.style.use("CMS")

    # --------------------------------------------------
    # Read histograms
    # --------------------------------------------------

    file = uproot.open(root_file)

    hists = {}

    for label, path in hist_names.items():
        h = file[path]

        values = h.values()
        edges = h.axes[0].edges()

        hists[label] = {
            "values": values,
            "edges": edges,
        }

    # --------------------------------------------------
    # Reference histogram (SM)
    # --------------------------------------------------

    sm = hists["SM"]["values"]
    edges = hists["SM"]["edges"]

    # Bin centers for ratio
    centers = 0.5 * (edges[:-1] + edges[1:])

    # --------------------------------------------------
    # Figure with ratio panel
    # --------------------------------------------------

    fig, (ax, rax) = plt.subplots(
        2,
        1,
        figsize=(8, 8),
        gridspec_kw={"height_ratios": [3, 1]},
        sharex=True,
    )

    # --------------------------------------------------
    # Upper pad: shapes
    # --------------------------------------------------

    colors = {
        "SM": "black",
        f"{op}": "orange",
        f"{op}_m1": "blue",
    }

    for label, hist in hists.items():
        hep.histplot(
            hist["values"],
            hist["edges"],
            label=label,
            ax=ax,
            histtype="step",
            linewidth=2,
            color=colors.get(label, None),
        )

    ax.legend()
    ax.set_ylabel("Events")
    ax.set_yscale("log")

    hep.cms.label(ax=ax, data=False)

    # --------------------------------------------------
    # Ratio pad
    # --------------------------------------------------

    for label, hist in hists.items():

        if label == "SM":
            continue

        ratio = np.divide(
            hist["values"],
            sm,
            out=np.zeros_like(hist["values"], dtype=float),
            where=sm != 0,
        )

        rax.step(
            centers,
            ratio,
            where="mid",
            linewidth=2,
            label=f"{label}/SM",
            color=colors.get(label, None),
        )

    rax.axhline(1.0, color="black", linestyle="--")

    rax.set_ylabel("Ratio")
    rax.set_xlabel(r"$m_{\ell\ell}$")
    rax.set_ylim(0.5, 1.5)

    # --------------------------------------------------
    # Save
    # --------------------------------------------------

    plt.tight_layout()
    plt.savefig(f"check/mll_shapes_{op}.png", dpi=300)
    plt.savefig(f"check/mll_shapes_{op}.pdf")

    print("Saved plots:")
    print("  mll_shapes.png")
    print("  mll_shapes.pdf")