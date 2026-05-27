#!/usr/bin/env python3
"""
make_restriction_table.py
=========================
Generate a PDF of the SMEFT operator restriction table.

Green  = operator included in the DY SMEFT scan (count > 0)
Red    = operator not included / set to zero

Usage:
    python3 make_restriction_table.py [--out operator_table.pdf]
"""

import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Full SMEFTsim operator list with event counts from the gridpack scan.
# count > 0  → green (included in analysis)
# count == 0 → red   (not included / restricted to zero)
operators = [
    ("cG",0),("cW",0),("cH",0),("cHbox",0),("cHDD",72),("cHG",0),
    ("cHW",0),("cHB",0),("cHWB",48),("cuHRe",0),("ctHRe",0),("cdHRe",0),
    ("cbHRe",0),("cuGRe",0),("ctGRe",0),("cuWRe",0),("ctWRe",0),("cuBRe",0),
    ("ctBRe",0),("cdGRe",0),("cbGRe",0),("cdWRe",0),("cbWRe",0),("cdBRe",0),
    ("cbBRe",0),("cHj1",12),("cHQ1",0),("cHj3",12),("cHQ3",0),("cHu",6),
    ("cHt",0),("cHd",6),("cHbq",0),("cHudRe",0),("cHtbRe",0),("cjj11",0),
    ("cjj18",0),("cjj31",0),("cjj38",0),("cQj11",0),("cQj18",0),("cQj31",0),
    ("cQj38",0),("cQQ1",0),("cQQ8",0),("cuu1",0),("cuu8",0),("ctt",0),
    ("ctu1",0),("ctu8",0),("cdd1",0),("cdd8",0),("cbb",0),("cbd1",0),
    ("cbd8",0),("cud1",0),("ctb1",0),("ctd1",0),("cbu1",0),("cud8",0),
    ("ctb8",0),("ctd8",0),("cbu8",0),("cutbd1Re",0),("cutbd8Re",0),("cju1",0),
    ("cQu1",0),("cju8",0),("cQu8",0),("ctj1",0),("ctj8",0),("cQt1",0),
    ("cQt8",0),("cjd1",0),("cjd8",0),("cQd1",0),("cQd8",0),("cbj1",0),
    ("cbj8",0),("cQb1",0),("cQb8",0),("cjQtu1Re",0),("cjQtu8Re",0),
    ("cjQbd1Re",0),("cjQbd8Re",0),("cjujd1Re",0),("cjujd8Re",0),
    ("cjujd11Re",0),("cjujd81Re",0),("cQtjd1Re",0),("cQtjd8Re",0),
    ("cjuQb1Re",0),("cjuQb8Re",0),("cQujb1Re",0),("cQujb8Re",0),
    ("cjtQd1Re",0),("cjtQd8Re",0),("cQtQb1Re",0),("cQtQb8Re",0),
    ("ceHRe",0),("ceWRe",0),("ceBRe",0),("cHl1",12),("cHl3",48),("cHe",12),
    ("cll",0),("cll1",60),("clj1",12),("clj3",12),("cQl1",12),("cQl3",12),
    ("cee",0),("ceu",6),("cte",0),("ced",6),("cbe",6),("cje",12),
    ("cQe",12),("clu",6),("ctl",0),("cld",6),("cbl",6),("cle",0),
    ("cledjRe",0),("clebQRe",0),("cleju1Re",0),("cleQt1Re",0),
    ("cleju3Re",0),("cleQt3Re",0),
]

N_COLS = 13
N_ROWS = (len(operators) + N_COLS - 1) // N_COLS

GREEN = "#c8f0d0"
RED   = "#f5c6cb"


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", default="operator_table.pdf",
                        help="Output filename (default: operator_table.pdf)")
    args = parser.parse_args()

    fig, ax = plt.subplots(figsize=(10.5, 7.5))
    fig.patch.set_facecolor("white")
    fig.patch.set_edgecolor("none")
    ax.set_facecolor("white")
    ax.set_xlim(0, N_COLS)
    ax.set_ylim(0, N_ROWS)
    ax.axis("off")

    for idx, (name, count) in enumerate(operators):
        col = idx % N_COLS
        row = N_ROWS - 1 - (idx // N_COLS)
        color = GREEN if count > 0 else RED

        ax.add_patch(plt.Rectangle((col, row), 1, 1,
                                   facecolor=color, edgecolor="white", linewidth=1.2))
        ax.text(col + 0.5, row + 0.5, name,
                ha="center", va="center",
                fontsize=7.5, fontfamily="monospace", color="#1a1a1a")

    # Legend
    legend_patches = [
        mpatches.Patch(facecolor=GREEN, edgecolor="#aaaaaa", label="Included in DY SMEFT scan"),
        mpatches.Patch(facecolor=RED,   edgecolor="#aaaaaa", label="Set to zero (not constrained)"),
    ]
    ax.legend(handles=legend_patches, loc="lower right", fontsize=8,
              frameon=True, framealpha=0.9, edgecolor="#cccccc")

    outpdf = args.out
    outpng = outpdf.replace(".pdf", ".png")
    fig.savefig(outpdf, format="pdf", bbox_inches="tight", edgecolor="none")
    fig.savefig(outpng, dpi=150, bbox_inches="tight")
    print(f"Saved -> {outpdf}")
    print(f"Saved -> {outpng}")


if __name__ == "__main__":
    main()
