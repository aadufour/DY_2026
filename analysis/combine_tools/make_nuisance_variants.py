#!/usr/bin/env python3
"""
make_nuisance_variants.py
==========================
Write a single datacard.txt variant with a given set of nuisances removed —
by default the theory uncertainties (EXCLUDE_LIST: QCDScale, PDFweight).

Does NOT touch shapes.root or re-run make_cards.py/spritz-cards-eft: combine
only reads a histo_{proc}_{syst}Up/Down pair if the datacard's systematics
block has a row for it, so the source dir's shapes.root works unchanged for
the variant. The output dir just symlinks it (and, by default, metadata.json
and jsonComb.json — see --sidecars) alongside the filtered datacard.txt.

Datacard layout this relies on (as written by analysis/spritz/make_cards.py):
    ...
    ----------------------------------------------------------------------------------------------------   <- header separator (1st dash-only line)
    bin         ...
    observation ...
    shapes  * ...
    shapes  data_obs ...

    bin        ...
    process    ...
    process    ...
    rate       ...
    ----------------------------------------------------------------------------------------------------   <- systematics separator (2nd dash-only line)
    QCDScale   shape   1.0  1.0  -    ...
    PDFweight  shape   1.0  -    1.0  ...
    ...
    lumi       lnN     1.0084 ...
    inc_mm_mll autoMCStats 10 0 1
    ...

Everything up to and including the 2nd dash-only line is copied verbatim.
After that, each line's first whitespace-separated token is its nuisance
name (e.g. "QCDScale", "lumi"); a line whose name is in the exclude set is
dropped, everything else is kept as-is.

Usage
-----
List the nuisances found in the real datacard (dry run, no files written):
    python3 make_nuisance_variants.py --datacard datacards/inc_mm/mll/datacard.txt --list

Default: drop EXCLUDE_LIST (QCDScale, PDFweight):
    python3 make_nuisance_variants.py \\
        --datacard datacards/inc_mm/mll/datacard.txt \\
        --outdir   nuisance_scan/inc_mm_mll_no_theory

Custom exclude list:
    python3 make_nuisance_variants.py \\
        --datacard datacards/inc_mm/mll/datacard.txt \\
        --outdir   nuisance_scan/inc_mm_mll_no_lumi \\
        --exclude  lumi
"""

import argparse
import os


EXCLUDE_LIST = ["QCDScale", "PDFweight"]


def is_dash_line(line):
    s = line.strip()
    return len(s) > 0 and set(s) == {"-"}


def split_datacard(lines):
    """Return (header_lines, syst_lines) split after the 2nd dash-only line."""
    dash_idx = [i for i, line in enumerate(lines) if is_dash_line(line)]
    if len(dash_idx) < 2:
        raise ValueError(
            f"Expected >= 2 dash-separator lines in the datacard, found {len(dash_idx)}. "
            "Datacard format may have changed — check split_datacard()."
        )
    split_at = dash_idx[1]
    return lines[: split_at + 1], lines[split_at + 1 :]


def classify_syst_lines(syst_lines):
    """
    Returns (nuisance_names_in_order, entries) where entries is a list of
    (name_or_None, raw_line) — name is None for lines that are always kept
    (blank, autoMCStats, or unrecognized).
    """
    names = []
    entries = []
    for line in syst_lines:
        stripped = line.strip()
        if not stripped:
            entries.append((None, line))
            continue
        if "autoMCStats" in stripped:
            entries.append((None, line))
            continue
        name = stripped.split()[0]
        entries.append((name, line))
        if name not in names:
            names.append(name)
    return names, entries


def build_variant_lines(header_lines, entries, drop_set):
    kept = []
    for name, line in entries:
        if name is not None and name in drop_set:
            continue
        kept.append(line)
    return header_lines + kept


def write_variant(outdir, source_dir, header_lines, entries, drop_set, sidecars):
    os.makedirs(outdir, exist_ok=True)

    lines = build_variant_lines(header_lines, entries, drop_set)
    with open(os.path.join(outdir, "datacard.txt"), "w") as f:
        f.writelines(lines)

    # Symlink only the requested sidecar files (default: shapes.root) — the source
    # dir may already contain a full prior run's output (model_*.root, scans,
    # plots, ...); we don't want any of that pulled into the variant folder.
    for fname in sidecars:
        src = os.path.abspath(os.path.join(source_dir, fname))
        if not os.path.exists(src):
            print(f"  WARNING: sidecar '{fname}' not found in {source_dir}, skipping")
            continue
        dst = os.path.join(outdir, fname)
        if os.path.islink(dst) or os.path.exists(dst):
            os.remove(dst)
        os.symlink(os.path.relpath(src, outdir), dst)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--datacard", required=True, help="Path to the full datacard.txt (from spritz-cards-eft)")
    parser.add_argument("--outdir", help="Directory to write the variant into (required unless --list)")
    parser.add_argument("--exclude", default=None,
                        help=f"Comma-separated nuisance names to drop (default: {','.join(EXCLUDE_LIST)})")
    parser.add_argument("--sidecars", default="shapes.root,metadata.json,jsonComb.json",
                        help="Comma-separated filenames from the source dir to symlink into the variant folder "
                             "(default: shapes.root, metadata.json, jsonComb.json — deliberately NOT everything in "
                             "the source dir, which may already hold a prior run's model_*.root/scans/plots)")
    parser.add_argument("--list", action="store_true", help="Just print the nuisance names found in the datacard and exit")
    args = parser.parse_args()

    with open(args.datacard) as f:
        lines = f.readlines()
    header_lines, syst_lines = split_datacard(lines)
    master_names, entries = classify_syst_lines(syst_lines)

    if args.list:
        print(f"{len(master_names)} nuisances found in {args.datacard}:")
        for n in master_names:
            print(f"  {n}")
        return

    if not args.outdir:
        parser.error("--outdir is required unless --list is given")

    drop_set = set(s.strip() for s in args.exclude.split(",") if s.strip()) if args.exclude else set(EXCLUDE_LIST)
    unknown = drop_set - set(master_names)
    if unknown:
        print(f"  WARNING: exclude name(s) {sorted(unknown)} not found in {args.datacard}, ignoring")

    source_dir = os.path.dirname(os.path.abspath(args.datacard))
    sidecars = [s.strip() for s in args.sidecars.split(",") if s.strip()]

    write_variant(args.outdir, source_dir, header_lines, entries, drop_set, sidecars)

    kept_names = [n for n in master_names if n not in drop_set]
    print(f"Wrote {os.path.join(args.outdir, 'datacard.txt')}")
    print(f"  kept={len(kept_names)}/{len(master_names)}  dropped={sorted(drop_set) or '-'}")


if __name__ == "__main__":
    main()
