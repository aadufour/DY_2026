#!/usr/bin/env python3
"""
make_nuisance_variants.py
==========================
Generate datacard.txt variants with different nuisance subsets, for isolating
which systematic causes a non-parabolic ("wiggly") likelihood scan.

Does NOT touch shapes.root or re-run make_cards.py/spritz-cards-eft: combine
only reads a histo_{proc}_{syst}Up/Down pair if the datacard's systematics
block has a row for it, so one shared shapes.root works for every variant.
Each variant folder just symlinks the source dir's shapes.root (and any other
sidecar files, e.g. metadata.json) and gets its own filtered datacard.txt.

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

Everything up to and including the 2nd dash-only line is copied verbatim into
every variant. After that, each line's first whitespace-separated token is
checked against THEORY_NUISANCES (QCDScale, PDFweight) — only those two are
candidates for removal/toggling in the variants below. Every other line
(lumi, other experimental systs, autoMCStats, blanks, anything unrecognized)
is always kept verbatim, since this tool isolates theoretical uncertainties
specifically, not the full nuisance set.

Usage
-----
List the nuisances found in the real datacard (dry run, no files written):
    python3 make_nuisance_variants.py --datacard datacards/inc_mm/mll/datacard.txt --list

Leave-one-out (start from the full set, drop one nuisance per variant):
    python3 make_nuisance_variants.py \\
        --datacard datacards/inc_mm/mll/datacard.txt \\
        --outdir   nuisance_scan/inc_mm_mll \\
        --mode     leave-one-out

Add-one-in (start stat-only, add one nuisance per variant):
    python3 make_nuisance_variants.py \\
        --datacard datacards/inc_mm/mll/datacard.txt \\
        --outdir   nuisance_scan/inc_mm_mll \\
        --mode     add-one-in

Incremental (start with 1 nuisance, add one more per variant up to the full
set — N variants for N nuisances found). Order defaults to datacard order;
pass --shuffle-seed for a reproducible random order instead:
    python3 make_nuisance_variants.py \\
        --datacard datacards/inc_mm/mll/datacard.txt \\
        --outdir   nuisance_scan/inc_mm_mll \\
        --mode     incremental [--shuffle-seed 42]

Custom variants from a JSON spec {variant_name: [nuisance_names_to_DROP]}:
    python3 make_nuisance_variants.py \\
        --datacard datacards/inc_mm/mll/datacard.txt \\
        --outdir   nuisance_scan/inc_mm_mll \\
        --mode     custom --spec my_variants.json
"""

import argparse
import json
import os
import random


THEORY_NUISANCES = {"QCDScale", "PDFweight"}


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
    (name_or_None, raw_line) — name is None for lines that are always kept:
    blank, autoMCStats, or any systematic not in THEORY_NUISANCES. Only
    QCDScale/PDFweight rows get a real name and are eligible to be dropped
    by the variant modes below.
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
        if name not in THEORY_NUISANCES:
            entries.append((None, line))
            continue
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


def write_variant(outdir, variant_name, source_dir, header_lines, entries, drop_set, kept_names, sidecars):
    variant_dir = os.path.join(outdir, variant_name)
    os.makedirs(variant_dir, exist_ok=True)

    lines = build_variant_lines(header_lines, entries, drop_set)
    with open(os.path.join(variant_dir, "datacard.txt"), "w") as f:
        f.writelines(lines)

    # Symlink only the requested sidecar files (default: shapes.root) — the source
    # dir may already contain a full prior run's output (model_*.root, scans,
    # plots, ...); we don't want any of that pulled into a fresh variant folder.
    for fname in sidecars:
        src = os.path.abspath(os.path.join(source_dir, fname))
        if not os.path.exists(src):
            print(f"  WARNING: sidecar '{fname}' not found in {source_dir}, skipping")
            continue
        dst = os.path.join(variant_dir, fname)
        if os.path.islink(dst) or os.path.exists(dst):
            os.remove(dst)
        os.symlink(os.path.relpath(src, variant_dir), dst)

    return variant_dir, kept_names


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--datacard", required=True, help="Path to the full datacard.txt (from spritz-cards-eft)")
    parser.add_argument("--outdir", help="Directory to write variant subfolders into (required unless --list)")
    parser.add_argument("--mode", choices=["leave-one-out", "add-one-in", "incremental", "custom"], help="Variant generation mode")
    parser.add_argument("--spec", help="JSON file {variant_name: [nuisance_names_to_drop]} — required for --mode custom")
    parser.add_argument("--shuffle-seed", type=int, default=None,
                        help="--mode incremental only: shuffle the add-order with this random seed instead of using datacard order")
    parser.add_argument("--sidecars", default="shapes.root",
                        help="Comma-separated filenames from the source dir to symlink into each variant folder "
                             "(default: shapes.root only — deliberately NOT everything in the source dir, which may "
                             "already hold a prior run's model_*.root/scans/plots)")
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

    if not args.outdir or not args.mode:
        parser.error("--outdir and --mode are required unless --list is given")

    source_dir = os.path.dirname(os.path.abspath(args.datacard))
    os.makedirs(args.outdir, exist_ok=True)
    sidecars = [s.strip() for s in args.sidecars.split(",") if s.strip()]

    variants = {}  # variant_name -> drop_set
    if args.mode == "leave-one-out":
        variants["full"] = set()
        for n in master_names:
            variants[f"no_{n}"] = {n}
    elif args.mode == "add-one-in":
        variants["stat_only"] = set(master_names)
        for n in master_names:
            variants[f"stat_plus_{n}"] = set(master_names) - {n}
    elif args.mode == "incremental":
        order = list(master_names)
        if args.shuffle_seed is not None:
            random.Random(args.shuffle_seed).shuffle(order)
        for i in range(1, len(order) + 1):
            variant_name = f"{i:02d}_plus_{order[i - 1]}"
            variants[variant_name] = set(master_names) - set(order[:i])
    elif args.mode == "custom":
        if not args.spec:
            parser.error("--mode custom requires --spec")
        with open(args.spec) as f:
            spec = json.load(f)
        for variant_name, drop_list in spec.items():
            unknown = set(drop_list) - set(master_names)
            if unknown:
                parser.error(f"variant '{variant_name}': unknown nuisance name(s) {unknown} not in {master_names}")
            variants[variant_name] = set(drop_list)

    manifest = {}
    if args.mode == "incremental":
        manifest["_meta"] = {"mode": args.mode, "order": order, "shuffle_seed": args.shuffle_seed}
    print(f"Writing {len(variants)} variant(s) into {args.outdir}/")
    for variant_name, drop_set in variants.items():
        kept_names = [n for n in master_names if n not in drop_set]
        variant_dir, kept = write_variant(
            args.outdir, variant_name, source_dir, header_lines, entries, drop_set, kept_names, sidecars
        )
        manifest[variant_name] = {"kept": kept_names, "dropped": sorted(drop_set)}
        print(f"  {variant_name:30s} kept={len(kept_names)}/{len(master_names)}  dropped={sorted(drop_set) or '-'}")

    with open(os.path.join(args.outdir, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"Wrote {os.path.join(args.outdir, 'manifest.json')}")


if __name__ == "__main__":
    main()
