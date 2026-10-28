#!/usr/bin/env python3

import argparse
import re
from pathlib import Path


def fasta_ids(path):
    ids = []

    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])

    return ids


def read_ids(path):
    with open(path) as f:
        return [
            line.strip()
            for line in f
            if line.strip()
        ]


def check_unique(ids, label):
    duplicates = sorted(
        {x for x in ids if ids.count(x) > 1}
    )

    if duplicates:
        raise SystemExit(
            f"Duplicate IDs in {label}: "
            + ", ".join(duplicates[:20])
        )


def check_tree(path, bait_ids):
    tree = Path(path).read_text()

    missing = []

    for bait in bait_ids:
        pattern = (
            r"(?<=[(,])"
            + re.escape(bait)
            + r"(?=[:),])"
        )

        if re.search(pattern, tree) is None:
            missing.append(bait)

    if missing:
        raise SystemExit(
            f"Baits missing from tree {path}: "
            + ", ".join(missing)
        )


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--all", required=True)
    parser.add_argument("--aligned", required=True)
    parser.add_argument("--trimmed", required=True)
    parser.add_argument("--bait-ids", required=True)
    parser.add_argument("--counts-out", required=True)

    parser.add_argument("--fasttree")
    parser.add_argument("--iqtree")

    return parser.parse_args()


def main():
    args = parse_args()

    all_file = Path(args.all)
    aligned_file = Path(args.aligned)
    trimmed_file = Path(args.trimmed)
    bait_file = Path(args.bait_ids)
    counts_file = Path(args.counts_out)

    for path in (
        all_file,
        aligned_file,
        trimmed_file,
        bait_file,
    ):
        if not path.is_file():
            raise SystemExit(f"File not found: {path}")

    all_ids = fasta_ids(all_file)
    aligned_ids = fasta_ids(aligned_file)
    trimmed_ids = fasta_ids(trimmed_file)
    bait_ids = read_ids(bait_file)

    check_unique(all_ids, "input dataset")
    check_unique(aligned_ids, "alignment")
    check_unique(trimmed_ids, "trimmed alignment")
    check_unique(bait_ids, "bait list")

    missing_aligned = sorted(
        set(all_ids) - set(aligned_ids)
    )

    if missing_aligned:
        raise SystemExit(
            "Sequences missing after MAFFT: "
            + ", ".join(missing_aligned[:20])
        )

    extra_aligned = sorted(
        set(aligned_ids) - set(all_ids)
    )

    if extra_aligned:
        raise SystemExit(
            "Unexpected sequences after MAFFT: "
            + ", ".join(extra_aligned[:20])
        )

    missing_trimmed = sorted(
        set(aligned_ids) - set(trimmed_ids)
    )

    missing_baits = sorted(
        set(bait_ids) - set(trimmed_ids)
    )

    if missing_baits:
        raise SystemExit(
            "Baits missing after trimming: "
            + ", ".join(missing_baits)
        )

    counts_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with open(counts_file, "w") as out:
        out.write("stage\tsequences\n")
        out.write(f"input\t{len(all_ids)}\n")
        out.write(f"aligned\t{len(aligned_ids)}\n")
        out.write(f"trimmed\t{len(trimmed_ids)}\n")

    if args.fasttree:
        fasttree_file = Path(args.fasttree)

        if not fasttree_file.is_file():
            raise SystemExit(
                f"FastTree file not found: {fasttree_file}"
            )

        check_tree(fasttree_file, bait_ids)

    if args.iqtree:
        iqtree_file = Path(args.iqtree)

        if not iqtree_file.is_file():
            raise SystemExit(
                f"IQ-TREE file not found: {iqtree_file}"
            )

        check_tree(iqtree_file, bait_ids)

    print(f"Input sequences: {len(all_ids)}")
    print(f"Aligned sequences: {len(aligned_ids)}")
    print(f"Trimmed sequences: {len(trimmed_ids)}")

    if missing_trimmed:
        print(
            f"Sequences removed during trimming: "
            f"{len(missing_trimmed)}"
        )

    print(f"Baits checked: {len(bait_ids)}")

    if args.fasttree:
        print("FastTree bait check: OK")

    if args.iqtree:
        print("IQ-TREE bait check: OK")


if __name__ == "__main__":
    main()
