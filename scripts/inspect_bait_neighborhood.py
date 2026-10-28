#!/usr/bin/env python3

import argparse
import csv
from Bio import Phylo


def load_references(path):
    references = set()

    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            tree_id = row.get("tree_id", "").strip()

            if tree_id:
                references.add(tree_id)

    return references


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--tree",
        required=True,
    )

    parser.add_argument(
        "--bait",
        required=True,
    )

    parser.add_argument(
        "--metadata",
        required=True,
    )

    parser.add_argument(
        "--top",
        type=int,
        default=50,
    )

    args = parser.parse_args()

    tree = Phylo.read(
        args.tree,
        "newick",
    )

    terminals = tree.get_terminals()

    references = load_references(
        args.metadata
    )

    bait = next(
        (
            leaf
            for leaf in terminals
            if leaf.name == args.bait
        ),
        None,
    )

    if bait is None:
        raise SystemExit(
            f"Bait not found: {args.bait}"
        )

    rows = []

    for leaf in terminals:
        if leaf is bait:
            continue

        name = leaf.name or ""

        if name in references:
            category = "REFERENCE"

        elif name.startswith(
            "Amborella_trichopoda"
        ):
            category = "AMBORELLA"

        else:
            category = "candidate"

        rows.append(
            (
                tree.distance(
                    bait,
                    leaf,
                ),
                category,
                name,
            )
        )

    rows.sort(
        key=lambda row: row[0]
    )

    print(
        "rank\tdistance\tcategory\tsequence"
    )

    for rank, (
        distance,
        category,
        name,
    ) in enumerate(
        rows[: args.top],
        start=1,
    ):
        print(
            f"{rank}\t"
            f"{distance:.6f}\t"
            f"{category}\t"
            f"{name}"
        )


if __name__ == "__main__":
    main()
