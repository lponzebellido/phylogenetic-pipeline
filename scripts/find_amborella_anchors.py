#!/usr/bin/env python3

import argparse
from Bio import Phylo


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Rank Amborella sequences according to the size of "
            "their MRCA clade with a functional bait."
        )
    )

    parser.add_argument(
        "--tree",
        required=True,
        help="Newick tree produced by IQ-TREE",
    )

    parser.add_argument(
        "--bait",
        required=True,
        help="Functional reference sequence ID",
    )

    parser.add_argument(
        "--top",
        type=int,
        default=10,
        help="Number of Amborella candidates to report",
    )

    args = parser.parse_args()

    tree = Phylo.read(args.tree, "newick")

    terminals = tree.get_terminals()

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
            f"Bait not found in tree: {args.bait}"
        )

    amborella = [
        leaf
        for leaf in terminals
        if leaf.name
        and leaf.name.startswith(
            "Amborella_trichopoda"
        )
    ]

    if not amborella:
        raise SystemExit(
            "No Amborella sequences found in tree"
        )

    rows = []

    for leaf in amborella:
        mrca = tree.common_ancestor(
            bait,
            leaf,
        )

        descendants = mrca.get_terminals()

        rows.append(
            {
                "amborella": leaf.name,
                "mrca_leaves": len(descendants),
                "distance": tree.distance(
                    bait,
                    leaf,
                ),
            }
        )

    rows.sort(
        key=lambda row: (
            row["mrca_leaves"],
            row["distance"],
        )
    )

    print(
        "rank\t"
        "mrca_leaves\t"
        "distance\t"
        "amborella_id"
    )

    for rank, row in enumerate(
        rows[: args.top],
        start=1,
    ):
        print(
            f"{rank}\t"
            f"{row['mrca_leaves']}\t"
            f"{row['distance']:.6f}\t"
            f"{row['amborella']}"
        )


if __name__ == "__main__":
    main()
