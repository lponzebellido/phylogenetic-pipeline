#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

from Bio import Phylo


def load_groups(path):
    groups = []

    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {"clade_id", "anchors", "color", "label"}
        missing = required - set(reader.fieldnames or [])

        if missing:
            raise SystemExit(
                "Missing columns: "
                + ", ".join(sorted(missing))
            )

        for row in reader:
            anchors = [
                x.strip()
                for x in row["anchors"].split(";")
                if x.strip()
            ]

            if len(anchors) < 2:
                raise SystemExit(
                    f"{row['clade_id']}: at least two anchors required"
                )

            row["anchor_list"] = anchors
            groups.append(row)

    return groups


def leaf_set(clade):
    return {
        leaf.name
        for leaf in clade.get_terminals()
        if leaf.name
    }


def find_smallest_split(tree, anchors):
    all_leaves = leaf_set(tree.root)
    anchor_set = set(anchors)

    missing = anchor_set - all_leaves

    if missing:
        raise SystemExit(
            "Anchors missing from tree: "
            + ", ".join(sorted(missing))
        )

    best = None

    for clade in tree.find_clades(order="preorder"):
        if clade is tree.root:
            continue

        side_a = leaf_set(clade)
        side_b = all_leaves - side_a

        for orientation, side in [
            ("descendant", side_a),
            ("complement", side_b),
        ]:
            if not anchor_set.issubset(side):
                continue

            candidate = {
                "size": len(side),
                "orientation": orientation,
                "members": side,
                "boundary_clade": clade,
            }

            if best is None or candidate["size"] < best["size"]:
                best = candidate

    return best


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--tree", required=True)
    parser.add_argument("--groups", required=True)
    parser.add_argument("--out-dir", required=True)

    args = parser.parse_args()

    tree = Phylo.read(args.tree, "newick")
    groups = load_groups(args.groups)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    total = len(tree.get_terminals())

    print(
        "clade_id\tanchors\t"
        "split_size\ttotal_leaves\tfraction\torientation"
    )

    summary_rows = []

    for group in groups:
        result = find_smallest_split(
            tree,
            group["anchor_list"],
        )

        if result is None:
            raise SystemExit(
                f"No split found for {group['clade_id']}"
            )

        member_file = (
            out_dir
            / f"{group['clade_id']}_members.txt"
        )

        with open(member_file, "w") as out:
            for member in sorted(result["members"]):
                out.write(member + "\n")

        fraction = result["size"] / total

        print(
            f"{group['clade_id']}\t"
            f"{';'.join(group['anchor_list'])}\t"
            f"{result['size']}\t"
            f"{total}\t"
            f"{fraction:.6f}\t"
            f"{result['orientation']}"
        )

        summary_rows.append(
            {
                "clade_id": group["clade_id"],
                "anchors": ";".join(group["anchor_list"]),
                "split_size": result["size"],
                "total_leaves": total,
                "fraction": f"{fraction:.6f}",
                "orientation": result["orientation"],
                "color": group["color"],
                "label": group["label"],
            }
        )

    summary_file = out_dir / "candidate_splits.tsv"

    with open(summary_file, "w", newline="") as out:
        fields = [
            "clade_id",
            "anchors",
            "split_size",
            "total_leaves",
            "fraction",
            "orientation",
            "color",
            "label",
        ]

        writer = csv.DictWriter(
            out,
            fieldnames=fields,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(summary_rows)


if __name__ == "__main__":
    main()
