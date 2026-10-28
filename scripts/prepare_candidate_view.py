#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

from Bio import Phylo


def read_groups(path):
    groups = []

    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {
            "clade_id",
            "color",
            "label",
        }

        missing = required - set(reader.fieldnames or [])

        if missing:
            raise SystemExit(
                "Missing group columns: "
                + ", ".join(sorted(missing))
            )

        for row in reader:
            groups.append(row)

    return groups


def read_members(path):
    with open(path) as handle:
        return {
            line.strip()
            for line in handle
            if line.strip()
        }


def terminal_names(clade):
    return {
        leaf.name
        for leaf in clade.get_terminals()
        if leaf.name
    }


def first_terminal(clade):
    leaves = [
        leaf
        for leaf in clade.get_terminals()
        if leaf.name
    ]

    if not leaves:
        raise SystemExit(
            "Encountered child clade without terminals"
        )

    return sorted(
        leaves,
        key=lambda leaf: leaf.name,
    )[0]


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--tree",
        required=True,
    )

    parser.add_argument(
        "--groups",
        required=True,
    )

    parser.add_argument(
        "--members-dir",
        required=True,
    )

    parser.add_argument(
        "--out-tree",
        required=True,
    )

    parser.add_argument(
        "--out-clades",
        required=True,
    )

    parser.add_argument(
        "--root-info",
        required=True,
    )

    args = parser.parse_args()

    tree = Phylo.read(
        args.tree,
        "newick",
    )

    groups = read_groups(
        args.groups
    )

    members_dir = Path(
        args.members_dir
    )

    group_members = {}

    for row in groups:
        clade_id = row["clade_id"]

        member_file = (
            members_dir
            / f"{clade_id}_members.txt"
        )

        if not member_file.is_file():
            raise SystemExit(
                f"Missing member file: {member_file}"
            )

        group_members[clade_id] = read_members(
            member_file
        )

    # Validate that candidate splits are disjoint.
    clade_ids = [
        row["clade_id"]
        for row in groups
    ]

    for i, a in enumerate(clade_ids):
        for b in clade_ids[i + 1:]:
            overlap = (
                group_members[a]
                & group_members[b]
            )

            if overlap:
                raise SystemExit(
                    f"{a} and {b} overlap by "
                    f"{len(overlap)} sequences"
                )

    all_leaves = terminal_names(
        tree.root
    )

    candidate_union = set().union(
        *group_members.values()
    )

    outside = sorted(
        all_leaves - candidate_union
    )

    if not outside:
        raise SystemExit(
            "No sequence exists outside candidate splits"
        )

    # Deterministic root used ONLY for visualization.
    root_name = outside[0]

    terminals = {
        leaf.name: leaf
        for leaf in tree.get_terminals()
        if leaf.name
    }

    root_leaf = terminals[
        root_name
    ]

    tree.root_with_outgroup(
        root_leaf
    )

    # Rebuild terminal lookup after rerooting.
    terminals = {
        leaf.name: leaf
        for leaf in tree.get_terminals()
        if leaf.name
    }

    candidate_rows = []

    for row in groups:
        clade_id = row["clade_id"]
        members = group_members[
            clade_id
        ]

        missing = (
            members
            - set(terminals)
        )

        if missing:
            raise SystemExit(
                f"{clade_id}: members missing "
                "after rerooting: "
                + ", ".join(
                    sorted(missing)
                )
            )

        targets = [
            terminals[name]
            for name in members
        ]

        clade = tree.common_ancestor(
            *targets
        )

        observed = terminal_names(
            clade
        )

        if observed != members:
            raise SystemExit(
                f"{clade_id}: rerooted MRCA "
                f"contains {len(observed)} leaves, "
                f"expected {len(members)}"
            )

        if len(clade.clades) < 2:
            raise SystemExit(
                f"{clade_id}: candidate clade "
                "does not have two child branches"
            )

        # Pick one terminal from two different
        # immediate children. Their MRCA is therefore
        # exactly the candidate-clade root.
        anchor_a_leaf = first_terminal(
            clade.clades[0]
        )

        anchor_b_leaf = first_terminal(
            clade.clades[1]
        )

        check = tree.common_ancestor(
            anchor_a_leaf,
            anchor_b_leaf,
        )

        if terminal_names(check) != members:
            raise SystemExit(
                f"{clade_id}: selected anchors "
                "do not recover full candidate clade"
            )

        candidate_rows.append(
            {
                "clade_id": clade_id,
                "anchor_a": anchor_a_leaf.name,
                "anchor_b": anchor_b_leaf.name,
                "color": row["color"],
                "label": row["label"],
            }
        )

        print(
            f"{clade_id}: "
            f"{len(members)} members | "
            f"{anchor_a_leaf.name} | "
            f"{anchor_b_leaf.name}"
        )

    out_tree = Path(
        args.out_tree
    )

    out_tree.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    Phylo.write(
        tree,
        out_tree,
        "newick",
    )

    out_clades = Path(
        args.out_clades
    )

    out_clades.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with open(
        out_clades,
        "w",
        newline="",
    ) as handle:
        fields = [
            "clade_id",
            "anchor_a",
            "anchor_b",
            "color",
            "label",
        ]

        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(
            candidate_rows
        )

    root_info = Path(
        args.root_info
    )

    with open(
        root_info,
        "w",
    ) as handle:
        handle.write(
            "IMPORTANT: visualization root only; "
            "not an inferred biological root.\n"
        )

        handle.write(
            f"root_leaf\t{root_name}\n"
        )

        handle.write(
            f"candidate_sequences\t"
            f"{len(candidate_union)}\n"
        )

        handle.write(
            f"outside_sequences\t"
            f"{len(outside)}\n"
        )

    print()
    print(
        f"Visualization root: {root_name}"
    )
    print(
        f"Tree: {out_tree}"
    )
    print(
        f"Candidate clades: {out_clades}"
    )


if __name__ == "__main__":
    main()
