#!/usr/bin/env python3

import argparse
import csv
import re
from pathlib import Path


EVIDENCE = {
    "direct": {
        "shape": 3,
        "color": "#d73027",
        "size": 12,
        "label": "Direct transport evidence",
    },
    "strong": {
        "shape": 2,
        "color": "#fc8d59",
        "size": 10,
        "label": "Strong functional/genetic evidence",
    },
    "supporting": {
        "shape": 4,
        "color": "#4575b4",
        "size": 10,
        "label": "Supporting evidence",
    },
    "uncertain": {
        "shape": 1,
        "color": "#969696",
        "size": 8,
        "label": "Uncertain",
    },
}


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--family", required=True)
    parser.add_argument("--tree", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--mapping")
    parser.add_argument("--out-dir", required=True)

    return parser.parse_args()


def load_metadata(path):
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)

    required = {
        "tree_id",
        "species",
        "family",
        "evidence_level",
        "use_as_bait",
    }

    missing = required - set(reader.fieldnames or [])

    if missing:
        raise SystemExit(
            "Missing metadata columns: "
            + ", ".join(sorted(missing))
        )

    return rows


#def tree_contains(tree, node_id):
#    pattern = (
#        r"(?:^|[(,])'?"
#        + re.escape(node_id)
#        + r"'?(?=[:),])"
#    )
#
#    return re.search(pattern, tree) is not None


def resolve_tree_id(tree, node_id):
    candidates = [node_id]

    sanitized = node_id.replace("@", "_")

    if sanitized != node_id:
        candidates.append(sanitized)

    for candidate in candidates:
        pattern = (
            r"(?:^|[(,])'?"
            + re.escape(candidate)
            + r"'?(?=[:),])"
        )

        if re.search(pattern, tree):
            return candidate

    return None

def load_amborella_ids(mapping_file):
    ids = []

    with open(mapping_file) as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            original = row.get("original_header", "")

            if "amborella" in original.lower():
                ids.append(row["normalized_id"])

    return ids


def write_reference_styles(path, family, references):
    with open(path, "w") as out:
        out.write("DATASET_STYLE\n")
        out.write("SEPARATOR TAB\n")
        out.write(
            f"DATASET_LABEL\t{family} reference styles\n"
        )
        out.write("COLOR\t#000000\n")
        out.write("DATA\n")

        for row in references:
            level = row["evidence_level"].lower()
            style = EVIDENCE.get(
                level,
                EVIDENCE["uncertain"],
            )

            node_id = row["tree_id"]
            color = style["color"]

            out.write(
                f"{node_id}\tbranch\tnode\t"
                f"{color}\t4\tnormal\n"
            )

            out.write(
                f"{node_id}\tlabel\tnode\t"
                f"{color}\t1.3\tbold\n"
            )


def write_reference_symbols(path, family, references):
    levels_present = []

    for row in references:
        level = row["evidence_level"].lower()

        if level not in EVIDENCE:
            level = "uncertain"

        if level not in levels_present:
            levels_present.append(level)

    with open(path, "w") as out:
        out.write("DATASET_SYMBOL\n")
        out.write("SEPARATOR TAB\n")
        out.write(
            f"DATASET_LABEL\t{family} reference evidence\n"
        )
        out.write("COLOR\t#000000\n")
        out.write("MAXIMUM_SIZE\t15\n")

        out.write("LEGEND_TITLE\tReference evidence\n")

        out.write(
            "LEGEND_SHAPES\t"
            + "\t".join(
                str(EVIDENCE[level]["shape"])
                for level in levels_present
            )
            + "\n"
        )

        out.write(
            "LEGEND_COLORS\t"
            + "\t".join(
                EVIDENCE[level]["color"]
                for level in levels_present
            )
            + "\n"
        )

        out.write(
            "LEGEND_LABELS\t"
            + "\t".join(
                EVIDENCE[level]["label"]
                for level in levels_present
            )
            + "\n"
        )

        out.write("DATA\n")

        for row in references:
            level = row["evidence_level"].lower()

            if level not in EVIDENCE:
                level = "uncertain"

            style = EVIDENCE[level]

            out.write(
                f"{row['tree_id']}\t"
                f"{style['shape']}\t"
                f"{style['size']}\t"
                f"{style['color']}\t"
                f"1\t1\t"
                f"{row['tree_id']}\n"
            )


def write_amborella_symbols(path, ids):
    with open(path, "w") as out:
        out.write("DATASET_SYMBOL\n")
        out.write("SEPARATOR TAB\n")
        out.write("DATASET_LABEL\tAmborella\n")
        out.write("COLOR\t#7B2CBF\n")
        out.write("MAXIMUM_SIZE\t15\n")
        out.write("LEGEND_TITLE\tTaxonomic reference\n")
        out.write("LEGEND_SHAPES\t2\n")
        out.write("LEGEND_COLORS\t#7B2CBF\n")
        out.write("LEGEND_LABELS\tAmborella trichopoda\n")
        out.write("DATA\n")

        for node_id in ids:
            out.write(
                f"{node_id}\t"
                "2\t"
                "10\t"
                "#7B2CBF\t"
                "1\t"
                "0\t"
                "Amborella\n"
            )


def load_candidate_clades(path):
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)

    required = {
        "clade_id",
        "anchor_a",
        "anchor_b",
        "color",
        "label",
    }

    missing = required - set(reader.fieldnames or [])

    if missing:
        raise SystemExit(
            "Missing candidate-clade columns: "
            + ", ".join(sorted(missing))
        )

    return rows


def write_candidate_clades(path, family, clades):
    with open(path, "w") as out:
        out.write("DATASET_STYLE\n")
        out.write("SEPARATOR TAB\n")
        out.write(f"DATASET_LABEL\t{family} candidate clades\n")
        out.write("COLOR\t#7B2CBF\n")
        out.write("DATA\n")

        for row in clades:
            node_id = f"{row['anchor_a']}|{row['anchor_b']}"
            color = row["color"].strip() or "#7B2CBF"

            out.write(
                f"{node_id}\t"
                "branch\t"
                "clade\t"
                f"{color}\t"
                "2\t"
                "normal\n"
            )


def main():
    args = parse_args()

    tree_file = Path(args.tree)
    metadata_file = Path(args.metadata)
    out_dir = Path(args.out_dir)

    if not tree_file.is_file():
        raise SystemExit(
            f"Tree file not found: {tree_file}"
        )

    if not metadata_file.is_file():
        raise SystemExit(
            f"Metadata file not found: {metadata_file}"
        )

    out_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    tree = tree_file.read_text()

    metadata = load_metadata(metadata_file)

    references = [
        row
        for row in metadata
        if row["family"] == args.family
        and row["use_as_bait"].lower() == "yes"
    ]

    if not references:
        raise SystemExit(
            f"No references found for family {args.family}"
        )

    #missing = [
    #    row["tree_id"]
    #    for row in references
    #    if not tree_contains(tree, row["tree_id"])
    #]

        resolved_references = []
        missing = []

        for row in references:
            resolved_id = resolve_tree_id(
                tree,
                row["tree_id"],
            )

            if resolved_id is None:
                missing.append(row["tree_id"])
                continue

            resolved_row = row.copy()
            resolved_row["tree_id"] = resolved_id
            resolved_references.append(resolved_row)

        references = resolved_references

        if missing:
            raise SystemExit(
                "Reference IDs missing from tree: "
                + ", ".join(missing)
            )

    styles_file = (
        out_dir
        / f"{args.family}_iTOL_reference_styles.txt"
    )

    symbols_file = (
        out_dir
        / f"{args.family}_iTOL_reference_symbols.txt"
    )

    write_reference_styles(
        styles_file,
        args.family,
        references,
    )

    write_reference_symbols(
        symbols_file,
        args.family,
        references,
    )

    amborella_count = 0

    if args.mapping:
        mapping_file = Path(args.mapping)

        if not mapping_file.is_file():
            raise SystemExit(
                f"Mapping file not found: {mapping_file}"
            )

        amborella_ids = load_amborella_ids(
            mapping_file
        )

     #   amborella_ids = [
     #       node_id
     #       for node_id in amborella_ids
     #       if tree_contains(tree, node_id)
     #   ]

        amborella_ids = [
            resolve_tree_id(tree, node_id)
            for node_id in amborella_ids
        ]

        amborella_ids = [
            node_id
            for node_id in amborella_ids
            if node_id is not None
        ]

        amborella_ids = list(
            dict.fromkeys(amborella_ids)
        )

        if amborella_ids:
            amborella_file = (
                out_dir
                / f"{args.family}_iTOL_amborella_symbols.txt"
            )

            write_amborella_symbols(
                amborella_file,
                amborella_ids,
            )

            amborella_count = len(
                amborella_ids
            )

    candidate_count = 0

    candidate_file = (
        metadata_file.parent
        / f"{args.family}_candidate_clades.tsv"
    )

    if candidate_file.is_file():
        candidate_clades = load_candidate_clades(
            candidate_file
        )

        resolved_clades = []
        missing_anchors = []

        for row in candidate_clades:
            anchor_a = resolve_tree_id(
                tree,
                row["anchor_a"],
            )

            anchor_b = resolve_tree_id(
                tree,
                row["anchor_b"],
            )

            if anchor_a is None:
                missing_anchors.append(
                    f"{row['clade_id']}:{row['anchor_a']}"
                )

            if anchor_b is None:
                missing_anchors.append(
                    f"{row['clade_id']}:{row['anchor_b']}"
                )

            if anchor_a is None or anchor_b is None:
                continue

            resolved = row.copy()
            resolved["anchor_a"] = anchor_a
            resolved["anchor_b"] = anchor_b
            resolved_clades.append(resolved)

        if missing_anchors:
            raise SystemExit(
                "Candidate-clade anchors missing from tree: "
                + ", ".join(missing_anchors)
            )

        if resolved_clades:
            candidate_out = (
                out_dir
                / f"{args.family}_iTOL_candidate_clades.txt"
            )

            write_candidate_clades(
                candidate_out,
                args.family,
                resolved_clades,
            )

            candidate_count = len(resolved_clades)


    print(
        f"References annotated: {len(references)}"
    )
    print(
        f"Amborella sequences annotated: {amborella_count}"
    )
    print(
        f"Candidate clades annotated: {candidate_count}"
    )
    print(
        f"Output directory: {out_dir}"
    )


if __name__ == "__main__":
    main()
