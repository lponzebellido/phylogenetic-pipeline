from __future__ import annotations

import itertools
from pathlib import Path

from Bio import Phylo

from .splits import find_smallest_split, leaf_set
from .utils import atomic_write_text, atomic_write_tsv, read_tsv, resolve_tree_id


def _leaf_distances(tree, start_name: str) -> dict[str, float]:
    """Return all patristic leaf distances with one linear tree traversal."""
    adjacency: dict[int, list[tuple[int, float]]] = {}
    nodes: dict[int, object] = {}
    start_id: int | None = None
    for parent in tree.find_clades(order="preorder"):
        parent_id = id(parent)
        nodes[parent_id] = parent
        adjacency.setdefault(parent_id, [])
        if not parent.clades and parent.name == start_name:
            start_id = parent_id
        for child in parent.clades:
            child_id = id(child)
            nodes[child_id] = child
            weight = float(child.branch_length or 0.0)
            adjacency.setdefault(parent_id, []).append((child_id, weight))
            adjacency.setdefault(child_id, []).append((parent_id, weight))
    if start_id is None:
        raise ValueError(f"Leaf missing from tree: {start_name}")
    result: dict[str, float] = {}
    stack = [(start_id, -1, 0.0)]
    while stack:
        node_id, parent_id, distance = stack.pop()
        node = nodes[node_id]
        if not node.clades and node.name:
            result[node.name] = distance
        for neighbor_id, weight in adjacency[node_id]:
            if neighbor_id != parent_id:
                stack.append((neighbor_id, node_id, distance + weight))
    return result


def analyze_neighborhoods(
    family: str,
    tree_path: Path,
    references_path: Path,
    output_dir: Path,
    top: int = 25,
) -> dict[str, int]:
    tree = Phylo.read(str(tree_path), "newick")
    leaves = leaf_set(tree.root)
    metadata = read_tsv(references_path, {"tree_id", "family", "use_as_bait"})
    references = []
    for row in metadata:
        if row["family"] != family or row["use_as_bait"].strip().lower() != "yes":
            continue
        node_id = resolve_tree_id(leaves, row["tree_id"].strip())
        if node_id is None:
            raise ValueError(f"Reference missing from tree: {row['tree_id']}")
        references.append(node_id)
    references = list(dict.fromkeys(references))
    amborella = sorted(name for name in leaves if name.lower().startswith("amborella_"))
    output_dir.mkdir(parents=True, exist_ok=True)
    distances = {reference: _leaf_distances(tree, reference) for reference in references}

    pair_rows = []
    for a, b in itertools.combinations(references, 2):
        split = find_smallest_split(tree, [a, b])
        pair_rows.append(
            {
                "reference_a": a,
                "reference_b": b,
                "patristic_distance": f"{distances[a][b]:.9f}",
                "smallest_split_size": split.size,
                "smallest_split_orientation": split.orientation,
            }
        )
    atomic_write_tsv(
        output_dir / f"{family}_reference_pairs.tsv",
        pair_rows,
        [
            "reference_a",
            "reference_b",
            "patristic_distance",
            "smallest_split_size",
            "smallest_split_orientation",
        ],
    )

    amborella_rows = []
    best_amborella: dict[str, list[dict[str, object]]] = {}
    for reference in references:
        candidates: list[dict[str, object]] = []
        for node_id in amborella:
            split = find_smallest_split(tree, [reference, node_id])
            candidates.append(
                {
                    "reference": reference,
                    "amborella_id": node_id,
                    "patristic_distance": distances[reference][node_id],
                    "smallest_split_size": split.size,
                    "smallest_split_orientation": split.orientation,
                }
            )
        candidates.sort(
            key=lambda row: (
                int(row["smallest_split_size"]),
                float(row["patristic_distance"]),
                str(row["amborella_id"]),
            )
        )
        best_amborella[reference] = candidates[:top]
        for rank, row in enumerate(candidates[:top], 1):
            amborella_rows.append(
                {
                    **row,
                    "rank": rank,
                    "patristic_distance": f"{float(row['patristic_distance']):.9f}",
                }
            )
    atomic_write_tsv(
        output_dir / f"{family}_reference_amborella.tsv",
        amborella_rows,
        [
            "reference",
            "rank",
            "amborella_id",
            "patristic_distance",
            "smallest_split_size",
            "smallest_split_orientation",
        ],
    )

    nearest_rows = []
    reference_set = set(references)
    for reference in references:
        ranked = []
        for node_id in leaves - {reference}:
            if node_id in reference_set:
                category = "REFERENCE"
            elif node_id in set(amborella):
                category = "AMBORELLA"
            else:
                category = "candidate"
            ranked.append((distances[reference][node_id], category, node_id))
        ranked.sort(key=lambda row: (row[0], row[2]))
        for rank, (distance, category, node_id) in enumerate(ranked[:top], 1):
            nearest_rows.append(
                {
                    "reference": reference,
                    "rank": rank,
                    "patristic_distance": f"{distance:.9f}",
                    "category": category,
                    "sequence": node_id,
                }
            )
    atomic_write_tsv(
        output_dir / f"{family}_reference_neighborhoods.tsv",
        nearest_rows,
        ["reference", "rank", "patristic_distance", "category", "sequence"],
    )

    lines = [
        f"# {family} candidate-group decision report",
        "",
        "This report supports, but does not replace, the biological decision recorded in ",
        f"`metadata/{family}_candidate_groups.tsv`.",
        "",
        "The split sizes below are root-independent edge-split sizes; the serialized Newick root is not interpreted biologically.",
        "",
        "## Best Amborella associations per reference",
        "",
        "| Reference | Amborella sequence | Split size | Patristic distance |",
        "|---|---|---:|---:|",
    ]
    for reference in references:
        if best_amborella[reference]:
            row = best_amborella[reference][0]
            lines.append(
                f"| {reference} | {row['amborella_id']} | {row['smallest_split_size']} | "
                f"{float(row['patristic_distance']):.6f} |"
            )
    lines.extend(
        [
            "",
            "Detailed tables:",
            "",
            f"- `{family}_reference_pairs.tsv`",
            f"- `{family}_reference_amborella.tsv`",
            f"- `{family}_reference_neighborhoods.tsv`",
            "",
        ]
    )
    atomic_write_text(
        output_dir / f"{family}_candidate_decision_report.md", "\n".join(lines)
    )
    return {
        "references": len(references),
        "amborella": len(amborella),
        "reference_pairs": len(pair_rows),
    }

