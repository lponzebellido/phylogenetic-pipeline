from __future__ import annotations

import itertools
import tempfile
from pathlib import Path

from Bio import Phylo

from .splits import leaf_set
from .utils import atomic_write_text, atomic_write_tsv, read_tsv, resolve_tree_id


def _read_members(path: Path) -> set[str]:
    return {line.strip() for line in path.read_text().splitlines() if line.strip()}


def _first_terminal(clade):
    leaves = sorted((leaf for leaf in clade.get_terminals() if leaf.name), key=lambda x: x.name)
    if not leaves:
        raise ValueError("Candidate child branch has no named terminal")
    return leaves[0]


def _reference_ids(path: Path | None, leaves: set[str]) -> list[str]:
    if path is None or not path.is_file():
        return sorted(leaves)[: min(12, len(leaves))]
    resolved = []
    for row in read_tsv(path, {"tree_id"}):
        node_id = resolve_tree_id(leaves, row["tree_id"].strip())
        if node_id:
            resolved.append(node_id)
    return sorted(set(resolved))


def create_candidate_view(
    tree_path: Path,
    groups_path: Path,
    members_dir: Path,
    output_tree: Path,
    output_clades: Path,
    root_info: Path,
    references_path: Path | None = None,
    tolerance: float = 1e-6,
) -> list[dict[str, str]]:
    original = Phylo.read(str(tree_path), "newick")
    original_leaves = leaf_set(original.root)
    groups = read_tsv(groups_path, {"clade_id", "color", "label"})
    if not groups:
        raise ValueError(f"No accepted candidate groups in {groups_path}")

    members: dict[str, set[str]] = {}
    for row in groups:
        member_file = members_dir / f"{row['clade_id']}_members.txt"
        if not member_file.is_file():
            raise FileNotFoundError(f"Missing candidate member file: {member_file}")
        values = _read_members(member_file)
        missing = values - original_leaves
        if missing:
            raise ValueError(
                f"{row['clade_id']}: members missing from original tree: "
                + ", ".join(sorted(missing)[:20])
            )
        members[row["clade_id"]] = values

    for a, b in itertools.combinations(members, 2):
        overlap = members[a] & members[b]
        if overlap:
            raise ValueError(f"{a} and {b} overlap by {len(overlap)} sequences")

    candidate_union = set().union(*members.values())
    outside = sorted(original_leaves - candidate_union)
    if not outside:
        raise ValueError("No sequence exists outside the accepted candidate splits")
    root_name = outside[0]

    references = _reference_ids(references_path, original_leaves)
    before_distances = {
        (a, b): original.distance(a, b) for a, b in itertools.combinations(references, 2)
    }

    terminals = {leaf.name: leaf for leaf in original.get_terminals() if leaf.name}
    original.root_with_outgroup(terminals[root_name])

    rows: list[dict[str, str]] = []
    terminals = {leaf.name: leaf for leaf in original.get_terminals() if leaf.name}
    for row in groups:
        expected = members[row["clade_id"]]
        clade = original.common_ancestor(*(terminals[name] for name in expected))
        observed = leaf_set(clade)
        if observed != expected:
            raise ValueError(
                f"{row['clade_id']}: visualization MRCA has {len(observed)} leaves; "
                f"expected {len(expected)}"
            )
        if len(clade.clades) < 2:
            raise ValueError(f"{row['clade_id']}: candidate root has fewer than two children")
        anchor_a = _first_terminal(clade.clades[0]).name
        anchor_b = _first_terminal(clade.clades[1]).name
        if leaf_set(original.common_ancestor(anchor_a, anchor_b)) != expected:
            raise ValueError(f"{row['clade_id']}: generated anchors do not recover the clade")
        rows.append(
            {
                "clade_id": row["clade_id"],
                "anchor_a": anchor_a,
                "anchor_b": anchor_b,
                "color": row["color"].strip() or "#7B2CBF",
                "label": row["label"],
            }
        )

    output_tree.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", suffix=".treefile", dir=output_tree.parent, delete=False
    ) as handle:
        temporary_tree = Path(handle.name)
    try:
        Phylo.write(original, str(temporary_tree), "newick")
        written = Phylo.read(str(temporary_tree), "newick")
        written_leaves = leaf_set(written.root)
        if written_leaves != original_leaves:
            raise ValueError("Candidate-view tree changed the leaf set")
        for row in rows:
            observed = leaf_set(written.common_ancestor(row["anchor_a"], row["anchor_b"]))
            if observed != members[row["clade_id"]]:
                raise ValueError(f"{row['clade_id']}: membership changed after Newick rewrite")
        max_delta = 0.0
        for pair, before in before_distances.items():
            delta = abs(written.distance(*pair) - before)
            max_delta = max(max_delta, delta)
            if delta > tolerance:
                raise ValueError(
                    f"Reference distance changed by {delta:.3g} for {pair}; "
                    f"tolerance is {tolerance}"
                )
        temporary_tree.replace(output_tree)
    finally:
        if temporary_tree.exists():
            temporary_tree.unlink()

    atomic_write_tsv(
        output_clades,
        rows,
        ["clade_id", "anchor_a", "anchor_b", "color", "label"],
    )
    info = (
        "key\tvalue\n"
        "warning\tVisualization root only; not an inferred biological root.\n"
        f"source_tree\t{tree_path}\n"
        f"root_leaf\t{root_name}\n"
        f"total_leaves\t{len(original_leaves)}\n"
        f"candidate_sequences\t{len(candidate_union)}\n"
        f"outside_sequences\t{len(outside)}\n"
        f"reference_pairs_checked\t{len(before_distances)}\n"
        f"maximum_distance_delta\t{max_delta:.12g}\n"
    )
    atomic_write_text(root_info, info)
    return rows

