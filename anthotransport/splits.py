from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from Bio import Phylo

from .utils import atomic_write_text, atomic_write_tsv, read_tsv, resolve_tree_id


@dataclass(frozen=True)
class Split:
    members: frozenset[str]
    orientation: str
    boundary_index: int

    @property
    def size(self) -> int:
        return len(self.members)


def leaf_set(clade) -> set[str]:
    return {leaf.name for leaf in clade.get_terminals() if leaf.name}


def find_smallest_split(tree, anchors: list[str]) -> Split:
    all_leaves = leaf_set(tree.root)
    anchor_set = set(anchors)
    missing = anchor_set - all_leaves
    if missing:
        raise ValueError(f"Anchors missing from tree: {', '.join(sorted(missing))}")

    # One postorder pass supplies subtree sizes and anchor counts. This avoids
    # repeatedly materializing every descendant set on trees with thousands of leaves.
    stats: dict[int, tuple[int, int]] = {}

    def visit(clade) -> tuple[int, int]:
        if not clade.clades:
            value = (1, int(bool(clade.name and clade.name in anchor_set)))
        else:
            children = [visit(child) for child in clade.clades]
            value = (sum(item[0] for item in children), sum(item[1] for item in children))
        stats[id(clade)] = value
        return value

    visit(tree.root)
    candidates: list[tuple[int, str, int, object]] = []
    for index, clade in enumerate(tree.find_clades(order="preorder")):
        if clade is tree.root:
            continue
        subtree_size, anchor_count = stats[id(clade)]
        if anchor_count == len(anchor_set):
            candidates.append((subtree_size, "descendant", index, clade))
        if anchor_count == 0:
            candidates.append((len(all_leaves) - subtree_size, "complement", index, clade))
    if not candidates:
        raise ValueError("No edge split contains every anchor")
    smallest = min(item[0] for item in candidates)
    finalists: list[Split] = []
    for size, orientation, index, clade in candidates:
        if size != smallest:
            continue
        descendant = leaf_set(clade)
        members = descendant if orientation == "descendant" else all_leaves - descendant
        finalists.append(Split(frozenset(members), orientation, index))
    return min(finalists, key=lambda item: (tuple(sorted(item.members)), item.orientation))


def load_groups(path: Path, leaf_names: set[str]) -> list[dict[str, object]]:
    rows = read_tsv(path, {"clade_id", "anchors", "color", "label"})
    if not rows:
        raise ValueError(
            f"{path} has no accepted groups. Fill this explicit scientific decision file first."
        )
    seen: set[str] = set()
    parsed: list[dict[str, object]] = []
    for row in rows:
        clade_id = row["clade_id"].strip()
        if not clade_id or clade_id in seen:
            raise ValueError(f"Missing or duplicate clade_id: {clade_id!r}")
        seen.add(clade_id)
        raw_anchors = [part.strip() for part in row["anchors"].split(";") if part.strip()]
        if len(raw_anchors) < 2:
            raise ValueError(f"{clade_id}: at least two anchors are required")
        resolved = []
        for anchor in raw_anchors:
            value = resolve_tree_id(leaf_names, anchor)
            if value is None:
                raise ValueError(f"{clade_id}: anchor missing from tree: {anchor}")
            resolved.append(value)
        parsed.append({**row, "clade_id": clade_id, "anchor_list": resolved})
    return parsed


def identify_candidate_splits(
    tree_path: Path,
    groups_path: Path,
    output_dir: Path,
    references_path: Path | None = None,
) -> list[dict[str, str]]:
    tree = Phylo.read(str(tree_path), "newick")
    all_leaves = leaf_set(tree.root)
    groups = load_groups(groups_path, all_leaves)
    output_dir.mkdir(parents=True, exist_ok=True)

    reference_ids: set[str] = set()
    if references_path and references_path.is_file():
        for row in read_tsv(references_path, {"tree_id"}):
            resolved = resolve_tree_id(all_leaves, row["tree_id"].strip())
            if resolved:
                reference_ids.add(resolved)
    amborella_ids = {name for name in all_leaves if name.lower().startswith("amborella_")}

    results: list[tuple[dict[str, object], Split, Path]] = []
    for group in groups:
        split = find_smallest_split(tree, group["anchor_list"])  # type: ignore[arg-type]
        member_file = output_dir / f"{group['clade_id']}_members.txt"
        atomic_write_text(
            member_file, "".join(f"{name}\n" for name in sorted(split.members))
        )
        results.append((group, split, member_file))

    overlap_rows: list[dict[str, object]] = []
    for index, (group_a, split_a, _) in enumerate(results):
        for group_b, split_b, _ in results[index + 1 :]:
            overlap = split_a.members & split_b.members
            overlap_rows.append(
                {
                    "clade_a": group_a["clade_id"],
                    "clade_b": group_b["clade_id"],
                    "overlap_count": len(overlap),
                    "overlap_examples": ";".join(sorted(overlap)[:10]),
                }
            )
    atomic_write_tsv(
        output_dir / "candidate_split_overlaps.tsv",
        overlap_rows,
        ["clade_a", "clade_b", "overlap_count", "overlap_examples"],
    )

    rows: list[dict[str, str]] = []
    for group, split, member_file in results:
        references = sorted(split.members & reference_ids)
        amborella = sorted(split.members & amborella_ids)
        rows.append(
            {
                "clade_id": str(group["clade_id"]),
                "anchors": ";".join(group["anchor_list"]),  # type: ignore[arg-type]
                "split_size": str(split.size),
                "total_leaves": str(len(all_leaves)),
                "fraction": f"{split.size / len(all_leaves):.6f}",
                "orientation": split.orientation,
                "member_file": str(member_file),
                "functional_references": ";".join(references),
                "amborella_sequences": ";".join(amborella),
                "color": str(group["color"]),
                "label": str(group["label"]),
            }
        )
    atomic_write_tsv(
        output_dir / "candidate_splits.tsv",
        rows,
        [
            "clade_id",
            "anchors",
            "split_size",
            "total_leaves",
            "fraction",
            "orientation",
            "member_file",
            "functional_references",
            "amborella_sequences",
            "color",
            "label",
        ],
    )

    overlapping = [row for row in overlap_rows if int(row["overlap_count"]) > 0]
    if overlapping:
        pairs = ", ".join(
            f"{row['clade_a']}/{row['clade_b']}={row['overlap_count']}"
            for row in overlapping
        )
        raise ValueError(
            "Accepted candidate splits overlap; downstream processing stopped. "
            f"Review {output_dir / 'candidate_split_overlaps.tsv'} ({pairs})"
        )
    return rows

