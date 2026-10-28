from __future__ import annotations

from pathlib import Path

from Bio import Phylo

from .splits import leaf_set
from .utils import atomic_write_text, read_tsv, resolve_tree_id


EVIDENCE = {
    "direct": {
        "shape": "3",
        "color": "#d73027",
        "size": "12",
        "label": "Direct transport evidence",
    },
    "strong": {
        "shape": "2",
        "color": "#fc8d59",
        "size": "10",
        "label": "Strong functional/genetic evidence",
    },
    "supporting": {
        "shape": "4",
        "color": "#4575b4",
        "size": "10",
        "label": "Supporting evidence",
    },
    "uncertain": {
        "shape": "1",
        "color": "#969696",
        "size": "8",
        "label": "Uncertain",
    },
}


def _style(level: str) -> dict[str, str]:
    return EVIDENCE.get(level.strip().lower(), EVIDENCE["uncertain"])


def _references(
    family: str, metadata: Path, leaves: set[str]
) -> list[dict[str, str]]:
    rows = read_tsv(
        metadata,
        {"tree_id", "species", "family", "evidence_level", "use_as_bait"},
    )
    selected = [
        row
        for row in rows
        if row["family"] == family and row["use_as_bait"].strip().lower() == "yes"
    ]
    if not selected:
        raise ValueError(f"No bait/reference rows found for family {family} in {metadata}")
    resolved: list[dict[str, str]] = []
    missing: list[str] = []
    for row in selected:
        node_id = resolve_tree_id(leaves, row["tree_id"].strip())
        if node_id is None:
            missing.append(row["tree_id"])
        else:
            resolved.append({**row, "tree_id": node_id})
    if missing:
        raise ValueError("Reference IDs missing from tree: " + ", ".join(missing))
    return resolved


def _amborella(mapping: Path | None, leaves: set[str]) -> list[str]:
    if mapping is None:
        return sorted(name for name in leaves if name.lower().startswith("amborella_"))
    rows = read_tsv(mapping, {"normalized_id", "original_header"})
    resolved: list[str] = []
    for row in rows:
        if "amborella" not in row["original_header"].lower():
            continue
        node_id = resolve_tree_id(leaves, row["normalized_id"])
        if node_id:
            resolved.append(node_id)
    return list(dict.fromkeys(resolved))


def _reference_style_lines(references: list[dict[str, str]]) -> list[str]:
    lines: list[str] = []
    for row in references:
        style = _style(row["evidence_level"])
        lines.extend(
            [
                f"{row['tree_id']}\tbranch\tnode\t{style['color']}\t4\tnormal",
                f"{row['tree_id']}\tlabel\tnode\t{style['color']}\t1.3\tbold",
            ]
        )
    return lines


def _candidate_style_lines(clades: list[dict[str, str]]) -> list[str]:
    return [
        f"{row['anchor_a']}|{row['anchor_b']}\tbranch\tclade\t"
        f"{row['color'].strip() or '#7B2CBF'}\t2\tnormal"
        for row in clades
    ]


def write_itol_annotations(
    family: str,
    tree_path: Path,
    metadata: Path,
    output_dir: Path,
    mapping: Path | None = None,
    candidate_clades: Path | None = None,
) -> dict[str, int]:
    tree = Phylo.read(str(tree_path), "newick")
    leaves = leaf_set(tree.root)
    references = _references(family, metadata, leaves)
    amborella = _amborella(mapping, leaves)
    output_dir.mkdir(parents=True, exist_ok=True)

    reference_lines = _reference_style_lines(references)
    style_header = (
        "DATASET_STYLE\n"
        "SEPARATOR TAB\n"
        f"DATASET_LABEL\t{family} reference styles\n"
        "COLOR\t#000000\n"
        "DATA\n"
    )
    atomic_write_text(
        output_dir / f"{family}_iTOL_reference_styles.txt",
        style_header + "\n".join(reference_lines) + "\n",
    )

    levels: list[str] = []
    for row in references:
        level = row["evidence_level"].strip().lower()
        if level not in EVIDENCE:
            level = "uncertain"
        if level not in levels:
            levels.append(level)
    symbol_lines = []
    for row in references:
        style = _style(row["evidence_level"])
        symbol_lines.append(
            f"{row['tree_id']}\t{style['shape']}\t{style['size']}\t"
            f"{style['color']}\t1\t1\t{row['tree_id']}"
        )
    symbol_text = (
        "DATASET_SYMBOL\n"
        "SEPARATOR TAB\n"
        f"DATASET_LABEL\t{family} reference evidence\n"
        "COLOR\t#000000\n"
        "MAXIMUM_SIZE\t15\n"
        "LEGEND_TITLE\tReference evidence\n"
        "LEGEND_SHAPES\t"
        + "\t".join(EVIDENCE[level]["shape"] for level in levels)
        + "\nLEGEND_COLORS\t"
        + "\t".join(EVIDENCE[level]["color"] for level in levels)
        + "\nLEGEND_LABELS\t"
        + "\t".join(EVIDENCE[level]["label"] for level in levels)
        + "\nDATA\n"
        + "\n".join(symbol_lines)
        + "\n"
    )
    atomic_write_text(output_dir / f"{family}_iTOL_reference_symbols.txt", symbol_text)

    amborella_lines = [
        f"{node_id}\t2\t10\t#7B2CBF\t1\t0\tAmborella" for node_id in amborella
    ]
    amborella_text = (
        "DATASET_SYMBOL\n"
        "SEPARATOR TAB\n"
        "DATASET_LABEL\tAmborella\n"
        "COLOR\t#7B2CBF\n"
        "MAXIMUM_SIZE\t15\n"
        "LEGEND_TITLE\tTaxonomic reference\n"
        "LEGEND_SHAPES\t2\n"
        "LEGEND_COLORS\t#7B2CBF\n"
        "LEGEND_LABELS\tAmborella trichopoda\n"
        "DATA\n"
        + "\n".join(amborella_lines)
        + "\n"
    )
    atomic_write_text(output_dir / f"{family}_iTOL_amborella_symbols.txt", amborella_text)

    candidate_count = 0
    if candidate_clades is not None and candidate_clades.is_file():
        clades = read_tsv(
            candidate_clades, {"clade_id", "anchor_a", "anchor_b", "color", "label"}
        )
        resolved_clades: list[dict[str, str]] = []
        for row in clades:
            anchor_a = resolve_tree_id(leaves, row["anchor_a"])
            anchor_b = resolve_tree_id(leaves, row["anchor_b"])
            if anchor_a is None or anchor_b is None:
                raise ValueError(
                    f"{row['clade_id']}: candidate-view anchor missing from tree"
                )
            resolved_clades.append({**row, "anchor_a": anchor_a, "anchor_b": anchor_b})
        candidate_lines = _candidate_style_lines(resolved_clades)
        candidate_text = (
            "DATASET_STYLE\n"
            "SEPARATOR TAB\n"
            f"DATASET_LABEL\t{family} candidate clades\n"
            "COLOR\t#7B2CBF\n"
            "DATA\n"
            + "\n".join(candidate_lines)
            + "\n"
        )
        atomic_write_text(
            output_dir / f"{family}_iTOL_candidate_clades.txt", candidate_text
        )

        # Candidate styles intentionally come first; reference styles then override them.
        combined = (
            "DATASET_STYLE\n"
            "SEPARATOR TAB\n"
            f"DATASET_LABEL\t{family} candidate and reference styles\n"
            "COLOR\t#000000\n"
            "DATA\n"
            + "\n".join(candidate_lines + reference_lines)
            + "\n"
        )
        atomic_write_text(output_dir / f"{family}_iTOL_combined_styles.txt", combined)
        candidate_count = len(resolved_clades)

    return {
        "references": len(references),
        "amborella": len(amborella),
        "candidate_clades": candidate_count,
    }

