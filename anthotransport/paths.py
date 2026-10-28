from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class FamilyPaths:
    family: str
    root: Path
    baits: Path
    references: Path
    candidate_groups: Path

    @property
    def results(self) -> Path:
        return self.root / "results"

    @property
    def hits(self) -> Path:
        return self.results / "hits"

    @property
    def dataset(self) -> Path:
        return self.results / "dataset"

    @property
    def alignment(self) -> Path:
        return self.results / "alignment"

    @property
    def trees(self) -> Path:
        return self.results / "trees"

    @property
    def qc(self) -> Path:
        return self.results / "qc"

    @property
    def itol(self) -> Path:
        return self.results / "itol"

    @property
    def clades(self) -> Path:
        return self.results / "clades"

    @property
    def logs(self) -> Path:
        return self.root / "logs"

    @property
    def candidates_fasta(self) -> Path:
        return self.dataset / f"{self.family}_candidates.fasta"

    @property
    def all_fasta(self) -> Path:
        return self.dataset / f"{self.family}_all.fasta"

    @property
    def aligned(self) -> Path:
        return self.alignment / f"{self.family}_aligned.fasta"

    @property
    def trimmed(self) -> Path:
        return self.alignment / f"{self.family}_aligned_trimmed.fasta"

    @property
    def fasttree(self) -> Path:
        return self.trees / f"{self.family}_FastTree.nwk"

    @property
    def iqtree(self) -> Path:
        return self.trees / f"{self.family}_IQTREE.treefile"

    @property
    def candidate_view(self) -> Path:
        return self.trees / f"{self.family}_IQTREE_candidate_view.treefile"

    @property
    def mapping(self) -> Path:
        return self.qc / f"{self.family}_header_mapping.tsv"

    @property
    def bait_ids(self) -> Path:
        return self.qc / f"{self.family}_bait_ids.txt"

    @property
    def selected_hits(self) -> Path:
        return self.qc / f"{self.family}_selected_hits.tsv"

    @property
    def counts(self) -> Path:
        return self.qc / f"{self.family}_sequence_counts.tsv"

    @property
    def candidate_clades(self) -> Path:
        return self.root / "metadata" / f"{self.family}_candidate_clades.tsv"

    @property
    def split_summary(self) -> Path:
        return self.clades / "candidate_splits.tsv"

    @property
    def overlap_report(self) -> Path:
        return self.clades / "candidate_split_overlaps.tsv"

    @property
    def root_info(self) -> Path:
        return self.qc / f"{self.family}_candidate_view_root.tsv"

    @property
    def manifest(self) -> Path:
        return self.qc / f"{self.family}_pipeline_manifest.json"

    @property
    def report(self) -> Path:
        return self.results / f"{self.family}_pipeline_report.md"

    @property
    def decision_report(self) -> Path:
        return self.qc / f"{self.family}_candidate_decision_report.md"

    def ensure_directories(self) -> None:
        for path in (
            self.root / "baits",
            self.root / "metadata",
            self.hits,
            self.dataset,
            self.alignment,
            self.trees,
            self.qc,
            self.itol,
            self.clades,
            self.logs,
        ):
            path.mkdir(parents=True, exist_ok=True)

