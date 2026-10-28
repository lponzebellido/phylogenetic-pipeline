from __future__ import annotations

import json
import os
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

from Bio import Phylo

from . import __version__
from .annotations import write_itol_annotations
from .candidate_view import create_candidate_view
from .neighborhoods import analyze_neighborhoods
from .paths import FamilyPaths
from .splits import identify_candidate_splits, leaf_set
from .utils import (
    CommandRunner,
    atomic_write_json,
    atomic_write_text,
    executable_version,
    fasta_ids,
    fingerprint_path,
    require_unique,
    signature,
    utc_now,
)


STAGES = [
    "init",
    "hits",
    "dataset",
    "alignment",
    "trim",
    "fasttree",
    "iqtree",
    "validate",
    "itol-base",
    "neighborhoods",
    "candidate-splits",
    "candidate-view",
    "itol-final",
    "report",
]

EXPENSIVE_STAGES = {"hits", "alignment", "fasttree", "iqtree"}


@dataclass(frozen=True)
class PipelineConfig:
    paths: FamilyPaths
    pipeline_root: Path
    proteomes: Path | None
    threads: int = 28
    blast_threads: int = 4
    hits_per_bait: int = 10
    similarity_cutoff: float = 30.0
    length_cutoff: int = 30
    occupancy: float = 0.1
    fasttree_bin: Path = Path("FastTree")
    iqtree_model: str = "LG+G"
    neighborhood_top: int = 25

    @property
    def shared_scripts(self) -> Path:
        return self.pipeline_root / "shared_scripts"


class Pipeline:
    def __init__(
        self,
        config: PipelineConfig,
        *,
        resume: bool = False,
        force: bool = False,
        adopt_existing: bool = False,
        dry_run: bool = False,
    ):
        self.config = config
        self.paths = config.paths
        self.resume = resume
        self.force = force
        self.adopt_existing = adopt_existing
        self.dry_run = dry_run
        self.runner = CommandRunner(
            self.paths.logs / f"{self.paths.family}_pipeline_commands.log", dry_run
        )
        self.manifest = self._load_manifest()

    def _load_manifest(self) -> dict[str, Any]:
        if self.paths.manifest.is_file():
            try:
                return json.loads(self.paths.manifest.read_text())
            except json.JSONDecodeError as exc:
                raise ValueError(f"Invalid pipeline manifest: {self.paths.manifest}") from exc
        return {
            "schema_version": 1,
            "pipeline_version": __version__,
            "family": self.paths.family,
            "family_dir": str(self.paths.root),
            "stages": {},
        }

    def _save_manifest(self) -> None:
        self.manifest["updated_at"] = utc_now()
        if not self.dry_run:
            atomic_write_json(self.paths.manifest, self.manifest)

    def stage_outputs(self, stage: str) -> list[Path]:
        p = self.paths
        outputs = {
            "init": [p.references, p.candidate_groups],
            "hits": [p.hits],
            "dataset": [
                p.candidates_fasta,
                p.all_fasta,
                p.mapping,
                p.qc / f"{p.family}_sequence_QC.tsv",
                p.bait_ids,
                p.selected_hits,
            ],
            "alignment": [p.aligned],
            "trim": [p.trimmed],
            "fasttree": [p.fasttree],
            "iqtree": [p.iqtree],
            "validate": [p.counts],
            "itol-base": [
                p.itol / f"{p.family}_iTOL_reference_styles.txt",
                p.itol / f"{p.family}_iTOL_reference_symbols.txt",
                p.itol / f"{p.family}_iTOL_amborella_symbols.txt",
            ],
            "neighborhoods": [
                p.qc / f"{p.family}_reference_pairs.tsv",
                p.qc / f"{p.family}_reference_amborella.tsv",
                p.qc / f"{p.family}_reference_neighborhoods.tsv",
                p.decision_report,
            ],
            "candidate-splits": [p.split_summary, p.overlap_report],
            "candidate-view": [p.candidate_view, p.candidate_clades, p.root_info],
            "itol-final": [
                p.itol / f"{p.family}_iTOL_reference_symbols.txt",
                p.itol / f"{p.family}_iTOL_amborella_symbols.txt",
                p.itol / f"{p.family}_iTOL_candidate_clades.txt",
                p.itol / f"{p.family}_iTOL_combined_styles.txt",
            ],
            "report": [p.report],
        }
        return outputs[stage]

    def stage_inputs(self, stage: str) -> list[Path]:
        c, p = self.config, self.paths
        scripts = c.shared_scripts
        values = {
            "init": [],
            "hits": [p.baits, c.proteomes, scripts / "collect_best_BLAST_hits.py"],
            "dataset": [
                p.baits,
                p.hits,
                scripts / "summarize_hits.py",
                scripts / "prepare_dataset.py",
            ],
            "alignment": [p.all_fasta],
            "trim": [p.aligned, scripts / "algntrim.py"],
            "fasttree": [p.trimmed, c.fasttree_bin],
            "iqtree": [p.trimmed],
            "validate": [
                p.all_fasta,
                p.aligned,
                p.trimmed,
                p.bait_ids,
                p.fasttree,
                p.iqtree,
                scripts / "validate_dataset.py",
            ],
            "itol-base": [p.iqtree, p.references, p.mapping],
            "neighborhoods": [p.iqtree, p.references],
            "candidate-splits": [p.iqtree, p.candidate_groups, p.references],
            "candidate-view": [
                p.iqtree,
                p.candidate_groups,
                p.split_summary,
                p.clades,
                p.references,
            ],
            "itol-final": [p.candidate_view, p.references, p.mapping, p.candidate_clades],
            "report": [
                p.iqtree,
                p.references,
                p.candidate_view,
                p.candidate_clades,
                p.itol,
            ],
        }
        return [path for path in values[stage] if path is not None]

    def stage_parameters(self, stage: str) -> dict[str, Any]:
        c = self.config
        common = {
            "family": self.paths.family,
            "family_dir": str(self.paths.root),
            "pipeline_version": __version__,
        }
        specific = {
            "hits": {
                "hits_per_bait": c.hits_per_bait,
                "similarity_cutoff": c.similarity_cutoff,
                "length_cutoff": c.length_cutoff,
                "blast_threads": c.blast_threads,
            },
            "alignment": {"threads": c.threads},
            "trim": {"occupancy": c.occupancy},
            "fasttree": {"model": "LG+GAMMA", "binary": str(c.fasttree_bin)},
            "iqtree": {"model": c.iqtree_model, "threads": c.threads},
            "neighborhoods": {"top": c.neighborhood_top},
        }
        return {**common, **specific.get(stage, {})}

    def _input_signature(self, stage: str) -> str:
        value = {
            "stage": stage,
            "parameters": self.stage_parameters(stage),
            "inputs": [fingerprint_path(path) for path in self.stage_inputs(stage)],
        }
        return signature(value)

    def _outputs_exist(self, stage: str) -> bool:
        outputs = self.stage_outputs(stage)
        if not outputs:
            return stage == "init"
        return all(path.exists() for path in outputs)

    def _is_current(self, stage: str, input_signature: str) -> bool:
        record = self.manifest.get("stages", {}).get(stage)
        return bool(
            record
            and record.get("status") in {"completed", "adopted"}
            and record.get("input_signature") == input_signature
            and self._outputs_exist(stage)
        )

    def _validate_existing(self, stage: str) -> None:
        if not self._outputs_exist(stage):
            raise ValueError(f"Cannot adopt {stage}: one or more expected outputs are missing")
        p = self.paths
        if stage == "hits":
            hits = list(p.hits.rglob("best_hit_sequences.fasta"))
            blasts = list(p.hits.rglob("blast_results.txt"))
            if not hits or len(hits) != len(blasts):
                raise ValueError("Cannot adopt hits: incomplete BLAST/best-hit file pairs")
        elif stage in {"dataset", "alignment", "trim"}:
            target = {
                "dataset": p.all_fasta,
                "alignment": p.aligned,
                "trim": p.trimmed,
            }[stage]
            ids = fasta_ids(target)
            if not ids:
                raise ValueError(f"Cannot adopt empty FASTA: {target}")
            require_unique(ids, str(target))
        elif stage in {"fasttree", "iqtree", "candidate-view"}:
            target = {
                "fasttree": p.fasttree,
                "iqtree": p.iqtree,
                "candidate-view": p.candidate_view,
            }[stage]
            tree = Phylo.read(str(target), "newick")
            if not tree.get_terminals():
                raise ValueError(f"Cannot adopt empty tree: {target}")

    def _record_stage(self, stage: str, input_signature: str, status: str) -> None:
        if "software_versions" not in self.manifest:
            self.manifest["software_versions"] = self.tool_versions()
        self.manifest.setdefault("stages", {})[stage] = {
            "status": status,
            "recorded_at": utc_now(),
            "input_signature": input_signature,
            "parameters": self.stage_parameters(stage),
            "inputs": [fingerprint_path(path) for path in self.stage_inputs(stage)],
            "outputs": [fingerprint_path(path) for path in self.stage_outputs(stage)],
        }
        self._save_manifest()

    def run(self, stages: list[str]) -> None:
        for stage in stages:
            self.run_stage(stage)

    def run_stage(self, stage: str) -> None:
        if stage not in STAGES:
            raise ValueError(f"Unknown stage: {stage}")
        self.paths.ensure_directories()
        input_signature = self._input_signature(stage)
        print(f"\n== {self.paths.family}: {stage} ==")

        if self.resume and self._is_current(stage, input_signature):
            print("current: skipped")
            return

        if self.adopt_existing and self._outputs_exist(stage):
            self._validate_existing(stage)
            print("existing outputs validated and adopted into the manifest")
            self._record_stage(stage, input_signature, "adopted")
            return

        if stage in EXPENSIVE_STAGES and self._outputs_exist(stage) and not self.force:
            raise ValueError(
                f"{stage} has existing outputs but no matching provenance record. "
                "Refusing an expensive rerun. Use --adopt-existing after checking the "
                "parameters, or --force to rerun it."
            )

        function: Callable[[], None] = getattr(self, f"stage_{stage.replace('-', '_')}")
        started = utc_now()
        try:
            function()
        except Exception:
            self.manifest.setdefault("stages", {})[stage] = {
                "status": "failed",
                "started_at": started,
                "failed_at": utc_now(),
                "input_signature": input_signature,
                "parameters": self.stage_parameters(stage),
            }
            self._save_manifest()
            raise
        if not self.dry_run and not self._outputs_exist(stage):
            raise RuntimeError(f"{stage} completed without all declared outputs")
        self._record_stage(stage, input_signature, "completed")

    def _require_proteomes(self) -> Path:
        if self.config.proteomes is None or not self.config.proteomes.is_dir():
            raise ValueError("--proteomes must name an existing directory for this stage")
        return self.config.proteomes

    def _archive(self, path: Path) -> None:
        if not path.exists():
            return
        timestamp = utc_now().replace(":", "").replace("+", "_")
        destination = self.paths.root / "archive" / timestamp / path.relative_to(self.paths.root)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(path), str(destination))

    def stage_init(self) -> None:
        self.paths.ensure_directories()
        if not self.paths.references.exists():
            template = self.paths.references
            header = (
                "tree_id\tspecies\tgene_id\tprotein_accession\tfamily\t"
                "functional_context\tevidence_type\tevidence_level\tuse_as_bait\t"
                "reference\tnotes\n"
            )
            if not self.dry_run:
                atomic_write_text(template, header)
            print(f"created reference metadata template: {template}")
        if not self.paths.candidate_groups.exists():
            template = self.paths.candidate_groups
            header = "clade_id\tanchors\tcolor\tlabel\n"
            if not self.dry_run:
                atomic_write_text(template, header)
            print(f"created candidate-group template: {template}")

    def stage_hits(self) -> None:
        proteomes = self._require_proteomes()
        if not self.paths.baits.is_file():
            raise FileNotFoundError(f"Bait FASTA not found: {self.paths.baits}")
        script = self.config.shared_scripts / "collect_best_BLAST_hits.py"
        with tempfile.TemporaryDirectory(dir=self.paths.results) as temporary:
            temp_hits = Path(temporary) / "hits"
            temp_hits.mkdir()
            self.runner.run(
                [
                    sys.executable,
                    str(script),
                    "--baits",
                    str(self.paths.baits),
                    "--out",
                    str(temp_hits) + "/",
                    "--subjectdir",
                    str(proteomes),
                    "--number",
                    str(self.config.hits_per_bait),
                    "--simcut",
                    str(self.config.similarity_cutoff),
                    "--lencut",
                    str(self.config.length_cutoff),
                    "--cpu",
                    str(self.config.blast_threads),
                ]
            )
            if self.dry_run:
                return
            hits = list(temp_hits.rglob("best_hit_sequences.fasta"))
            blasts = list(temp_hits.rglob("blast_results.txt"))
            if not hits or len(hits) != len(blasts):
                raise RuntimeError("Homolog collection produced incomplete per-proteome outputs")
            self._archive(self.paths.hits)
            shutil.move(str(temp_hits), str(self.paths.hits))

    def stage_dataset(self) -> None:
        proteomes = self._require_proteomes()
        scripts = self.config.shared_scripts
        with tempfile.TemporaryDirectory(dir=self.paths.results) as temporary:
            base = Path(temporary)
            dataset = base / "dataset"
            qc = base / "qc"
            dataset.mkdir()
            qc.mkdir()
            self.runner.run(
                [
                    sys.executable,
                    str(scripts / "summarize_hits.py"),
                    "--baits",
                    str(self.paths.baits),
                    "--subject-dir",
                    str(proteomes),
                    "--hits-dir",
                    str(self.paths.hits),
                    "--out",
                    str(qc / self.paths.selected_hits.name),
                    "--number",
                    str(self.config.hits_per_bait),
                    "--simcut",
                    str(self.config.similarity_cutoff),
                    "--lencut",
                    str(self.config.length_cutoff),
                ]
            )
            self.runner.run(
                [
                    sys.executable,
                    str(scripts / "prepare_dataset.py"),
                    "--family",
                    self.paths.family,
                    "--baits",
                    str(self.paths.baits),
                    "--hits-dir",
                    str(self.paths.hits),
                    "--dataset-dir",
                    str(dataset),
                    "--qc-dir",
                    str(qc),
                ]
            )
            if self.dry_run:
                return
            for source in [*dataset.iterdir(), *qc.iterdir()]:
                destination_dir = self.paths.dataset if source.parent == dataset else self.paths.qc
                destination = destination_dir / source.name
                self._archive(destination)
                shutil.move(str(source), str(destination))

    def stage_alignment(self) -> None:
        temporary = self.paths.alignment / f".{self.paths.family}_aligned.tmp"
        log = self.paths.alignment / f"{self.paths.family}_mafft.log"
        self.runner.run(
            ["mafft", "--thread", str(self.config.threads), str(self.paths.all_fasta)],
            stdout_path=temporary,
            stderr_path=log,
        )
        if not self.dry_run:
            self._archive(self.paths.aligned)
            os.replace(temporary, self.paths.aligned)

    def stage_trim(self) -> None:
        temporary = self.paths.alignment / f".{self.paths.family}_trimmed.tmp"
        log = self.paths.alignment / f"{self.paths.family}_trim.log"
        self.runner.run(
            [
                sys.executable,
                str(self.config.shared_scripts / "algntrim.py"),
                "--in",
                str(self.paths.aligned),
                "--out",
                str(temporary),
                "--occ",
                str(self.config.occupancy),
            ],
            stdout_path=log,
            stderr_path=log,
        )
        if not self.dry_run:
            self._archive(self.paths.trimmed)
            os.replace(temporary, self.paths.trimmed)

    def stage_fasttree(self) -> None:
        temporary = self.paths.trees / f".{self.paths.family}_FastTree.tmp"
        log = self.paths.trees / f"{self.paths.family}_FastTree.log"
        self.runner.run(
            [
                str(self.config.fasttree_bin),
                "-lg",
                "-gamma",
                "-log",
                str(log),
                str(self.paths.trimmed),
            ],
            stdout_path=temporary,
        )
        if not self.dry_run:
            self._archive(self.paths.fasttree)
            os.replace(temporary, self.paths.fasttree)

    def stage_iqtree(self) -> None:
        with tempfile.TemporaryDirectory(dir=self.paths.trees) as temporary:
            prefix = Path(temporary) / f"{self.paths.family}_IQTREE"
            self.runner.run(
                [
                    "iqtree2",
                    "-s",
                    str(self.paths.trimmed),
                    "-m",
                    self.config.iqtree_model,
                    "-T",
                    str(self.config.threads),
                    "--prefix",
                    str(prefix),
                ]
            )
            if self.dry_run:
                return
            if not prefix.with_suffix(".treefile").is_file():
                raise RuntimeError("IQ-TREE did not produce a .treefile")
            for existing in self.paths.trees.glob(f"{self.paths.family}_IQTREE.*"):
                self._archive(existing)
            for source in Path(temporary).iterdir():
                shutil.move(str(source), str(self.paths.trees / source.name))

    def stage_validate(self) -> None:
        self.runner.run(
            [
                sys.executable,
                str(self.config.shared_scripts / "validate_dataset.py"),
                "--all",
                str(self.paths.all_fasta),
                "--aligned",
                str(self.paths.aligned),
                "--trimmed",
                str(self.paths.trimmed),
                "--bait-ids",
                str(self.paths.bait_ids),
                "--counts-out",
                str(self.paths.counts),
                "--fasttree",
                str(self.paths.fasttree),
                "--iqtree",
                str(self.paths.iqtree),
            ]
        )

    def stage_itol_base(self) -> None:
        counts = write_itol_annotations(
            self.paths.family,
            self.paths.iqtree,
            self.paths.references,
            self.paths.itol,
            self.paths.mapping,
        )
        print(f"annotations: {counts}")

    def stage_neighborhoods(self) -> None:
        counts = analyze_neighborhoods(
            self.paths.family,
            self.paths.iqtree,
            self.paths.references,
            self.paths.qc,
            self.config.neighborhood_top,
        )
        print(f"neighborhood analysis: {counts}")

    def stage_candidate_splits(self) -> None:
        rows = identify_candidate_splits(
            self.paths.iqtree,
            self.paths.candidate_groups,
            self.paths.clades,
            self.paths.references,
        )
        print("candidate split sizes: " + ", ".join(f"{r['clade_id']}={r['split_size']}" for r in rows))

    def stage_candidate_view(self) -> None:
        rows = create_candidate_view(
            self.paths.iqtree,
            self.paths.candidate_groups,
            self.paths.clades,
            self.paths.candidate_view,
            self.paths.candidate_clades,
            self.paths.root_info,
            self.paths.references,
        )
        print(f"candidate clades validated: {len(rows)}")

    def stage_itol_final(self) -> None:
        counts = write_itol_annotations(
            self.paths.family,
            self.paths.candidate_view,
            self.paths.references,
            self.paths.itol,
            self.paths.mapping,
            self.paths.candidate_clades,
        )
        print(f"final annotations: {counts}")

    def stage_report(self) -> None:
        tree = Phylo.read(str(self.paths.iqtree), "newick")
        total_leaves = len(leaf_set(tree.root))
        stages = self.manifest.get("stages", {})
        stage_lines = [
            f"| {name} | {stages.get(name, {}).get('status', 'not recorded')} |"
            for name in STAGES
        ]
        candidate_ready = self.paths.candidate_view.is_file()
        tree_to_upload = self.paths.candidate_view if candidate_ready else self.paths.iqtree
        itol_files = [
            self.paths.itol / f"{self.paths.family}_iTOL_reference_symbols.txt",
            self.paths.itol / f"{self.paths.family}_iTOL_amborella_symbols.txt",
        ]
        if candidate_ready:
            itol_files.insert(0, self.paths.itol / f"{self.paths.family}_iTOL_combined_styles.txt")
        text = "\n".join(
            [
                f"# {self.paths.family} pipeline report",
                "",
                f"Generated: {utc_now()}",
                f"Total IQ-TREE leaves: {total_leaves}",
                "",
                "## Stage status",
                "",
                "| Stage | Status |",
                "|---|---|",
                *stage_lines,
                "",
                "## iTOL upload set",
                "",
                f"Tree: `{tree_to_upload}`",
                *[f"- `{path}`" for path in itol_files],
                "",
                "The candidate-view tree is a visualization-only rerooting. Its root has no inferred biological meaning; the original IQ-TREE tree remains the primary ML result.",
                "",
                f"Manifest: `{self.paths.manifest}`",
                f"Candidate decision support: `{self.paths.decision_report}`",
                "",
            ]
        )
        atomic_write_text(self.paths.report, text)

    def tool_versions(self) -> dict[str, str]:
        return {
            "python": executable_version([sys.executable, "--version"]),
            "mafft": executable_version(["mafft", "--version"]),
            "iqtree2": executable_version(["iqtree2", "--version"]),
            "fasttree": executable_version([str(self.config.fasttree_bin)]),
        }


def selected_stages(
    stage: str | None, from_stage: str | None, to_stage: str | None
) -> list[str]:
    if stage:
        return [stage]
    start = STAGES.index(from_stage) if from_stage else 0
    end = STAGES.index(to_stage) + 1 if to_stage else len(STAGES)
    if start >= end:
        raise ValueError("--from-stage must not come after --to-stage")
    return STAGES[start:end]

