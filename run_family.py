#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path

from anthotransport.paths import FamilyPaths
from anthotransport.pipeline import Pipeline, PipelineConfig, STAGES, selected_stages


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Run a reproducible, resumable transporter-family phylogeny pipeline."
    )
    result.add_argument("--family", required=True, help="Short family name, e.g. ABCC")
    location = result.add_mutually_exclusive_group(required=True)
    location.add_argument("--family-dir", type=Path, help="Existing family directory")
    location.add_argument("--workdir", type=Path, help="Directory for a new family")
    result.add_argument("--baits", type=Path, help="Bait FASTA; inferred when omitted")
    result.add_argument(
        "--references", type=Path, help="Reference metadata TSV; inferred when omitted"
    )
    result.add_argument(
        "--candidate-groups",
        type=Path,
        help="Accepted candidate-groups TSV; inferred when omitted",
    )
    result.add_argument("--proteomes", type=Path, help="Proteome FASTA directory")
    result.add_argument("--threads", type=int, default=28)
    result.add_argument("--blast-threads", type=int, default=4)
    result.add_argument("--hits-per-bait", type=int, default=10)
    result.add_argument("--similarity-cutoff", type=float, default=30.0)
    result.add_argument("--length-cutoff", type=int, default=30)
    result.add_argument("--occupancy", type=float, default=0.1)
    result.add_argument(
        "--fasttree-bin",
        type=Path,
        default=Path(os.environ.get("FASTTREE_BIN", "FastTree")),
        help="FastTree executable or path (default: FASTTREE_BIN or FastTree on PATH)",
    )
    result.add_argument("--iqtree-model", default="LG+G")
    result.add_argument("--neighborhood-top", type=int, default=25)
    selection = result.add_mutually_exclusive_group()
    selection.add_argument("--stage", choices=STAGES)
    selection.add_argument("--from-stage", choices=STAGES)
    result.add_argument("--to-stage", choices=STAGES)
    result.add_argument("--resume", action="store_true")
    result.add_argument(
        "--adopt-existing",
        action="store_true",
        help="Validate existing outputs and record them instead of rerunning",
    )
    result.add_argument(
        "--force",
        action="store_true",
        help="Allow an expensive stage to rerun despite existing untracked/stale output",
    )
    result.add_argument("--dry-run", action="store_true")
    return result


def main() -> int:
    args = parser().parse_args()
    if args.to_stage and args.stage:
        raise SystemExit("--to-stage cannot be combined with --stage")
    root = (args.family_dir or args.workdir).resolve()
    family = args.family.strip()
    paths = FamilyPaths(
        family=family,
        root=root,
        baits=(args.baits or root / "baits" / f"{family}_baits.fasta").resolve(),
        references=(
            args.references or root / "metadata" / f"{family}_references.tsv"
        ).resolve(),
        candidate_groups=(
            args.candidate_groups
            or root / "metadata" / f"{family}_candidate_groups.tsv"
        ).resolve(),
    )
    if min(args.threads, args.blast_threads, args.hits_per_bait, args.length_cutoff) < 1:
        raise SystemExit("Thread, hit, and length values must be positive")
    fasttree_bin = args.fasttree_bin.expanduser()
    if not fasttree_bin.is_absolute() and fasttree_bin.parent == Path("."):
        located = shutil.which(str(fasttree_bin))
        if located:
            fasttree_bin = Path(located).resolve()
    else:
        fasttree_bin = fasttree_bin.resolve()
    config = PipelineConfig(
        paths=paths,
        pipeline_root=Path(__file__).resolve().parent,
        proteomes=args.proteomes.resolve() if args.proteomes else None,
        threads=args.threads,
        blast_threads=args.blast_threads,
        hits_per_bait=args.hits_per_bait,
        similarity_cutoff=args.similarity_cutoff,
        length_cutoff=args.length_cutoff,
        occupancy=args.occupancy,
        fasttree_bin=fasttree_bin,
        iqtree_model=args.iqtree_model,
        neighborhood_top=args.neighborhood_top,
    )
    pipeline = Pipeline(
        config,
        resume=args.resume,
        force=args.force,
        adopt_existing=args.adopt_existing,
        dry_run=args.dry_run,
    )
    try:
        pipeline.run(selected_stages(args.stage, args.from_stage, args.to_stage))
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

