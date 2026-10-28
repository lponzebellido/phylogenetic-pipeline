# Phylogenetic pipeline

This pipeline searches a plant proteome collection with characterized reference proteins, builds broad phylogenetic trees and prepares evidence for manually defining candidate lineages.

This repository is not presented as production-ready software. It is a traceable snapshot for lab review and for deciding whether the workflow should remain standalone or whether generic parts should later be integrated into KIPEs.

The active ABCC, MATE and P3A analysis directories are not part of this repository and were not modified while preparing it.

## What the pipeline does

For one transporter family, the workflow:

1. checks the family directory and creates blank metadata files;
2. searches each plant proteome with the reference proteins using BLASTp;
3. retains up to 10 hits per bait and species after the current identity and alignment-length filters;
4. combines and normalizes the selected proteins and reference sequences;
5. aligns the proteins with MAFFT;
6. trims low-occupancy alignment columns;
7. builds FastTree and IQ-TREE trees;
8. validates sequence counts, identifiers and tree membership;
9. creates reference and Amborella annotations for iTOL;
10. reports the phylogenetic neighborhoods around characterized references;
11. pauses for a manual candidate-group decision;
12. extracts the accepted edge splits, creates a visualization-only tree and prepares final iTOL files.

The current stage order is:

```text
init -> hits -> dataset -> alignment -> trim -> fasttree -> iqtree -> validate
     -> itol-base -> neighborhoods -> candidate-splits -> candidate-view
     -> itol-final -> report
```

## Current scientific settings

The review copy preserves the settings used in the current analyses:

- BLASTp e-value: `0.001`
- minimum identity: `30%`
- minimum alignment length: `30` amino acids
- selected hits: up to `10` per bait and species
- BLAST does not use `--max_target_seqs`
- trimming occupancy: `0.1`, meaning columns with more than 90% gaps are removed
- FastTree: `-lg -gamma`
- IQ-TREE model: `LG+G`

The pipeline does not currently request IQ-TREE bootstrap or SH-aLRT support. Changing models, thresholds, hit selection or support calculations requires a separate scientific decision and was deliberately not done while preparing this repository.

## Required software

The Python code was reviewed with Python 3.12 and requires:

```text
biopython==1.88
```

The full workflow also needs these programs on `PATH`:

- NCBI BLAST+ (`blastp` and `makeblastdb`)
- MAFFT
- FastTree
- IQ-TREE 2 (`iqtree2`)

Create a local environment with:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

If the FastTree executable is not named `FastTree` or is not on `PATH`, pass `--fasttree-bin /path/to/FastTree` or set `FASTTREE_BIN`.

## Required inputs

### Reference FASTA

Each analysis needs one protein FASTA containing the characterized bait sequences:

```text
FAMILY/baits/FAMILY_baits.fasta
```

FASTA identifiers must be unique. The actual ABCC, MATE and P3A sequences are intentionally not included in this review repository; see [references/README.md](references/README.md).

### Reference metadata

The required table is:

```text
FAMILY/metadata/FAMILY_references.tsv
```

It has these columns:

```text
tree_id species gene_id protein_accession family functional_context evidence_type evidence_level use_as_bait reference notes
```

`tree_id` must match the bait identifier as it appears in the normalized tree. `family` must match the value passed to `--family`. Rows used as functional references must have `use_as_bait` set to `yes`.

### Plant proteomes

`--proteomes` must point to a directory containing one protein FASTA per species. Supported filename extensions are `.fasta`, `.fa` and `.faa`, including uppercase variants.

The legacy homolog-collection script uses the part of each filename before the first period as the species/dataset directory name. Those prefixes therefore need to be unique.

The 260-species plant proteome collection used for ABCC, MATE and P3A is external and is not included in Git.

### Candidate groups

After reviewing the neighborhood reports, create or edit:

```text
FAMILY/metadata/FAMILY_candidate_groups.tsv
```

The columns are:

```text
clade_id anchors color label
```

`anchors` is a semicolon-separated list of reference and/or candidate tree IDs. This file is a manual, auditable scientific decision. The software checks the chosen splits for missing anchors and overlap, but it does not decide which lineages are biologically valid.

## Start a new family

The commands below use `XYZ` as a placeholder.

Create the blank family structure:

```bash
python run_family.py --family XYZ --workdir /data/analyses/XYZ --stage init
```

This creates the result directories and these two header-only files:

```text
/data/analyses/XYZ/metadata/XYZ_references.tsv
/data/analyses/XYZ/metadata/XYZ_candidate_groups.tsv
```

Next:

1. copy the reference proteins to `/data/analyses/XYZ/baits/XYZ_baits.fasta`;
2. complete `XYZ_references.tsv`;
3. keep `XYZ_candidate_groups.tsv` header-only for now;
4. check that the proteome directory contains the expected FASTA files.

Run through the neighborhood decision point:

```bash
python run_family.py --family XYZ --family-dir /data/analyses/XYZ --proteomes /data/proteomes --threads 28 --to-stage neighborhoods --resume
```

Inspect the files in `results/qc/`, especially:

```text
XYZ_reference_pairs.tsv
XYZ_reference_amborella.tsv
XYZ_reference_neighborhoods.tsv
XYZ_candidate_decision_report.md
```

After manually completing `metadata/XYZ_candidate_groups.tsv`, continue with:

```bash
python run_family.py --family XYZ --family-dir /data/analyses/XYZ --proteomes /data/proteomes --threads 28 --from-stage candidate-splits --resume
```

For long VM jobs, the same command can be started with `nohup` and monitored through its log. Do not launch a second copy while the first is still running.

## Resume and provenance behavior

`--resume` skips a stage only when its manifest record, input fingerprint, parameters and expected outputs still match.

Existing expensive output is not silently overwritten. For a legacy analysis, `--adopt-existing` validates compatible output and records it in the manifest. Use `--force` only for an intentional rerun.

The manifest is written to:

```text
FAMILY/results/qc/FAMILY_pipeline_manifest.json
```

External commands are recorded in:

```text
FAMILY/logs/FAMILY_pipeline_commands.log
```

## Family directory layout

```text
FAMILY/
|-- baits/
|   `-- FAMILY_baits.fasta
|-- metadata/
|   |-- FAMILY_references.tsv
|   |-- FAMILY_candidate_groups.tsv
|   `-- FAMILY_candidate_clades.tsv              generated later
|-- logs/                                        generated
|-- archive/                                     generated when output is replaced
`-- results/                                     generated
    |-- hits/
    |-- dataset/
    |-- alignment/
    |-- trees/
    |-- qc/
    |-- itol/
    `-- clades/
```

Family work directories can live outside this repository, which is recommended for large analyses.

## Main outputs

- selected homologues and normalized identifier mapping;
- combined candidate/reference FASTA;
- raw and trimmed alignments;
- FastTree and IQ-TREE trees;
- dataset and tree validation summaries;
- reference-pair, Amborella and nearest-neighbor tables;
- accepted candidate-split membership lists;
- a visualization-only candidate-view tree;
- iTOL reference, Amborella and candidate-lineage annotations;
- a stage manifest and short pipeline report.

The original IQ-TREE tree remains the primary inferred tree. The candidate-view tree is rerooted only to display the accepted splits clearly; its root has no biological interpretation.

## Manual steps that remain

The workflow does not automate:

- choosing and checking characterized reference proteins;
- curating evidence fields in the reference metadata;
- deciding whether references define one or several candidate lineages;
- completing `FAMILY_candidate_groups.tsv`;
- interpreting taxonomic composition or functional relevance;
- selecting focused follow-up trees or wet-lab candidates;
- uploading and styling files in iTOL.

These are documented as manual steps because they require biological judgment.


## Examples and tests

`configs/template/` contains blank metadata schemas. `configs/examples/` contains small metadata examples from the initial ABCC, MATE and P3A use cases. No protein sequences or generated trees are included.

Run the standalone tests with:

```bash
python -m unittest discover -s tests -v
```

The root-independent split tests and initialization test need no external dataset or phylogeny programs. The ABCC integration test is skipped unless an external analysis root is supplied:

```bash
PHYLOGENY_PROJECT_ROOT=/path/to/analysis/root python -m unittest tests.test_abcc_integration -v
```

A full end-to-end example is not bundled because the real workflow requires external proteomes and all phylogeny executables. See [examples/README.md](examples/README.md) for the safe review example.

## Current project-specific assumptions

Some parts are reusable but not yet fully generic:

- Amborella sequences receive a dedicated taxonomic annotation;
- evidence levels and iTOL colors follow the current transporter project;
- the input is a directory of one proteome FASTA per species;
- the legacy BLAST script derives a species label from the filename prefix;
- candidate groups are represented as tree edge splits defined by anchors;
- default models, filters and hit counts come from ABCC, MATE and P3A;
- the ABCC regression test needs external ABCC results.

These assumptions are candidates for later review, not silent refactoring.
