# Transporter phylogeny workflow

This is a clean copy of the workflow that I am currently using for the ABCC, MATE and P3A transporter analyses.

It is being shared so the current scripts can be inspected.

## What is automated

`run_phylogeny.sh` runs these steps:

1. search each plant proteome with the bait proteins using BLASTp;
2. retain the best homologues that pass the current filters;
3. combine the selected proteins and bait sequences;
4. normalize FASTA identifiers and write QC tables;
5. align the proteins with MAFFT;
6. trim the alignment;
7. build a FastTree tree;
8. build an IQ-TREE tree;
9. validate sequence and tree contents;
10. prepare the basic reference and Amborella files for iTOL.

The current parameters are kept in `run_phylogeny.sh`:

- up to 10 hits per bait and species;
- minimum identity of 30%;
- minimum alignment length of 30 amino acids;
- BLASTp e-value of 0.001;
- FastTree with `-lg -gamma`;
- IQ-TREE with `LG+G`.

## What is still manual

After the main tree is built, the characterized references must be inspected manually.

The additional scripts can be used to:

- inspect the closest tree sequences around each reference;
- compare references with Amborella sequences;
- define candidate groups in `metadata/FAMILY_candidate_groups.tsv`;
- extract the selected candidate splits;
- check that the selected groups do not overlap;
- prepare a rerooted tree for visualization.

These steps require scientific judgment. The scripts do not decide which candidate groups are biologically correct.

The four scripts used for this part were developed during the ABCC analysis:

```text
scripts/inspect_bait_neighborhood.py
scripts/find_amborella_anchors.py
scripts/find_candidate_splits.py
scripts/prepare_candidate_view.py
```

MATE and P3A had not yet received accepted candidate-group tables when this clean copy was prepared.

## Required software

- Python 3
- Biopython
- NCBI BLAST+ (`blastp` and `makeblastdb`)
- MAFFT
- FastTree
- IQ-TREE 2 (`iqtree2`)

The current shell script uses `/vol/tests/test1/FastTree` as its default FastTree location. This can be replaced without editing the script by setting `FASTTREE_BIN` when the analysis is started.

## Required input

For a family called `FAMILY`, provide:

```text
baits/FAMILY_baits.fasta
metadata/FAMILY_references.tsv
/path/to/plant_proteomes/
```

The proteome directory must contain one protein FASTA per species. The accepted extensions are `.fasta`, `.fa` and `.faa`, including uppercase variants.

`metadata/FAMILY_references.tsv` describes the characterized proteins and their evidence. A header-only template is included. Examples from ABCC, MATE and P3A are available in `metadata_examples/`.

## Starting another family

Copy this folder to the location where the new analysis should run. Replace `FAMILY` with the real family name.

Prepare the bait FASTA:

```text
baits/FAMILY_baits.fasta
```

Copy and rename the metadata template:

```bash
cp metadata/FAMILY_references.tsv metadata/NEWFAMILY_references.tsv
```

Complete the reference rows, then start the analysis from the project folder:

```bash
nohup ./run_phylogeny.sh NEWFAMILY baits/NEWFAMILY_baits.fasta /path/to/plant_proteomes > NEWFAMILY_full_run.log 2>&1 &
```

If FastTree is somewhere else:

```bash
FASTTREE_BIN=/path/to/FastTree nohup ./run_phylogeny.sh NEWFAMILY baits/NEWFAMILY_baits.fasta /path/to/plant_proteomes > NEWFAMILY_full_run.log 2>&1 &
```

Follow the run with:

```bash
tail -f NEWFAMILY_full_run.log
```

The script recreates `results/hits/` unless `REUSE_HITS=1` is set. Check that another run is not active before starting it again.

When the tree has been inspected and candidate groups are ready to be defined, copy the second template:

```bash
cp metadata/FAMILY_candidate_groups.tsv metadata/NEWFAMILY_candidate_groups.tsv
```

## Generated output

The main generated directories are:

```text
results/hits/        BLAST databases, raw searches and selected sequences
results/dataset/     combined candidate and reference FASTA files
results/alignment/   MAFFT and trimmed alignments
results/trees/       FastTree and IQ-TREE files
results/qc/          mappings, selected-hit tables and validation output
results/itol/        iTOL annotation files
results/clades/      candidate members after manual group definition
```

Generated results, comprehensive proteomes, alignments, trees, logs and BLAST databases are intentionally not included in this clean copy.

## Files in this copy

- `run_phylogeny.sh` and the six common scripts were identical in ABCC, MATE and P3A.
- The four candidate-inspection scripts were copied from ABCC without modification.
- `metadata_examples/` contains only the small reference tables and the accepted ABCC candidate-group table.
- No sequence dataset or completed analysis output is included.
