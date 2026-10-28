#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 FAMILY BAITS SUBJECTDIR" >&2
    exit 1
fi

FAMILY="$1"

if [[ ! -f "$2" ]]; then
    echo "Bait file not found: $2" >&2
    exit 1
fi

if [[ ! -d "$3" ]]; then
    echo "Subject directory not found: $3" >&2
    exit 1
fi

BAITS="$(realpath "$2")"
SUBJECTDIR="$(realpath "$3")"

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_DIR="$ROOT_DIR/scripts"
RESULTS_DIR="$ROOT_DIR/results"

COLLECT_SCRIPT="$SCRIPT_DIR/collect_best_BLAST_hits.py"
SUMMARIZE_SCRIPT="$SCRIPT_DIR/summarize_hits.py"
PREPARE_SCRIPT="$SCRIPT_DIR/prepare_dataset.py"
TRIM_SCRIPT="$SCRIPT_DIR/algntrim.py"
VALIDATE_SCRIPT="$SCRIPT_DIR/validate_dataset.py"
ITOL_SCRIPT="$SCRIPT_DIR/make_itol_annotations.py"

METADATA="$ROOT_DIR/metadata/${FAMILY}_references.tsv"

FASTTREE_BIN="${FASTTREE_BIN:-/vol/tests/test1/FastTree}"
IQTREE_MODEL="${IQTREE_MODEL:-LG+G}"

HITS_PER_BAIT="${HITS_PER_BAIT:-10}"
SIMCUT="${SIMCUT:-30}"
LENCUT="${LENCUT:-30}"
COLLECT_CPU="${COLLECT_CPU:-4}"
REUSE_HITS="${REUSE_HITS:-0}"

HITS_DIR="$RESULTS_DIR/hits"
DATASET_DIR="$RESULTS_DIR/dataset"
ALIGNMENT_DIR="$RESULTS_DIR/alignment"
TREE_DIR="$RESULTS_DIR/trees"
QC_DIR="$RESULTS_DIR/qc"
ITOL_DIR="$RESULTS_DIR/itol"

CANDIDATES="$DATASET_DIR/${FAMILY}_candidates.fasta"
ALL="$DATASET_DIR/${FAMILY}_all.fasta"

SELECTED_HITS="$QC_DIR/${FAMILY}_selected_hits.tsv"
MAPPING="$QC_DIR/${FAMILY}_header_mapping.tsv"
BAIT_IDS="$QC_DIR/${FAMILY}_bait_ids.txt"
COUNTS="$QC_DIR/${FAMILY}_sequence_counts.tsv"

ALIGNED="$ALIGNMENT_DIR/${FAMILY}_aligned.fasta"
TRIMMED="$ALIGNMENT_DIR/${FAMILY}_aligned_trimmed.fasta"

FASTTREE_OUT="$TREE_DIR/${FAMILY}_FastTree.nwk"
FASTTREE_LOG="$TREE_DIR/${FAMILY}_FastTree.log"

IQTREE_PREFIX="$TREE_DIR/${FAMILY}_IQTREE"
IQTREE_TREE="${IQTREE_PREFIX}.treefile"

for file in \
    "$COLLECT_SCRIPT" \
    "$SUMMARIZE_SCRIPT" \
    "$PREPARE_SCRIPT" \
    "$TRIM_SCRIPT" \
    "$VALIDATE_SCRIPT" \
    "$ITOL_SCRIPT" \
    "$METADATA"
do
    if [[ ! -f "$file" ]]; then
        echo "Required file not found: $file" >&2
        exit 1
    fi
done

for program in python3 mafft iqtree2
do
    if ! command -v "$program" >/dev/null 2>&1; then
        echo "Required program not found: $program" >&2
        exit 1
    fi
done

if [[ ! -x "$FASTTREE_BIN" ]]; then
    echo "FastTree executable not found: $FASTTREE_BIN" >&2
    exit 1
fi

mkdir -p \
    "$HITS_DIR" \
    "$DATASET_DIR" \
    "$ALIGNMENT_DIR" \
    "$TREE_DIR" \
    "$QC_DIR" \
    "$ITOL_DIR"

if [[ "$REUSE_HITS" != "1" ]]; then
    rm -rf "$HITS_DIR"
    mkdir -p "$HITS_DIR"
fi

echo
echo " ${FAMILY}: homolog collection "

python3 "$COLLECT_SCRIPT" \
    --baits "$BAITS" \
    --out "$HITS_DIR/" \
    --subjectdir "$SUBJECTDIR" \
    --number "$HITS_PER_BAIT" \
    --simcut "$SIMCUT" \
    --lencut "$LENCUT" \
    --cpu "$COLLECT_CPU"

echo
echo " ${FAMILY}: BLAST hit summary "

python3 "$SUMMARIZE_SCRIPT" \
    --baits "$BAITS" \
    --subject-dir "$SUBJECTDIR" \
    --hits-dir "$HITS_DIR" \
    --out "$SELECTED_HITS" \
    --number "$HITS_PER_BAIT" \
    --simcut "$SIMCUT" \
    --lencut "$LENCUT"

echo
echo " ${FAMILY}: dataset preparation "

python3 "$PREPARE_SCRIPT" \
    --family "$FAMILY" \
    --baits "$BAITS" \
    --hits-dir "$HITS_DIR" \
    --dataset-dir "$DATASET_DIR" \
    --qc-dir "$QC_DIR"

echo
echo " ${FAMILY}: MAFFT alignment "

mafft \
    --thread -1 \
    "$ALL" \
    > "$ALIGNED" \
    2> "$ALIGNMENT_DIR/${FAMILY}_mafft.log"

echo
echo " ${FAMILY}: alignment trimming "

python3 "$TRIM_SCRIPT" \
    --in "$ALIGNED" \
    --out "$TRIMMED" \
    > "$ALIGNMENT_DIR/${FAMILY}_trim.log" \
    2>&1

echo
echo " ${FAMILY}: FastTree "

"$FASTTREE_BIN" \
    -lg \
    -gamma \
    -log "$FASTTREE_LOG" \
    "$TRIMMED" \
    > "$FASTTREE_OUT"

echo
echo " ${FAMILY}: IQ-TREE "

rm -f "${IQTREE_PREFIX}".*

iqtree2 \
    -s "$TRIMMED" \
    -m "$IQTREE_MODEL" \
    -T AUTO \
    --prefix "$IQTREE_PREFIX"

echo
echo " ${FAMILY}: validation "

python3 "$VALIDATE_SCRIPT" \
    --all "$ALL" \
    --aligned "$ALIGNED" \
    --trimmed "$TRIMMED" \
    --bait-ids "$BAIT_IDS" \
    --counts-out "$COUNTS" \
    --fasttree "$FASTTREE_OUT" \
    --iqtree "$IQTREE_TREE"

echo
echo "${FAMILY}: iTOL annotations"

python3 "$ITOL_SCRIPT" \
    --family "$FAMILY" \
    --tree "$IQTREE_TREE" \
    --metadata "$METADATA" \
    --mapping "$MAPPING" \
    --out-dir "$ITOL_DIR"

echo
echo "Analysis complete"
echo "Family: $FAMILY"
echo "Baits: $(grep -c '^>' "$BAITS")"
echo "Candidates: $(grep -c '^>' "$CANDIDATES")"
echo "Total sequences: $(grep -c '^>' "$ALL")"
echo "Alignment: $ALIGNED"
echo "Trimmed alignment: $TRIMMED"
echo "FastTree: $FASTTREE_OUT"
echo "IQ-TREE: $IQTREE_TREE"
echo "Selected hits: $SELECTED_HITS"
echo "QC: $QC_DIR"
echo "iTOL annotations: $ITOL_DIR"
