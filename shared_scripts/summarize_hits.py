#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict
from pathlib import Path

def read_fasta_lengths(path):
    lengths = {}
    header = None
    seq = []

    with open(path) as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    sequence = "".join(seq).upper().rstrip("*")
                    lengths[header] = len(sequence)

                header = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)

    if header is not None:
        sequence = "".join(seq).upper().rstrip("*")
        lengths[header] = len(sequence)

    return lengths

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--baits", required=True)
    parser.add_argument("--subject-dir", required=True)
    parser.add_argument("--hits-dir", required=True)
    parser.add_argument("--out", required=True)

    parser.add_argument("--number", type=int, default=10)
    parser.add_argument("--simcut", type=float, default=30.0)
    parser.add_argument("--lencut", type=int, default=30)

    return parser.parse_args()


def main():
    args = parse_args()

    baits_file = Path(args.baits)
    subject_dir = Path(args.subject_dir)
    hits_dir = Path(args.hits_dir)
    output_file = Path(args.out)

    if not baits_file.is_file():
        raise SystemExit(f"Bait file not found: {baits_file}")

    if not subject_dir.is_dir():
        raise SystemExit(f"Subject directory not found: {subject_dir}")

    if not hits_dir.is_dir():
        raise SystemExit(f"Hits directory not found: {hits_dir}")

    bait_lengths = read_fasta_lengths(baits_file)

    subject_files = []

    for pattern in (
        "*.fasta",
        "*.fa",
        "*.faa",
        "*.FASTA",
        "*.FA",
        "*.FAA",
    ):
        subject_files.extend(subject_dir.glob(pattern))

    subject_map = {}

    for path in subject_files:
        species = path.name.split(".")[0]

        if species in subject_map:
            raise SystemExit(
                f"Duplicate subject basename: {species}"
            )

        subject_map[species] = path

    result_files = sorted(
        hits_dir.glob("*/blast_results.txt")
    )

    if not result_files:
        raise SystemExit(
            f"No blast_results.txt files found under {hits_dir}"
        )

    output_file.parent.mkdir(parents=True, exist_ok=True)

    rows = []

    for result_file in result_files:
        species = result_file.parent.name

        subject_file = subject_map.get(species)

        if subject_file is None:
            raise SystemExit(
                f"No subject FASTA found for species: {species}"
            )

        subject_lengths = read_fasta_lengths(subject_file)

        hits = defaultdict(list)

        with open(result_file) as f:
            for line in f:
                parts = line.rstrip().split("\t")

                if len(parts) < 12:
                    continue

                query = parts[0]
                subject = parts[1]

                pident = float(parts[2])
                alignment_length = int(parts[3])

                if pident < args.simcut:
                    continue

                if alignment_length < args.lencut:
                    continue

                qstart = int(parts[6])
                qend = int(parts[7])
                sstart = int(parts[8])
                send = int(parts[9])

                evalue = float(parts[10])
                bitscore = float(parts[11])

                qlen = bait_lengths.get(query)
                slen = subject_lengths.get(subject)

                qcov = (
                    (abs(qend - qstart) + 1) / qlen * 100.0
                    if qlen
                    else None
                )

                scov = (
                    (abs(send - sstart) + 1) / slen * 100.0
                    if slen
                    else None
                )

                hits[query].append(
                    {
                        "subject": subject,
                        "pident": pident,
                        "alignment_length": alignment_length,
                        "query_coverage": qcov,
                        "subject_coverage": scov,
                        "evalue": evalue,
                        "bitscore": bitscore,
                    }
                )

        for bait, bait_hits in hits.items():
            bait_hits.sort(
                key=lambda x: x["bitscore"],
                reverse=True,
            )

            for rank, hit in enumerate(
                bait_hits[:args.number],
                start=1,
            ):
                rows.append(
                    {
                        "species": species,
                        "bait": bait,
                        "hit": hit["subject"],
                        "rank": rank,
                        "pident": f"{hit['pident']:.2f}",
                        "alignment_length": hit["alignment_length"],
                        "query_coverage": (
                            f"{hit['query_coverage']:.2f}"
                            if hit["query_coverage"] is not None
                            else "NA"
                        ),
                        "subject_coverage": (
                            f"{hit['subject_coverage']:.2f}"
                            if hit["subject_coverage"] is not None
                            else "NA"
                        ),
                        "evalue": hit["evalue"],
                        "bitscore": hit["bitscore"],
                    }
                )

    fields = [
        "species",
        "bait",
        "hit",
        "rank",
        "pident",
        "alignment_length",
        "query_coverage",
        "subject_coverage",
        "evalue",
        "bitscore",
    ]

    with open(output_file, "w", newline="") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=fields,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(rows)

    print(f"BLAST result files: {len(result_files)}")
    print(f"Selected hit records: {len(rows)}")
    print(f"Output: {output_file}")


if __name__ == "__main__":
    main()
