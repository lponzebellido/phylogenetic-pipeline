#!/usr/bin/env python3

import argparse
import hashlib
import re
from pathlib import Path

def clean_sequence(sequence):
    return sequence.upper().rstrip("*")

def read_fasta(path):
    header = None
    seq = []

    with open(path) as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    yield header, clean_sequence("".join(seq))

                header = line[1:].strip()
                seq = []
            else:
                seq.append(line)

    if header is not None:
        yield header, clean_sequence("".join(seq))


def normalize_header(header):
    name = re.sub(r"\s+", "_", header)
    name = re.sub(r"[^A-Za-z0-9_.@-]", "_", name)
    return name


def species_from_header(header):
    if "@" in header:
        return header.split("@", 1)[0]
    return ""


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--family", required=True)
    parser.add_argument("--baits", required=True)
    parser.add_argument("--hits-dir", required=True)
    parser.add_argument("--dataset-dir", required=True)
    parser.add_argument("--qc-dir", required=True)

    return parser.parse_args()


def main():
    args = parse_args()

    family = args.family

    baits_file = Path(args.baits)
    hits_dir = Path(args.hits_dir)
    dataset_dir = Path(args.dataset_dir)
    qc_dir = Path(args.qc_dir)

    if not baits_file.is_file():
        raise SystemExit(f"Bait file not found: {baits_file}")

    if not hits_dir.is_dir():
        raise SystemExit(f"Hits directory not found: {hits_dir}")

    dataset_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)

    candidates_file = dataset_dir / f"{family}_candidates.fasta"
    all_file = dataset_dir / f"{family}_all.fasta"

    mapping_file = qc_dir / f"{family}_header_mapping.tsv"
    qc_file = qc_dir / f"{family}_sequence_QC.tsv"
    bait_ids_file = qc_dir / f"{family}_bait_ids.txt"

    hit_files = sorted(hits_dir.rglob("best_hit_sequences.fasta"))

    if not hit_files:
        raise SystemExit(
            f"No best_hit_sequences.fasta files found under {hits_dir}"
        )

    bait_records = list(read_fasta(baits_file))

    if not bait_records:
        raise SystemExit("No bait sequences found")

    bait_headers = [header for header, _ in bait_records]

    if len(bait_headers) != len(set(bait_headers)):
        raise SystemExit("Duplicate bait headers detected")

    candidate_records = []

    for hit_file in hit_files:
        candidate_records.extend(read_fasta(hit_file))

    if not candidate_records:
        raise SystemExit("No candidate sequences found")

    with open(candidates_file, "w") as out:
        for header, sequence in candidate_records:
            out.write(f">{header}\n{sequence}\n")

    records = []

    for source, source_records in (
        ("bait", bait_records),
        ("candidate", candidate_records),
    ):
        for header, sequence in source_records:
            if not sequence:
                raise SystemExit(f"Empty sequence: {header}")

            records.append(
                {
                    "source": source,
                    "original": header,
                    "sequence": sequence,
                }
            )

    seen = {}

    for record in records:
        base = normalize_header(record["original"])

        seen[base] = seen.get(base, 0) + 1

        if seen[base] == 1:
            normalized = base
        else:
            normalized = f"{base}_dup{seen[base]}"

        record["base"] = base
        record["normalized"] = normalized

    bait_hashes = {
        hashlib.sha256(sequence.encode()).hexdigest()
        for _, sequence in bait_records
    }

    bait_headers_first_token = {
        header.split()[0]
        for header, _ in bait_records
    }

    with open(all_file, "w") as out:
        for record in records:
            out.write(
                f">{record['normalized']}\n"
                f"{record['sequence']}\n"
            )

    with open(mapping_file, "w") as out:
        out.write(
            "normalized_id\tsource\toriginal_header\n"
        )

        for record in records:
            out.write(
                f"{record['normalized']}\t"
                f"{record['source']}\t"
                f"{record['original']}\n"
            )

    with open(bait_ids_file, "w") as out:
        for record in records:
            if record["source"] == "bait":
                out.write(record["normalized"] + "\n")

    with open(qc_file, "w") as out:
        out.write(
            "normalized_id\t"
            "source\t"
            "species\t"
            "length\t"
            "sha256\t"
            "same_header_as_bait\t"
            "same_sequence_as_bait\t"
            "renamed_duplicate\n"
        )

        for record in records:
            digest = hashlib.sha256(
                record["sequence"].encode()
            ).hexdigest()

            first_token = record["original"].split()[0]

            same_header = (
                record["source"] == "candidate"
                and first_token in bait_headers_first_token
            )

            same_sequence = (
                record["source"] == "candidate"
                and digest in bait_hashes
            )

            renamed_duplicate = (
                record["normalized"] != record["base"]
            )

            out.write(
                f"{record['normalized']}\t"
                f"{record['source']}\t"
                f"{species_from_header(record['original'])}\t"
                f"{len(record['sequence'])}\t"
                f"{digest}\t"
                f"{int(same_header)}\t"
                f"{int(same_sequence)}\t"
                f"{int(renamed_duplicate)}\n"
            )

    print(f"Hit files: {len(hit_files)}")
    print(f"Baits: {len(bait_records)}")
    print(f"Candidates: {len(candidate_records)}")
    print(f"Total: {len(records)}")


if __name__ == "__main__":
    main()
