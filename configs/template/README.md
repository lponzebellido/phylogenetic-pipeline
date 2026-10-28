# Blank family metadata

These are header-only copies of the two manually maintained TSV files.

Normally you do not need to copy them yourself. Run:

```bash
python run_family.py --family XYZ --workdir /data/analyses/XYZ --stage init
```

The command creates correctly named files in `/data/analyses/XYZ/metadata/`.

Complete `XYZ_references.tsv` before running the phylogeny. Keep `XYZ_candidate_groups.tsv` header-only until you have reviewed the neighborhood reports.
