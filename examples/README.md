# Safe review example

This repository provides a structural example rather than a bundled full run.

To create a blank example family:

```bash
python run_family.py --family EXAMPLE --workdir /tmp/EXAMPLE --stage init
```

That command creates the expected directory layout and header-only metadata files. It does not need BLAST, MAFFT, FastTree, IQ-TREE or a proteome dataset.

A complete run is not bundled because it would require protein sequences, a per-species proteome collection and all external phylogeny programs. The standalone unit tests use small synthetic trees and temporary files instead.

Run them with:

```bash
python -m unittest discover -s tests -v
```
