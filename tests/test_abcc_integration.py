from __future__ import annotations

import itertools
import os
import tempfile
import unittest
from pathlib import Path

from Bio import Phylo

from anthotransport.annotations import write_itol_annotations
from anthotransport.candidate_view import create_candidate_view
from anthotransport.splits import identify_candidate_splits, leaf_set


class AbccIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        project_root = os.environ.get("PHYLOGENY_PROJECT_ROOT")
        if not project_root:
            raise unittest.SkipTest(
                "set PHYLOGENY_PROJECT_ROOT to run the external ABCC integration test"
            )
        cls.project = Path(project_root)
        cls.abcc = cls.project / "ABCC"
        required = [
            cls.abcc / "results/trees/ABCC_IQTREE.treefile",
            cls.abcc / "metadata/ABCC_candidate_groups.tsv",
            cls.abcc / "metadata/ABCC_references.tsv",
            cls.abcc / "results/qc/ABCC_header_mapping.tsv",
        ]
        if not all(path.is_file() for path in required):
            raise unittest.SkipTest("ABCC reference outputs are not available")

    def test_validated_abcc_downstream_results(self):
        original_path = self.abcc / "results/trees/ABCC_IQTREE.treefile"
        groups = self.abcc / "metadata/ABCC_candidate_groups.tsv"
        references = self.abcc / "metadata/ABCC_references.tsv"
        mapping = self.abcc / "results/qc/ABCC_header_mapping.tsv"
        original = Phylo.read(str(original_path), "newick")
        original_leaves = leaf_set(original.root)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            clades = root / "clades"
            rows = identify_candidate_splits(original_path, groups, clades, references)
            sizes = {row["clade_id"]: int(row["split_size"]) for row in rows}
            self.assertEqual(
                sizes,
                {
                    "ABCC_C1": 1129,
                    "ABCC_C2": 1186,
                    "ABCC_C3": 515,
                    "ABCC_C4": 1037,
                    "ABCC_C5": 881,
                },
            )
            member_sets = {
                clade_id: set((clades / f"{clade_id}_members.txt").read_text().splitlines())
                for clade_id in sizes
            }
            for a, b in itertools.combinations(member_sets, 2):
                self.assertFalse(member_sets[a] & member_sets[b])
            union = set().union(*member_sets.values())
            self.assertEqual(len(union), 4748)
            self.assertEqual(len(original_leaves - union), 1393)

            view = root / "ABCC_IQTREE_candidate_view.treefile"
            candidate_clades = root / "ABCC_candidate_clades.tsv"
            root_info = root / "ABCC_candidate_view_root.tsv"
            created = create_candidate_view(
                original_path,
                groups,
                clades,
                view,
                candidate_clades,
                root_info,
                references,
            )
            self.assertEqual(len(created), 5)
            self.assertEqual(leaf_set(Phylo.read(str(view), "newick").root), original_leaves)

            counts = write_itol_annotations(
                "ABCC", view, references, root / "itol", mapping, candidate_clades
            )
            self.assertEqual(
                counts, {"references": 9, "amborella": 12, "candidate_clades": 5}
            )
            combined = (root / "itol/ABCC_iTOL_combined_styles.txt").read_text()
            first_candidate = combined.index("branch\tclade")
            first_reference = combined.index("branch\tnode")
            self.assertLess(first_candidate, first_reference)


if __name__ == "__main__":
    unittest.main()

