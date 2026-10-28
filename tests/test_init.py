from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from anthotransport.paths import FamilyPaths
from anthotransport.pipeline import Pipeline, PipelineConfig


class InitStageTests(unittest.TestCase):
    def test_init_creates_blank_family_metadata(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "XYZ"
            paths = FamilyPaths(
                family="XYZ",
                root=root,
                baits=root / "baits/XYZ_baits.fasta",
                references=root / "metadata/XYZ_references.tsv",
                candidate_groups=root / "metadata/XYZ_candidate_groups.tsv",
            )
            config = PipelineConfig(
                paths=paths,
                pipeline_root=Path(__file__).resolve().parents[1],
                proteomes=None,
            )
            Pipeline(config).stage_init()

            self.assertEqual(
                paths.references.read_text().splitlines()[0],
                "tree_id\tspecies\tgene_id\tprotein_accession\tfamily\t"
                "functional_context\tevidence_type\tevidence_level\tuse_as_bait\t"
                "reference\tnotes",
            )
            self.assertEqual(
                paths.candidate_groups.read_text(),
                "clade_id\tanchors\tcolor\tlabel\n",
            )
            self.assertTrue((root / "results").is_dir())
            self.assertTrue((root / "logs").is_dir())


if __name__ == "__main__":
    unittest.main()
