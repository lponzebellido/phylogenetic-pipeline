from __future__ import annotations

import copy
import unittest
from io import StringIO

from Bio import Phylo

from anthotransport.neighborhoods import _leaf_distances
from anthotransport.splits import find_smallest_split


class RootIndependentSplitTests(unittest.TestCase):
    def setUp(self):
        self.tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);"), "newick")

    def test_smallest_split_is_independent_of_serialized_root(self):
        first = find_smallest_split(self.tree, ["D", "E"])
        rerooted = copy.deepcopy(self.tree)
        rerooted.root_with_outgroup("A")
        second = find_smallest_split(rerooted, ["D", "E"])
        self.assertEqual(first.members, frozenset({"D", "E"}))
        self.assertEqual(first.members, second.members)

    def test_complement_side_is_considered(self):
        rerooted = copy.deepcopy(self.tree)
        rerooted.root_with_outgroup("A")
        result = find_smallest_split(rerooted, ["D", "E"])
        self.assertEqual(result.members, frozenset({"D", "E"}))

    def test_missing_anchor_fails(self):
        with self.assertRaisesRegex(ValueError, "missing"):
            find_smallest_split(self.tree, ["A", "NOT_PRESENT"])

    def test_linear_distance_traversal_matches_biopython(self):
        observed = _leaf_distances(self.tree, "A")
        for leaf in ("B", "C", "D", "E"):
            self.assertAlmostEqual(observed[leaf], self.tree.distance("A", leaf))


if __name__ == "__main__":
    unittest.main()

