import importlib.util
import pathlib
import unittest

import networkx as nx

MODULE_PATH = pathlib.Path(__file__).parents[1] / "02_match_combinations.py"
SPEC = importlib.util.spec_from_file_location("matcher", MODULE_PATH)
MATCHER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MATCHER)


class MatcherTests(unittest.TestCase):
    def test_induced_label_distinguishes_edge_counts(self):
        graph = nx.Graph()
        graph.add_nodes_from(["a", "b", "c"])
        self.assertEqual(MATCHER.canonical_induced_label(graph, ("a", "b", "c")), "000")
        graph.add_edge("a", "b")
        self.assertEqual(MATCHER.canonical_induced_label(graph, ("a", "b", "c")), "100")
        graph.add_edge("a", "c")
        self.assertEqual(MATCHER.canonical_induced_label(graph, ("a", "b", "c")), "110")

    def test_community_pattern_is_label_invariant(self):
        communities = {"a": 7, "b": 7, "c": 9}
        self.assertEqual(MATCHER.community_pattern(communities, ("a", "b", "c")), "2+1")

    def test_pair_metrics(self):
        graph = nx.path_graph(["a", "b", "c"])
        neighbors = {n: set(graph.neighbors(n)) for n in graph.nodes}
        common, jaccard = MATCHER.pair_metrics(graph, neighbors, "a", "c")
        self.assertEqual(common, 1)
        self.assertEqual(jaccard, 1.0)

    def test_bin_edges_collapse_ties(self):
        edges = MATCHER.quantile_edges([0, 0, 0, 1, 2])
        self.assertEqual(edges, sorted(set(edges)))
        self.assertGreaterEqual(MATCHER.bin_value(2, edges), MATCHER.bin_value(0, edges))


if __name__ == "__main__":
    unittest.main()
