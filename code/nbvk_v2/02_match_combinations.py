from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
import random
import sys
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

SEED = 20260726
ANCHORS = ("ARG1", "LCN2", "LTF")


def canonical_induced_label(graph: nx.Graph, nodes: tuple[str, ...]) -> str:
    values = [
        int(graph.has_edge(nodes[i], nodes[j]))
        for i, j in ((0, 1), (0, 2), (1, 2))
    ]
    return "".join(map(str, sorted(values, reverse=True)))


def community_pattern(communities: dict[str, int], nodes: tuple[str, ...]) -> str:
    counts: dict[int, int] = {}
    for node in nodes:
        counts[communities[node]] = counts.get(communities[node], 0) + 1
    return "+".join(map(str, sorted(counts.values(), reverse=True)))


def quantile_edges(values: list[float]) -> list[float]:
    edges = np.quantile(np.asarray(values, dtype=float), [0, .2, .4, .6, .8, 1])
    return sorted(set(float(x) for x in edges))


def bin_value(value: float, edges: list[float]) -> int:
    if len(edges) <= 1:
        return 0
    return int(np.searchsorted(np.asarray(edges[1:-1]), value, side="right"))


def pair_metrics(
    graph: nx.Graph, neighbor_sets: dict[str, set[str]], a: str, b: str
) -> tuple[int, float]:
    common = len(neighbor_sets[a] & neighbor_sets[b])
    union = len(neighbor_sets[a] | neighbor_sets[b])
    return common, common / union if union else 0.0


def make_pools(features: pd.DataFrame, pool_size: int = 100) -> dict[str, list[str]]:
    cols = ["log1p_degree", "betweenness_percentile", "participation"]
    z = features[cols].copy()
    z = (z - z.mean()) / z.std(ddof=1)
    z.index = features["gene"]
    pools: dict[str, list[str]] = {}
    for anchor in ANCHORS:
        candidates = z.drop(index=list(ANCHORS))
        dist = np.sqrt(((candidates - z.loc[anchor]) ** 2).sum(axis=1))
        ranked = pd.DataFrame({"gene": dist.index, "distance": dist.values})
        ranked = ranked.sort_values(["distance", "gene"], kind="mergesort")
        pools[anchor] = ranked.head(pool_size)["gene"].tolist()
    return pools


def main(root: Path) -> None:
    derived = root / "results/nbvk_v2/derived"
    derived.mkdir(parents=True, exist_ok=True)
    features = pd.read_csv(derived / "node_features_score700.csv")
    edges = pd.read_csv(derived / "lcc_edges_score700.csv")
    graph = nx.from_pandas_edgelist(edges, "protein_a", "protein_b")
    communities = dict(zip(features["gene"], features["community"]))
    pools = make_pools(features)

    pool_rows = []
    for anchor, genes in pools.items():
        for rank, gene in enumerate(genes, 1):
            pool_rows.append({"anchor": anchor, "rank": rank, "gene": gene})
    pd.DataFrame(pool_rows).to_csv(derived / "anchor_pools.csv", index=False)

    neighbor_sets = {n: set(graph.neighbors(n)) for n in graph.nodes}
    nodes = sorted(graph.nodes)
    common_values: list[int] = []
    jaccard_values: list[float] = []
    for a_i, a in enumerate(nodes):
        for b in nodes[a_i + 1 :]:
            common, jaccard = pair_metrics(graph, neighbor_sets, a, b)
            common_values.append(common)
            jaccard_values.append(jaccard)
    common_edges = quantile_edges(common_values)
    jaccard_edges = quantile_edges(jaccard_values)

    raw_triples = set()
    for triple in itertools.product(*(pools[a] for a in ANCHORS)):
        if len(set(triple)) == 3:
            raw_triples.add(tuple(sorted(triple)))
    union_values = [
        len(neighbor_sets[a] | neighbor_sets[b] | neighbor_sets[c])
        for a, b, c in raw_triples
    ]
    union_edges = quantile_edges(union_values)

    all_lengths = dict(nx.all_pairs_shortest_path_length(graph))

    def signature(triple: tuple[str, str, str]) -> tuple:
        a, b, c = triple
        pairs = ((a, b), (a, c), (b, c))
        dcat = sorted(
            1 if all_lengths[x][y] == 1 else 2 if all_lengths[x][y] == 2 else 3
            for x, y in pairs
        )
        commons, jaccards = zip(
            *(pair_metrics(graph, neighbor_sets, x, y) for x, y in pairs)
        )
        common_bins = sorted(bin_value(x, common_edges) for x in commons)
        jaccard_bins = sorted(bin_value(x, jaccard_edges) for x in jaccards)
        union_size = len(neighbor_sets[a] | neighbor_sets[b] | neighbor_sets[c])
        return (
            canonical_induced_label(graph, triple),
            tuple(dcat),
            community_pattern(communities, triple),
            tuple(common_bins),
            tuple(jaccard_bins),
            bin_value(union_size, union_edges),
        )

    observed = signature(tuple(ANCHORS))
    signatures = {triple: signature(triple) for triple in sorted(raw_triples)}
    legal = [triple for triple, sig in signatures.items() if sig == observed]
    # Prespecified relaxation 1: drop community occupation only.
    relax1 = [
        triple
        for triple, sig in signatures.items()
        if sig[:2] == observed[:2] and sig[3:] == observed[3:]
    ]
    # Prespecified relaxation 2: retain induced structure and pairwise distances.
    relax2 = [
        triple
        for triple, sig in signatures.items()
        if sig[:2] == observed[:2]
    ]
    pd.DataFrame(legal, columns=["gene1", "gene2", "gene3"]).to_csv(
        derived / "legal_strict_triplets.csv", index=False
    )
    pd.DataFrame(relax1, columns=["gene1", "gene2", "gene3"]).to_csv(
        derived / "legal_relax1_triplets.csv", index=False
    )
    pd.DataFrame(relax2, columns=["gene1", "gene2", "gene3"]).to_csv(
        derived / "legal_relax2_triplets.csv", index=False
    )
    metadata = {
        "seed": SEED,
        "anchors": ANCHORS,
        "raw_unique_triplets": len(raw_triples),
        "strict_legal_triplets": len(legal),
        "relax1_legal_triplets": len(relax1),
        "relax2_legal_triplets": len(relax2),
        "acceptance_rate": len(legal) / len(raw_triples),
        "observed_signature": repr(observed),
        "common_neighbor_bin_edges": common_edges,
        "jaccard_bin_edges": jaccard_edges,
        "triplet_union_bin_edges": union_edges,
        "confirmatory_support": len(legal) >= 1000,
    }
    (derived / "matching_metadata.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(
        f"MATCHING_READY raw={len(raw_triples)} strict={len(legal)} "
        f"relax1={len(relax1)} relax2={len(relax2)} "
        f"support={metadata['confirmatory_support']}"
    )


if __name__ == "__main__":
    repository_root = (
        Path(sys.argv[1]).resolve()
        if len(sys.argv) > 1
        else Path(__file__).resolve().parents[2]
    )
    main(repository_root)
