from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
PRIMARY = ROOT / "results/nbvk_v1/primary"
OUT = ROOT / "results/nbvk_v1/exploratory_posthoc"
OUT.mkdir(parents=True, exist_ok=True)

edges_frame = pd.read_csv(
    ROOT / "data/network_inputs/string_edges_score700.csv"
)
edges = {
    tuple(sorted(pair))
    for pair in edges_frame[["protein_a", "protein_b"]].itertuples(
        index=False, name=None
    )
}
adjacency: dict[str, set[str]] = {}
for left, right in edges:
    adjacency.setdefault(left, set()).add(right)
    adjacency.setdefault(right, set()).add(left)

observed = pd.read_csv(PRIMARY / "observed_nbvk.csv").set_index("set_name")
combinations = [
    "ARG1_LCN2",
    "ARG1_LTF",
    "LCN2_LTF",
    "ARG1_LCN2_LTF",
]


def internal_edge_count(nodes: list[str]) -> int:
    return sum(
        tuple(sorted((nodes[i], nodes[j]))) in edges
        for i in range(len(nodes))
        for j in range(i + 1, len(nodes))
    )


rows = []
for set_name in combinations:
    observed_nodes = set_name.split("_")
    observed_internal_edges = internal_edge_count(observed_nodes)
    null = pd.read_csv(
        PRIMARY / f"null_distributions/topology_{set_name}.csv"
    )
    sampled_nodes = null["sampled_set"].str.split("+", regex=False)
    null["internal_edge_count"] = sampled_nodes.map(internal_edge_count)
    conditioned = null[
        null["internal_edge_count"].eq(observed_internal_edges)
    ]
    exceedances = int(
        (
            conditioned["Syn_Bliss"]
            >= observed.loc[set_name, "Syn_Bliss"]
        ).sum()
    )
    p_value = (1 + exceedances) / (len(conditioned) + 1)
    rows.append(
        {
            "set_name": set_name,
            "observed_internal_edge_count": observed_internal_edges,
            "conditioned_draws": len(conditioned),
            "exceedances": exceedances,
            "p_empirical_right_plus1": p_value,
            "unique_sampled_sets_unconditioned": null[
                "sampled_set"
            ].nunique(),
            "analysis_status": (
                "post hoc exploratory stress test; not preregistered"
            ),
        }
    )

result = pd.DataFrame(rows)
result["FDR_BH"] = np.nan
order = np.argsort(result["p_empirical_right_plus1"].to_numpy())
ranked = (
    result["p_empirical_right_plus1"].to_numpy()[order]
    * len(result)
    / np.arange(1, len(result) + 1)
)
ranked = np.minimum.accumulate(ranked[::-1])[::-1]
adjusted = np.empty(len(result))
adjusted[order] = np.minimum(ranked, 1)
result["FDR_BH"] = adjusted
result.to_csv(
    OUT / "internal_edge_conditioned_synergy.csv",
    index=False,
    encoding="utf-8-sig",
)

pair_context = pd.DataFrame(
    [
        {
            "pair": "LCN2+LTF",
            "direct_edge": tuple(sorted(("LCN2", "LTF"))) in edges,
            "shared_neighbor_count": len(
                adjacency["LCN2"] & adjacency["LTF"]
            ),
            "shared_neighbors": "|".join(
                sorted(adjacency["LCN2"] & adjacency["LTF"])
            ),
            "triple_synergy_fraction_from_pair": (
                observed.loc["LCN2_LTF", "Syn_Bliss"]
                / observed.loc["ARG1_LCN2_LTF", "Syn_Bliss"]
            ),
        }
    ]
)
pair_context.to_csv(
    OUT / "LCN2_LTF_pair_context.csv",
    index=False,
    encoding="utf-8-sig",
)
print("POSTHOC_COMBINATION_STRESS_TEST_COMPLETE")
