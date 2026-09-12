# NBVK implementation audit

## 2026-07-26 pre-run checks

1. Python/Scipy all-pairs shortest paths required approximately 0.301 seconds per 1,226-node graph; existing R igraph 2.3.3 required approximately 0.037 seconds. The implementation therefore uses R/igraph with the identical frozen graph and distance definitions.
2. Singleton knockout damage is precomputed once for every eligible node and reused in synergy calculations. This is an exact computational cache and does not change any statistic.
3. The first primary-script launch stopped before permutations because:
   - the sandbox returned `NA` for `parallel::detectCores()`;
   - a pre-run assertion used an estimated LCC edge count of 13,412.
4. Read-only component verification established 1,226 nodes and 13,416 edges in the largest component; the remaining seven connected nodes contain four edges. The assertion was corrected to the observed frozen-input count.
5. Worker count was fixed at eight. These are implementation corrections made before any null distributions or inferential results were produced. There is no analytical deviation from the pre-registration.
6. The second launch calculated the seven observed knockouts, then stopped before null generation because pandas-written Boolean values (`True`/`False`) were read by R as character strings. `make_topology_pools()` was corrected to parse this field explicitly. No formula, pool definition, observed set or decision rule was changed; the complete analysis is rerun from the beginning.
7. Independent post-run verification reproduced every empirical P value from the saved 10,000-row null files but detected incorrect multiplicity code: `ave(..., method="BH")` treated `method` as an additional grouping argument. The code was replaced by explicit within-baseline `p.adjust(..., method="BH")`, a regression test was added, and inference/decision tables were regenerated from the unchanged complete null distributions. Observed effects, sampled sets, null values and empirical P values were not modified.
8. During the weighted-network sensitivity run, repeated Monte Carlo draws were being recomputed despite representing identical node sets. The incomplete weighted group was stopped; no weighted null file or inference table had been completed. The evaluator was changed to compute each unique sampled node set once and map its exact result back to all original 10,000 draws. This exact cache preserves the sampled rows, empirical distribution and registered statistics. Completed score>=0.400 outputs were retained.
