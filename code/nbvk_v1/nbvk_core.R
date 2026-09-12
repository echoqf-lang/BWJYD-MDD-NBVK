suppressPackageStartupMessages(library(igraph))

build_graph <- function(edge_file, threshold = 0.700, weighted = FALSE) {
  edges <- read.csv(edge_file, check.names = FALSE)
  edges <- edges[edges$combined_score >= threshold, , drop = FALSE]
  edges <- edges[edges$protein_a != edges$protein_b, , drop = FALSE]
  pair_key <- apply(
    cbind(pmin(edges$protein_a, edges$protein_b),
          pmax(edges$protein_a, edges$protein_b)),
    1,
    paste,
    collapse = "|"
  )
  edges <- edges[!duplicated(pair_key), , drop = FALSE]
  graph <- graph_from_data_frame(
    edges[, c("protein_a", "protein_b", "combined_score")],
    directed = FALSE
  )
  E(graph)$distance <- if (weighted) 1 / E(graph)$combined_score else 1
  graph
}

largest_component <- function(graph) {
  membership <- components(graph, mode = "weak")
  selected <- which(membership$membership == which.max(membership$csize))
  induced_subgraph(graph, selected)
}

reciprocal_efficiency_sum <- function(distance_matrix) {
  if (nrow(distance_matrix) < 2L) return(0)
  upper <- distance_matrix[upper.tri(distance_matrix)]
  sum(ifelse(is.finite(upper) & upper > 0, 1 / upper, 0))
}

nbvk_damage <- function(graph, original_distances, knockout,
                        weighted = FALSE) {
  knockout <- unique(as.character(knockout))
  retained <- setdiff(V(graph)$name, knockout)
  if (length(retained) < 2L) stop("At least two retained nodes are required")

  denominator <- reciprocal_efficiency_sum(
    original_distances[retained, retained, drop = FALSE]
  )
  reduced <- induced_subgraph(graph, vids = retained)
  reduced_distances <- distances(
    reduced,
    weights = if (weighted) E(reduced)$distance else NA,
    algorithm = if (weighted) "dijkstra" else "unweighted"
  )
  numerator <- reciprocal_efficiency_sum(reduced_distances)
  comp <- components(reduced, mode = "weak")
  retained_pairs <- choose(length(retained), 2)
  connected_pairs <- sum(choose(comp$csize, 2))

  data.frame(
    knockout = paste(sort(knockout), collapse = "+"),
    k = length(knockout),
    retained_nodes = length(retained),
    efficiency_sum_before = denominator,
    efficiency_sum_after = numerator,
    D_E = 1 - numerator / denominator,
    LCC_after = max(comp$csize),
    D_LCC = 1 - max(comp$csize) / length(retained),
    components_after = comp$no,
    disconnected_pair_fraction = 1 - connected_pairs / retained_pairs,
    stringsAsFactors = FALSE
  )
}

bliss_expected <- function(single_damages) {
  1 - prod(1 - single_damages)
}

synergy_metrics <- function(combination_damage, single_damages) {
  data.frame(
    bliss_expected = bliss_expected(single_damages),
    Syn_Bliss = combination_damage - bliss_expected(single_damages),
    Syn_add = combination_damage - sum(single_damages)
  )
}

empirical_summary <- function(observed, null_values) {
  null_values <- null_values[is.finite(null_values)]
  null_sd <- sd(null_values)
  data.frame(
    observed = observed,
    null_mean = mean(null_values),
    null_sd = null_sd,
    null_q025 = unname(quantile(null_values, 0.025)),
    null_median = unname(quantile(null_values, 0.5)),
    null_q975 = unname(quantile(null_values, 0.975)),
    z = if (null_sd > 0) (observed - mean(null_values)) / null_sd else NA_real_,
    empirical_percentile = mean(null_values <= observed),
    p_empirical_right = (1 + sum(null_values >= observed)) /
      (length(null_values) + 1)
  )
}

bh_within_group <- function(p_values, groups) {
  adjusted <- rep(NA_real_, length(p_values))
  for (group in unique(groups)) {
    selected <- which(groups == group)
    adjusted[selected] <- p.adjust(p_values[selected], method = "BH")
  }
  adjusted
}

make_topology_pools <- function(node_metrics, anchors, pool_size = 50L) {
  if (!is.logical(node_metrics$is_isolate)) {
    node_metrics$is_isolate <- tolower(node_metrics$is_isolate) == "true"
  }
  metrics <- node_metrics[
    !node_metrics$gene %in% anchors & !node_metrics$is_isolate,
    c("gene", "degree", "betweenness_centrality"),
    drop = FALSE
  ]
  all_metrics <- node_metrics[
    !node_metrics$is_isolate,
    c("gene", "degree", "betweenness_centrality"),
    drop = FALSE
  ]
  all_metrics$log_degree <- log1p(all_metrics$degree)
  all_metrics$betweenness_percentile <- rank(
    all_metrics$betweenness_centrality,
    ties.method = "average"
  ) / nrow(all_metrics)
  means <- colMeans(all_metrics[, c("log_degree", "betweenness_percentile")])
  sds <- apply(
    all_metrics[, c("log_degree", "betweenness_percentile")],
    2,
    sd
  )
  all_metrics$z_degree <- (all_metrics$log_degree - means[["log_degree"]]) /
    sds[["log_degree"]]
  all_metrics$z_betweenness <- (
    all_metrics$betweenness_percentile -
      means[["betweenness_percentile"]]
  ) / sds[["betweenness_percentile"]]

  result <- list()
  for (anchor in anchors) {
    anchor_row <- all_metrics[all_metrics$gene == anchor, , drop = FALSE]
    candidates <- all_metrics[
      !all_metrics$gene %in% anchors,
      ,
      drop = FALSE
    ]
    candidates$match_distance <- sqrt(
      (candidates$z_degree - anchor_row$z_degree)^2 +
        (candidates$z_betweenness - anchor_row$z_betweenness)^2
    )
    candidates <- candidates[
      order(candidates$match_distance, candidates$gene),
      ,
      drop = FALSE
    ]
    result[[anchor]] <- head(candidates, pool_size)
  }
  result
}
