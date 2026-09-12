suppressPackageStartupMessages({
  library(igraph)
  library(parallel)
})
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args)) normalizePath(args[[1]]) else normalizePath(file.path(script_dir, "../.."))
source(file.path(ROOT, "code/nbvk_v1/nbvk_core.R"))

PPI <- file.path(ROOT, "data/network_inputs")
PRIMARY <- file.path(ROOT, "results/nbvk_v1/primary")
RESULTS <- file.path(ROOT, "results/nbvk_v1/sensitivity")
NULL_DIR <- file.path(RESULTS, "null_distributions")
META <- file.path(ROOT, "environment/nbvk")
dir.create(NULL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(META, recursive = TRUE, showWarnings = FALSE)

B <- 10000L
WORKERS <- 8L
ANCHORS <- c("ARG1", "LCN2", "LTF")
SETS <- list(
  ARG1 = "ARG1",
  LCN2 = "LCN2",
  LTF = "LTF",
  ARG1_LCN2 = c("ARG1", "LCN2"),
  ARG1_LTF = c("ARG1", "LTF"),
  LCN2_LTF = c("LCN2", "LTF"),
  ARG1_LCN2_LTF = c("ARG1", "LCN2", "LTF")
)

calculate_observed <- function(graph, weighted, label) {
  distance_matrix <- distances(
    graph,
    weights = if (weighted) E(graph)$distance else NA,
    algorithm = if (weighted) "dijkstra" else "unweighted"
  )
  dimnames(distance_matrix) <- list(V(graph)$name, V(graph)$name)
  rows <- lapply(names(SETS), function(set_name) {
    result <- nbvk_damage(
      graph,
      distance_matrix,
      SETS[[set_name]],
      weighted = weighted
    )
    result$set_name <- set_name
    if (length(SETS[[set_name]]) > 1L) {
      singles <- vapply(
        SETS[[set_name]],
        function(gene) nbvk_damage(
          graph,
          distance_matrix,
          gene,
          weighted = weighted
        )$D_E,
        numeric(1)
      )
      result <- cbind(result, synergy_metrics(result$D_E, singles))
    } else {
      result$bliss_expected <- NA_real_
      result$Syn_Bliss <- NA_real_
      result$Syn_add <- NA_real_
    }
    result$network_variant <- label
    result
  })
  do.call(rbind, rows)
}

edge700 <- file.path(PPI, "string_edges_score700.csv")
edge400 <- file.path(PPI, "string_edges_score400.csv")
metrics_all <- read.csv(file.path(PPI, "full_ppi_node_metrics_score700.csv"))

g700_full_edges <- build_graph(edge700, 0.700, FALSE)
missing_nodes <- setdiff(metrics_all$gene, V(g700_full_edges)$name)
g700_all <- add_vertices(g700_full_edges, length(missing_nodes), name = missing_nodes)
stopifnot(vcount(g700_all) == 1291L, ecount(g700_all) == 13420L)

g400 <- largest_component(build_graph(edge400, 0.400, FALSE))
g700_weighted <- largest_component(build_graph(edge700, 0.700, TRUE))
g700_lcc <- largest_component(g700_full_edges)

ordered <- metrics_all[metrics_all$gene %in% V(g700_lcc)$name, , drop = FALSE]
ordered <- ordered[
  order(-ordered$degree, -ordered$betweenness_centrality, ordered$gene),
  ,
  drop = FALSE
]
top_n <- ceiling(0.01 * vcount(g700_lcc))
top1_genes <- head(ordered$gene, top_n)
stopifnot(!any(ANCHORS %in% top1_genes))
g700_without_top1 <- delete_vertices(g700_lcc, top1_genes)

variant_results <- do.call(
  rbind,
  list(
    calculate_observed(g700_all, FALSE, "full_1291_nodes"),
    calculate_observed(g400, FALSE, "score400_LCC"),
    calculate_observed(g700_weighted, TRUE, "score700_weighted_LCC"),
    calculate_observed(g700_without_top1, FALSE, "score700_LCC_without_top1pct_degree")
  )
)
write.csv(
  variant_results,
  file.path(RESULTS, "network_variant_observed.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(rank = seq_along(top1_genes), gene = top1_genes),
  file.path(RESULTS, "excluded_top1pct_degree_nodes.csv"),
  row.names = FALSE
)
cat("Network-variant observed analyses complete.\n")

# Primary graph resources reused for degree-only and alternate-seed nulls.
graph <- g700_lcc
original_distances <- distances(graph, weights = NA, algorithm = "unweighted")
dimnames(original_distances) <- list(V(graph)$name, V(graph)$name)
eligible <- setdiff(V(graph)$name, ANCHORS)
singleton <- read.csv(file.path(PRIMARY, "all_eligible_singleton_damage.csv"))
singleton_lookup <- setNames(singleton$D_E, singleton$gene)

damage_fast <- function(nodes) {
  if (length(nodes) == 1L) return(unname(singleton_lookup[nodes]))
  nbvk_damage(graph, original_distances, nodes)$D_E
}

evaluate_sample <- function(nodes) {
  combined <- damage_fast(nodes)
  singles <- unname(singleton_lookup[nodes])
  syn <- if (length(nodes) > 1L) {
    synergy_metrics(combined, singles)
  } else {
    data.frame(bliss_expected = NA_real_, Syn_Bliss = NA_real_, Syn_add = NA_real_)
  }
  c(
    sampled_set = paste(sort(nodes), collapse = "+"),
    D_E = combined,
    bliss_expected = syn$bliss_expected,
    Syn_Bliss = syn$Syn_Bliss,
    Syn_add = syn$Syn_add
  )
}

evaluate_samples <- function(samples, label, output_file) {
  cat(sprintf("Running %s: %d samples...\n", label, length(samples)))
  values <- mclapply(
    samples,
    evaluate_sample,
    mc.cores = WORKERS,
    mc.preschedule = TRUE
  )
  frame <- as.data.frame(do.call(rbind, values), stringsAsFactors = FALSE)
  for (column in c("D_E", "bliss_expected", "Syn_Bliss", "Syn_add")) {
    frame[[column]] <- as.numeric(frame[[column]])
  }
  frame$iteration <- seq_len(nrow(frame))
  frame$null_label <- label
  frame <- frame[, c(
    "null_label", "iteration", "sampled_set", "D_E",
    "bliss_expected", "Syn_Bliss", "Syn_add"
  )]
  write.csv(frame, output_file, row.names = FALSE)
  cat(sprintf("Completed %s.\n", label))
  frame
}

make_degree_pools <- function(node_metrics, anchors, pool_size = 50L) {
  x <- node_metrics[node_metrics$gene %in% V(graph)$name, , drop = FALSE]
  x$z_log_degree <- as.numeric(scale(log1p(x$degree)))
  pools <- list()
  for (anchor in anchors) {
    anchor_z <- x$z_log_degree[x$gene == anchor]
    candidates <- x[!x$gene %in% anchors, , drop = FALSE]
    candidates$match_distance <- abs(candidates$z_log_degree - anchor_z)
    candidates <- candidates[
      order(candidates$match_distance, candidates$gene),
      ,
      drop = FALSE
    ]
    pools[[anchor]] <- head(candidates, pool_size)
  }
  pools
}

sample_from_pools <- function(anchor_set, pools) {
  repeat {
    sampled <- vapply(
      anchor_set,
      function(anchor) sample(pools[[anchor]]$gene, 1L),
      character(1)
    )
    if (length(unique(sampled)) == length(sampled)) return(sampled)
  }
}

observed <- read.csv(file.path(PRIMARY, "observed_nbvk.csv"))
run_inference <- function(nulls, baseline, output_prefix) {
  damage_rows <- list()
  synergy_rows <- list()
  for (set_name in names(SETS)) {
    obs <- observed[observed$set_name == set_name, ]
    null <- nulls[[set_name]]
    damage <- empirical_summary(obs$D_E, null$D_E)
    damage$set_name <- set_name
    damage$baseline <- baseline
    damage_rows[[length(damage_rows) + 1L]] <- damage
    if (length(SETS[[set_name]]) > 1L) {
      synergy <- empirical_summary(obs$Syn_Bliss, null$Syn_Bliss)
      synergy$set_name <- set_name
      synergy$baseline <- baseline
      synergy_rows[[length(synergy_rows) + 1L]] <- synergy
    }
  }
  damage <- do.call(rbind, damage_rows)
  damage$FDR_BH <- p.adjust(damage$p_empirical_right, method = "BH")
  synergy <- do.call(rbind, synergy_rows)
  synergy$FDR_BH <- p.adjust(synergy$p_empirical_right, method = "BH")
  write.csv(
    damage,
    file.path(RESULTS, paste0(output_prefix, "_damage_inference.csv")),
    row.names = FALSE
  )
  write.csv(
    synergy,
    file.path(RESULTS, paste0(output_prefix, "_synergy_inference.csv")),
    row.names = FALSE
  )
}

# Degree-only matching sensitivity.
degree_pools <- make_degree_pools(metrics_all, ANCHORS, 50L)
write.csv(
  do.call(rbind, lapply(names(degree_pools), function(anchor) {
    x <- degree_pools[[anchor]]
    data.frame(
      anchor = anchor,
      pool_rank = seq_len(nrow(x)),
      gene = x$gene,
      degree = x$degree,
      match_distance = x$match_distance
    )
  })),
  file.path(RESULTS, "degree_only_matching_pools.csv"),
  row.names = FALSE
)
set.seed(20260726L)
degree_nulls <- list()
for (set_name in names(SETS)) {
  samples <- replicate(
    B,
    sample_from_pools(SETS[[set_name]], degree_pools),
    simplify = FALSE
  )
  degree_nulls[[set_name]] <- evaluate_samples(
    samples,
    sprintf("degree_only_%s", set_name),
    file.path(NULL_DIR, sprintf("degree_only_%s.csv", set_name))
  )
}
run_inference(degree_nulls, "degree_only_matched", "degree_only")

# Alternate seed Monte Carlo stability for both registered baselines.
topology_pools_table <- read.csv(
  file.path(PRIMARY, "topology_matching_pools.csv")
)
topology_pools <- split(topology_pools_table, topology_pools_table$anchor)
set.seed(20260727L)
uniform_by_k <- list()
for (k in 1:3) {
  samples <- replicate(B, sample(eligible, k, replace = FALSE), simplify = FALSE)
  uniform_by_k[[as.character(k)]] <- evaluate_samples(
    samples,
    sprintf("run02_uniform_k%d", k),
    file.path(NULL_DIR, sprintf("run02_uniform_k%d.csv", k))
  )
}
topology_alt <- list()
for (set_name in names(SETS)) {
  samples <- replicate(
    B,
    sample_from_pools(SETS[[set_name]], topology_pools),
    simplify = FALSE
  )
  topology_alt[[set_name]] <- evaluate_samples(
    samples,
    sprintf("run02_topology_%s", set_name),
    file.path(NULL_DIR, sprintf("run02_topology_%s.csv", set_name))
  )
}
run_inference(topology_alt, "topology_matched_run02", "run02_topology")

uniform_alt_by_set <- lapply(
  names(SETS),
  function(set_name) uniform_by_k[[as.character(length(SETS[[set_name]]))]]
)
names(uniform_alt_by_set) <- names(SETS)
run_inference(
  uniform_alt_by_set,
  "uniform_random_run02",
  "run02_uniform"
)

writeLines(
  capture.output(sessionInfo()),
  file.path(META, "R_sessionInfo_sensitivity.txt")
)
writeLines(
  c(
    "permutations_per_set=10000",
    "workers=8",
    "primary_seed=20260726",
    "alternate_seed=20260727",
    "analytical_deviations=none"
  ),
  file.path(META, "sensitivity_run_parameters.txt")
)
cat("SENSITIVITY_NBVK_COMPLETE\n")
