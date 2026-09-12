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
RESULTS <- file.path(ROOT, "results/nbvk_v1/primary")
NULL_DIR <- file.path(RESULTS, "null_distributions")
META <- file.path(ROOT, "environment/nbvk")
LOGS <- file.path(ROOT, "docs/nbvk")
dir.create(NULL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(META, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGS, recursive = TRUE, showWarnings = FALSE)

SEED <- 20260726L
B <- 10000L
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
# The sandbox cannot query macOS CPU topology (`detectCores()` returns NA).
# Use a fixed worker count so the execution environment is reproducible.
workers <- 8L
cat(sprintf("NBVK primary run: seed=%d, B=%d, workers=%d\n", SEED, B, workers))

edge_file <- file.path(PPI, "string_edges_score700.csv")
graph_full <- build_graph(edge_file, threshold = 0.700, weighted = FALSE)
graph <- largest_component(graph_full)
stopifnot(vcount(graph) == 1226L, ecount(graph) == 13416L)
stopifnot(all(ANCHORS %in% V(graph)$name))

original_distances <- distances(graph, weights = NA, algorithm = "unweighted")
dimnames(original_distances) <- list(V(graph)$name, V(graph)$name)

observed <- do.call(
  rbind,
  lapply(names(SETS), function(set_name) {
    row <- nbvk_damage(graph, original_distances, SETS[[set_name]])
    row$set_name <- set_name
    if (length(SETS[[set_name]]) > 1L) {
      singles <- vapply(
        SETS[[set_name]],
        function(gene) nbvk_damage(graph, original_distances, gene)$D_E,
        numeric(1)
      )
      row <- cbind(row, synergy_metrics(row$D_E, singles))
    } else {
      row$bliss_expected <- NA_real_
      row$Syn_Bliss <- NA_real_
      row$Syn_add <- NA_real_
    }
    row
  })
)
observed <- observed[, c("set_name", setdiff(names(observed), "set_name"))]
write.csv(observed, file.path(RESULTS, "observed_nbvk.csv"), row.names = FALSE)
cat("Observed seven prespecified perturbations saved.\n")

node_metrics <- read.csv(file.path(PPI, "full_ppi_node_metrics_score700.csv"))
node_metrics <- node_metrics[node_metrics$gene %in% V(graph)$name, , drop = FALSE]
pools <- make_topology_pools(node_metrics, ANCHORS, pool_size = 50L)
pool_table <- do.call(
  rbind,
  lapply(names(pools), function(anchor) {
    x <- pools[[anchor]]
    x$anchor <- anchor
    x$pool_rank <- seq_len(nrow(x))
    x[, c(
      "anchor", "pool_rank", "gene", "degree",
      "betweenness_centrality", "match_distance"
    )]
  })
)
write.csv(pool_table, file.path(RESULTS, "topology_matching_pools.csv"), row.names = FALSE)

eligible <- setdiff(V(graph)$name, ANCHORS)

damage_one <- function(nodes) {
  nbvk_damage(graph, original_distances, nodes)$D_E
}

cat(sprintf("Precomputing %d singleton knockouts...\n", length(eligible)))
singleton_values <- unlist(
  mclapply(eligible, damage_one, mc.cores = workers, mc.preschedule = TRUE),
  use.names = FALSE
)
singleton_damage <- data.frame(gene = eligible, D_E = singleton_values)
write.csv(
  singleton_damage,
  file.path(RESULTS, "all_eligible_singleton_damage.csv"),
  row.names = FALSE
)
singleton_lookup <- setNames(singleton_damage$D_E, singleton_damage$gene)

evaluate_sample <- function(nodes) {
  combined <- damage_one(nodes)
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
    mc.cores = workers,
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

set.seed(SEED)
uniform_nulls <- list()
for (k in 1:3) {
  samples <- replicate(B, sample(eligible, k, replace = FALSE), simplify = FALSE)
  uniform_nulls[[as.character(k)]] <- evaluate_samples(
    samples,
    sprintf("uniform_k%d", k),
    file.path(NULL_DIR, sprintf("uniform_k%d.csv", k))
  )
}

sample_matched_set <- function(anchor_set) {
  repeat {
    sampled <- vapply(
      anchor_set,
      function(anchor) sample(pools[[anchor]]$gene, 1L),
      character(1)
    )
    if (length(unique(sampled)) == length(sampled)) return(sampled)
  }
}

topology_nulls <- list()
for (set_name in names(SETS)) {
  anchor_set <- SETS[[set_name]]
  samples <- replicate(
    B,
    sample_matched_set(anchor_set),
    simplify = FALSE
  )
  topology_nulls[[set_name]] <- evaluate_samples(
    samples,
    sprintf("topology_%s", set_name),
    file.path(NULL_DIR, sprintf("topology_%s.csv", set_name))
  )
}

summaries <- list()
for (set_name in names(SETS)) {
  obs <- observed[observed$set_name == set_name, ]
  k <- length(SETS[[set_name]])
  topo <- topology_nulls[[set_name]]
  uni <- uniform_nulls[[as.character(k)]]
  topo_summary <- empirical_summary(obs$D_E, topo$D_E)
  topo_summary$set_name <- set_name
  topo_summary$baseline <- "topology_matched"
  uni_summary <- empirical_summary(obs$D_E, uni$D_E)
  uni_summary$set_name <- set_name
  uni_summary$baseline <- "uniform_random"
  summaries[[length(summaries) + 1L]] <- topo_summary
  summaries[[length(summaries) + 1L]] <- uni_summary
}
damage_inference <- do.call(rbind, summaries)
damage_inference$FDR_BH <- bh_within_group(
  damage_inference$p_empirical_right,
  damage_inference$baseline
)
damage_inference <- damage_inference[, c(
  "set_name", "baseline", "observed", "null_mean", "null_sd",
  "null_q025", "null_median", "null_q975", "z",
  "empirical_percentile", "p_empirical_right", "FDR_BH"
)]
write.csv(
  damage_inference,
  file.path(RESULTS, "damage_inference.csv"),
  row.names = FALSE
)

synergy_summaries <- list()
combo_names <- names(SETS)[lengths(SETS) > 1L]
for (set_name in combo_names) {
  obs <- observed[observed$set_name == set_name, ]
  topo <- topology_nulls[[set_name]]
  k <- length(SETS[[set_name]])
  uni <- uniform_nulls[[as.character(k)]]
  topo_summary <- empirical_summary(obs$Syn_Bliss, topo$Syn_Bliss)
  topo_summary$set_name <- set_name
  topo_summary$baseline <- "topology_matched"
  uni_summary <- empirical_summary(obs$Syn_Bliss, uni$Syn_Bliss)
  uni_summary$set_name <- set_name
  uni_summary$baseline <- "uniform_random"
  synergy_summaries[[length(synergy_summaries) + 1L]] <- topo_summary
  synergy_summaries[[length(synergy_summaries) + 1L]] <- uni_summary
}
synergy_inference <- do.call(rbind, synergy_summaries)
synergy_inference$FDR_BH <- bh_within_group(
  synergy_inference$p_empirical_right,
  synergy_inference$baseline
)
synergy_inference <- synergy_inference[, c(
  "set_name", "baseline", "observed", "null_mean", "null_sd",
  "null_q025", "null_median", "null_q975", "z",
  "empirical_percentile", "p_empirical_right", "FDR_BH"
)]
write.csv(
  synergy_inference,
  file.path(RESULTS, "synergy_inference.csv"),
  row.names = FALSE
)

ordered_hubs <- node_metrics[
  order(
    -node_metrics$degree,
    -node_metrics$betweenness_centrality,
    node_metrics$gene
  ),
  ,
  drop = FALSE
]
ordered_hubs <- ordered_hubs[!ordered_hubs$gene %in% ANCHORS, , drop = FALSE]
hub_rows <- do.call(
  rbind,
  lapply(1:3, function(k) {
    selected <- head(ordered_hubs$gene, k)
    result <- nbvk_damage(graph, original_distances, selected)
    result$hub_genes <- paste(selected, collapse = "+")
    result
  })
)
write.csv(hub_rows, file.path(RESULTS, "hub_reference.csv"), row.names = FALSE)

triple_damage <- damage_inference[
  damage_inference$set_name == "ARG1_LCN2_LTF" &
    damage_inference$baseline == "topology_matched",
]
triple_synergy <- synergy_inference[
  synergy_inference$set_name == "ARG1_LCN2_LTF" &
    synergy_inference$baseline == "topology_matched",
]
decision <- data.frame(
  hypothesis = c("triple_damage", "triple_synergy"),
  supported = c(
    observed$D_E[observed$set_name == "ARG1_LCN2_LTF"] > 0 &&
      triple_damage$FDR_BH < 0.05 &&
      triple_damage$observed > triple_damage$null_q975,
    observed$Syn_Bliss[observed$set_name == "ARG1_LCN2_LTF"] > 0 &&
      triple_synergy$FDR_BH < 0.05
  ),
  rule = c(
    "D_E>0; topology-matched BH FDR<0.05; observed>null_q975",
    "Syn_Bliss>0; topology-matched synergy-family BH FDR<0.05"
  )
)
write.csv(decision, file.path(RESULTS, "preregistered_decisions.csv"), row.names = FALSE)

writeLines(
  capture.output(sessionInfo()),
  file.path(META, "R_sessionInfo_primary.txt")
)
writeLines(
  c(
    sprintf("seed=%d", SEED),
    sprintf("permutations_per_set=%d", B),
    sprintf("workers=%d", workers),
    "deviations=none"
  ),
  file.path(META, "primary_run_parameters.txt")
)
cat("PRIMARY_NBVK_COMPLETE\n")
