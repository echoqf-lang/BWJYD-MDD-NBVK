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
RESULTS <- file.path(ROOT, "results/nbvk_v1/sensitivity")
NULL_DIR <- file.path(RESULTS, "variant_topology_nulls")
dir.create(NULL_DIR, recursive = TRUE, showWarnings = FALSE)

B <- 10000L
WORKERS <- 8L
SEED <- 20260726L
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
observed_variants <- read.csv(
  file.path(RESULTS, "network_variant_observed.csv")
)

run_variant <- function(label, graph, weighted) {
  cat(sprintf("Preparing %s (%d nodes, %d edges)...\n",
              label, vcount(graph), ecount(graph)))
  original <- distances(
    graph,
    weights = if (weighted) E(graph)$distance else NA,
    algorithm = if (weighted) "dijkstra" else "unweighted"
  )
  dimnames(original) <- list(V(graph)$name, V(graph)$name)

  node_metrics <- data.frame(
    gene = V(graph)$name,
    degree = degree(graph),
    betweenness_centrality = betweenness(
      graph,
      directed = FALSE,
      weights = if (weighted) E(graph)$distance else NA,
      normalized = TRUE
    ),
    is_isolate = degree(graph) == 0
  )
  pools <- make_topology_pools(node_metrics, ANCHORS, 50L)
  pool_table <- do.call(rbind, lapply(names(pools), function(anchor) {
    x <- pools[[anchor]]
    data.frame(
      network_variant = label,
      anchor = anchor,
      pool_rank = seq_len(nrow(x)),
      gene = x$gene,
      degree = x$degree,
      betweenness_centrality = x$betweenness_centrality,
      match_distance = x$match_distance
    )
  }))
  write.csv(
    pool_table,
    file.path(RESULTS, sprintf("%s_topology_matching_pools.csv", label)),
    row.names = FALSE
  )

  eligible <- setdiff(V(graph)$name, ANCHORS)
  damage <- function(nodes) {
    nbvk_damage(graph, original, nodes, weighted = weighted)$D_E
  }
  cat(sprintf("Precomputing %s singleton damage...\n", label))
  singleton_values <- unlist(
    mclapply(eligible, damage, mc.cores = WORKERS, mc.preschedule = TRUE),
    use.names = FALSE
  )
  singleton_lookup <- setNames(singleton_values, eligible)
  write.csv(
    data.frame(gene = eligible, D_E = singleton_values),
    file.path(RESULTS, sprintf("%s_singleton_damage.csv", label)),
    row.names = FALSE
  )

  evaluate <- function(nodes) {
    combined <- if (length(nodes) == 1L) {
      unname(singleton_lookup[nodes])
    } else {
      damage(nodes)
    }
    singles <- unname(singleton_lookup[nodes])
    syn <- if (length(nodes) > 1L) {
      synergy_metrics(combined, singles)
    } else {
      data.frame(bliss_expected = NA_real_, Syn_Bliss = NA_real_, Syn_add = NA_real_)
    }
    c(
      sampled_set = paste(sort(nodes), collapse = "+"),
      D_E = combined,
      Syn_Bliss = syn$Syn_Bliss,
      Syn_add = syn$Syn_add
    )
  }

  sample_set <- function(anchor_set) {
    repeat {
      selected <- vapply(
        anchor_set,
        function(anchor) sample(pools[[anchor]]$gene, 1L),
        character(1)
      )
      if (length(unique(selected)) == length(selected)) return(selected)
    }
  }

  set.seed(SEED)
  nulls <- list()
  for (set_name in names(SETS)) {
    samples <- replicate(B, sample_set(SETS[[set_name]]), simplify = FALSE)
    cat(sprintf("Running %s %s...\n", label, set_name))
    sample_keys <- vapply(
      samples,
      function(nodes) paste(sort(nodes), collapse = "+"),
      character(1)
    )
    unique_indices <- !duplicated(sample_keys)
    unique_samples <- samples[unique_indices]
    unique_keys <- sample_keys[unique_indices]
    cat(sprintf("  %d unique sets among %d Monte Carlo draws.\n",
                length(unique_samples), length(samples)))
    unique_values <- mclapply(
      unique_samples,
      evaluate,
      mc.cores = WORKERS,
      mc.preschedule = TRUE
    )
    unique_frame <- as.data.frame(
      do.call(rbind, unique_values),
      stringsAsFactors = FALSE
    )
    frame <- unique_frame[match(sample_keys, unique_keys), , drop = FALSE]
    for (column in c("D_E", "Syn_Bliss", "Syn_add")) {
      frame[[column]] <- as.numeric(frame[[column]])
    }
    frame$iteration <- seq_len(nrow(frame))
    write.csv(
      frame,
      file.path(NULL_DIR, sprintf("%s_topology_%s.csv", label, set_name)),
      row.names = FALSE
    )
    nulls[[set_name]] <- frame
  }

  damage_rows <- list()
  synergy_rows <- list()
  for (set_name in names(SETS)) {
    obs <- observed_variants[
      observed_variants$network_variant == label &
        observed_variants$set_name == set_name,
    ]
    ds <- empirical_summary(obs$D_E, nulls[[set_name]]$D_E)
    ds$network_variant <- label
    ds$set_name <- set_name
    damage_rows[[length(damage_rows) + 1L]] <- ds
    if (length(SETS[[set_name]]) > 1L) {
      ss <- empirical_summary(obs$Syn_Bliss, nulls[[set_name]]$Syn_Bliss)
      ss$network_variant <- label
      ss$set_name <- set_name
      synergy_rows[[length(synergy_rows) + 1L]] <- ss
    }
  }
  damage_result <- do.call(rbind, damage_rows)
  damage_result$FDR_BH <- p.adjust(
    damage_result$p_empirical_right,
    method = "BH"
  )
  synergy_result <- do.call(rbind, synergy_rows)
  synergy_result$FDR_BH <- p.adjust(
    synergy_result$p_empirical_right,
    method = "BH"
  )
  write.csv(
    damage_result,
    file.path(RESULTS, sprintf("%s_topology_damage_inference.csv", label)),
    row.names = FALSE
  )
  write.csv(
    synergy_result,
    file.path(RESULTS, sprintf("%s_topology_synergy_inference.csv", label)),
    row.names = FALSE
  )
  cat(sprintf("Completed %s.\n", label))
}

score400 <- largest_component(
  build_graph(file.path(PPI, "string_edges_score400.csv"), 0.400, FALSE)
)
score700_weighted <- largest_component(
  build_graph(file.path(PPI, "string_edges_score700.csv"), 0.700, TRUE)
)
if (!file.exists(file.path(
  RESULTS,
  "score400_LCC_topology_damage_inference.csv"
))) {
  run_variant("score400_LCC", score400, FALSE)
} else {
  cat("score400_LCC inference already complete; preserving existing outputs.\n")
}
run_variant("score700_weighted_LCC", score700_weighted, TRUE)
cat("VARIANT_TOPOLOGY_INFERENCE_COMPLETE\n")
