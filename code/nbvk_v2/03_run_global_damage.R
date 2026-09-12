suppressPackageStartupMessages({
  library(igraph)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(file.path(script_dir, "../.."))
source(file.path(root, "code/nbvk_v1/nbvk_core.R"))

derived <- file.path(root, "results/nbvk_v2/derived")
results_dir <- file.path(root, "results/nbvk_v2")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

edges <- read.csv(file.path(derived, "lcc_edges_score700.csv"), check.names = FALSE)
g <- graph_from_data_frame(
  edges[, c("protein_a", "protein_b", "combined_score")],
  directed = FALSE
)
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
E(g)$distance <- 1
original_distances <- distances(g, weights = NA, algorithm = "unweighted")

anchors <- c("ARG1", "LCN2", "LTF")
observed <- nbvk_damage(g, original_distances, anchors, weighted = FALSE)

triplets <- read.csv(file.path(derived, "legal_relax2_triplets.csv"))
strict <- read.csv(file.path(derived, "legal_strict_triplets.csv"))
relax1 <- read.csv(file.path(derived, "legal_relax1_triplets.csv"))
key <- function(x) paste(sort(as.character(x)), collapse = "+")
strict_keys <- apply(strict, 1, key)
relax1_keys <- apply(relax1, 1, key)

available_cores <- detectCores(logical = FALSE)
if (is.na(available_cores)) available_cores <- 2L
n_cores <- min(4L, max(1L, available_cores))
damage_one <- function(i) {
  genes <- as.character(unlist(triplets[i, ], use.names = FALSE))
  out <- nbvk_damage(g, original_distances, genes, weighted = FALSE)
  out$in_strict <- out$knockout %in% strict_keys
  out$in_relax1 <- out$knockout %in% relax1_keys
  out
}
null_list <- mclapply(seq_len(nrow(triplets)), damage_one, mc.cores = n_cores)
null <- do.call(rbind, null_list)

summarize_level <- function(selected, level) {
  values <- null$D_E[selected]
  b <- sum(values >= observed$D_E)
  data.frame(
    level = level,
    n = length(values),
    observed_D_E = observed$D_E,
    null_mean = mean(values),
    null_median = median(values),
    null_sd = sd(values),
    null_q025 = unname(quantile(values, 0.025)),
    null_q975 = unname(quantile(values, 0.975)),
    percentile = mean(values <= observed$D_E),
    tail_b = b,
    p_tail_raw = b / length(values),
    p_tail_plus1 = (b + 1) / (length(values) + 1),
    confirmatory_allowed = level == "strict" && length(values) >= 1000,
    stringsAsFactors = FALSE
  )
}
summary <- rbind(
  summarize_level(null$in_strict, "strict"),
  summarize_level(null$in_relax1, "relax1_drop_community"),
  summarize_level(rep(TRUE, nrow(null)), "relax2_induced_and_distance_only")
)

write.csv(observed, file.path(results_dir, "observed_global_damage.csv"), row.names = FALSE)
write.csv(null, file.path(results_dir, "global_damage_null_relax2.csv"), row.names = FALSE)
write.csv(summary, file.path(results_dir, "global_damage_inference.csv"), row.names = FALSE)
cat("GLOBAL_DAMAGE_READY\n")
print(summary)
