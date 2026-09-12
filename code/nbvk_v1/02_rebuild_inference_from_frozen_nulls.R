script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args)) normalizePath(args[[1]]) else normalizePath(file.path(script_dir, "../.."))
source(file.path(ROOT, "code/nbvk_v1/nbvk_core.R"))

RESULTS <- file.path(ROOT, "results/nbvk_v1/primary")
NULL_DIR <- file.path(RESULTS, "null_distributions")

SETS <- list(
  ARG1 = "ARG1",
  LCN2 = "LCN2",
  LTF = "LTF",
  ARG1_LCN2 = c("ARG1", "LCN2"),
  ARG1_LTF = c("ARG1", "LTF"),
  LCN2_LTF = c("LCN2", "LTF"),
  ARG1_LCN2_LTF = c("ARG1", "LCN2", "LTF")
)

observed <- read.csv(file.path(RESULTS, "observed_nbvk.csv"))

damage_rows <- list()
for (set_name in names(SETS)) {
  obs <- observed[observed$set_name == set_name, ]
  k <- length(SETS[[set_name]])
  null_specs <- list(
    topology_matched = file.path(
      NULL_DIR,
      sprintf("topology_%s.csv", set_name)
    ),
    uniform_random = file.path(NULL_DIR, sprintf("uniform_k%d.csv", k))
  )
  for (baseline in names(null_specs)) {
    null <- read.csv(null_specs[[baseline]])
    stopifnot(nrow(null) == 10000L, all(is.finite(null$D_E)))
    summary <- empirical_summary(obs$D_E, null$D_E)
    summary$set_name <- set_name
    summary$baseline <- baseline
    damage_rows[[length(damage_rows) + 1L]] <- summary
  }
}
damage <- do.call(rbind, damage_rows)
damage$FDR_BH <- bh_within_group(damage$p_empirical_right, damage$baseline)
damage <- damage[, c(
  "set_name", "baseline", "observed", "null_mean", "null_sd",
  "null_q025", "null_median", "null_q975", "z",
  "empirical_percentile", "p_empirical_right", "FDR_BH"
)]
write.csv(damage, file.path(RESULTS, "damage_inference.csv"), row.names = FALSE)

synergy_rows <- list()
for (set_name in names(SETS)[lengths(SETS) > 1L]) {
  obs <- observed[observed$set_name == set_name, ]
  k <- length(SETS[[set_name]])
  null_specs <- list(
    topology_matched = file.path(
      NULL_DIR,
      sprintf("topology_%s.csv", set_name)
    ),
    uniform_random = file.path(NULL_DIR, sprintf("uniform_k%d.csv", k))
  )
  for (baseline in names(null_specs)) {
    null <- read.csv(null_specs[[baseline]])
    stopifnot(nrow(null) == 10000L, all(is.finite(null$Syn_Bliss)))
    summary <- empirical_summary(obs$Syn_Bliss, null$Syn_Bliss)
    summary$set_name <- set_name
    summary$baseline <- baseline
    synergy_rows[[length(synergy_rows) + 1L]] <- summary
  }
}
synergy <- do.call(rbind, synergy_rows)
synergy$FDR_BH <- bh_within_group(synergy$p_empirical_right, synergy$baseline)
synergy <- synergy[, c(
  "set_name", "baseline", "observed", "null_mean", "null_sd",
  "null_q025", "null_median", "null_q975", "z",
  "empirical_percentile", "p_empirical_right", "FDR_BH"
)]
write.csv(
  synergy,
  file.path(RESULTS, "synergy_inference.csv"),
  row.names = FALSE
)

triple_damage <- damage[
  damage$set_name == "ARG1_LCN2_LTF" &
    damage$baseline == "topology_matched",
]
triple_synergy <- synergy[
  synergy$set_name == "ARG1_LCN2_LTF" &
    synergy$baseline == "topology_matched",
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
write.csv(
  decision,
  file.path(RESULTS, "preregistered_decisions.csv"),
  row.names = FALSE
)
cat("INFERENCE_REBUILD_COMPLETE\n")
