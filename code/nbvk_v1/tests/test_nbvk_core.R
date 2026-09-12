script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
source(file.path(script_dir, "../nbvk_core.R"))

# Path A-B-C-D: deleting B disconnects A from C-D.
graph <- graph_from_edgelist(
  matrix(c("A", "B", "B", "C", "C", "D"), ncol = 2, byrow = TRUE),
  directed = FALSE
)
E(graph)$distance <- 1
original <- distances(graph, weights = NA, algorithm = "unweighted")
damage_b <- nbvk_damage(graph, original, "B")

# Retained pairs: A-C, A-D, C-D. Before reciprocal sum = 1/2 + 1/3 + 1.
# After deletion only C-D remains = 1.
expected <- 1 - 1 / (1 / 2 + 1 / 3 + 1)
stopifnot(abs(damage_b$D_E - expected) < 1e-12)
stopifnot(damage_b$components_after == 2)
stopifnot(abs(damage_b$disconnected_pair_fraction - 2 / 3) < 1e-12)

# Triangle: deleting one node does not damage communication among retained nodes.
triangle <- make_ring(3)
V(triangle)$name <- c("A", "B", "C")
E(triangle)$distance <- 1
triangle_dist <- distances(triangle, weights = NA, algorithm = "unweighted")
stopifnot(abs(nbvk_damage(triangle, triangle_dist, "A")$D_E) < 1e-12)

# Bliss and additive synergy formulas.
syn <- synergy_metrics(0.30, c(0.10, 0.10))
stopifnot(abs(syn$bliss_expected - 0.19) < 1e-12)
stopifnot(abs(syn$Syn_Bliss - 0.11) < 1e-12)
stopifnot(abs(syn$Syn_add - 0.10) < 1e-12)

# Empirical right-tail p uses the registered +1 correction.
summary <- empirical_summary(3, c(1, 2, 3, 4))
stopifnot(abs(summary$p_empirical_right - 3 / 5) < 1e-12)

# BH correction must occur within, not across, baseline families.
p <- c(0.01, 0.04, 0.02, 0.20)
groups <- c("A", "A", "B", "B")
expected_bh <- c(0.02, 0.04, 0.04, 0.20)
stopifnot(all.equal(bh_within_group(p, groups), expected_bh))

cat("ALL_NBVK_CORE_TESTS_PASS\n")
