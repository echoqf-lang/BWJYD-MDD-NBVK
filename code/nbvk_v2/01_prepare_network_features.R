suppressPackageStartupMessages(library(igraph))

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(file.path(script_dir, "../.."))
edge_file <- file.path(root, "data/network_inputs/string_edges_score700.csv")
out_dir <- file.path(root, "results/nbvk_v2/derived")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

edges <- read.csv(edge_file, check.names = FALSE)
edges <- edges[
  edges$combined_score >= 0.700 & edges$protein_a != edges$protein_b,
  ,
  drop = FALSE
]
edges$key <- paste(
  pmin(edges$protein_a, edges$protein_b),
  pmax(edges$protein_a, edges$protein_b),
  sep = "|"
)
edges <- edges[order(edges$key, -edges$combined_score), , drop = FALSE]
edges <- edges[!duplicated(edges$key), , drop = FALSE]

g <- graph_from_data_frame(
  edges[, c("protein_a", "protein_b", "combined_score")],
  directed = FALSE
)
anchors <- c("ARG1", "LCN2", "LTF")
if (!all(anchors %in% V(g)$name)) stop("Anchor missing from graph")
cc <- components(g)
anchor_components <- unique(cc$membership[match(anchors, V(g)$name)])
if (length(anchor_components) != 1L) stop("Anchors are not in one component")
g <- induced_subgraph(g, which(cc$membership == anchor_components[[1]]))

set.seed(20260726)
comm <- cluster_leiden(
  g,
  objective_function = "modularity",
  resolution = 1.0,
  n_iterations = -1
)
membership <- membership(comm)
names(membership) <- V(g)$name

deg <- degree(g)
btw <- betweenness(g, directed = FALSE, normalized = TRUE)
btw_pct <- rank(btw, ties.method = "average") / length(btw)
part <- numeric(vcount(g))
names(part) <- V(g)$name
for (v in V(g)$name) {
  if (deg[[v]] == 0) {
    part[[v]] <- 0
  } else {
    nbr <- neighbors(g, v)$name
    counts <- table(membership[nbr])
    part[[v]] <- 1 - sum((counts / deg[[v]])^2)
  }
}

node_features <- data.frame(
  gene = V(g)$name,
  degree = as.numeric(deg[V(g)$name]),
  log1p_degree = log1p(as.numeric(deg[V(g)$name])),
  betweenness = as.numeric(btw[V(g)$name]),
  betweenness_percentile = as.numeric(btw_pct[V(g)$name]),
  participation = as.numeric(part[V(g)$name]),
  community = as.integer(membership[V(g)$name]),
  stringsAsFactors = FALSE
)
write.csv(
  node_features,
  file.path(out_dir, "node_features_score700.csv"),
  row.names = FALSE
)
write.csv(
  edges[edges$protein_a %in% V(g)$name & edges$protein_b %in% V(g)$name, ],
  file.path(out_dir, "lcc_edges_score700.csv"),
  row.names = FALSE
)

cat(
  sprintf(
    "FEATURES_READY nodes=%d edges=%d communities=%d\n",
    vcount(g), ecount(g), length(unique(membership))
  )
)
