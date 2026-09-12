#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
check_only <- "--check" %in% args
simulate <- "--simulate" %in% args

root <- normalizePath(".", mustWork = TRUE)
project_lib <- file.path(root, ".r-lib")
if (dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))
manifest_path <- file.path(root, "data", "derived", "geo_sample_manifest.tsv")
raw_root <- file.path(root, "data", "raw_geo", "expression")
clean_root <- Sys.getenv("PREPROCESS_OUTPUT_ROOT", unset = file.path(root, "data", "clean"))

stopf <- function(...) stop(sprintf(...), call. = FALSE)

sha256 <- function(path) {
  out <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  if (!length(out)) stopf("SHA256 failed: %s", path)
  strsplit(out[[1]], "[[:space:]]+")[[1]][1]
}

read_soft_descriptions <- function(path) {
  x <- readLines(path, warn = FALSE)
  gsm <- sub("^!Sample_geo_accession = ", "", x[grepl("^!Sample_geo_accession = ", x)])
  desc <- sub("^!Sample_description = ", "", x[grepl("^!Sample_description = ", x)])
  if (length(gsm) != length(desc)) {
    stopf("SOFT GSM/description count mismatch in %s: %d vs %d", path, length(gsm), length(desc))
  }
  setNames(gsm, desc)
}

read_soft_titles <- function(path) {
  x <- readLines(path, warn = FALSE)
  gsm <- sub("^!Sample_geo_accession = ", "", x[grepl("^!Sample_geo_accession = ", x)])
  title <- sub("^!Sample_title = ", "", x[grepl("^!Sample_title = ", x)])
  if (length(gsm) != length(title)) {
    stopf("SOFT GSM/title count mismatch in %s: %d vs %d", path, length(gsm), length(title))
  }
  setNames(gsm, title)
}

write_tsv_gz <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

validate_output <- function(dataset, expected_gsm, raw_path, expected_sha = NULL) {
  out_dir <- file.path(clean_root, dataset)
  expr_path <- file.path(out_dir, "expression_gene.tsv.gz")
  meta_path <- file.path(out_dir, "sample_metadata.tsv")
  qc_path <- file.path(out_dir, "qc_metrics.tsv")
  needed <- c(expr_path, meta_path, qc_path)
  if (dataset == "GSE251778") {
    needed <- c(needed, file.path(out_dir, "raw_counts_gene.tsv.gz"),
                file.path(out_dir, "normalization_factors.tsv"))
  }
  missing <- needed[!file.exists(needed)]
  if (length(missing)) stopf("%s missing outputs: %s", dataset, paste(missing, collapse = ", "))
  expr <- read.delim(expr_path, check.names = FALSE)
  meta <- read.delim(meta_path, check.names = FALSE)
  qc <- read.delim(qc_path, check.names = FALSE)
  if (!identical(names(expr)[1], "Gene Symbol")) stopf("%s first expression column is not Gene Symbol", dataset)
  if (anyDuplicated(expr[["Gene Symbol"]])) stopf("%s has duplicated Gene Symbol", dataset)
  if (any(!nzchar(expr[["Gene Symbol"]]))) stopf("%s has empty Gene Symbol", dataset)
  if (!identical(names(expr)[-1], expected_gsm)) stopf("%s expression column order differs from manifest", dataset)
  if (!identical(meta$GSM, expected_gsm)) stopf("%s metadata row order differs from manifest", dataset)
  if (!identical(qc$GSM, expected_gsm)) stopf("%s QC row order differs from manifest", dataset)
  m <- as.matrix(expr[-1])
  storage.mode(m) <- "double"
  if (any(!is.finite(m))) stopf("%s expression contains non-finite values", dataset)
  actual_sha <- sha256(raw_path)
  if (!is.null(expected_sha) && !identical(actual_sha, expected_sha)) {
    stopf("%s raw SHA mismatch: %s != %s", dataset, actual_sha, expected_sha)
  }
  message(sprintf("CHECK_OK %s N=%d genes=%d raw_sha256=%s", dataset, length(expected_gsm), nrow(expr), actual_sha))
}

compare_clean_dirs <- function(reference, candidate, datasets, tolerance = 1e-10) {
  for (dataset in datasets) {
    ref_dir <- file.path(reference, dataset)
    new_dir <- file.path(candidate, dataset)
    files <- sort(list.files(ref_dir))
    if (!identical(files, sort(list.files(new_dir)))) stopf("%s output file set differs", dataset)
    for (filename in files) {
      a <- read.delim(file.path(ref_dir, filename), check.names = FALSE)
      b <- read.delim(file.path(new_dir, filename), check.names = FALSE)
      if (!identical(names(a), names(b)) || !identical(dim(a), dim(b))) stopf("%s/%s schema differs", dataset, filename)
      for (j in seq_along(a)) {
        if (is.numeric(a[[j]])) {
          if (!isTRUE(all.equal(a[[j]], b[[j]], tolerance = tolerance, check.attributes = FALSE))) {
            stopf("%s/%s numeric content differs in %s", dataset, filename, names(a)[j])
          }
        } else if (!identical(a[[j]], b[[j]])) {
          stopf("%s/%s text content differs in %s", dataset, filename, names(a)[j])
        }
      }
    }
  }
}

technical_qc <- function(expr, library_size = NULL, zero_fraction = NULL) {
  n <- ncol(expr)
  pc <- prcomp(t(expr), center = TRUE, scale. = FALSE)$x
  pc5 <- matrix(NA_real_, n, 5L, dimnames = list(colnames(expr), paste0("PC", 1:5)))
  pc5[, seq_len(min(5L, ncol(pc)))] <- pc[, seq_len(min(5L, ncol(pc))), drop = FALSE]
  corr <- apply(cor(expr, method = "spearman"), 2, median)
  low_corr <- median(corr) - 3 * mad(corr, constant = 1)
  flags <- rep("", n)
  add_flag <- function(hit, label) {
    flags[hit] <<- ifelse(nzchar(flags[hit]), paste(flags[hit], label, sep = ";"), label)
  }
  add_flag(corr < low_corr, "low_median_sample_correlation")
  thresholds <- paste0("corr_low=", signif(low_corr, 6))
  if (!is.null(library_size)) {
    log_lib <- log2(library_size)
    low_lib <- median(log_lib) - 3 * mad(log_lib, constant = 1)
    add_flag(log_lib < low_lib, "low_log_library_size")
    thresholds <- paste(thresholds, paste0("log2_library_low=", signif(low_lib, 6)), sep = ";")
  }
  if (!is.null(zero_fraction)) {
    high_zero <- median(zero_fraction) + 3 * mad(zero_fraction, constant = 1)
    add_flag(zero_fraction > high_zero, "high_zero_fraction")
    thresholds <- paste(thresholds, paste0("zero_fraction_high=", signif(high_zero, 6)), sep = ";")
  }
  data.frame(sample_correlation_median = corr, pc5, qc_flag = ifelse(nzchar(flags), flags, "none"),
             threshold = thresholds,
             review_status = ifelse(nzchar(flags), "retain_soft_flag", "retain"), check.names = FALSE)
}

simulate_pipeline <- function() {
  if (!requireNamespace("edgeR", quietly = TRUE) || !requireNamespace("limma", quietly = TRUE)) {
    stopf("Simulation requires edgeR and limma")
  }
  counts <- matrix(c(100, 120, 90, 110, 0, 0, 0, 0, 20, 25, 0, 0), nrow = 3, byrow = TRUE)
  rownames(counts) <- c("KEEP1", "DROP", "KEEP2")
  colnames(counts) <- paste0("S", 1:4)
  group <- factor(c("control", "control", "case", "case"))
  cpm <- edgeR::cpm(counts)
  keep <- rowSums(cpm >= 1) >= ceiling(min(table(group)) * 0.2)
  if (!identical(names(which(keep)), c("KEEP1", "KEEP2"))) stopf("Simulation CPM filter failed")
  dge <- edgeR::DGEList(counts[keep, , drop = FALSE], group = group)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  design <- model.matrix(~ group)
  v <- limma::voom(dge, design, plot = FALSE)
  if (!all(is.finite(v$E)) || !identical(dim(v$E), c(2L, 4L))) stopf("Simulation TMM+voom failed")
  message("SIMULATION_OK CPM filter and TMM+voom")
}

if (simulate) {
  simulate_pipeline()
  quit(save = "no", status = 0)
}

manifest <- read.delim(manifest_path, check.names = FALSE)

dataset <- "GSE251778"
raw_path <- file.path(raw_root, dataset, "GSE251778_Mokhtari_rna_counts.csv.gz")
soft_path <- file.path(root, "data", "raw_geo", paste0(dataset, "_samples_brief.soft"))
expected_sha <- "c78ac40912c8077400c238f95ca3e1d53dbf12912b5d0b02b90fbf2fcdc967f7"
meta <- manifest[manifest$dataset_id == dataset & manifest$include_primary == TRUE, , drop = FALSE]
expected_gsm <- meta$GSM

if (check_only) {
  tmp <- tempfile("rnaseq-check-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  script <- normalizePath("analysis/05_preprocess_rnaseq.R")
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)),
                    env = c(paste0("PREPROCESS_OUTPUT_ROOT=", shQuote(tmp)),
                            paste0("R_LIBS_USER=", shQuote(project_lib))))
  if (status != 0L) stopf("temporary RNA-seq rebuild failed")
  compare_clean_dirs(file.path(root, "data", "clean"), tmp, c("GSE251778", "GSE101521"))
  validate_output(dataset, expected_gsm, raw_path, expected_sha)
  brain_meta <- manifest[manifest$dataset_id == "GSE101521" & manifest$include_primary == TRUE, , drop = FALSE]
  validate_output("GSE101521", brain_meta$GSM,
                  file.path(raw_root, "GSE101521", "GSE101521_totalRNA_counts.csv.gz"),
                  "041dc245eda3effef1664da0ff9341b5a77a608104af27f9b4fea578b6f1eb7c")
  message("READ_ONLY_CHECK_OK RNA-seq formal outputs equal temporary rebuild")
  quit(save = "no", status = 0)
}

for (pkg in c("edgeR", "limma", "AnnotationDbi", "org.Hs.eg.db")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stopf("Missing required package: %s", pkg)
}
if (!file.exists(raw_path)) stopf("Missing official raw counts: %s", raw_path)
if (!identical(sha256(raw_path), expected_sha)) stopf("%s raw SHA mismatch before preprocessing", dataset)

raw <- read.csv(raw_path, check.names = FALSE)
if (!"gene_ID" %in% names(raw)) stopf("%s lacks gene_ID", dataset)
count_names <- setdiff(names(raw), c("", "X", "gene_ID"))
counts <- as.matrix(raw[, count_names, drop = FALSE])
storage.mode(counts) <- "double"
if (any(!is.finite(counts)) || any(counts < 0) || any(counts != floor(counts))) {
  stopf("%s is not a finite non-negative integer raw-count matrix", dataset)
}
rownames(counts) <- sub("\\..*$", "", raw$gene_ID)

desc_to_gsm <- read_soft_descriptions(soft_path)
gsm <- unname(desc_to_gsm[count_names])
if (anyNA(gsm) || any(!nzchar(gsm))) stopf("%s count columns do not all map to GSM via official SOFT descriptions", dataset)
if (anyDuplicated(gsm)) stopf("%s count columns map to duplicated GSM", dataset)
if (!setequal(gsm, expected_gsm)) {
  stopf("%s raw/manifest GSM mismatch; missing=%s extra=%s", dataset,
        paste(setdiff(expected_gsm, gsm), collapse = ","),
        paste(setdiff(gsm, expected_gsm), collapse = ","))
}
counts <- counts[, match(expected_gsm, gsm), drop = FALSE]
colnames(counts) <- expected_gsm

unique_ensembl <- unique(rownames(counts))
symbol_by_ensembl <- AnnotationDbi::mapIds(
  org.Hs.eg.db::org.Hs.eg.db, keys = unique_ensembl,
  keytype = "ENSEMBL", column = "SYMBOL", multiVals = "filter"
)
symbols <- symbol_by_ensembl[rownames(counts)]
reliable <- !is.na(symbols) & nzchar(symbols)
counts <- counts[reliable, , drop = FALSE]
symbols <- unname(symbols[reliable])
counts <- rowsum(counts, group = symbols, reorder = FALSE)

group <- factor(meta$diagnosis, levels = c("control", "MDD"))
if (anyNA(group)) stopf("%s contains unexpected diagnosis labels", dataset)
min_group <- min(table(group))
min_samples <- ceiling(as.numeric(min_group) * 0.2)
keep <- rowSums(edgeR::cpm(counts) >= 1) >= min_samples
filtered <- counts[keep, , drop = FALSE]
if (!nrow(filtered)) stopf("%s CPM filter removed all genes", dataset)

dge <- edgeR::DGEList(filtered, group = group)
dge <- edgeR::calcNormFactors(dge, method = "TMM")
design <- model.matrix(~ group)
v <- limma::voom(dge, design, plot = FALSE)

effective_lib <- dge$samples$lib.size * dge$samples$norm.factors
reconstructed <- log2(sweep(filtered + 0.5, 2, effective_lib + 1, "/") * 1e6)
if (!isTRUE(all.equal(unname(reconstructed), unname(v$E), tolerance = 1e-10))) {
  stopf("%s saved counts/norm factors cannot reconstruct voom QC logCPM", dataset)
}

expr <- data.frame("Gene Symbol" = rownames(v$E), v$E, check.names = FALSE)
out_dir <- file.path(clean_root, dataset)
write_tsv_gz(expr, file.path(out_dir, "expression_gene.tsv.gz"))
write_tsv_gz(data.frame("Gene Symbol" = rownames(filtered), filtered, check.names = FALSE),
             file.path(out_dir, "raw_counts_gene.tsv.gz"))
norm <- data.frame(
  GSM = expected_gsm,
  library_size = dge$samples$lib.size,
  tmm_norm_factor = dge$samples$norm.factors,
  effective_library_size = effective_lib,
  stringsAsFactors = FALSE
)
write.table(norm, file.path(out_dir, "normalization_factors.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(meta, file.path(out_dir, "sample_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
zero_fraction <- colMeans(counts == 0)
tech <- technical_qc(v$E, colSums(counts), zero_fraction)
qc <- data.frame(
  GSM = expected_gsm,
  library_size = colSums(counts),
  tmm_norm_factor = dge$samples$norm.factors,
  effective_library_size = norm$effective_library_size,
  zero_fraction = zero_fraction,
  tech,
  expression_gene_usage = "QC_only_recompute_voom_after_final_design",
  stringsAsFactors = FALSE
)
write.table(qc, file.path(out_dir, "qc_metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

validate_output(dataset, expected_gsm, raw_path, expected_sha)

# User-authorized deviation: the only official GSE101521 individual-level
# total-RNA matrix contains DESeq2-normalized, non-integer counts.  Treat it as
# continuous processed expression (log2(x+1)); never use TMM/voom on this file.
dataset <- "GSE101521"
raw_path <- file.path(raw_root, dataset, "GSE101521_totalRNA_counts.csv.gz")
soft_path <- file.path(root, "data", "raw_geo", paste0(dataset, "_samples_brief.soft"))
expected_sha <- "041dc245eda3effef1664da0ff9341b5a77a608104af27f9b4fea578b6f1eb7c"
meta <- manifest[manifest$dataset_id == dataset & manifest$include_primary == TRUE, , drop = FALSE]
expected_gsm <- meta$GSM
if (!file.exists(raw_path) || !identical(sha256(raw_path), expected_sha)) {
  stopf("%s missing or SHA mismatch", dataset)
}
raw <- read.csv(raw_path, check.names = FALSE)
feature_id <- sub("\\..*$", "", raw[[1]])
processed_names <- names(raw)[-1]
processed <- as.matrix(raw[-1])
storage.mode(processed) <- "double"
if (any(!is.finite(processed)) || any(processed < 0) || all(processed == floor(processed))) {
  stopf("%s is not the expected finite non-negative normalized-count matrix", dataset)
}
title_to_gsm <- read_soft_titles(soft_path)
gsm <- unname(title_to_gsm[processed_names])
if (anyNA(gsm) || anyDuplicated(gsm) || !setequal(gsm, expected_gsm)) {
  stopf("%s official column titles do not map exactly to manifest GSM", dataset)
}
processed <- processed[, match(expected_gsm, gsm), drop = FALSE]
colnames(processed) <- expected_gsm

unique_ensembl <- unique(feature_id)
symbol_by_ensembl <- AnnotationDbi::mapIds(
  org.Hs.eg.db::org.Hs.eg.db, keys = unique_ensembl,
  keytype = "ENSEMBL", column = "SYMBOL", multiVals = "filter"
)
symbols <- symbol_by_ensembl[feature_id]
reliable <- !is.na(symbols) & nzchar(symbols)
processed <- processed[reliable, , drop = FALSE]
symbols <- unname(symbols[reliable])
processed <- rowsum(processed, group = symbols, reorder = FALSE)
min_samples <- ceiling(min(table(meta$diagnosis)) * 0.2)
keep <- rowSums(processed > 0) >= min_samples
processed <- processed[keep, , drop = FALSE]
q99 <- unname(quantile(processed, 0.99))
if (q99 <= 50) stopf("%s distribution does not meet frozen log2 rule (q99=%g)", dataset, q99)
log_expr <- log2(processed + 1)

out_dir <- file.path(clean_root, dataset)
expr <- data.frame("Gene Symbol" = rownames(log_expr), log_expr, check.names = FALSE)
write_tsv_gz(expr, file.path(out_dir, "expression_gene.tsv.gz"))
write.table(meta, file.path(out_dir, "sample_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
qc <- data.frame(
  GSM = expected_gsm,
  processed_library_sum = colSums(processed),
  zero_fraction = colMeans(processed == 0),
  technical_qc(log_expr, colSums(processed), colMeans(processed == 0)),
  transform = "log2(x+1)",
  data_level = "official_normalized_counts",
  stringsAsFactors = FALSE
)
write.table(qc, file.path(out_dir, "qc_metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
validate_output(dataset, expected_gsm, raw_path, expected_sha)
