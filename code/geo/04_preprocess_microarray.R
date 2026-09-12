#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
check_only <- "--check" %in% args
simulate <- "--simulate" %in% args
root <- normalizePath(".", mustWork = TRUE)
project_lib <- file.path(root, ".r-lib")
if (dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))
clean_root <- Sys.getenv("PREPROCESS_OUTPUT_ROOT", unset = file.path(root, "data", "clean"))

stopf <- function(...) stop(sprintf(...), call. = FALSE)
sha256 <- function(path) {
  out <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  if (!length(out)) stopf("SHA256 failed: %s", path)
  strsplit(out[[1]], "[[:space:]]+")[[1]][1]
}
write_tsv_gz <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}
validate_output <- function(dataset, expected_gsm, raw_paths, expected_sha) {
  out <- file.path(clean_root, dataset)
  paths <- file.path(out, c("expression_gene.tsv.gz", "sample_metadata.tsv",
                            "qc_metrics.tsv", "probe_gene_map.tsv"))
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stopf("%s missing outputs: %s", dataset, paste(missing, collapse = ", "))
  expr <- read.delim(paths[1], check.names = FALSE)
  meta <- read.delim(paths[2], check.names = FALSE)
  qc <- read.delim(paths[3], check.names = FALSE)
  map <- read.delim(paths[4], check.names = FALSE)
  if (!identical(names(expr)[1], "Gene Symbol") || anyDuplicated(expr[["Gene Symbol"]])) {
    stopf("%s Gene Symbol is absent or non-unique", dataset)
  }
  if (!identical(names(expr)[-1], expected_gsm)) stopf("%s expression columns differ from manifest", dataset)
  if (!identical(meta$GSM, expected_gsm) || !identical(qc$GSM, expected_gsm)) {
    stopf("%s metadata/QC order differs from manifest", dataset)
  }
  if (any(!is.finite(as.matrix(expr[-1])))) stopf("%s expression contains non-finite values", dataset)
  if (anyDuplicated(map$Gene.Symbol[map$selected_primary])) stopf("%s selected map has duplicated genes", dataset)
  actual <- vapply(raw_paths, sha256, character(1))
  if (!identical(unname(actual), unname(expected_sha))) stopf("%s raw SHA mismatch", dataset)
  message(sprintf("CHECK_OK %s N=%d genes=%d raw_sha256=%s",
                  dataset, length(expected_gsm), nrow(expr), paste(actual, collapse = ",")))
}
compare_clean_dirs <- function(reference, candidate, datasets, tolerance = 1e-10) {
  for (dataset in datasets) {
    files <- sort(list.files(file.path(reference, dataset)))
    if (!identical(files, sort(list.files(file.path(candidate, dataset))))) stopf("%s output file set differs", dataset)
    for (filename in files) {
      a <- read.delim(file.path(reference, dataset, filename), check.names = FALSE)
      b <- read.delim(file.path(candidate, dataset, filename), check.names = FALSE)
      if (!identical(names(a), names(b)) || !identical(dim(a), dim(b))) stopf("%s/%s schema differs", dataset, filename)
      for (j in seq_along(a)) {
        ok <- if (is.numeric(a[[j]])) isTRUE(all.equal(a[[j]], b[[j]], tolerance = tolerance, check.attributes = FALSE)) else identical(a[[j]], b[[j]])
        if (!ok) stopf("%s/%s differs in %s", dataset, filename, names(a)[j])
      }
    }
  }
}
technical_qc_array <- function(expr) {
  pc <- prcomp(t(expr), center = TRUE, scale. = FALSE)$x
  pc5 <- matrix(NA_real_, ncol(expr), 5L, dimnames = list(colnames(expr), paste0("PC", 1:5)))
  pc5[, seq_len(min(5L, ncol(pc)))] <- pc[, seq_len(min(5L, ncol(pc))), drop = FALSE]
  corr <- apply(cor(expr, method = "spearman"), 2, median)
  low <- median(corr) - 3 * mad(corr, constant = 1)
  flag <- corr < low
  data.frame(sample_correlation_median = corr, pc5,
             qc_flag = ifelse(flag, "low_median_sample_correlation", "none"),
             threshold = paste0("corr_low=", signif(low, 6)),
             review_status = ifelse(flag, "retain_soft_flag", "retain"), check.names = FALSE)
}
select_primary_probe <- function(expr, symbols) {
  means <- rowMeans(expr)
  idx <- split(seq_along(symbols), symbols)
  chosen <- vapply(idx, function(i) i[which.max(means[i])], integer(1))
  chosen[order(symbols[chosen])]
}

read_gpl570_symbols <- function(path) {
  head <- readLines(gzfile(path), n = 100L, warn = FALSE)
  header_line <- match(TRUE, grepl("^ID\\t", head))
  if (is.na(header_line)) stopf("GPL570 annotation header not found")
  x <- read.delim(gzfile(path), skip = header_line - 1L, check.names = FALSE, quote = "")
  setNames(x[["Gene symbol"]], x[["ID"]])
}

read_gpl6848_symbols <- function(path) {
  lines <- readLines(path, warn = FALSE)
  begin <- match("!platform_table_begin", lines) + 1L
  end <- match("!platform_table_end", lines) - 1L
  if (is.na(begin) || is.na(end)) stopf("GPL6848 data table is incomplete")
  x <- read.delim(text = lines[begin:end], check.names = FALSE, quote = "")
  setNames(x[["GENE_SYMBOL"]], x[["ID"]])
}
read_gpl32193_symbols <- function(path) {
  head <- readLines(gzfile(path), n = 100L, warn = FALSE)
  header_line <- match(TRUE, grepl("^ID\\t", head))
  if (is.na(header_line)) stopf("GPL32193 annotation header not found")
  x <- read.delim(gzfile(path), skip = header_line - 1L, check.names = FALSE, quote = "")
  setNames(x[["GENE_SYMBOL"]], as.character(x[["ID"]]))
}

process_official_processed <- function(dataset, matrix_path, platform_path, platform_reader,
                                       expected_matrix_sha, expected_platform_sha,
                                       log2_if_q99_gt50 = FALSE) {
  suppressPackageStartupMessages(requireNamespace("GEOquery"))
  suppressPackageStartupMessages(requireNamespace("Biobase"))
  manifest <- read.delim(file.path(root, "data", "derived", "geo_sample_manifest.tsv"), check.names = FALSE)
  meta <- manifest[manifest$dataset_id == dataset & manifest$include_primary == TRUE, , drop = FALSE]
  expected_gsm <- meta$GSM
  if (!identical(sha256(matrix_path), expected_matrix_sha) ||
      !identical(sha256(platform_path), expected_platform_sha)) stopf("%s input SHA mismatch", dataset)
  gse <- GEOquery::getGEO(filename = matrix_path, getGPL = FALSE)
  E <- Biobase::exprs(gse)
  missing_gsm <- setdiff(expected_gsm, colnames(E))
  if (length(missing_gsm) || anyDuplicated(colnames(E))) {
    stopf("%s processed matrix/manifest GSM mismatch; missing=%s duplicated_columns=%d", dataset,
          paste(missing_gsm, collapse = ","), anyDuplicated(colnames(E)))
  }
  E <- E[, match(expected_gsm, colnames(E)), drop = FALSE]
  if (log2_if_q99_gt50) {
    q99 <- unname(quantile(E, 0.99, na.rm = TRUE))
    if (q99 <= 50) stopf("%s does not meet fixed processed-intensity log2 rule", dataset)
    E <- log2(E + 1)
  }
  symbols <- unname(platform_reader(platform_path)[rownames(E)])
  reliable <- !is.na(symbols) & nzchar(symbols) & !grepl("///|;|,", symbols)
  finite_all_samples <- rowSums(!is.finite(E)) == 0L
  variable <- rep(FALSE, nrow(E))
  variable[finite_all_samples] <- apply(E[finite_all_samples, , drop = FALSE], 1, stats::sd) > 0
  keep <- reliable & finite_all_samples & variable
  chosen <- select_primary_probe(E[keep, , drop = FALSE], symbols[keep])
  geneE <- E[keep, , drop = FALSE][chosen, , drop = FALSE]
  rownames(geneE) <- symbols[keep][chosen]
  probe_ids <- rownames(E)[keep]
  probe_map <- data.frame(
    probe_id = rownames(E),
    Gene.Symbol = symbols,
    reliable_annotation = reliable,
    finite_all_included_samples = finite_all_samples,
    nonzero_variance = variable,
    mean_expression_all_included_samples = rowMeans(E, na.rm = TRUE),
    selected_primary = rownames(E) %in% probe_ids[chosen],
    stringsAsFactors = FALSE
  )
  out <- file.path(clean_root, dataset)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write_tsv_gz(data.frame("Gene Symbol" = rownames(geneE), geneE, check.names = FALSE),
               file.path(out, "expression_gene.tsv.gz"))
  write.table(meta, file.path(out, "sample_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  qc <- data.frame(
    GSM = expected_gsm,
    expression_median = apply(geneE, 2, median),
    expression_IQR = apply(geneE, 2, IQR),
    technical_qc_array(geneE),
    data_level = "official_processed_expression",
    deviation = "deviation_requested_by_user",
    stringsAsFactors = FALSE
  )
  write.table(qc, file.path(out, "qc_metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(probe_map, file.path(out, "probe_gene_map.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  validate_output(dataset, expected_gsm, c(matrix_path, platform_path),
                  c(expected_matrix_sha, expected_platform_sha))
}

if (simulate) {
  x <- matrix(c(8, 8, 8, 9, 9, 9, 4, 4, 4), nrow = 3, byrow = TRUE)
  rownames(x) <- c("P1", "P2", "P3")
  chosen <- select_primary_probe(x, c("A", "A", "B"))
  if (!identical(rownames(x)[chosen], c("P2", "P3"))) stopf("Simulation primary-probe rule failed")
  message("SIMULATION_OK highest-across-all-samples primary-probe rule")
  quit(save = "no", status = 0)
}

if (check_only) {
  tmp_check <- tempfile("microarray-check-")
  dir.create(tmp_check)
  on.exit(unlink(tmp_check, recursive = TRUE), add = TRUE)
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    c("--vanilla", shQuote(normalizePath("analysis/04_preprocess_microarray.R"))),
                    env = c(paste0("PREPROCESS_OUTPUT_ROOT=", shQuote(tmp_check)),
                            paste0("R_LIBS_USER=", shQuote(project_lib))))
  if (status != 0L) stopf("temporary microarray rebuild failed")
  compare_clean_dirs(file.path(root, "data", "clean"), tmp_check,
                     c("GSE98793", "GSE19738", "GSE39653"))
  message("READ_ONLY_CHECK_OK microarray formal outputs equal temporary rebuild")
  quit(save = "no", status = 0)
}

gse98793_matrix <- file.path(root, "data", "raw_geo", "expression", "GSE98793",
                             "GSE98793_series_matrix.txt.gz")
gpl570_annot <- file.path(root, "data", "raw_geo", "platforms", "GPL570.annot.gz")
process_official_processed(
  "GSE98793", gse98793_matrix, gpl570_annot, read_gpl570_symbols,
  "ea7e32f3a38b17f282450aae16c47fea53f3e36aeecca6c7d5f52b0d0489fa91",
  "d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394"
)


gse19738_matrix <- file.path(root, "data", "raw_geo", "expression", "GSE19738",
                             "GSE19738_series_matrix.txt.gz")
gpl6848_annot <- file.path(root, "data", "raw_geo", "platforms", "GPL6848_data.soft")
process_official_processed(
  "GSE19738", gse19738_matrix, gpl6848_annot, read_gpl6848_symbols,
  "8285b10b917939acc5c4f36f7a69a1b56c6a024ec591f7b6fd5e8a96f49e5811",
  "13b4004c352b81009e6bbd192cfc607860a366e93c19288f190ac193e6c332dd"
)

dataset <- "GSE39653"
manifest <- read.delim(file.path(root, "data", "derived", "geo_sample_manifest.tsv"), check.names = FALSE)
meta <- manifest[manifest$dataset_id == dataset & manifest$include_primary == TRUE, , drop = FALSE]
expected_gsm <- meta$GSM
raw_dir <- file.path(root, "data", "raw_geo", "expression", dataset)
raw_matrix <- file.path(raw_dir, "GSE39653_non-normalized.txt.gz")
series_matrix <- file.path(raw_dir, "GSE39653_series_matrix.txt.gz")
raw_tar <- file.path(raw_dir, "GSE39653_RAW.tar")
raw_paths <- c(raw_matrix, series_matrix, raw_tar)
expected_sha <- c(
  "9ecc9fd34ed73797311fd134cb16a8f773387743bee76a3d69145bb32fc07d4c",
  "947b26302824d9ed95437003063199b8d0e0414b6c5f9c21877967bfd89b53e1",
  "19dbc61b4da9ac807799a96081f4eae405079dad108b82a9d2bdb74678cdb948"
)

for (pkg in c("limma", "GEOquery", "Biobase")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stopf("Missing required package: %s", pkg)
}
if (any(!file.exists(raw_paths))) stopf("%s missing official raw/mapping inputs", dataset)

raw <- read.delim(raw_matrix, check.names = FALSE)
intensity_cols <- seq(2, ncol(raw), 2)
detection_cols <- seq(3, ncol(raw), 2)
E <- as.matrix(raw[, intensity_cols, drop = FALSE])
D <- as.matrix(raw[, detection_cols, drop = FALSE])
storage.mode(E) <- storage.mode(D) <- "double"
rownames(E) <- rownames(D) <- raw[[1]]
if (ncol(E) != 53L || any(!is.finite(E)) || any(!is.finite(D))) stopf("%s malformed non-normalized matrix", dataset)

# The non-normalized file uses anonymous SAMPLE n labels.  The official series
# matrix is used only to prove the GSM mapping, never as the analysis input.
gse <- GEOquery::getGEO(filename = series_matrix, getGPL = FALSE)
published <- Biobase::exprs(gse)
common <- intersect(rownames(E), rownames(published))
if (length(common) < 10000L) stopf("%s insufficient probes for label verification", dataset)
C <- cor(E[common, , drop = FALSE], published[common, , drop = FALSE], method = "spearman")
best <- apply(C, 1, which.max)
best_cor <- apply(C, 1, max)
second_cor <- apply(C, 1, function(z) sort(z, decreasing = TRUE)[2])
gsm_all <- colnames(published)[best]
if (length(unique(gsm_all)) != ncol(E) || any(best_cor < 0.999999) || any(best_cor - second_cor < 0.05)) {
  stopf("%s anonymous columns do not have unique near-perfect GSM mappings", dataset)
}
if (!all(expected_gsm %in% gsm_all)) stopf("%s manifest GSM absent from verified mapping", dataset)
sel <- match(expected_gsm, gsm_all)
E <- E[, sel, drop = FALSE]
D <- D[, sel, drop = FALSE]
colnames(E) <- colnames(D) <- expected_gsm

tmp <- tempfile("gpl10558-")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
member <- "GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt.gz"
utils::untar(raw_tar, files = member, exdir = tmp)
lines <- readLines(gzfile(file.path(tmp, member)), warn = FALSE)
begin <- match("[Probes]", lines) + 1L
end <- match("[Controls]", lines) - 1L
if (is.na(begin) || is.na(end)) stopf("%s platform manifest lacks probe/control sections", dataset)
annot <- read.delim(text = lines[begin:end], check.names = FALSE, quote = "")
annot <- annot[!duplicated(annot$Probe_Id), , drop = FALSE]
symbol <- annot$Symbol[match(rownames(E), annot$Probe_Id)]
reliable <- !is.na(symbol) & nzchar(symbol) & symbol != "NA"
detected <- rowSums(D <= 0.05) >= ceiling(ncol(D) * 0.2)
keep <- reliable & detected
if (!any(keep)) stopf("%s probe QC removed all probes", dataset)

elist <- methods::new("EListRaw")
elist$E <- E
elist$other$Detection <- D
normalized <- limma::neqc(elist, offset = 16, robust = TRUE)
normE <- normalized$E[keep, , drop = FALSE]
symbols <- symbol[keep]
chosen <- select_primary_probe(normE, symbols)
geneE <- normE[chosen, , drop = FALSE]
rownames(geneE) <- symbols[chosen]

probe_map <- data.frame(
  probe_id = rownames(E),
  Gene.Symbol = symbol,
  reliable_annotation = reliable,
  detection_qc_pass = detected,
  mean_expression_all_included_samples = rowMeans(normalized$E),
  selected_primary = FALSE,
  stringsAsFactors = FALSE
)
probe_map$selected_primary[match(rownames(normE)[chosen], probe_map$probe_id)] <- TRUE

out <- file.path(clean_root, dataset)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
expr <- data.frame("Gene Symbol" = rownames(geneE), geneE, check.names = FALSE)
write_tsv_gz(expr, file.path(out, "expression_gene.tsv.gz"))
write.table(meta, file.path(out, "sample_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
qc <- data.frame(
  GSM = expected_gsm,
  raw_median_intensity = apply(E, 2, median),
  detected_probe_fraction = colMeans(D <= 0.05),
  technical_qc_array(geneE),
  anonymous_label_mapping_correlation = best_cor[sel],
  anonymous_label_mapping_margin = best_cor[sel] - second_cor[sel],
  data_level = "official_non_normalized_intensity_with_detection_p",
  deviation = "pre_registered_platform_preprocessing",
  series_matrix_usage = "anonymous_column_to_GSM_mapping_only",
  stringsAsFactors = FALSE
)
write.table(qc, file.path(out, "qc_metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(probe_map, file.path(out, "probe_gene_map.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
validate_output(dataset, expected_gsm, raw_paths, expected_sha)
