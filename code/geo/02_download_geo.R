#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

datasets <- c(
  "GSE98793", "GSE251778", "GSE201332", "GSE19738",
  "GSE39653", "GSE101521", "GSE144136",
  "GSE38206", "GSE80655", "GSE102556"
)

args <- commandArgs(trailingOnly = TRUE)
check_only <- "--check" %in% args
self_test <- "--self-test" %in% args
project_root <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(project_root, "code", "geo", "README.md"))) {
  stop("Run from the project root.")
}

raw_dir <- file.path(project_root, "data", "raw_geo")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
download_log_path <- file.path(raw_dir, "download_manifest.tsv")
checksums_path <- file.path(raw_dir, "CHECKSUMS.sha256")

geo_bucket <- function(gse) {
  number <- as.integer(sub("^GSE", "", gse))
  sprintf("GSE%dnnn", number %/% 1000L)
}

sources <- do.call(rbind, lapply(datasets, function(gse) {
  bucket <- geo_bucket(gse)
  data.frame(
    dataset_id = gse,
    kind = c("series_brief_soft", "samples_brief_soft", "supplementary_index"),
    url = c(
      sprintf("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&targ=self&form=text&view=brief", gse),
      sprintf("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&targ=gsm&form=text&view=brief", gse),
      sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/", bucket, gse)
    ),
    filename = c(
      sprintf("%s_series_brief.soft", gse),
      sprintf("%s_samples_brief.soft", gse),
      sprintf("%s_supplementary_index.html", gse)
    ),
    stringsAsFactors = FALSE
  )
}))

sha256_file <- function(path) {
  out <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  if (length(out) != 1L) stop("Unable to compute SHA256 for ", path)
  strsplit(out, "[[:space:]]+")[[1L]][1L]
}

assert_metadata_only <- function(path, kind) {
  if (!grepl("soft$", path)) return(invisible(TRUE))
  lines <- readLines(path, warn = FALSE)
  forbidden <- c("!sample_table_begin", "!series_matrix_table_begin")
  hits <- vapply(forbidden, function(x) any(grepl(x, lines, fixed = TRUE)), logical(1))
  if (any(hits)) {
    stop("Expression table marker found in metadata-only file: ", basename(path))
  }
  if (kind == "samples_brief_soft" && !any(grepl("^\\^SAMPLE = GSM", lines))) {
    stop("No GEO sample records found in ", basename(path))
  }
  invisible(TRUE)
}

verify_snapshot <- function(log_path, checksum_path, data_dir = raw_dir) {
  problems <- character()
  if (!file.exists(log_path) || !file.exists(checksum_path)) {
    return("raw_manifest_or_checksums_missing")
  }
  log <- read.delim(log_path, check.names = FALSE, colClasses = "character")
  required_columns <- c("dataset_id", "kind", "filename", "url", "accessed_at",
                        "file_size_bytes", "sha256")
  if (!identical(names(log), required_columns)) return("download_manifest_schema_mismatch")
  expected <- sources$filename
  if (!setequal(log$filename, expected) || anyDuplicated(log$filename)) {
    problems <- c(problems, "download_manifest_file_set_mismatch")
  }
  expected_order <- match(log$filename, sources$filename)
  if (anyNA(expected_order) || !identical(log$url, sources$url[expected_order])) {
    problems <- c(problems, "download_manifest_URL_mismatch")
  }
  iso8601_tz <- "^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}(Z|[+-][0-9]{2}:[0-9]{2})$"
  if (any(!grepl(iso8601_tz, log$accessed_at)) || length(unique(log$accessed_at)) != 1L) {
    problems <- c(problems, "download_manifest_first_accessed_at_mismatch")
  }
  checksum_lines <- readLines(checksum_path, warn = FALSE)
  checksum_match <- regexec("^([0-9a-f]{64})  ([^/]+)$", checksum_lines)
  checksum_parts <- regmatches(checksum_lines, checksum_match)
  if (length(checksum_parts) != nrow(log) || any(lengths(checksum_parts) != 3L)) {
    problems <- c(problems, "checksums_schema_mismatch")
    checksums <- data.frame(filename = character(), sha256 = character())
  } else {
    checksums <- data.frame(
      filename = vapply(checksum_parts, `[`, character(1), 3L),
      sha256 = vapply(checksum_parts, `[`, character(1), 2L),
      stringsAsFactors = FALSE
    )
    if (anyDuplicated(checksums$filename) || !identical(checksums$filename, log$filename) ||
        !identical(checksums$sha256, log$sha256)) {
      problems <- c(problems, "checksums_download_manifest_mismatch")
    }
  }
  for (i in seq_len(nrow(log))) {
    path <- file.path(data_dir, log$filename[i])
    if (!file.exists(path)) {
      problems <- c(problems, "raw_metadata_missing")
      next
    }
    actual_size <- file.info(path)$size
    actual_sha <- sha256_file(path)
    if (!identical(as.numeric(log$file_size_bytes[i]), as.numeric(actual_size))) {
      problems <- c(problems, "raw_metadata_size_mismatch")
    }
    if (!identical(log$sha256[i], actual_sha)) problems <- c(problems, "raw_metadata_SHA_mismatch")
    tryCatch(assert_metadata_only(path, log$kind[i]),
             error = function(e) problems <<- c(problems, "expression_table_marker_found"))
  }
  unique(problems)
}

if (self_test) {
  baseline <- verify_snapshot(download_log_path, checksums_path)
  if (length(baseline)) stop("Baseline snapshot failed: ", paste(baseline, collapse = ","))
  td <- tempfile("geo-download-self-test-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  checksum_copy <- file.path(td, "CHECKSUMS.sha256")
  file.copy(checksums_path, checksum_copy)
  Sys.chmod(checksum_copy, "0644")
  x <- readLines(checksum_copy)
  x[1L] <- sub("^[0-9a-f]", ifelse(substr(x[1L], 1L, 1L) == "0", "1", "0"), x[1L])
  writeLines(x, checksum_copy)
  detected <- verify_snapshot(download_log_path, checksum_copy)
  if (!"checksums_download_manifest_mismatch" %in% detected) {
    stop("Checksum tamper self-test failed; detected: ", paste(detected, collapse = ","))
  }
  manifest_copy <- file.path(td, "download_manifest.tsv")
  file.copy(download_log_path, manifest_copy)
  Sys.chmod(manifest_copy, "0644")
  log <- read.delim(manifest_copy, check.names = FALSE, colClasses = "character")
  log$accessed_at[1L] <- "2026/07/24 14:30:03"
  write.table(log, manifest_copy, sep = "\t", quote = FALSE, row.names = FALSE)
  detected_time <- verify_snapshot(manifest_copy, checksums_path)
  if (!"download_manifest_first_accessed_at_mismatch" %in% detected_time) {
    stop("Non-ISO accessed_at self-test failed; detected: ", paste(detected_time, collapse = ","))
  }
  message("PASS: real-file CHECKSUMS tampering and non-ISO snapshot time detected.")
  quit(save = "no", status = 0L)
}

if (check_only) {
  problems <- verify_snapshot(download_log_path, checksums_path)
  if (length(problems)) stop("Raw snapshot check failed: ", paste(problems, collapse = ", "))
  message("PASS: CHECKSUMS, download manifest, 30 raw files, and metadata-only guards agree.")
  quit(save = "no", status = 0L)
}

download_one <- function(url, destination) {
  tmp <- paste0(destination, ".part")
  on.exit(unlink(tmp), add = TRUE)
  status <- system2(
    "curl",
    c("-L", "--fail", "--http1.1", "--retry", "5", "--retry-all-errors",
      "--retry-delay", "2", "--connect-timeout", "30",
      "--max-time", "300", "-sS", shQuote(url), "-o", shQuote(tmp))
  )
  if (!identical(status, 0L)) stop("Download failed: ", url)
  if (file.info(tmp)$size == 0) stop("Empty response: ", url)
  if (!isTRUE(file.rename(tmp, destination))) stop("Atomic rename failed: ", destination)
}

first_acquired_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
first_acquired_at <- sub("([+-][0-9]{2})([0-9]{2})$", "\\1:\\2", first_acquired_at)
records <- vector("list", nrow(sources))
for (i in seq_len(nrow(sources))) {
  destination <- file.path(raw_dir, sources$filename[i])
  message(sprintf("[%d/%d] %s", i, nrow(sources), sources$filename[i]))
  download_one(sources$url[i], destination)
  assert_metadata_only(destination, sources$kind[i])
  records[[i]] <- data.frame(
    dataset_id = sources$dataset_id[i],
    kind = sources$kind[i],
    filename = sources$filename[i],
    url = sources$url[i],
    accessed_at = first_acquired_at,
    file_size_bytes = file.info(destination)$size,
    sha256 = sha256_file(destination),
    stringsAsFactors = FALSE
  )
}

download_log <- do.call(rbind, records)
write.table(download_log, download_log_path, sep = "\t", quote = FALSE, row.names = FALSE)

checksum_lines <- sprintf("%s  %s", download_log$sha256, download_log$filename)
writeLines(checksum_lines, checksums_path, useBytes = TRUE)

Sys.chmod(list.files(raw_dir, full.names = TRUE), mode = "0444")
message("DONE: downloaded ", nrow(download_log), " metadata/index files; no expression tables.")
