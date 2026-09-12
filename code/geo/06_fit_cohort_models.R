#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(20260724)

args <- commandArgs(trailingOnly = TRUE)
check_only <- "--check" %in% args
simulate_only <- "--simulate" %in% args
root <- normalizePath(".", mustWork = TRUE)
project_lib <- file.path(root, ".r-lib")
if (dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))
output_root <- Sys.getenv("COHORT_MODEL_OUTPUT_ROOT", unset = root)

stopf <- function(...) stop(sprintf(...), call. = FALSE)
for (pkg in c("limma", "edgeR")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stopf("Missing required package: %s", pkg)
}

write_tsv <- function(x, path, gzip = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- if (gzip) gzfile(path, "wt") else file(path, "wt")
  on.exit(close(con), add = TRUE)
  write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

read_expression <- function(dataset) {
  x <- read.delim(file.path(root, "data", "clean", dataset, "expression_gene.tsv.gz"),
                  check.names = FALSE)
  gene <- x[[1]]
  m <- as.matrix(x[-1])
  storage.mode(m) <- "double"
  rownames(m) <- gene
  if (anyDuplicated(gene) || any(!is.finite(m))) stopf("%s invalid clean expression", dataset)
  m
}

read_metadata <- function(dataset) {
  x <- read.delim(file.path(root, "data", "clean", dataset, "sample_metadata.tsv"),
                  check.names = FALSE, na.strings = c("", "NA"))
  if (anyDuplicated(x$GSM)) stopf("%s duplicated GSM metadata", dataset)
  x
}

normalize_missing <- function(x) {
  x <- trimws(as.character(x))
  x[x == "" | toupper(x) == "NR"] <- NA_character_
  x
}

prepare_metadata <- function(meta, requested) {
  meta$diagnosis <- factor(meta$diagnosis, levels = c("control", "MDD"))
  if (anyNA(meta$diagnosis)) stopf("%s unexpected diagnosis; control reference required",
                                   unique(meta$dataset_id))
  used <- character()
  records <- list()
  for (v in requested) {
    raw <- normalize_missing(meta[[v]])
    status <- reason <- ""
    if (all(is.na(raw))) {
      status <- "excluded"; reason <- "NR"
    } else if (anyNA(raw)) {
      status <- "excluded"; reason <- "partial_missing_no_imputation"
    } else {
      value <- if (v %in% c("age", "RIN", "PMI")) suppressWarnings(as.numeric(raw)) else factor(raw)
      if (anyNA(value)) {
        status <- "excluded"; reason <- "non_numeric_or_invalid"
      } else if (length(unique(value)) < 2L) {
        status <- "excluded"; reason <- "constant"
      } else {
        trial <- meta
        trial[[v]] <- value
        terms <- c("diagnosis", used, v)
        mm <- model.matrix(reformulate(terms), trial)
        if (qr(mm)$rank < ncol(mm)) {
          status <- "excluded"; reason <- "not_identifiable_or_collinear"
        } else {
          meta[[v]] <- value
          used <- c(used, v)
          status <- "included"; reason <- "pre_specified_available"
        }
      }
    }
    records[[length(records) + 1L]] <- data.frame(
      covariate = v, status = status, reason = reason, stringsAsFactors = FALSE)
  }
  design <- model.matrix(reformulate(c("diagnosis", used)), meta)
  if (qr(design)$rank < ncol(design)) stopf("%s BLOCKED: rank-deficient final design",
                                            unique(meta$dataset_id))
  coef_name <- grep("^diagnosisMDD$", colnames(design), value = TRUE)
  if (length(coef_name) != 1L) stopf("%s MDD coefficient not uniquely encoded",
                                     unique(meta$dataset_id))
  list(meta = meta, used = used, audit = do.call(rbind, records),
       design = design, coef = coef_name)
}

effect_table <- function(fit, coef_name, ncase, ncontrol, covariates, genes = rownames(fit$coefficients)) {
  j <- match(coef_name, colnames(fit$coefficients))
  beta <- fit$coefficients[, j]
  tt <- fit$t[, j]
  se_unmoderated <- fit$stdev.unscaled[, j] * fit$sigma
  se_moderated <- fit$stdev.unscaled[, j] * sqrt(fit$s2.post)
  p <- fit$p.value[, j]
  if (any(!is.finite(beta)) ||
      any(!is.finite(se_unmoderated)) || any(se_unmoderated <= 0) ||
      any(!is.finite(se_moderated)) || any(se_moderated <= 0) ||
      any(!is.finite(tt)) || any(!is.finite(p))) {
    stopf("Non-finite cohort effect statistics")
  }
  data.frame(
    gene = genes, beta = beta, logFC = beta,
    SE = se_moderated, SE_unmoderated = se_unmoderated,
    SE_moderated = se_moderated, moderated_t = tt,
    residual_scale = fit$sigma, df = fit$df.residual,
    df_residual = fit$df.residual, df_total = fit$df.total,
    p = p, FDR = p.adjust(p, method = "BH"),
    N = ncase + ncontrol, Ncase = ncase, Ncontrol = ncontrol,
    covariates = if (length(covariates)) paste(covariates, collapse = "+") else "none",
    direction = ifelse(beta > 0, "MDD_higher", ifelse(beta < 0, "MDD_lower", "zero")),
    stringsAsFactors = FALSE, check.names = FALSE
  )
}

fit_limma <- function(expr, design_info) {
  fit <- limma::lmFit(expr, design_info$design)
  fit <- limma::eBayes(fit) # frozen standard default: robust=FALSE
  tab <- effect_table(fit, design_info$coef,
                      sum(design_info$meta$diagnosis == "MDD"),
                      sum(design_info$meta$diagnosis == "control"),
                      design_info$used)
  list(table = tab, fit = fit)
}

validate_alignment <- function(expr, meta) {
  if (!identical(colnames(expr), meta$GSM)) {
    stopf("Sample mismatch: expression columns are not identical to metadata GSM order")
  }
}

fit_microarray <- function(dataset, requested) {
  expr <- read_expression(dataset)
  meta <- read_metadata(dataset)
  validate_alignment(expr, meta)
  di <- prepare_metadata(meta, requested)
  ans <- fit_limma(expr, di)
  list(table = ans$table, audit = di$audit, used = di$used,
       design = colnames(di$design), ncase = sum(meta$diagnosis == "MDD"),
       ncontrol = sum(meta$diagnosis == "control"))
}

fit_rnaseq_stratum <- function(counts, meta, norm, sex_value) {
  take <- meta$sex == sex_value
  smeta <- meta[take, , drop = FALSE]
  scounts <- counts[, smeta$GSM, drop = FALSE]
  snorm <- norm[match(smeta$GSM, norm$GSM), , drop = FALSE]
  if (anyNA(snorm$GSM)) stopf("GSE251778 norm-factor sample mismatch")
  di <- prepare_metadata(smeta, c("age", "medication", "batch"))
  dge <- edgeR::DGEList(counts = scounts)
  dge$samples$lib.size <- snorm$library_size
  dge$samples$norm.factors <- snorm$tmm_norm_factor
  v <- limma::voom(dge, di$design, plot = FALSE)
  fit <- limma::eBayes(limma::lmFit(v, di$design))
  tab <- effect_table(fit, di$coef, sum(smeta$diagnosis == "MDD"),
                      sum(smeta$diagnosis == "control"), di$used)
  list(table = tab, audit = di$audit, design = colnames(di$design),
       weights = v$weights)
}

combine_strata <- function(a, b) {
  if (!identical(a$gene, b$gene)) stopf("Sex-stratum genes differ")
  wa_mod <- 1 / a$SE_moderated^2
  wb_mod <- 1 / b$SE_moderated^2
  se_mod <- sqrt(1 / (wa_mod + wb_mod))
  beta_mod <- (wa_mod * a$beta + wb_mod * b$beta) / (wa_mod + wb_mod)
  wa_unmod <- 1 / a$SE_unmoderated^2
  wb_unmod <- 1 / b$SE_unmoderated^2
  se_unmod <- sqrt(1 / (wa_unmod + wb_unmod))
  beta_unmod <- (wa_unmod * a$beta + wb_unmod * b$beta) / (wa_unmod + wb_unmod)
  z <- beta_mod / se_mod
  p <- 2 * pnorm(-abs(z))
  q <- wa_mod * (a$beta - beta_mod)^2 + wb_mod * (b$beta - beta_mod)^2
  q_p <- pchisq(q, df = 1, lower.tail = FALSE)
  i2 <- ifelse(q > 0, pmax(0, (q - 1) / q) * 100, 0)
  same_direction <- sign(a$beta) == sign(b$beta)
  heterogeneity_flag <- q_p < 0.05 | i2 > 75 | !same_direction
  data.frame(
    gene = a$gene, beta = beta_mod, logFC = beta_mod,
    SE = se_mod, SE_unmoderated = se_unmod, SE_moderated = se_mod,
    beta_unmoderated_combined = beta_unmod,
    beta_moderated_combined = beta_mod, moderated_z = z,
    beta_male = a$beta, beta_female = b$beta,
    SE_unmoderated_male = a$SE_unmoderated,
    SE_unmoderated_female = b$SE_unmoderated,
    SE_moderated_male = a$SE_moderated,
    SE_moderated_female = b$SE_moderated,
    same_direction = same_direction, Cochran_Q = q, Q_p = q_p, I2 = i2,
    heterogeneity_flag = heterogeneity_flag,
    residual_scale = NA_real_, df = a$df_residual + b$df_residual,
    df_residual = a$df_residual + b$df_residual, df_total = Inf,
    p = p, FDR = p.adjust(p, "BH"),
    N = a$Ncase + b$Ncase + a$Ncontrol + b$Ncontrol,
    Ncase = a$Ncase + b$Ncase, Ncontrol = a$Ncontrol + b$Ncontrol,
    covariates = "sex_stratified_inverse_variance_fixed_effect",
    direction = ifelse(beta_mod > 0, "MDD_higher", ifelse(beta_mod < 0, "MDD_lower", "zero")),
    stringsAsFactors = FALSE, check.names = FALSE
  )
}

fit_gse251778 <- function() {
  dataset <- "GSE251778"
  raw <- read.delim(file.path(root, "data", "clean", dataset, "raw_counts_gene.tsv.gz"),
                    check.names = FALSE)
  counts <- as.matrix(raw[-1]); storage.mode(counts) <- "double"; rownames(counts) <- raw[[1]]
  meta <- read_metadata(dataset)
  validate_alignment(counts, meta)
  norm <- read.delim(file.path(root, "data", "clean", dataset, "normalization_factors.tsv"),
                     check.names = FALSE)
  required_norm <- c("GSM", "library_size", "tmm_norm_factor", "effective_library_size")
  if (!all(required_norm %in% names(norm)) || anyDuplicated(norm$GSM) ||
      !setequal(norm$GSM, colnames(counts))) {
    stopf("%s invalid normalization-factor GSM schema, duplicates, missing, or extra samples", dataset)
  }
  norm <- norm[match(colnames(counts), norm$GSM), , drop = FALSE]
  if (!isTRUE(all.equal(as.numeric(norm$library_size), as.numeric(colSums(counts)),
                        tolerance = 1e-10, check.attributes = FALSE))) {
    stopf("%s library_size does not equal raw-count column sum", dataset)
  }
  if (any(!is.finite(norm$library_size)) || any(norm$library_size <= 0) ||
      any(!is.finite(norm$tmm_norm_factor)) || any(norm$tmm_norm_factor <= 0) ||
      any(!is.finite(norm$effective_library_size)) || any(norm$effective_library_size <= 0) ||
      !isTRUE(all.equal(norm$effective_library_size,
                        norm$library_size * norm$tmm_norm_factor,
                        tolerance = 1e-8, check.attributes = FALSE))) {
    stopf("%s invalid normalization factors or effective library sizes", dataset)
  }
  male <- fit_rnaseq_stratum(counts, meta, norm, "Male")
  female <- fit_rnaseq_stratum(counts, meta, norm, "Female")
  combined <- combine_strata(male$table, female$table)
  audit <- rbind(
    transform(male$audit, stratum = "Male"),
    transform(female$audit, stratum = "Female"),
    data.frame(covariate = c("sex", "platform_id"), status = c("stratified", "excluded"),
               reason = c("pre_specified_sex_stratification",
                          "perfectly_coupled_with_sex_not_jointly_identifiable"),
               stratum = "combined")
  )
  list(table = combined, male = male$table, female = female$table,
       audit = audit, design = paste("Male:", paste(male$design, collapse = ","),
                                    "| Female:", paste(female$design, collapse = ",")),
       ncase = sum(meta$diagnosis == "MDD"), ncontrol = sum(meta$diagnosis == "control"))
}

expect_error <- function(expr, pattern) {
  ok <- FALSE
  tryCatch(force(expr), error = function(e) ok <<- grepl(pattern, conditionMessage(e)))
  if (!ok) stopf("Negative test did not detect: %s", pattern)
}

simulate_pipeline <- function() {
  set.seed(20260724)
  n <- 24L
  meta <- data.frame(dataset_id = "SIM", GSM = paste0("S", seq_len(n)),
                     diagnosis = rep(c("control", "MDD"), each = n / 2),
                     age = seq(20, 66, length.out = n), sex = rep(c("F", "M"), n / 2),
                     medication = "NR", batch = rep(c("B1", "B2"), each = n / 2))
  design_info <- prepare_metadata(meta, c("age", "sex", "medication"))
  expr <- matrix(rnorm(80 * n, sd = 0.2), 80, n,
                 dimnames = list(paste0("G", seq_len(80)), meta$GSM))
  expr[1:5, meta$diagnosis == "MDD"] <- expr[1:5, meta$diagnosis == "MDD"] + 1.2
  ans <- fit_limma(expr, design_info)$table
  if (any(ans$beta[1:5] <= 0.8)) stopf("Simulation direction/effect recovery failed")
  if (any(!is.finite(ans$SE_unmoderated)) || any(ans$SE_unmoderated <= 0) ||
      any(!is.finite(ans$SE_moderated)) || any(ans$SE_moderated <= 0)) {
    stopf("Simulation dual-SE recovery failed")
  }

  counts <- matrix(rnbinom(80 * n, mu = 100, size = 10), 80, n,
                   dimnames = list(rownames(expr), meta$GSM))
  counts[1:5, meta$diagnosis == "MDD"] <- counts[1:5, meta$diagnosis == "MDD"] * 4L
  dge <- edgeR::calcNormFactors(edgeR::DGEList(counts))
  v <- limma::voom(dge, design_info$design, plot = FALSE)
  if (!all(is.finite(v$weights)) || length(unique(signif(c(v$weights), 6))) < 2L) {
    stopf("Simulation voom weights invalid")
  }
  weighted_fit <- limma::eBayes(limma::lmFit(v, design_info$design))
  vf <- effect_table(weighted_fit, design_info$coef,
                     sum(design_info$meta$diagnosis == "MDD"),
                     sum(design_info$meta$diagnosis == "control"),
                     design_info$used)
  if (any(vf$beta[1:5] <= 0)) stopf("Simulation voom direction recovery failed")

  s1 <- ans; s2 <- ans
  s1$beta <- s1$logFC <- ans$beta + 0.05
  s2$beta <- s2$logFC <- ans$beta - 0.05
  cm <- combine_strata(s1, s2)
  if (any(cm$beta[1:5] <= 0.8) ||
      any(cm$SE_moderated >= pmax(s1$SE_moderated, s2$SE_moderated)) ||
      any(cm$SE_unmoderated >= pmax(s1$SE_unmoderated, s2$SE_unmoderated))) {
    stopf("Simulation sex-stratum fixed-effect combination failed")
  }

  reversed <- meta
  reversed$diagnosis <- ifelse(reversed$diagnosis == "MDD", "control", "MDD")
  original_case_ids <- meta$GSM[meta$diagnosis == "MDD"]
  expect_error({
    if (!identical(reversed$GSM[reversed$diagnosis == "MDD"], original_case_ids))
      stopf("diagnosis reversal detected")
  }, "diagnosis reversal")
  expect_error(validate_alignment(expr[, rev(seq_len(n))], meta), "Sample mismatch")
  tampered <- ans; tampered$beta[1] <- tampered$beta[1] + 1
  expect_error({
    if (!isTRUE(all.equal(ans, tampered, tolerance = 1e-12))) stopf("effect table tampering detected")
  }, "tampering")
  message("SIMULATION_OK direction dual_SE weighted_voom_fit sex_fixed_effect negative_tests_label_and_integrity")
}

compare_tables <- function(reference_root, candidate_root) {
  rel <- c(
    file.path("data", "derived", "cohort_effects",
              c("GSE98793_effects.tsv.gz", "GSE19738_effects.tsv.gz",
                "GSE39653_effects.tsv.gz", "GSE251778_effects.tsv.gz",
                "GSE251778_Male_effects.tsv.gz", "GSE251778_Female_effects.tsv.gz",
                "GSE101521_effects.tsv.gz")),
    file.path("results", "tables",
              c("cohort_deg_summary.tsv", "cohort_model_audit.tsv",
                "cohort_model_session.tsv",
                "GSE98793_DEG.tsv.gz", "GSE19738_DEG.tsv.gz",
                "GSE39653_DEG.tsv.gz", "GSE251778_DEG.tsv.gz",
                "GSE101521_DEG.tsv.gz"))
  )
  for (f in rel) {
    a <- read.delim(file.path(reference_root, f), check.names = FALSE)
    b <- read.delim(file.path(candidate_root, f), check.names = FALSE)
    if (!isTRUE(all.equal(a, b, tolerance = 1e-10, check.attributes = FALSE))) {
      stopf("CHECK table differs: %s", f)
    }
  }
}

run_analysis <- function() {
  specs <- list(
    GSE98793 = c("age", "sex", "medication", "batch"),
    GSE19738 = c("age", "sex", "medication", "batch"),
    GSE39653 = c("age", "sex", "medication", "batch"),
    GSE101521 = c("age", "sex", "medication", "RIN", "PMI", "batch")
  )
  fits <- list()
  for (dataset in names(specs)) fits[[dataset]] <- fit_microarray(dataset, specs[[dataset]])
  fits$GSE251778 <- fit_gse251778()
  order <- c("GSE98793", "GSE19738", "GSE39653", "GSE251778", "GSE101521")
  effect_dir <- file.path(output_root, "data", "derived", "cohort_effects")
  table_dir <- file.path(output_root, "results", "tables")
  audit <- summary <- list()
  for (dataset in order) {
    z <- fits[[dataset]]
    write_tsv(z$table, file.path(effect_dir, paste0(dataset, "_effects.tsv.gz")), TRUE)
    deg <- z$table[z$table$FDR < 0.05, , drop = FALSE]
    write_tsv(deg, file.path(table_dir, paste0(dataset, "_DEG.tsv.gz")), TRUE)
    summary[[dataset]] <- data.frame(
      dataset = dataset, analysis_role = if (dataset == "GSE101521") "exploratory" else
        if (dataset == "GSE39653") "PBMC_validation" else "primary_cohort",
      Ncase = z$ncase, Ncontrol = z$ncontrol, genes_tested = nrow(z$table),
      DEG_FDR_lt_0.05 = nrow(deg),
      DEG_MDD_higher = sum(deg$beta > 0), DEG_MDD_lower = sum(deg$beta < 0),
      sex_opposite_direction = if (dataset == "GSE251778") sum(!z$table$same_direction) else NA_integer_,
      heterogeneity_flagged = if (dataset == "GSE251778") sum(z$table$heterogeneity_flag) else NA_integer_,
      DEG_heterogeneity_flagged = if (dataset == "GSE251778") sum(deg$heterogeneity_flag) else NA_integer_,
      design = paste(z$design, collapse = ","), eBayes_robust = FALSE,
      stringsAsFactors = FALSE)
    a <- z$audit
    a$dataset <- dataset
    if (!"stratum" %in% names(a)) a$stratum <- "all"
    contract_reason <- if (dataset == "GSE251778") {
      paste0("Task5: standardize within sex strata using Hedges effect and unmoderated ",
             "sampling SE/residual scale, then combine; do not meta-analyze raw log2 beta ",
             "or moderated SE")
    } else {
      paste0("Task5: use preregistered Hedges standardized effect with unmoderated sampling ",
             "SE/residual scale; do not meta-analyze raw log2 beta or moderated SE")
    }
    a <- rbind(a, data.frame(dataset = dataset, stratum = "all",
                             covariate = "meta_effect_contract", status = "deferred_to_Task5",
                             reason = contract_reason, stringsAsFactors = FALSE))
    audit[[dataset]] <- a[, c("dataset", "stratum", "covariate", "status", "reason")]
  }
  write_tsv(fits$GSE251778$male,
            file.path(effect_dir, "GSE251778_Male_effects.tsv.gz"), TRUE)
  write_tsv(fits$GSE251778$female,
            file.path(effect_dir, "GSE251778_Female_effects.tsv.gz"), TRUE)
  write_tsv(do.call(rbind, summary), file.path(table_dir, "cohort_deg_summary.tsv"))
  write_tsv(do.call(rbind, audit), file.path(table_dir, "cohort_model_audit.tsv"))
  session <- c(sprintf("seed\t%d", 20260724),
               sprintf("R\t%s", R.version.string),
               sprintf("limma\t%s", as.character(packageVersion("limma"))),
               sprintf("edgeR\t%s", as.character(packageVersion("edgeR"))),
               "eBayes_robust\tFALSE",
               "SE_alias\tSE equals SE_moderated; Task5 must use SE_unmoderated",
               "single_cohort_p_FDR\tmoderated_t with df_total; BH over all tested genes",
               paste0("GSE251778_combined_p_FDR\tmoderated-SE inverse-variance fixed-effect ",
                      "normal z; BH over all tested genes"))
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(session, file.path(table_dir, "cohort_model_session.tsv"))
}

if (simulate_only) {
  simulate_pipeline()
  quit(save = "no", status = 0)
}

if (check_only) {
  simulate_pipeline()
  tmp <- tempfile("cohort-model-check-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  script <- normalizePath("analysis/06_fit_cohort_models.R")
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    c("--vanilla", shQuote(script)),
                    env = c(paste0("COHORT_MODEL_OUTPUT_ROOT=", shQuote(tmp)),
                            paste0("R_LIBS_USER=", shQuote(project_lib))))
  if (status != 0L) stopf("Temporary cohort model rebuild failed")
  compare_tables(root, tmp)
  message("READ_ONLY_CHECK_OK all cohort tables equal temporary rebuild")
  quit(save = "no", status = 0)
}

run_analysis()
message("COHORT_MODELS_OK")
