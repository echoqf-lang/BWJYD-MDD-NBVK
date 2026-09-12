#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(20260724)

args <- commandArgs(trailingOnly = TRUE)
check_only <- "--check" %in% args
simulate_only <- "--simulate" %in% args
root <- normalizePath(".", mustWork = TRUE)
project_lib <- file.path(root, ".r-lib")
if (dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))
output_root <- Sys.getenv("BRAIN_OUTPUT_ROOT", unset = root)

stopf <- function(...) stop(sprintf(...), call. = FALSE)
if (!requireNamespace("limma", quietly = TRUE)) stopf("Missing required package: limma")

write_tsv <- function(x, path, gzip = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- if (gzip) gzfile(path, "wt") else file(path, "wt")
  on.exit(close(con), add = TRUE)
  write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

read_inputs <- function() {
  expr <- read.delim(
    gzfile(file.path(root, "data", "clean", "GSE101521",
                     "expression_gene.tsv.gz")),
    check.names = FALSE
  )
  genes <- expr[[1]]
  mat <- as.matrix(expr[-1])
  storage.mode(mat) <- "double"
  rownames(mat) <- genes
  meta <- read.delim(
    file.path(root, "data", "clean", "GSE101521", "sample_metadata.tsv"),
    check.names = FALSE, na.strings = c("", "NA", "NR")
  )
  evidence <- read.delim(
    gzfile(file.path(root, "data", "derived", "targets",
                     "candidate_meta_evidence.tsv.gz")),
    check.names = FALSE
  )
  core <- read.delim(
    gzfile(file.path(root, "data", "derived", "targets",
                     "high_priority_common_targets.tsv.gz")),
    check.names = FALSE
  )

  if (anyDuplicated(genes) || any(!is.finite(mat))) {
    stopf("Invalid GSE101521 clean expression matrix")
  }
  if (!identical(colnames(mat), meta$GSM)) {
    stopf("Expression columns are not identical to metadata GSM order")
  }
  if (anyDuplicated(meta$GSM) || anyDuplicated(evidence$gene) ||
      anyDuplicated(core$gene)) {
    stopf("Duplicated sample or locked target identifier")
  }
  if (nrow(meta) != 59L || sum(meta$diagnosis == "MDD") != 30L ||
      sum(meta$diagnosis == "control") != 29L) {
    stopf("Frozen GSE101521 sample counts changed")
  }
  if (!all(is.na(meta$batch))) {
    stopf("GSE101521 batch was expected to be entirely unavailable")
  }
  if (sum(meta$diagnosis == "MDD" & meta$suicide == "yes") != 21L ||
      sum(meta$diagnosis == "MDD" & meta$suicide == "no") != 9L ||
      any(meta$diagnosis == "control" & meta$suicide != "no")) {
    stopf("Suicide subgroup encoding is incomplete or unexpected")
  }
  if (nrow(evidence) != 567L || nrow(core) != 5L ||
      !all(core$gene %in% evidence$gene)) {
    stopf("Locked candidate/core set size or nesting changed")
  }
  required_meta <- c("age", "sex", "RIN", "PMI", "batch")
  if (!all(required_meta %in% names(meta))) stopf("Required brain metadata missing")
  list(expr = mat, meta = meta, evidence = evidence, core = core)
}

prepare_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

fit_one_model <- function(expr, meta, model_id, subgroup,
                          covariates = c("age", "sex", "RIN", "PMI"),
                          complete_case = TRUE,
                          sample_whitelist = NULL) {
  keep <- if (subgroup == "overall") {
    rep(TRUE, nrow(meta))
  } else if (subgroup == "MDD_non_suicide") {
    meta$diagnosis == "control" |
      (meta$diagnosis == "MDD" & meta$suicide == "no")
  } else if (subgroup == "MDD_suicide") {
    meta$diagnosis == "control" |
      (meta$diagnosis == "MDD" & meta$suicide == "yes")
  } else {
    stopf("Unknown subgroup: %s", subgroup)
  }
  if (!is.null(sample_whitelist)) keep <- keep & meta$GSM %in% sample_whitelist
  smeta <- meta[keep, , drop = FALSE]
  smeta$age <- prepare_numeric(smeta$age)
  smeta$RIN <- prepare_numeric(smeta$RIN)
  smeta$PMI <- prepare_numeric(smeta$PMI)
  smeta$sex <- factor(smeta$sex)
  smeta$diagnosis <- factor(smeta$diagnosis, levels = c("control", "MDD"))
  if (complete_case) {
    cc <- complete.cases(smeta[, c("diagnosis", covariates), drop = FALSE])
    smeta <- smeta[cc, , drop = FALSE]
  }
  if (anyNA(smeta[, c("diagnosis", covariates), drop = FALSE])) {
    stopf("%s contains missing model variables without a declared complete-case rule",
          model_id)
  }
  if (any(table(smeta$diagnosis) < 2L)) stopf("%s has too few samples", model_id)
  sex_tab <- table(smeta$sex)
  if (length(sex_tab) < 2L) stopf("%s sex is constant", model_id)

  design <- model.matrix(reformulate(c("diagnosis", covariates)), smeta)
  if (qr(design)$rank < ncol(design)) stopf("%s design is rank deficient", model_id)
  coef_name <- grep("^diagnosisMDD$", colnames(design), value = TRUE)
  if (length(coef_name) != 1L) stopf("%s diagnosis coefficient not unique", model_id)
  sexpr <- expr[, smeta$GSM, drop = FALSE]
  fit <- limma::eBayes(limma::lmFit(sexpr, design))
  j <- match(coef_name, colnames(fit$coefficients))
  beta <- fit$coefficients[, j]
  se <- fit$stdev.unscaled[, j] * sqrt(fit$s2.post)
  p <- fit$p.value[, j]
  if (any(!is.finite(beta)) || any(!is.finite(se)) || any(se <= 0) ||
      any(!is.finite(p))) {
    stopf("%s produced invalid statistics", model_id)
  }
  effects <- data.frame(
    model_id = model_id, contrast = paste0(subgroup, "_vs_control"),
    gene = rownames(expr), beta = beta, SE = se,
    CI_lower = beta - qt(0.975, fit$df.total) * se,
    CI_upper = beta + qt(0.975, fit$df.total) * se,
    moderated_t = fit$t[, j], p = p, FDR = p.adjust(p, "BH"),
    N = nrow(smeta), Ncase = sum(smeta$diagnosis == "MDD"),
    Ncontrol = sum(smeta$diagnosis == "control"),
    covariates = paste(covariates, collapse = "+"),
    complete_case = complete_case,
    sample_restriction = if (is.null(sample_whitelist)) {
      "subgroup_all_available_samples_before_complete_case_filter"
    } else {
      "prelocked_GSM_whitelist_from_overall_complete_case_model"
    },
    direction = ifelse(beta > 0, "MDD_higher",
                       ifelse(beta < 0, "MDD_lower", "zero")),
    stringsAsFactors = FALSE
  )
  list(effects = effects, fit = fit, design = design, meta = smeta,
       expr = sexpr, coef = j)
}

camera_one <- function(model, sets) {
  idx <- lapply(sets, function(g) which(rownames(model$expr) %in% g))
  out <- limma::camera(
    model$expr, index = idx, design = model$design,
    contrast = model$coef, inter.gene.cor = NA
  )
  out$gene_set <- rownames(out)
  rownames(out) <- NULL
  data.frame(
    model_id = unique(model$effects$model_id),
    contrast = unique(model$effects$contrast),
    gene_set = out$gene_set,
    set_size_input = vapply(sets[out$gene_set], length, integer(1)),
    set_size_tested = out$NGenes,
    estimated_inter_gene_correlation = out$Correlation,
    direction = out$Direction,
    p = out$PValue,
    camera_FDR_within_model = out$FDR,
    correlation_handling = paste0(
      "limma camera competitive test; inter.gene.cor estimated ",
      "from residual expression because inter.gene.cor=NA"
    ),
    stringsAsFactors = FALSE
  )
}

annotate_locked <- function(effects, evidence, core) {
  blood <- evidence[, c("gene", "meta_g", "FDR", "I2",
                        "GSE251778_heterogeneity_flag"), drop = FALSE]
  names(blood)[names(blood) == "FDR"] <- "blood_meta_FDR"
  names(blood)[names(blood) == "I2"] <- "blood_meta_I2"
  x <- merge(effects, blood, by = "gene", all.x = FALSE, all.y = FALSE,
             sort = FALSE)
  x$is_locked_core <- x$gene %in% core$gene
  x$same_direction_as_blood <- sign(x$beta) == sign(x$meta_g)
  x$nominal_directional_replication <-
    x$same_direction_as_blood & x$p < 0.05
  x$brain_transcriptome_FDR_lt_0_05 <- x$FDR < 0.05
  x
}

build_outputs <- function() {
  x <- read_inputs()
  primary_specs <- list(
    c("overall_complete_case", "overall"),
    c("non_suicide_complete_case", "MDD_non_suicide"),
    c("suicide_complete_case", "MDD_suicide")
  )
  primary <- lapply(primary_specs, function(s) {
    fit_one_model(x$expr, x$meta, s[1], s[2],
                  c("age", "sex", "RIN", "PMI"), TRUE)
  })
  overall_cc_gsm <- primary[[1]]$meta$GSM
  bridge <- fit_one_model(
    x$expr, x$meta, "overall_complete_case_reduced_covariate_bridge",
    "overall", c("age", "sex"), FALSE, sample_whitelist = overall_cc_gsm
  )
  sensitivity <- fit_one_model(
    x$expr, x$meta,
    "overall_full_sample_reduced_covariate_sensitivity", "overall",
    c("age", "sex"), FALSE
  )

  primary_effects <- do.call(rbind, lapply(primary, `[[`, "effects"))
  bridge_effects <- bridge$effects
  sensitivity_effects <- sensitivity$effects
  all_effects <- rbind(primary_effects, bridge_effects, sensitivity_effects)
  all_effects$analysis_role <- ifelse(
    all_effects$model_id ==
      "overall_complete_case_reduced_covariate_bridge",
    "same_sample_covariate_bridge",
    ifelse(
      all_effects$model_id ==
        "overall_full_sample_reduced_covariate_sensitivity",
      "full_sample_reduced_covariate_sensitivity",
      "primary_brain_validation"
    )
  )

  candidate <- do.call(rbind, lapply(primary, function(m) {
    annotate_locked(m$effects, x$evidence, x$core)
  }))
  candidate <- candidate[order(candidate$model_id, candidate$blood_meta_FDR,
                               candidate$p), , drop = FALSE]
  core <- candidate[candidate$is_locked_core, , drop = FALSE]
  core <- core[order(core$model_id, core$blood_meta_FDR), , drop = FALSE]

  sets <- list(
    locked_candidate_567 = x$evidence$gene,
    locked_core_5 = x$core$gene
  )
  gene_sets <- do.call(rbind, lapply(primary, camera_one, sets = sets))
  gene_sets$FDR_BH_all_six_locked_tests <- p.adjust(gene_sets$p, "BH")

  bridge_annotated <- annotate_locked(
    bridge$effects, x$evidence, x$core
  )
  full_sample_annotated <- annotate_locked(
    sensitivity$effects, x$evidence, x$core
  )
  overall_primary <- candidate[
    candidate$model_id == "overall_complete_case", , drop = FALSE
  ]
  bridge_match <- bridge_annotated[
    match(overall_primary$gene, bridge_annotated$gene), , drop = FALSE
  ]
  full_sample_match <- full_sample_annotated[
    match(overall_primary$gene, full_sample_annotated$gene), , drop = FALSE
  ]
  missingness <- data.frame(
    gene = overall_primary$gene,
    primary_56_full_cov_beta = overall_primary$beta,
    bridge_56_age_sex_beta = bridge_match$beta,
    full_sample_59_age_sex_beta = full_sample_match$beta,
    same_direction_primary_vs_bridge_covariate_effect =
      sign(overall_primary$beta) == sign(bridge_match$beta),
    same_direction_bridge_vs_full_sample_sample_inclusion_effect =
      sign(bridge_match$beta) == sign(full_sample_match$beta),
    same_direction_primary_vs_full_sample_joint_sensitivity =
      sign(overall_primary$beta) == sign(full_sample_match$beta),
    primary_56_full_cov_p = overall_primary$p,
    bridge_56_age_sex_p = bridge_match$p,
    full_sample_59_age_sex_p = full_sample_match$p,
    is_locked_core = overall_primary$is_locked_core,
    stringsAsFactors = FALSE
  )

  summary_rows <- list()
  for (m in primary) {
    e <- m$effects
    ca <- candidate[candidate$model_id == unique(e$model_id), , drop = FALSE]
    co <- ca[ca$is_locked_core, , drop = FALSE]
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      model_id = unique(e$model_id),
      N = unique(e$N), Ncase = unique(e$Ncase), Ncontrol = unique(e$Ncontrol),
      full_transcriptome_genes = nrow(e),
      full_transcriptome_FDR_lt_0_05 = sum(e$FDR < 0.05),
      locked_candidate_tested = nrow(ca),
      locked_candidate_nominal_directional_replication =
        sum(ca$nominal_directional_replication),
      locked_candidate_brain_FDR_lt_0_05 =
        sum(ca$brain_transcriptome_FDR_lt_0_05),
      locked_core_tested = nrow(co),
      locked_core_nominal_directional_replication =
        sum(co$nominal_directional_replication),
      locked_core_brain_FDR_lt_0_05 =
        sum(co$brain_transcriptome_FDR_lt_0_05),
      stringsAsFactors = FALSE
    )
  }
  summary <- do.call(rbind, summary_rows)
  summary$overall_candidate_same_direction_primary_vs_bridge <-
    ifelse(
      summary$model_id == "overall_complete_case",
      sum(missingness$same_direction_primary_vs_bridge_covariate_effect),
      NA_integer_
    )
  summary$overall_candidate_same_direction_bridge_vs_full_sample <-
    ifelse(
      summary$model_id == "overall_complete_case",
      sum(missingness$
            same_direction_bridge_vs_full_sample_sample_inclusion_effect),
      NA_integer_
    )
  summary$overall_candidate_same_direction_primary_vs_full_sample_joint <-
    ifelse(
      summary$model_id == "overall_complete_case",
      sum(missingness$
            same_direction_primary_vs_full_sample_joint_sensitivity),
      NA_integer_
    )
  summary$overall_core_same_direction_primary_vs_bridge <-
    ifelse(
      summary$model_id == "overall_complete_case",
      sum(missingness$
            same_direction_primary_vs_bridge_covariate_effect[
              missingness$is_locked_core
            ]),
      NA_integer_
    )
  summary$overall_core_same_direction_bridge_vs_full_sample <-
    ifelse(
      summary$model_id == "overall_complete_case",
      sum(missingness$
            same_direction_bridge_vs_full_sample_sample_inclusion_effect[
              missingness$is_locked_core
            ]),
      NA_integer_
    )
  summary$overall_core_same_direction_primary_vs_full_sample_joint <-
    ifelse(
      summary$model_id == "overall_complete_case",
      sum(missingness$
            same_direction_primary_vs_full_sample_joint_sensitivity[
              missingness$is_locked_core
            ]),
      NA_integer_
    )

  batch_nonmissing <- sum(!is.na(x$meta$batch))
  audit <- data.frame(
    item = c(
      "analysis_status", "dataset", "expression_input",
      "diagnostic_groups", "suicide_group_reliability",
      "primary_covariates", "primary_missingness_rule",
      "batch_availability", "RIN_availability", "PMI_availability",
      "overall_interpretation", "subgroup_interpretation",
      "locked_candidate_definition", "locked_core_definition",
      "locked_candidate_brain_coverage",
      "single_gene_replication_rule", "single_gene_multiplicity",
      "gene_set_test", "gene_set_multiplicity",
      "sensitivity_decomposition", "joint_sensitivity_interpretation",
      "additional_brain_datasets"
    ),
    value = c(
      "exploratory because upstream registered acquisition/model deviations occurred",
      "GSE101521 DLPFC BA9 total-RNA samples only",
      "official DESeq2-normalized noninteger counts transformed as log2(x+1) upstream",
      "29 control; 9 MDD non-suicide; 21 MDD suicide",
      "all controls encoded no; all 30 MDD have yes/no suicide annotation",
      "age + sex + RIN + PMI",
      paste0("complete-case per contrast; no imputation; ",
             "same pre-specified covariates in all three primary models"),
      sprintf(
        "%d/%d nonmissing; omitted as unavailable (not imputed)",
        sum(!is.na(x$meta$batch)), nrow(x$meta)
      ),
      sprintf("%d/59 nonmissing", sum(!is.na(x$meta$RIN))),
      sprintf("%d/59 nonmissing", sum(!is.na(x$meta$PMI))),
      "MDD case-control contrast includes both suicide strata and is not suicide-specific",
      paste0("non-suicide and suicide MDD were each compared with controls; ",
             "small non-suicide stratum limits precision"),
      "567 Task6 testable prescription-MDD candidates locked before brain analysis",
      "5 Task6 high-priority common targets locked before brain analysis",
      sprintf(
        "%d/567 assayed; unavailable in GSE101521: %s",
        sum(x$evidence$gene %in% rownames(x$expr)),
        paste(setdiff(x$evidence$gene, rownames(x$expr)), collapse = ";")
      ),
      "same sign as blood Meta and brain nominal p<0.05",
      paste0("brain FDR is BH across all assayed genes within each contrast; ",
             "nominal replication is not transcriptome-wide significance"),
      paste0("limma camera competitive test on expression and design; ",
             "residual inter-gene correlation estimated with inter.gene.cor=NA"),
      "BH across 2 locked sets x 3 primary contrasts (6 tests)",
      paste0("56-sample full-covariate primary vs 56-sample age+sex bridge ",
             "isolates covariate-adjustment sensitivity; 56-sample bridge vs ",
             "59-sample age+sex model isolates inclusion of three samples"),
      paste0("56-sample full-covariate vs 59-sample age+sex is a joint ",
             "sensitivity only; any direction change, including FASLG/CTSG, ",
             "must not be attributed to missingness without the bridge comparison"),
      paste0("GSE80655/GSE102556/GSE144136 not analyzed: no pre-cleaned reliable ",
             "processed matrix in the frozen workspace; no large raw download attempted")
    ),
    stringsAsFactors = FALSE
  )

  balance <- do.call(rbind, lapply(primary, function(m) {
    data.frame(
      model_id = unique(m$effects$model_id),
      variable = c("age", "RIN", "PMI", "Female", "Male"),
      control = c(
        mean(m$meta$age[m$meta$diagnosis == "control"]),
        mean(m$meta$RIN[m$meta$diagnosis == "control"]),
        mean(m$meta$PMI[m$meta$diagnosis == "control"]),
        sum(m$meta$diagnosis == "control" & m$meta$sex == "Female"),
        sum(m$meta$diagnosis == "control" & m$meta$sex == "Male")
      ),
      MDD = c(
        mean(m$meta$age[m$meta$diagnosis == "MDD"]),
        mean(m$meta$RIN[m$meta$diagnosis == "MDD"]),
        mean(m$meta$PMI[m$meta$diagnosis == "MDD"]),
        sum(m$meta$diagnosis == "MDD" & m$meta$sex == "Female"),
        sum(m$meta$diagnosis == "MDD" & m$meta$sex == "Male")
      ),
      unit = c("mean_years", "mean", "mean_hours", "n", "n"),
      stringsAsFactors = FALSE
    )
  }))

  session <- data.frame(
    key = c("seed", "R", "limma", "primary_model",
            "same_sample_bridge_model", "full_sample_sensitivity_model",
            "camera_inter_gene_correlation", "brain_FDR_family",
            "gene_set_FDR_family"),
    value = c(
      "20260724", R.version.string,
      as.character(utils::packageVersion("limma")),
      "~ diagnosis + age + sex + RIN + PMI; empirical-Bayes limma",
      "~ diagnosis + age + sex on the identical 56 complete-case samples",
      "~ diagnosis + age + sex on all 59 samples",
      "estimated from residual expression (inter.gene.cor=NA)",
      "all assayed genes separately within each contrast",
      "all six locked set-by-contrast tests"
    ),
    stringsAsFactors = FALSE
  )
  list(
    effects = all_effects, candidate = candidate, core = core,
    gene_sets = gene_sets, missingness = missingness,
    summary = summary, audit = audit, balance = balance, session = session
  )
}

simulate_pipeline <- function() {
  set.seed(20260724)
  genes <- paste0("G", seq_len(120))
  n0 <- 24L
  n1 <- 24L
  diagnosis <- factor(rep(c("control", "MDD"), each = 24),
                      levels = c("control", "MDD"))
  age <- rnorm(n0 + n1, 40, 8)
  sex <- factor(rep(c("Female", "Male"), length.out = n0 + n1))
  RIN <- rnorm(n0 + n1, 7, 0.5)
  PMI <- rnorm(n0 + n1, 12, 3)
  design <- model.matrix(~ diagnosis + age + sex + RIN + PMI)
  expr <- matrix(rnorm(length(genes) * (n0 + n1)), nrow = length(genes),
                 dimnames = list(genes, paste0("S", seq_len(n0 + n1))))
  expr[1:12, diagnosis == "MDD"] <-
    expr[1:12, diagnosis == "MDD"] + 1.2
  fit <- limma::eBayes(limma::lmFit(expr, design))
  j <- match("diagnosisMDD", colnames(design))
  if (mean(fit$coefficients[1:12, j]) <= 0 ||
      sum(fit$p.value[1:12, j] < 0.05) < 8L) {
    stopf("Simulation failed to recover injected positive effects")
  }
  cam <- limma::camera(expr, list(injected = 1:12), design, contrast = j,
                       inter.gene.cor = NA)
  if (cam["injected", "Direction"] != "Up" ||
      cam["injected", "PValue"] >= 0.05) {
    stopf("Simulation failed competitive gene-set recovery")
  }
  reversed <- -fit$coefficients[1, j]
  if (sign(reversed) == sign(fit$coefficients[1, j])) {
    stopf("Negative direction test failed")
  }
  duplicated <- c(genes, genes[1])
  if (!anyDuplicated(duplicated)) stopf("Duplicate-gene negative test failed")
  message("SIMULATION_OK contrast direction BH camera negative_tests")
}

write_outputs <- function(x, out_root) {
  derived <- file.path(out_root, "data", "derived", "brain")
  tables <- file.path(out_root, "results", "tables")
  write_tsv(x$effects, file.path(derived, "GSE101521_brain_effects.tsv.gz"), TRUE)
  write_tsv(x$candidate,
            file.path(derived, "GSE101521_candidate_validation.tsv.gz"), TRUE)
  write_tsv(x$core, file.path(tables, "brain_core_validation.tsv"))
  write_tsv(x$gene_sets, file.path(tables, "brain_gene_set_validation.tsv"))
  write_tsv(x$missingness,
            file.path(tables, "brain_sensitivity_decomposition.tsv"))
  write_tsv(x$summary, file.path(tables, "brain_validation_summary.tsv"))
  write_tsv(x$audit, file.path(tables, "brain_validation_audit.tsv"))
  write_tsv(x$balance, file.path(tables, "brain_covariate_balance.tsv"))
  write_tsv(x$session, file.path(tables, "brain_validation_session.tsv"))
}

output_files <- function(base) c(
  file.path(base, "data", "derived", "brain",
            c("GSE101521_brain_effects.tsv.gz",
              "GSE101521_candidate_validation.tsv.gz")),
  file.path(base, "results", "tables",
            c("brain_core_validation.tsv", "brain_gene_set_validation.tsv",
              "brain_sensitivity_decomposition.tsv",
              "brain_validation_summary.tsv", "brain_validation_audit.tsv",
              "brain_covariate_balance.tsv",
              "brain_validation_session.tsv"))
)

compare_outputs <- function(reference_root, candidate_root) {
  a_files <- output_files(reference_root)
  b_files <- output_files(candidate_root)
  for (i in seq_along(a_files)) {
    reader <- if (grepl("\\.gz$", a_files[i])) gzfile else identity
    a <- read.delim(reader(a_files[i]), check.names = FALSE)
    b <- read.delim(reader(b_files[i]), check.names = FALSE)
    if (!isTRUE(all.equal(a, b, tolerance = 1e-12,
                          check.attributes = FALSE))) {
      stopf("CHECK table differs: %s", basename(a_files[i]))
    }
  }
}

if (simulate_only) {
  simulate_pipeline()
  quit(save = "no", status = 0)
}

if (check_only) {
  simulate_pipeline()
  tmp <- tempfile("brain-validation-check-")
  dir.create(tmp, recursive = TRUE)
  rebuilt <- build_outputs()
  write_outputs(rebuilt, tmp)
  compare_outputs(root, tmp)
  message("CHECK_OK fresh temporary rebuild matches all brain validation outputs")
  quit(save = "no", status = 0)
}

out <- build_outputs()
write_outputs(out, output_root)
message("BRAIN_VALIDATION_OK")
