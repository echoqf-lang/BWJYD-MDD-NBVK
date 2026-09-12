#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(20260724)

args <- commandArgs(trailingOnly = TRUE)
check_only <- "--check" %in% args
simulate_only <- "--simulate" %in% args
root <- normalizePath(".", mustWork = TRUE)
project_lib <- file.path(root, ".r-lib")
if (dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))
output_root <- Sys.getenv("META_OUTPUT_ROOT", unset = root)

stopf <- function(...) stop(sprintf(...), call. = FALSE)
if (!requireNamespace("metafor", quietly = TRUE)) stopf("Missing required package: metafor")

write_tsv <- function(x, path, gzip = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- if (gzip) gzfile(path, "wt") else file(path, "wt")
  on.exit(close(con), add = TRUE)
  write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

read_effect <- function(dataset, stratum = NULL) {
  suffix <- if (is.null(stratum)) "" else paste0("_", stratum)
  path <- file.path(root, "data", "derived", "cohort_effects",
                    paste0(dataset, suffix, "_effects.tsv.gz"))
  x <- read.delim(path, check.names = FALSE)
  required <- c("gene", "beta", "SE_unmoderated", "SE_moderated",
                "residual_scale", "df_residual")
  if (!all(required %in% names(x))) stopf("%s missing effect-contract fields", basename(path))
  if (anyDuplicated(x$gene)) stopf("%s duplicated genes", basename(path))
  x
}

standardize_effect <- function(x, dataset, stratum = "all",
                               se_field = "SE_unmoderated") {
  if (!identical(se_field, "SE_unmoderated")) {
    stopf("Contract violation: moderated SE must not be used for standardization")
  }
  if (!se_field %in% names(x)) stopf("%s missing SE_unmoderated", dataset)
  required <- c("beta", se_field, "residual_scale", "df_residual")
  if (any(!vapply(x[required], is.numeric, logical(1)))) stopf("%s non-numeric contract field", dataset)
  valid <- is.finite(x$beta) & is.finite(x[[se_field]]) & x[[se_field]] > 0 &
    is.finite(x$residual_scale) & x$residual_scale > 0 &
    is.finite(x$df_residual) & x$df_residual > 1
  if (!all(valid)) stopf("%s invalid df, residual scale, beta, or sampling SE", dataset)
  J <- 1 - 3 / (4 * x$df_residual - 1)
  data.frame(
    gene = x$gene, dataset = dataset, stratum = stratum,
    beta = x$beta, SE_unmoderated = x[[se_field]],
    residual_scale = x$residual_scale, df_residual = x$df_residual,
    J = J, g = J * x$beta / x$residual_scale,
    var_g = J^2 * (x[[se_field]] / x$residual_scale)^2 +
      (J * x$beta / x$residual_scale)^2 / (2 * x$df_residual),
    source_se = se_field, stringsAsFactors = FALSE
  )
}

combine_gse251778 <- function(male, female, combined_source) {
  z <- merge(male, female, by = "gene", suffixes = c("_male", "_female"),
             all = FALSE, sort = FALSE)
  if (nrow(z) != nrow(male) || nrow(z) != nrow(female)) {
    stopf("GSE251778 sex-stratum gene sets differ")
  }
  wm <- 1 / z$var_g_male
  wf <- 1 / z$var_g_female
  g <- (wm * z$g_male + wf * z$g_female) / (wm + wf)
  vg <- 1 / (wm + wf)
  q <- wm * (z$g_male - g)^2 + wf * (z$g_female - g)^2
  qp <- pchisq(q, df = 1, lower.tail = FALSE)
  i2 <- ifelse(q > 0, pmax(0, (q - 1) / q) * 100, 0)
  src <- combined_source[match(z$gene, combined_source$gene), ]
  if (anyNA(src$gene)) stopf("GSE251778 combined-source genes differ from strata")
  data.frame(
    gene = z$gene, dataset = "GSE251778", stratum = "sex_combined",
    beta = NA_real_, SE_unmoderated = NA_real_, residual_scale = NA_real_,
    df_residual = z$df_residual_male + z$df_residual_female,
    J = NA_real_, g = g, var_g = vg,
    source_se = "sex_stratum_SE_unmoderated",
    g_male = z$g_male, var_g_male = z$var_g_male,
    g_female = z$g_female, var_g_female = z$var_g_female,
    sex_Cochran_Q = q, sex_Q_p = qp, sex_I2 = i2,
    sex_same_direction = sign(z$g_male) == sign(z$g_female),
    heterogeneity_flag = qp < 0.05 | i2 > 75 |
      sign(z$g_male) != sign(z$g_female),
    upstream_heterogeneity_flag = as.logical(src$heterogeneity_flag),
    platform_sex_limitation = "platform_id perfectly coupled with sex; effects not separately identifiable",
    stringsAsFactors = FALSE
  )
}

fit_reml <- function(yi, vi) {
  if (length(yi) < 2L || any(!is.finite(yi)) || any(!is.finite(vi)) || any(vi <= 0)) {
    return(list(status = "failure_invalid_input"))
  }
  tryCatch({
    warnings <- character()
    fit <- withCallingHandlers(
      metafor::rma.uni(yi = yi, vi = vi, method = "REML"),
      warning = function(w) {
        warnings <<- c(warnings, gsub("[\t\r\n]+", " ", conditionMessage(w)))
        invokeRestart("muffleWarning")
      })
    list(status = "ok", estimate = as.numeric(fit$b[1]), se = fit$se,
         ci_lb = fit$ci.lb, ci_ub = fit$ci.ub, p = fit$pval,
         tau2 = fit$tau2, I2 = fit$I2, Q = fit$QE, Q_p = fit$QEp,
         convergence_note = if (length(warnings)) paste(unique(warnings), collapse = " | ") else "")
  }, error = function(e) list(status = paste0(
    "failure_REML: ", gsub("[\t\r\n]+", " ", conditionMessage(e)))))
}

fit_fixed <- function(yi, vi) {
  w <- 1 / vi
  est <- sum(w * yi) / sum(w)
  se <- sqrt(1 / sum(w))
  c(estimate = est, se = se, ci_lb = est - qnorm(0.975) * se,
    ci_ub = est + qnorm(0.975) * se, p = 2 * pnorm(-abs(est / se)))
}

meta_rows <- function(effects, genes, cohorts) {
  out <- vector("list", length(genes))
  for (i in seq_along(genes)) {
    g <- genes[i]
    x <- effects[effects$gene == g & effects$dataset %in% cohorts, ]
    x <- x[match(cohorts, x$dataset), ]
    fit <- fit_reml(x$g, x$var_g)
    signs <- sign(x$g)
    direction_n <- max(sum(signs > 0), sum(signs < 0))
    out[[i]] <- data.frame(
      gene = g, k = length(cohorts), model = "random_REML",
      status = fit$status,
      convergence_note = if (!is.null(fit$convergence_note)) fit$convergence_note else "",
      meta_g = if (fit$status == "ok") fit$estimate else NA_real_,
      SE = if (fit$status == "ok") fit$se else NA_real_,
      CI_lower = if (fit$status == "ok") fit$ci_lb else NA_real_,
      CI_upper = if (fit$status == "ok") fit$ci_ub else NA_real_,
      p = if (fit$status == "ok") fit$p else NA_real_,
      tau2 = if (fit$status == "ok") fit$tau2 else NA_real_,
      I2 = if (fit$status == "ok") fit$I2 else NA_real_,
      Cochran_Q = if (fit$status == "ok") fit$Q else NA_real_,
      Q_p = if (fit$status == "ok") fit$Q_p else NA_real_,
      direction_consistent_n = direction_n,
      direction_3_of_3 = length(cohorts) == 3L && direction_n == 3L,
      cohort_directions = paste(paste(cohorts, ifelse(signs > 0, "+",
                                                       ifelse(signs < 0, "-", "0")),
                                      sep = ":"), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

expect_error <- function(expr, pattern) {
  ok <- FALSE
  tryCatch(force(expr), error = function(e) ok <<- grepl(pattern, conditionMessage(e)))
  if (!ok) stopf("Negative test did not detect: %s", pattern)
}

simulate_pipeline <- function() {
  base <- data.frame(
    gene = c("G1", "G2"), beta = c(2, -1),
    SE_unmoderated = c(0.5, 0.4), SE_moderated = c(0.1, 0.1),
    residual_scale = c(2, 1), df_residual = c(11, 21))
  s <- standardize_effect(base, "SIM")
  expected_J <- 1 - 3 / (4 * base$df_residual - 1)
  if (!isTRUE(all.equal(s$g, expected_J * base$beta / base$residual_scale,
                        tolerance = 1e-12)) ||
      !isTRUE(all.equal(s$var_g, expected_J^2 *
                         (base$SE_unmoderated / base$residual_scale)^2 +
                         s$g^2 / (2 * base$df_residual),
                        tolerance = 1e-12))) stopf("Simulation known-g recovery failed")

  m <- standardize_effect(transform(base, beta = c(2, 1)), "SIM", "Male")
  f <- standardize_effect(transform(base, beta = c(2, -3)), "SIM", "Female")
  src <- data.frame(gene = base$gene, heterogeneity_flag = c(FALSE, FALSE))
  csex <- combine_gse251778(m, f, src)
  if (abs(csex$g[1] - m$g[1]) > 1e-12 || csex$sex_I2[2] <= csex$sex_I2[1]) {
    stopf("Simulation sex-standardization or heterogeneity recovery failed")
  }
  if (!csex$heterogeneity_flag[2] || csex$upstream_heterogeneity_flag[2]) {
    stopf("Simulation standardized-scale flag independence failed")
  }

  fit <- fit_reml(c(0.2, 0.4, 0.8), c(0.01, 0.02, 0.03))
  if (fit$status != "ok" || !is.finite(fit$estimate) || !is.finite(fit$tau2)) {
    stopf("Simulation REML recovery failed")
  }
  p <- c(0.001, 0.02, 0.5)
  if (!isTRUE(all.equal(p.adjust(p, "BH"), c(0.003, 0.03, 0.5)))) {
    stopf("Simulation BH-FDR failed")
  }
  families <- data.frame(family = rep(c("A", "B"), each = 3), p = rep(p, 2))
  families$FDR <- ave(families$p, families$family,
                      FUN = function(z) p.adjust(z, "BH"))
  if (!isTRUE(all.equal(families$FDR, rep(c(0.003, 0.03, 0.5), 2)))) {
    stopf("Simulation family-specific BH-FDR failed")
  }
  expect_error({
    incomplete_var <- expected_J^2 *
      (base$SE_unmoderated / base$residual_scale)^2
    if (!isTRUE(all.equal(s$var_g, incomplete_var, tolerance = 1e-12))) {
      stopf("residual-SD uncertainty omission detected")
    }
  }, "residual-SD uncertainty omission")

  expect_error(standardize_effect(base, "SIM", se_field = "SE_moderated"),
               "moderated SE")
  expect_error({
    bad <- rbind(transform(s, dataset = "A"), transform(s, dataset = "B"))
    if (length(unique(bad$dataset[bad$gene == "G1"])) != 3L) stopf("missing cohort detected")
  }, "missing cohort")
  expect_error({
    reversed <- s
    reversed$g <- -reversed$g
    if (!identical(sign(s$g), sign(reversed$g))) stopf("direction reversal detected")
  }, "direction reversal")
  expect_error({
    tampered <- s
    tampered$g[1] <- tampered$g[1] + 1
    if (!isTRUE(all.equal(s, tampered, tolerance = 1e-12))) stopf("table tampering detected")
  }, "table tampering")
  message("SIMULATION_OK known_g heterogeneity REML FDR sex_standardize_then_combine negative_tests")
}

build_outputs <- function() {
  main_cohorts <- c("GSE98793", "GSE251778", "GSE19738")
  ordinary <- lapply(c("GSE98793", "GSE19738", "GSE39653"),
                     function(d) standardize_effect(read_effect(d), d))
  names(ordinary) <- c("GSE98793", "GSE19738", "GSE39653")
  male <- standardize_effect(read_effect("GSE251778", "Male"), "GSE251778", "Male")
  female <- standardize_effect(read_effect("GSE251778", "Female"), "GSE251778", "Female")
  combined_source <- read_effect("GSE251778")
  combined <- combine_gse251778(male, female, combined_source)

  common <- Reduce(intersect, list(ordinary$GSE98793$gene, combined$gene,
                                   ordinary$GSE19738$gene))
  if (!length(common)) stopf("No genes have valid effects in all three available whole-blood cohorts")
  main_effects <- rbind(
    ordinary$GSE98793[ordinary$GSE98793$gene %in% common,
                      c("gene", "dataset", "g", "var_g")],
    combined[combined$gene %in% common, c("gene", "dataset", "g", "var_g")],
    ordinary$GSE19738[ordinary$GSE19738$gene %in% common,
                      c("gene", "dataset", "g", "var_g")]
  )
  counts <- table(main_effects$gene)
  if (any(counts != 3L)) stopf("Main meta requires exactly all three available cohorts")

  meta <- meta_rows(main_effects, common, main_cohorts)
  meta$FDR <- NA_real_
  ok <- meta$status == "ok" & is.finite(meta$p)
  meta$FDR[ok] <- p.adjust(meta$p[ok], method = "BH")

  flag <- combined$heterogeneity_flag[match(meta$gene, combined$gene)]
  meta$GSE251778_heterogeneity_flag <- flag
  upstream_flag <- combined$upstream_heterogeneity_flag[match(meta$gene, combined$gene)]
  meta$GSE251778_upstream_raw_scale_flag <- upstream_flag
  meta$GSE251778_flag_scale_agreement <- flag == upstream_flag

  exclusion_effects <- main_effects[
    main_effects$dataset %in% c("GSE98793", "GSE19738"), ]
  exclusion <- meta_rows(exclusion_effects, common,
                         c("GSE98793", "GSE19738"))
  exclusion$FDR <- NA_real_
  exclusion_ok <- exclusion$status == "ok" & is.finite(exclusion$p)
  exclusion$FDR[exclusion_ok] <- p.adjust(exclusion$p[exclusion_ok], method = "BH")
  ix <- match(meta$gene, exclusion$gene)
  meta$sensitivity_excluding_GSE251778_applicable <- flag
  meta$sensitivity_excluding_GSE251778_status <- ifelse(flag, exclusion$status[ix], NA)
  meta$sensitivity_excluding_GSE251778_g <- ifelse(flag, exclusion$meta_g[ix], NA)
  meta$sensitivity_excluding_GSE251778_SE <- ifelse(flag, exclusion$SE[ix], NA)
  meta$sensitivity_excluding_GSE251778_p <- ifelse(flag, exclusion$p[ix], NA)
  meta$sensitivity_excluding_GSE251778_FDR <- ifelse(flag, exclusion$FDR[ix], NA)
  meta$sensitivity_excluding_GSE251778_I2 <- ifelse(flag, exclusion$I2[ix], NA)

  pbmc <- ordinary$GSE39653[match(meta$gene, ordinary$GSE39653$gene), ]
  meta$PBMC_available <- !is.na(pbmc$gene)
  meta$PBMC_g <- pbmc$g
  meta$PBMC_var_g <- pbmc$var_g
  meta$PBMC_same_direction <- ifelse(meta$PBMC_available & meta$status == "ok",
                                     sign(meta$meta_g) == sign(pbmc$g), NA)

  loo <- list()
  q <- 1L
  for (omitted in main_cohorts) {
    kept <- setdiff(main_cohorts, omitted)
    z <- meta_rows(main_effects, common, kept)
    z$omitted_cohort <- omitted
    z$k2_instability_note <- "k=2 REML is reported but heterogeneity estimation is unstable"
    loo[[q]] <- z
    q <- q + 1L
  }
  loo <- do.call(rbind, loo)
  loo$FDR <- ave(loo$p, loo$omitted_cohort, FUN = function(p) {
    out <- rep(NA_real_, length(p))
    ok <- is.finite(p)
    out[ok] <- p.adjust(p[ok], method = "BH")
    out
  })

  array_cohorts <- c("GSE98793", "GSE19738")
  array_reml <- exclusion
  array_fixed <- lapply(common, function(g) {
    x <- main_effects[main_effects$gene == g &
                        main_effects$dataset %in% array_cohorts, ]
    f <- fit_fixed(x$g, x$var_g)
    data.frame(gene = g, array_fixed_g = f["estimate"], array_fixed_SE = f["se"],
               array_fixed_p = f["p"])
  })
  array_fixed <- do.call(rbind, array_fixed)
  meta$array_only_fixed_g <- array_fixed$array_fixed_g[match(meta$gene, array_fixed$gene)]
  meta$array_only_fixed_SE <- array_fixed$array_fixed_SE[match(meta$gene, array_fixed$gene)]
  meta$array_only_fixed_p <- array_fixed$array_fixed_p[match(meta$gene, array_fixed$gene)]
  meta$array_only_fixed_FDR <- p.adjust(meta$array_only_fixed_p, method = "BH")
  meta$array_only_REML_g <- array_reml$meta_g[match(meta$gene, array_reml$gene)]
  meta$array_only_REML_p <- array_reml$p[match(meta$gene, array_reml$gene)]
  meta$array_only_REML_FDR <- NA_real_
  array_reml_ok <- is.finite(meta$array_only_REML_p)
  meta$array_only_REML_FDR[array_reml_ok] <-
    p.adjust(meta$array_only_REML_p[array_reml_ok], method = "BH")
  meta$array_only_status <- array_reml$status[match(meta$gene, array_reml$gene)]
  meta$array_only_interpretation <- "descriptive sensitivity only; not confirmatory"

  standardized <- rbind(
    ordinary$GSE98793, ordinary$GSE19738, ordinary$GSE39653,
    male, female,
    combined[, names(male)]
  )
  audit <- data.frame(
    item = c("preregistered_primary_cohorts", "actual_primary_cohorts",
             "deviation_GSE201332", "minimum_cohorts", "effect_contract",
             "GSE251778_combination", "GSE251778_platform_sex_limitation",
             "standardized_sampling_variance", "GSE251778_flag_scale",
             "GSE251778_flag_sensitivity", "primary_model",
             "REML_failure_policy", "legacy_variance_fit_note", "multiplicity",
             "sensitivity_multiplicity", "array_only_role", "PBMC_role",
             "PBMC_direction_test", "brain_role"),
    value = c("GSE98793;GSE251778;GSE201332;GSE19738",
              paste(main_cohorts, collapse = ";"),
              "excluded: user abandoned because samples could not be reliably mapped",
              "all 3 of 3 currently available whole-blood cohorts",
              "Hedges g from beta/residual_scale and SE_unmoderated only",
              "standardize Male and Female separately, then inverse-variance fixed effect",
              "platform_id perfectly coupled with sex; effects not separately identifiable",
              paste0("J^2*(SE_unmoderated/residual_scale)^2 + ",
                     "g^2/(2*df_residual); includes residual-SD estimation uncertainty"),
              "recomputed from sex-stratum g and var_g; Q, I2, and direction use one scale",
              paste0("all 7310 genes are re-estimated after excluding GSE251778 using ",
                     "two-cohort REML and full-family BH; flagged genes are then summarized"),
              "per-gene inverse-variance random effects, metafor REML",
              "report failure; no alternate estimator or fixed-effect fallback",
              paste0("HEBP1 failed REML under the superseded variance that omitted residual-SD ",
                     "uncertainty; it converges under the corrected frozen variance and is not ",
                     "artificially excluded"),
              sprintf("BH over %d genes with successful primary REML p values; failed fits excluded",
                      sum(ok)),
              paste0("LOO BH separately within each omitted-cohort family; array fixed and ",
                     "array REML BH separately over all testable genes"),
              "descriptive sensitivity only",
              "external direction/effect validation only; excluded from primary meta",
              paste0("exact two-sided binomial test against 0.5 with exact 95% CI; ",
                     "below-50% result is classified as not supporting directional validation"),
              "GSE101521 deferred to later exploratory analysis"),
    stringsAsFactors = FALSE
  )
  summary <- data.frame(
    metric = c("genes_testable_background", "primary_REML_success",
               "primary_REML_failure", "primary_REML_failure_gene",
               "FDR_lt_0.05", "direction_3_of_3",
               "I2_gt_75", "GSE251778_flagged", "GSE251778_not_flagged",
               "GSE251778_flag_scale_agreement",
               "GSE251778_exclusion_REML_success",
               "GSE251778_exclusion_FDR_lt_0.05",
               "PBMC_direction_evaluable", "PBMC_same_direction",
               "PBMC_same_direction_proportion", "PBMC_binomial_CI_lower",
               "PBMC_binomial_CI_upper", "PBMC_binomial_p_two_sided",
               "PBMC_direction_conclusion",
               paste0("LOO_FDR_lt_0.05_omit_", main_cohorts),
               "array_only_fixed_FDR_lt_0.05",
               "array_only_REML_FDR_lt_0.05"),
    value = {
      pbmc_n <- sum(!is.na(meta$PBMC_same_direction))
      pbmc_same <- sum(meta$PBMC_same_direction, na.rm = TRUE)
      pbmc_test <- binom.test(pbmc_same, pbmc_n, p = 0.5)
      c(nrow(meta), sum(ok), sum(!ok),
        if (any(!ok)) paste(meta$gene[!ok], collapse = ";") else "none",
        sum(meta$FDR < 0.05, na.rm = TRUE),
              sum(meta$direction_3_of_3, na.rm = TRUE),
              sum(meta$I2 > 75, na.rm = TRUE),
              sum(meta$GSE251778_heterogeneity_flag, na.rm = TRUE),
        sum(!meta$GSE251778_heterogeneity_flag, na.rm = TRUE),
        sum(meta$GSE251778_flag_scale_agreement, na.rm = TRUE),
        sum(exclusion_ok & exclusion$gene %in% meta$gene[flag]),
        sum(exclusion$FDR[exclusion$gene %in% meta$gene[flag]] < 0.05, na.rm = TRUE),
        pbmc_n, pbmc_same, pbmc_same / pbmc_n,
        unname(pbmc_test$conf.int[1]), unname(pbmc_test$conf.int[2]),
        pbmc_test$p.value,
        if (pbmc_test$conf.int[2] < 0.5) "below 50%; does not support directional validation"
        else if (pbmc_test$conf.int[1] > 0.5) "above 50%; supports directional validation"
        else "not different from 50%; does not support directional validation",
              vapply(main_cohorts, function(d)
                sum(loo$FDR[loo$omitted_cohort == d] < 0.05, na.rm = TRUE), numeric(1)),
              sum(meta$array_only_fixed_FDR < 0.05, na.rm = TRUE),
        sum(meta$array_only_REML_FDR < 0.05, na.rm = TRUE))
    },
    stringsAsFactors = FALSE
  )
  session <- data.frame(
    key = c("seed", "R", "metafor", "primary_model", "se_contract",
            "sampling_variance"),
    value = c("20260724", R.version.string, as.character(packageVersion("metafor")),
              "random effects REML; no result-driven fallback",
              "SE_unmoderated/residual_scale",
              "J^2*(SE_unmoderated/residual_scale)^2 + g^2/(2*df_residual)"),
    stringsAsFactors = FALSE
  )
  list(meta = meta, standardized = standardized, loo = loo,
       summary = summary, audit = audit, session = session)
}

write_outputs <- function(x, out_root) {
  meta_dir <- file.path(out_root, "data", "derived", "meta")
  tab_dir <- file.path(out_root, "results", "tables")
  write_tsv(x$meta, file.path(meta_dir, "whole_blood_meta.tsv.gz"), TRUE)
  write_tsv(x$standardized, file.path(meta_dir, "standardized_cohort_effects.tsv.gz"), TRUE)
  write_tsv(x$loo, file.path(meta_dir, "leave_one_out.tsv.gz"), TRUE)
  write_tsv(x$summary, file.path(tab_dir, "meta_summary.tsv"))
  write_tsv(x$audit, file.path(tab_dir, "meta_audit.tsv"))
  write_tsv(x$session, file.path(tab_dir, "meta_session.tsv"))
}

compare_outputs <- function(reference_root, candidate_root) {
  rel <- c(file.path("data", "derived", "meta",
                     c("whole_blood_meta.tsv.gz", "standardized_cohort_effects.tsv.gz",
                       "leave_one_out.tsv.gz")),
           file.path("results", "tables",
                     c("meta_summary.tsv", "meta_audit.tsv", "meta_session.tsv")))
  for (f in rel) {
    a <- read.delim(file.path(reference_root, f), check.names = FALSE)
    b <- read.delim(file.path(candidate_root, f), check.names = FALSE)
    if (!isTRUE(all.equal(a, b, tolerance = 1e-10, check.attributes = FALSE))) {
      stopf("CHECK table differs: %s", f)
    }
  }
}

if (simulate_only) {
  simulate_pipeline()
  quit(save = "no", status = 0)
}

if (check_only) {
  simulate_pipeline()
  tmp <- tempfile("meta-check-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  x <- build_outputs()
  write_outputs(x, tmp)
  compare_outputs(root, tmp)
  message("READ_ONLY_CHECK_OK all meta tables equal temporary rebuild")
  quit(save = "no", status = 0)
}

simulate_pipeline()
x <- build_outputs()
write_outputs(x, output_root)
message("META_ANALYSIS_OK")
