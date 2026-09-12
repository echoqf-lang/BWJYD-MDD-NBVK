#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
check_only <- "--check" %in% args
self_test <- any(c("--self-test", "--negative-test") %in% args)
root <- normalizePath(getwd(), mustWork = TRUE)
raw_dir <- file.path(root, "data", "raw_geo")
derived_dir <- file.path(root, "data", "derived")
results_dir <- file.path(root, "results", "tables")
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

dataset_spec <- data.frame(
  dataset_id = c("GSE98793", "GSE251778", "GSE201332", "GSE19738",
                 "GSE39653", "GSE101521", "GSE144136",
                 "GSE38206", "GSE80655", "GSE102556"),
  analysis_role = c(rep("primary_whole_blood", 4), "pbmc_validation",
                    "primary_brain_validation", "cell_type_localization",
                    rep("exploratory_only", 3)),
  expected_design = c(
    "128 MDD + 64 control whole-blood donors",
    "80 MDD + 89 control whole-blood donors; RNA-seq subseries only",
    "20 untreated MDD + 20 control whole-blood donors",
    "33 MDD + 34 control donors; basal only; 65 paired LPS samples excluded",
    "21 MDD + 24 control PBMC donors; 8 bipolar samples excluded",
    "59 total-RNA DLPFC donors (21 MDD-suicide, 9 MDD-nonsuicide, 29 control); 27 miRNA libraries excluded",
    "17 MDD-suicide + 17 control male DLPFC donors for snRNA-seq localization",
    "18 donors (9 MDE + 9 control) at baseline and week 8; exploratory baseline only",
    "281 brain-region samples from control/MDD/BD/SCZ; exploratory only",
    "Human and mouse, multiple brain regions; exploratory human records only"
  ),
  stringsAsFactors = FALSE
)

parse_soft_records <- function(path) {
  x <- readLines(path, warn = FALSE)
  starts <- grep("^\\^SAMPLE = ", x)
  if (!length(starts)) stop("No SAMPLE records: ", path)
  ends <- c(starts[-1L] - 1L, length(x))
  lapply(seq_along(starts), function(i) {
    block <- x[starts[i]:ends[i]]
    keys <- sub(" = .*", "", block)
    vals <- sub("^[^=]+ = ?", "", block)
    split(vals, keys)
  })
}

first_value <- function(rec, key) {
  value <- rec[[key]]
  if (is.null(value) || !length(value) || !nzchar(value[1L])) "NR" else value[1L]
}

characteristics <- function(rec) {
  values <- rec[["!Sample_characteristics_ch1"]]
  if (is.null(values)) return(setNames(character(), character()))
  pos <- regexpr(":", values, fixed = TRUE)
  keep <- pos > 0L
  keys <- tolower(trimws(substr(values[keep], 1L, pos[keep] - 1L)))
  vals <- trimws(substr(values[keep], pos[keep] + 1L, nchar(values[keep])))
  setNames(vals, keys)
}

char_value <- function(ch, keys) {
  hit <- keys[keys %in% names(ch)]
  if (!length(hit)) "NR" else unname(ch[hit[1L]])
}

clean_nr <- function(x) {
  ifelse(is.na(x) | !nzchar(trimws(x)) | toupper(trimws(x)) %in% c("NA", "N/A", "UNKNOWN"), "NR", trimws(x))
}

derive_donor <- function(gse, title, description) {
  if (gse == "GSE19738") {
    z <- sub(".*-(CON[0-9]+|MDD[0-9]+)$", "\\1", title, ignore.case = TRUE)
  } else if (gse == "GSE38206") {
    z <- sub("^PBMC_([PC][0-9]+)-.*$", "\\1", title, ignore.case = TRUE)
  } else if (gse == "GSE80655") {
    z <- sub("_.*$", "", title)
  } else if (gse %in% c("GSE102556", "GSE144136")) {
    z <- sub(":.*$", "", title)
  } else if (gse == "GSE101521") {
    z <- description
  } else {
    z <- title
  }
  if (!nzchar(z) || identical(z, title) && gse %in% c("GSE19738", "GSE38206")) "NR" else paste(gse, z, sep = ":")
}

derive_one <- function(gse, rec) {
  ch <- characteristics(rec)
  gsm <- sub("^SAMPLE = ", "", first_value(rec, "^SAMPLE"))
  title <- first_value(rec, "!Sample_title")
  description <- first_value(rec, "!Sample_description")
  platform <- first_value(rec, "!Sample_platform_id")
  source <- first_value(rec, "!Sample_source_name_ch1")
  organism <- first_value(rec, "!Sample_organism_ch1")
  taxid <- first_value(rec, "!Sample_taxid_ch1")
  molecule <- first_value(rec, "!Sample_molecule_ch1")
  library_strategy <- first_value(rec, "!Sample_library_strategy")
  group_raw <- char_value(ch, c("subject group", "group", "subject status", "disease",
                                "diagnosis", "status", "clinical diagnosis", "phenotype"))
  diagnosis <- "other"
  if (grepl("control|healthy|cntl|con\\b|ctrl", group_raw, ignore.case = TRUE)) diagnosis <- "control"
  if (grepl("major depress|mdd|mde", group_raw, ignore.case = TRUE)) diagnosis <- "MDD"
  if (grepl("bipolar", group_raw, ignore.case = TRUE)) diagnosis <- "BD"
  if (grepl("schizo", group_raw, ignore.case = TRUE)) diagnosis <- "SCZ"
  if (grepl("stress", group_raw, ignore.case = TRUE)) diagnosis <- "stress"

  tissue <- char_value(ch, "tissue")
  if (tissue == "NR") tissue <- source
  condition <- "NR"
  timepoint <- "NR"
  stimulation <- "NR"
  condition_source <- "NR"
  timepoint_source <- "NR"
  stimulation_source <- "NR"
  if (gse %in% c("GSE98793", "GSE251778", "GSE201332", "GSE39653")) {
    condition <- "basal"
    timepoint <- "baseline"
    stimulation <- "none"
    condition_source <- "Series_overall_design + frozen cross-sectional eligibility"
    timepoint_source <- condition_source
    stimulation_source <- condition_source
  }
  if (gse == "GSE19738" && grepl("LPS", char_value(ch, "treatment"), ignore.case = TRUE)) {
    condition <- "LPS"; timepoint <- "5-6h"; stimulation <- "LPS_5-6h"
    condition_source <- "Sample_characteristics_ch1:treatment"
    timepoint_source <- "Sample_source_name_ch1"
    stimulation_source <- "Sample_characteristics_ch1:treatment"
  } else if (gse == "GSE19738") {
    condition <- "basal"; timepoint <- "baseline"; stimulation <- "none"
    condition_source <- "Sample_characteristics_ch1:treatment"
    timepoint_source <- "Sample_source_name_ch1"
    stimulation_source <- "Sample_characteristics_ch1:treatment"
  }
  if (gse == "GSE38206") {
    collection <- char_value(ch, "sample collection")
    if (grepl("8 week", collection, ignore.case = TRUE)) {
      condition <- if (diagnosis == "MDD") "clinical_remission" else "matched_followup"
      timepoint <- "week_8"
    } else {
      condition <- if (diagnosis == "MDD") "severe_episode" else "matched_baseline"
      timepoint <- "baseline"
    }
    condition_source <- "Sample_characteristics_ch1:sample collection + diagnosis-derived label"
    timepoint_source <- "Sample_characteristics_ch1:sample collection"
  }
  brain_region <- char_value(ch, "brain region")
  if (brain_region == "NR" && grepl("brain|cortex|accumbens|insula|subiculum", tissue, ignore.case = TRUE)) {
    brain_region <- tissue
  }
  if (gse == "GSE144136") brain_region <- "DLPFC_BA9"

  include <- TRUE
  reasons <- character()
  if (gse == "GSE19738" && stimulation != "none") {
    include <- FALSE; reasons <- c(reasons, "LPS_stimulated")
  }
  if (gse == "GSE39653" && diagnosis == "BD") {
    include <- FALSE; reasons <- c(reasons, "bipolar_disorder")
  }
  if (gse == "GSE101521" && platform != "GPL16791") {
    include <- FALSE; reasons <- c(reasons, "miRNA_library_not_total_RNA")
  }
  if (gse %in% c("GSE38206", "GSE80655", "GSE102556")) {
    include <- FALSE; reasons <- c(reasons, "exploratory_not_primary")
  }
  if (gse == "GSE38206" && timepoint != "baseline") reasons <- c(reasons, "post_treatment_week_8")
  if (gse == "GSE80655" && diagnosis %in% c("BD", "SCZ")) reasons <- c(reasons, "non_MDD_psychiatric_diagnosis")
  if (gse == "GSE102556" && !grepl("postmortem brain", source, ignore.case = TRUE)) reasons <- c(reasons, "non_human_arm")

  suicide <- "NR"
  if (gse == "GSE101521") {
    suicide <- if (grepl("suicides \\(MDD-S\\)", group_raw, ignore.case = TRUE)) "yes" else
      if (diagnosis %in% c("MDD", "control")) "no" else "NR"
  } else if (gse == "GSE144136") {
    suicide <- if (diagnosis == "MDD") "yes" else if (diagnosis == "control") "NR" else "NR"
  } else if (gse == "GSE102556") {
    suicide <- if (grepl("suicide", char_value(ch, "cause of death"), ignore.case = TRUE)) "yes" else
      if (char_value(ch, "cause of death") != "NR") "no" else "NR"
  }

  medication <- char_value(ch, c("medication", "treatment"))
  medication_source <- if (medication == "NR") "NR" else "Sample_characteristics_ch1"
  if (gse == "GSE201332") {
    medication <- "never_treated"; medication_source <- "Series_summary"
  }
  if (gse == "GSE101521") {
    medication <- "medication_free"; medication_source <- "Series_summary"
  }

  age <- char_value(ch, c("age", "age (yrs)", "age at death"))
  sex <- char_value(ch, c("sex", "gender"))
  gad <- if (gse == "GSE98793") char_value(ch, "anxiety") else "NR"
  rin <- char_value(ch, "rin")
  pmi <- char_value(ch, c("pmi", "post-mortem interval"))
  batch <- char_value(ch, "batch")
  donor_id <- derive_donor(gse, title, description)
  donor_source <- if (donor_id == "NR") "NR" else
    if (gse == "GSE101521") "Sample_description" else "Sample_title"
  if (gse == "GSE101521") {
    # The series design states that GPL16791 contains one total-RNA library per
    # donor (59 donors). GEO does not provide a uniquely verifiable crosswalk
    # from the 27 miRNA libraries to those donors; do not match on RIN or other
    # biological covariates.
    if (platform == "GPL16791") {
      donor_id <- paste(gse, gsm, sep = ":")
      donor_source <- "Series_overall_design + one total-RNA GSM per donor"
    } else {
      donor_id <- "NR"
      donor_source <- "NR"
    }
  }
  donor_mapping_evidence <- if (donor_id == "NR") "unresolved" else
    if (gse == "GSE101521") "verified_series_design_unique_total_RNA_GSM" else
      "verified_explicit_sample_identifier"

  values <- c(age = age, sex = sex, medication = medication, GAD = gad,
              suicide = suicide, brain_region = brain_region, RIN = rin,
              PMI = pmi, batch = batch)
  values <- clean_nr(values)
  sources <- ifelse(values == "NR", "NR", "Sample_characteristics_ch1")
  names(sources) <- paste0("source_", names(values))
  sources["source_medication"] <- medication_source
  if (gse == "GSE101521") sources["source_suicide"] <- "Sample_characteristics_ch1:diagnosis"
  if (gse == "GSE144136") sources["source_suicide"] <- "Series_overall_design"
  if (gse == "GSE102556") sources["source_suicide"] <- "Sample_characteristics_ch1:cause of death"
  if (brain_region != "NR") sources["source_brain_region"] <- if (gse == "GSE144136") "Series_summary" else "Sample_characteristics_ch1"

  data.frame(
    dataset_id = gse, GSM = gsm,
    donor_id = donor_id,
    diagnosis = diagnosis, include_primary = include,
    exclusion_reason = if (length(reasons)) paste(unique(reasons), collapse = ";") else "included",
    tissue = clean_nr(tissue), condition = condition, timepoint = timepoint,
    stimulation = stimulation, age = values["age"], sex = values["sex"],
    medication = values["medication"], GAD = values["GAD"],
    suicide = values["suicide"], brain_region = values["brain_region"],
    RIN = values["RIN"], PMI = values["PMI"], batch = values["batch"],
    platform_id = platform, organism = organism, taxid = taxid,
    molecule = molecule, library_strategy = library_strategy, sample_title = title,
    metadata_source = sprintf("data/raw_geo/%s_samples_brief.soft", gse),
    source_donor_id = donor_source,
    donor_mapping_evidence = donor_mapping_evidence,
    source_diagnosis = "Sample_characteristics_ch1",
    source_tissue = if (char_value(ch, "tissue") == "NR") "Sample_source_name_ch1" else "Sample_characteristics_ch1",
    source_condition = condition_source,
    source_timepoint = timepoint_source,
    source_stimulation = stimulation_source,
    as.list(sources), check.names = FALSE, stringsAsFactors = FALSE
  )
}

sample_manifest <- do.call(rbind, lapply(dataset_spec$dataset_id, function(gse) {
  path <- file.path(raw_dir, sprintf("%s_samples_brief.soft", gse))
  records <- parse_soft_records(path)
  gsm_caret <- vapply(records, first_value, character(1), "^SAMPLE")
  gsm_field <- vapply(records, first_value, character(1), "!Sample_geo_accession")
  if (!identical(gsm_caret, gsm_field)) stop(gse, ": ^SAMPLE and !Sample_geo_accession differ.")
  series_lines <- readLines(file.path(raw_dir, sprintf("%s_series_brief.soft", gse)), warn = FALSE)
  series_gsm <- sub("^!Series_sample_id = ", "", grep("^!Series_sample_id = ", series_lines, value = TRUE))
  if (!setequal(gsm_caret, series_gsm) || anyDuplicated(series_gsm)) {
    stop(gse, ": series and sample SOFT GSM sets differ.")
  }
  do.call(rbind, lapply(records, function(rec) derive_one(gse, rec)))
}))
rownames(sample_manifest) <- NULL

repeat_type <- rep("none", nrow(sample_manifest))
for (gse in unique(sample_manifest$dataset_id)) {
  idx_g <- which(sample_manifest$dataset_id == gse & sample_manifest$donor_id != "NR")
  for (donor in unique(sample_manifest$donor_id[idx_g])) {
    idx <- idx_g[sample_manifest$donor_id[idx_g] == donor]
    if (length(idx) < 2L) next
    if (length(unique(sample_manifest$brain_region[idx])) > 1L) {
      repeat_type[idx] <- "cross_brain_region"
    } else if (length(unique(paste(sample_manifest$condition[idx], sample_manifest$timepoint[idx],
                                   sample_manifest$stimulation[idx]))) > 1L) {
      repeat_type[idx] <- "longitudinal_or_stimulation"
    } else {
      repeat_type[idx] <- "technical_or_modality"
    }
  }
}
sample_manifest$repeat_type <- repeat_type

count_donors <- function(x) {
  x <- unique(x[x != "NR"])
  if (length(x)) length(x) else NA_integer_
}

repeat_excess <- function(d, type) {
  z <- d[d$repeat_type == type & d$donor_id != "NR", , drop = FALSE]
  if (!nrow(z)) return(0L)
  sum(vapply(split(z$GSM, z$donor_id), function(x) max(length(x) - 1L, 0L), integer(1)))
}

audit_rows <- lapply(seq_len(nrow(dataset_spec)), function(i) {
  gse <- dataset_spec$dataset_id[i]
  d <- sample_manifest[sample_manifest$dataset_id == gse, ]
  inc <- d[d$include_primary, ]
  data.frame(
    dataset_id = gse,
    gsm_count = nrow(d),
    independent_donor_count = count_donors(d$donor_id),
    technical_or_modality_repeat_samples = sum(d$repeat_type == "technical_or_modality"),
    longitudinal_or_stimulation_repeat_samples = sum(d$repeat_type == "longitudinal_or_stimulation"),
    cross_brain_region_repeat_samples = sum(d$repeat_type == "cross_brain_region"),
    technical_or_modality_excess_observations = repeat_excess(d, "technical_or_modality"),
    longitudinal_or_stimulation_excess_observations = repeat_excess(d, "longitudinal_or_stimulation"),
    cross_brain_region_excess_observations = repeat_excess(d, "cross_brain_region"),
    included_case_gsm = sum(inc$diagnosis == "MDD"),
    included_control_gsm = sum(inc$diagnosis == "control"),
    included_unique_case_donors = count_donors(inc$donor_id[inc$diagnosis == "MDD"]),
    included_unique_control_donors = count_donors(inc$donor_id[inc$diagnosis == "control"]),
    excluded_gsm = sum(!d$include_primary),
    stringsAsFactors = FALSE
  )
})
audit <- do.call(rbind, audit_rows)

dataset_manifest <- merge(dataset_spec, audit, by = "dataset_id", sort = FALSE)
dataset_manifest <- dataset_manifest[match(dataset_spec$dataset_id, dataset_manifest$dataset_id), ]
dataset_manifest$expected_case <- c(128, 80, 20, 33, 21, 30, 17, 0, 0, 0)
dataset_manifest$expected_control <- c(64, 89, 20, 34, 24, 29, 17, 0, 0, 0)
dataset_manifest$expected_excluded <- c(0, 0, 0, 65, 8, 27, 0, 36, 281, 341)

frozen_design_problems <- function(samples, dataset_filter = dataset_spec$dataset_id) {
  problems <- character()
  included <- if (is.logical(samples$include_primary)) samples$include_primary else samples$include_primary == "TRUE"
  for (i in seq_len(nrow(dataset_manifest))) {
    gse <- dataset_manifest$dataset_id[i]
    if (!gse %in% dataset_filter) next
    d <- samples[samples$dataset_id == gse, , drop = FALSE]
    inc <- included[samples$dataset_id == gse]
    observed <- c(sum(inc & d$diagnosis == "MDD"), sum(inc & d$diagnosis == "control"), sum(!inc))
    expected <- unlist(dataset_manifest[i, c("expected_case", "expected_control", "expected_excluded")],
                       use.names = FALSE)
    if (!identical(as.numeric(observed), as.numeric(expected))) {
      problems <- c(problems, paste0(gse, "_case_control_exclusion_count_mismatch"))
    }
  }
  count <- function(gse, expr) sum(samples$dataset_id == gse & expr)
  if (count("GSE19738", samples$condition == "basal") != 67L ||
      count("GSE19738", samples$condition == "LPS") != 65L ||
      count("GSE19738", grepl("LPS_stimulated", samples$exclusion_reason, fixed = TRUE)) != 65L) {
    problems <- c(problems, "GSE19738_condition_or_reason_mismatch")
  }
  if (count("GSE39653", grepl("bipolar_disorder", samples$exclusion_reason, fixed = TRUE)) != 8L) {
    problems <- c(problems, "GSE39653_BD_reason_mismatch")
  }
  if (count("GSE101521", samples$platform_id == "GPL16791") != 59L ||
      count("GSE101521", samples$platform_id == "GPL15520") != 27L ||
      count("GSE101521", grepl("miRNA_library_not_total_RNA", samples$exclusion_reason, fixed = TRUE)) != 27L) {
    problems <- c(problems, "GSE101521_platform_or_reason_mismatch")
  }
  diagnosis_levels <- c("MDD", "control")
  condition_levels <- c("severe_episode", "clinical_remission",
                        "matched_baseline", "matched_followup")
  expected_38206 <- matrix(
    c(9L, 9L, 0L, 0L,
      0L, 0L, 9L, 9L),
    nrow = 2L, byrow = TRUE,
    dimnames = list(diagnosis_levels, condition_levels)
  )
  d38206 <- samples[samples$dataset_id == "GSE38206", , drop = FALSE]
  actual_38206 <- table(
    factor(d38206$diagnosis, levels = diagnosis_levels),
    factor(d38206$condition, levels = condition_levels)
  )
  known_cells <- d38206$diagnosis %in% diagnosis_levels & d38206$condition %in% condition_levels
  if (!all(known_cells) ||
      !identical(as.integer(actual_38206), as.integer(expected_38206))) {
    problems <- c(problems, "GSE38206_diagnosis_specific_condition_mismatch")
  }
  platform_expected <- list(
    GSE98793 = c(GPL570 = 192L), GSE251778 = c(GPL20301 = 64L, GPL24676 = 105L),
    GSE201332 = c(GPL32193 = 40L), GSE19738 = c(GPL6848 = 132L),
    GSE39653 = c(GPL10558 = 53L), GSE101521 = c(GPL15520 = 27L, GPL16791 = 59L),
    GSE144136 = c(GPL20301 = 34L), GSE38206 = c(GPL13607 = 36L),
    GSE80655 = c(GPL11154 = 281L), GSE102556 = c(GPL11154 = 263L, GPL13112 = 78L)
  )
  for (gse in names(platform_expected)) {
    observed <- sort(table(samples$platform_id[samples$dataset_id == gse]))
    expected <- sort(platform_expected[[gse]])
    if (!identical(names(observed), names(expected)) ||
        !identical(as.integer(observed), as.integer(expected))) {
      problems <- c(problems, paste0(gse, "_platform_count_mismatch"))
    }
  }
  species_102556 <- table(samples$organism[samples$dataset_id == "GSE102556"])
  if (!identical(as.integer(species_102556[c("Homo sapiens", "Mus musculus")]), c(263L, 78L))) {
    problems <- c(problems, "GSE102556_species_count_mismatch")
  }
  unique(problems)
}

dataset_manifest$design_verification <- vapply(dataset_manifest$dataset_id, function(gse) {
  if (length(frozen_design_problems(sample_manifest, gse))) "DIFF_REQUIRES_GSM_REVIEW" else "matches_frozen_design"
}, character(1))
dataset_manifest$cohort_reuse_assessment <- "not_assessed"
dataset_manifest$cohort_reuse_assessment[dataset_manifest$dataset_id == "GSE251778"] <-
  "GSE38206 comparison: unverified independence; different institute/city, study era, platform, sample codes; no direct subject-link evidence"
dataset_manifest$cohort_reuse_assessment[dataset_manifest$dataset_id == "GSE38206"] <-
  "GSE251778 comparison: unverified independence; different institute/city, study era, platform, sample codes; no direct subject-link evidence"

validate_content <- function(samples, datasets) {
  problems <- character()
  if (anyDuplicated(samples$GSM)) problems <- c(problems, "duplicate_GSM")
  if (any(samples$include_primary & duplicated(paste(samples$dataset_id, samples$donor_id)) &
          samples$donor_id != "NR")) problems <- c(problems, "same_donor_included_more_than_once")
  if (any(samples$include_primary & grepl("LPS", samples$stimulation, ignore.case = TRUE))) problems <- c(problems, "LPS_included")
  if (any(samples$include_primary & samples$diagnosis == "BD")) problems <- c(problems, "BD_included")
  if (any(samples$include_primary & samples$dataset_id == "GSE38206" & samples$timepoint != "baseline")) {
    problems <- c(problems, "post_treatment_included")
  }
  if (any(is.na(samples$exclusion_reason) | !nzchar(samples$exclusion_reason) |
          samples$exclusion_reason == "NR")) problems <- c(problems, "missing_exclusion_reason")
  if (any(samples$include_primary & samples$exclusion_reason != "included") ||
      any(!samples$include_primary & samples$exclusion_reason == "included")) {
    problems <- c(problems, "inconsistent_exclusion_reason")
  }
  for (field in c("condition", "timepoint", "stimulation")) {
    source_field <- paste0("source_", field)
    if (any((samples[[field]] == "NR") != (samples[[source_field]] == "NR"))) {
      problems <- c(problems, paste0(field, "_source_inconsistent"))
    }
  }
  if (any(datasets$design_verification != "matches_frozen_design")) problems <- c(problems, "design_count_mismatch")
  problems <- c(problems, frozen_design_problems(samples))
  unique(problems)
}

read_tsv_character <- function(path) {
  if (!file.exists(path)) stop("Missing formal output: ", path)
  x <- read.delim(path, check.names = FALSE, colClasses = "character",
                  na.strings = character(), quote = "", comment.char = "")
  x[is.na(x)] <- "NR"
  x
}

as_disk_character <- function(x) {
  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  write.table(x, tmp, sep = "\t", quote = FALSE, row.names = FALSE, na = "NR")
  read_tsv_character(tmp)
}

geo_bucket <- function(gse) sprintf("GSE%dnnn", as.integer(sub("^GSE", "", gse)) %/% 1000L)
expected_download_rows <- do.call(rbind, lapply(dataset_spec$dataset_id, function(gse) {
  data.frame(
    dataset_id = gse,
    kind = c("series_brief_soft", "samples_brief_soft", "supplementary_index"),
    filename = c(sprintf("%s_series_brief.soft", gse), sprintf("%s_samples_brief.soft", gse),
                 sprintf("%s_supplementary_index.html", gse)),
    url = c(
      sprintf("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&targ=self&form=text&view=brief", gse),
      sprintf("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&targ=gsm&form=text&view=brief", gse),
      sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/", geo_bucket(gse), gse)
    ), stringsAsFactors = FALSE
  )
}))

disk_check <- function(sample_path, dataset_path, audit_path, download_path) {
  problems <- validate_content(sample_manifest, dataset_manifest)
  actual_sample <- tryCatch(read_tsv_character(sample_path), error = function(e) NULL)
  actual_dataset <- tryCatch(read_tsv_character(dataset_path), error = function(e) NULL)
  actual_audit <- tryCatch(read_tsv_character(audit_path), error = function(e) NULL)
  if (!is.null(actual_sample) &&
      any(is.na(actual_sample$exclusion_reason) | !nzchar(actual_sample$exclusion_reason) |
          actual_sample$exclusion_reason == "NR")) {
    problems <- c(problems, "formal_sample_exclusion_reason_missing")
  }
  if (!is.null(actual_sample) && length(frozen_design_problems(actual_sample))) {
    problems <- c(problems, "formal_case_control_or_frozen_design_mismatch")
  }
  if (is.null(actual_sample) || !identical(actual_sample, as_disk_character(sample_manifest))) {
    problems <- c(problems, "formal_sample_manifest_mismatch")
  }
  if (is.null(actual_dataset) || !identical(actual_dataset, as_disk_character(dataset_manifest))) {
    problems <- c(problems, "formal_dataset_manifest_mismatch")
  }
  if (is.null(actual_audit) || !identical(actual_audit, as_disk_character(audit))) {
    problems <- c(problems, "formal_independence_audit_mismatch")
  }
  log <- tryCatch(read_tsv_character(download_path), error = function(e) NULL)
  required <- c("dataset_id", "kind", "filename", "url", "accessed_at", "file_size_bytes", "sha256")
  if (is.null(log) || !identical(names(log), required)) {
    problems <- c(problems, "download_manifest_schema_mismatch")
  } else {
    key <- paste(log$dataset_id, log$kind, log$filename, sep = "|")
    expected_key <- paste(expected_download_rows$dataset_id, expected_download_rows$kind,
                          expected_download_rows$filename, sep = "|")
    ord <- match(key, expected_key)
    if (anyNA(ord) || anyDuplicated(key) || nrow(log) != nrow(expected_download_rows) ||
        !identical(log$url, expected_download_rows$url[ord])) {
      problems <- c(problems, "download_manifest_URL_mismatch")
    }
    if (any(!grepl("^2026-07-24T[0-9]{2}:[0-9]{2}:[0-9]{2}\\+08:00$", log$accessed_at))) {
      problems <- c(problems, "download_manifest_accessed_at_mismatch")
    }
    for (i in seq_len(nrow(log))) {
      path <- file.path(raw_dir, log$filename[i])
      if (!file.exists(path)) {
        problems <- c(problems, "raw_metadata_missing")
      } else {
        actual_sha <- strsplit(system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE),
                               "[[:space:]]+")[[1L]][1L]
        if (!identical(actual_sha, log$sha256[i])) problems <- c(problems, "download_manifest_SHA_mismatch")
        if (!identical(as.character(file.info(path)$size), log$file_size_bytes[i])) {
          problems <- c(problems, "download_manifest_size_mismatch")
        }
      }
    }
  }
  unique(problems)
}

formal_paths <- list(
  sample = file.path(derived_dir, "geo_sample_manifest.tsv"),
  dataset = file.path(derived_dir, "geo_dataset_manifest.tsv"),
  audit = file.path(results_dir, "geo_sample_independence_audit.tsv"),
  download = file.path(raw_dir, "download_manifest.tsv")
)

if (self_test) {
  baseline <- disk_check(formal_paths$sample, formal_paths$dataset, formal_paths$audit, formal_paths$download)
  if (length(baseline)) stop("Baseline check failed before self-test: ", paste(baseline, collapse = ", "))
  test_case <- function(name, target, mutate, expected_problem) {
    td <- tempfile(pattern = "geo-self-test-")
    dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)
    paths <- lapply(formal_paths, function(x) {
      y <- file.path(td, basename(x)); file.copy(x, y); y
    })
    x <- read_tsv_character(paths[[target]])
    x <- mutate(x)
    Sys.chmod(paths[[target]], mode = "0644")
    write.table(x, paths[[target]], sep = "\t", quote = FALSE, row.names = FALSE, na = "NR")
    found <- disk_check(paths$sample, paths$dataset, paths$audit, paths$download)
    if (!expected_problem %in% found) stop("Self-test did not detect ", name, "; got: ", paste(found, collapse = ","))
    TRUE
  }
  tests <- c(
    sample_tamper = test_case("formal sample manifest tampering", "sample",
      function(x) { x$age[1] <- "999"; x }, "formal_sample_manifest_mismatch"),
    dataset_role_tamper = test_case("dataset role tampering", "dataset",
      function(x) { x$analysis_role[1] <- "exploratory_only"; x }, "formal_dataset_manifest_mismatch"),
    audit_tamper = test_case("audit tampering", "audit",
      function(x) { x$gsm_count[1] <- "1"; x }, "formal_independence_audit_mismatch"),
    exclusion_reason_missing = test_case("per-GSM exclusion reason missing", "sample",
      function(x) { x$exclusion_reason[1] <- "NR"; x }, "formal_sample_exclusion_reason_missing"),
    case_control_swap = test_case("case/control swap", "sample",
      function(x) {
        a <- which(x$dataset_id == "GSE98793" & x$diagnosis == "MDD")[1]
        x$diagnosis[a] <- "control"
        x
      }, "formal_case_control_or_frozen_design_mismatch"),
    followup_condition_swap = test_case("GSE38206 diagnosis-specific follow-up condition swap", "sample",
      function(x) {
        a <- which(x$dataset_id == "GSE38206" & x$diagnosis == "MDD" &
                     x$timepoint == "week_8")[1]
        b <- which(x$dataset_id == "GSE38206" & x$diagnosis == "control" &
                     x$timepoint == "week_8")[1]
        x$condition[c(a, b)] <- x$condition[c(b, a)]
        x
      }, "formal_case_control_or_frozen_design_mismatch"),
    download_url_tamper = test_case("download URL tampering", "download",
      function(x) { x$url[1] <- "https://example.invalid/"; x }, "download_manifest_URL_mismatch"),
    download_accessed_tamper = test_case("download accessed_at tampering", "download",
      function(x) { x$accessed_at[1] <- "NR"; x }, "download_manifest_accessed_at_mismatch"),
    download_sha_tamper = test_case("download SHA tampering", "download",
      function(x) {
        last <- substr(x$sha256[1], 64, 64)
        x$sha256[1] <- paste0(substr(x$sha256[1], 1, 63), ifelse(last == "0", "1", "0"))
        x
      },
      "download_manifest_SHA_mismatch")
  )
  if (!all(tests)) stop("Self-test failed.")
  message("PASS: 9 real-path tamper tests detected formal outputs, exclusion reason, case/control and follow-up condition errors, and download URL/accessed_at/SHA changes.")
  quit(save = "no", status = 0L)
}

problems <- validate_content(sample_manifest, dataset_manifest)
if (length(problems)) stop("Validation failed: ", paste(problems, collapse = ", "))

if (check_only) {
  problems <- disk_check(formal_paths$sample, formal_paths$dataset, formal_paths$audit, formal_paths$download)
  if (length(problems)) stop("Formal-output check failed: ", paste(problems, collapse = ", "))
  message("PASS: raw rebuild exactly matches all three formal outputs; inclusion, independence, URL/time/SHA checks passed.")
  print(dataset_manifest[, c("dataset_id", "gsm_count", "independent_donor_count",
                             "included_case_gsm", "included_control_gsm", "excluded_gsm",
                             "design_verification")], row.names = FALSE)
  quit(save = "no", status = 0L)
}

write.table(sample_manifest, file.path(derived_dir, "geo_sample_manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NR")
write.table(dataset_manifest, file.path(derived_dir, "geo_dataset_manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NR")
write.table(audit, file.path(results_dir, "geo_sample_independence_audit.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NR")

message("DONE: wrote ", nrow(sample_manifest), " sample records across ", nrow(dataset_manifest), " GEO datasets.")
print(dataset_manifest[, c("dataset_id", "gsm_count", "independent_donor_count",
                           "included_case_gsm", "included_control_gsm", "excluded_gsm",
                           "design_verification")], row.names = FALSE)
