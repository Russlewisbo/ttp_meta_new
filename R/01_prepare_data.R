#' Prepare analysis-ready effect-size table from the v4 workbook.
#'
#' Reads TTP_MetaAnalysis_Extraction_Complete_v4.xlsx, joins the four extraction
#' tables, computes per-study log-OR + variance via metafor::escalc(), and
#' returns one tidy tibble per outcome family (mortality, persistent_bacteremia,
#' microbiological_clearance, relapse).
#'
#' Run from the project root:
#'   source("R/01_prepare_data.R")
#'   tbl_es_all <- prepare_effect_sizes(
#'     "data/TTP_MetaAnalysis_Extraction_Complete_v4.xlsx",
#'     include_inclusive_pool = FALSE
#'   )
#'
#' @param xlsx_path Path to v4 workbook.
#' @param include_inclusive_pool If FALSE (default), drops rows whose
#'   study.notes contain "[INCLUSIVE_ONLY]" — yielding the n=100 strict pool.
#'   If TRUE, keeps all 140 studies for sensitivity analysis.
#' @param drop_contaminated_ttp If TRUE, drops rows whose notes contain
#'   "[CONTAMINATED_TTP]" — TTP used as part of the CRBSI case definition
#'   rather than a freely measured prognostic variable.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(stringr)
  library(metafor)
  library(purrr)
})

prepare_effect_sizes <- function(xlsx_path,
                                 include_inclusive_pool = FALSE,
                                 drop_contaminated_ttp = TRUE) {
  stopifnot(file.exists(xlsx_path))

  tbl_study      <- readxl::read_excel(xlsx_path, sheet = "tbl_study")
  tbl_population <- readxl::read_excel(xlsx_path, sheet = "tbl_population")
  tbl_ttp        <- readxl::read_excel(xlsx_path, sheet = "tbl_ttp")
  tbl_outcomes   <- readxl::read_excel(xlsx_path, sheet = "tbl_outcomes")

  # ---- Filter to chosen analysis pool ----
  inclusive_ids <- tbl_study |>
    filter(str_detect(coalesce(notes, ""), fixed("[INCLUSIVE_ONLY]"))) |>
    pull(study_id)

  contaminated_ids <- tbl_study |>
    filter(str_detect(coalesce(notes, ""), fixed("[CONTAMINATED_TTP]"))) |>
    pull(study_id)

  if (!include_inclusive_pool) {
    keep_ids <- setdiff(tbl_study$study_id, inclusive_ids)
  } else {
    keep_ids <- tbl_study$study_id
  }
  if (drop_contaminated_ttp) {
    keep_ids <- setdiff(keep_ids, contaminated_ids)
  }

  tbl_study      <- tbl_study      |> filter(study_id %in% keep_ids)
  tbl_population <- tbl_population |> filter(study_id %in% keep_ids)
  tbl_ttp        <- tbl_ttp        |> filter(study_id %in% keep_ids)
  tbl_outcomes   <- tbl_outcomes   |> filter(study_id %in% keep_ids)

  message(sprintf("Pool: %d studies (inclusive_pool=%s, drop_contaminated=%s)",
                  length(keep_ids), include_inclusive_pool, drop_contaminated_ttp))

  # ---- Effect-size computation ----
  # Strategy by row, in order of preference:
  #   (a) If 2x2 counts available (events_short_ttp, n_short_ttp, events_long_ttp, n_long_ttp) → escalc OR
  #   (b) Else if effect_estimate + ci_lower + ci_upper reported → back-calculate log-OR/HR/RR and SE
  #   (c) Else if or_per_hour + CI reported (continuous) → per-hour log-OR; not pooled with (a)/(b) directly
  #   (d) Else if beta_per_hour + se_beta → continuous

  tbl_es <- tbl_outcomes |>
    left_join(select(tbl_population, study_id, arm_id, pathogen_class,
                     pathogen_species, infection_source, pct_appropriate_empiric,
                     blood_culture_system, source_control),
              by = c("study_id", "arm_id")) |>
    left_join(select(tbl_ttp, study_id, arm_id, ttp_reporting, ttp_cutpoint_hours,
                     ttp_cutpoint_basis),
              by = c("study_id", "arm_id")) |>
    left_join(select(tbl_study, study_id, pub_year, study_design, rob_overall),
              by = "study_id")

  # ---- (a) 2x2 escalc -------------------------------------------------------
  es_2x2 <- tbl_es |>
    filter(!is.na(events_short_ttp), !is.na(n_short_ttp),
           !is.na(events_long_ttp), !is.na(n_long_ttp),
           n_short_ttp > 0, n_long_ttp > 0) |>
    mutate(
      ai = events_short_ttp, n1i = n_short_ttp,
      ci = events_long_ttp,  n2i = n_long_ttp
    ) |>
    {\(df) {
      e <- metafor::escalc(measure = "OR", ai = df$ai, bi = df$n1i - df$ai,
                            ci = df$ci, di = df$n2i - df$ci,
                            add = 0.5, to = "only0")
      df$yi  <- as.numeric(e$yi)
      df$vi  <- as.numeric(e$vi)
      df$sei <- sqrt(df$vi)
      df$es_source <- "two_by_two"
      df
    }}()

  # ---- (b) Reported OR/HR/RR with CI ---------------------------------------
  es_reported <- tbl_es |>
    anti_join(es_2x2, by = c("study_id", "arm_id", "outcome_type", "outcome_timepoint_days")) |>
    filter(!is.na(effect_estimate), !is.na(ci_lower), !is.na(ci_upper),
           effect_estimate > 0, ci_lower > 0, ci_upper > 0,
           effect_measure %in% c("OR", "HR", "RR", "IRR")) |>
    mutate(
      yi = log(effect_estimate),
      sei = (log(ci_upper) - log(ci_lower)) / (2 * qnorm(0.975)),
      vi = sei^2,
      es_source = if_else(adjusted == TRUE, "reported_adjusted", "reported_unadjusted")
    )

  # ---- (c) Continuous TTP: per-hour OR or beta -----------------------------
  es_continuous <- tbl_es |>
    anti_join(es_2x2,      by = c("study_id", "arm_id", "outcome_type", "outcome_timepoint_days")) |>
    anti_join(es_reported, by = c("study_id", "arm_id", "outcome_type", "outcome_timepoint_days")) |>
    filter((!is.na(beta_per_hour) & !is.na(se_beta)) |
           (!is.na(or_per_hour)  & !is.na(or_per_hour_ci_low) & !is.na(or_per_hour_ci_high))) |>
    mutate(
      yi = if_else(!is.na(beta_per_hour),
                   beta_per_hour,
                   log(or_per_hour)),
      sei = if_else(!is.na(se_beta),
                    se_beta,
                    (log(or_per_hour_ci_high) - log(or_per_hour_ci_low)) / (2 * qnorm(0.975))),
      vi = sei^2,
      es_source = "continuous_per_hour"
    )

  out <- bind_rows(es_2x2, es_reported, es_continuous) |>
    mutate(
      organism_saureus = str_detect(coalesce(pathogen_species, ""),
                                    regex("S\\. aureus|staphylococcus aureus", ignore_case = TRUE)),
      es_id = row_number()
    ) |>
    relocate(es_id, study_id, arm_id, outcome_type, outcome_timepoint_days,
             yi, vi, sei, es_source)

  attr(out, "n_studies") <- length(keep_ids)
  attr(out, "tbl_study") <- tbl_study
  attr(out, "tbl_population") <- tbl_population
  attr(out, "tbl_ttp") <- tbl_ttp
  attr(out, "tbl_outcomes") <- tbl_outcomes

  out
}

# Convenience: filter for one outcome family
es_for_outcome <- function(es, outcome) {
  es |> dplyr::filter(outcome_type == outcome)
}
