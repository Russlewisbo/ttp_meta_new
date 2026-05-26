#' Fit the Bayesian random-effects meta-analysis suite.
#'
#' Mirrors the original v1 analysis (nine models): pooled, gram_pos, gram_neg,
#' pathogen meta-regression, cutpoint meta-regression, 2x2-only, reported-only,
#' adjustment-status, publication-year. Fits with brms (Stan backend).
#'
#' Usage:
#'   source("R/01_prepare_data.R")
#'   source("R/02_fit_models.R")
#'   es_mort <- prepare_effect_sizes("data/TTP_MetaAnalysis_Extraction_Complete_v4.xlsx") |>
#'     es_for_outcome("mortality")
#'   fits <- fit_all_models(es_mort, cache_dir = "_brms_cache")
#'
#' This will compile each Stan model once (~30-60s each first time, then cached).

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(tidyr)
})

# Priors used throughout (same as original analysis):
.prior_intercept <- prior(normal(0, 1),   class = "Intercept")
.prior_tau       <- prior(cauchy(0, 0.5), class = "sd")
.prior_mod       <- prior(normal(0, 0.5), class = "b")

.brms_defaults <- function(seed) list(
  chains    = 4,
  iter      = 4000,
  warmup    = 1000,
  cores     = parallel::detectCores(),
  control   = list(adapt_delta = 0.99, max_treedepth = 12),
  seed      = seed,
  silent    = 2,
  refresh   = 0
)

fit_one <- function(es, formula, priors, name, cache_dir, seed = 1) {
  args <- c(
    list(formula = formula,
         data    = es,
         prior   = priors,
         file    = file.path(cache_dir, name),
         file_refit = "on_change"),
    .brms_defaults(seed))
  do.call(brm, args)
}

fit_all_models <- function(es_mort, cache_dir = "_brms_cache") {
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

  fits <- list()

  # 1. Overall pooled effect
  fits$fit_base <- fit_one(
    es_mort,
    formula = bf(yi | se(sei) ~ 1 + (1 | study_id)),
    priors  = c(.prior_intercept, .prior_tau),
    name = "fit_base", cache_dir = cache_dir
  )

  # 2-3. Pathogen subgroups
  es_gpos <- es_mort |> filter(pathogen_class == "gram_positive")
  es_gneg <- es_mort |> filter(pathogen_class == "gram_negative")
  if (nrow(es_gpos) >= 3)
    fits$fit_gram_pos <- fit_one(
      es_gpos,
      formula = bf(yi | se(sei) ~ 1 + (1 | study_id)),
      priors  = c(.prior_intercept, .prior_tau),
      name = "fit_gram_pos", cache_dir = cache_dir)
  if (nrow(es_gneg) >= 3)
    fits$fit_gram_neg <- fit_one(
      es_gneg,
      formula = bf(yi | se(sei) ~ 1 + (1 | study_id)),
      priors  = c(.prior_intercept, .prior_tau),
      name = "fit_gram_neg", cache_dir = cache_dir)

  es_fungal <- es_mort |> filter(pathogen_class == "fungal")
  if (nrow(es_fungal) >= 3)
    fits$fit_fungal <- fit_one(
      es_fungal,
      formula = bf(yi | se(sei) ~ 1 + (1 | study_id)),
      priors  = c(.prior_intercept, .prior_tau),
      name = "fit_fungal", cache_dir = cache_dir)

  # 4. Pathogen meta-regression
  es_metareg <- es_mort |>
    filter(pathogen_class %in% c("gram_positive", "gram_negative")) |>
    mutate(pathogen_class = factor(pathogen_class))
  fits$fit_metareg <- fit_one(
    es_metareg,
    formula = bf(yi | se(sei) ~ 1 + pathogen_class + (1 | study_id)),
    priors  = c(.prior_intercept, .prior_tau, .prior_mod),
    name = "fit_metareg", cache_dir = cache_dir)

  # 5. ⭐ TTP cutpoint meta-regression
  es_cut <- es_mort |>
    filter(es_source == "two_by_two", !is.na(ttp_cutpoint_hours)) |>
    mutate(cutpoint_c = ttp_cutpoint_hours - median(ttp_cutpoint_hours, na.rm = TRUE))
  if (nrow(es_cut) >= 5)
    fits$fit_cutpoint <- fit_one(
      es_cut,
      formula = bf(yi | se(sei) ~ 1 + cutpoint_c + (1 | study_id)),
      priors  = c(.prior_intercept, .prior_tau, .prior_mod),
      name = "fit_cutpoint", cache_dir = cache_dir)

  # 6. 2x2 only
  es_2x2 <- es_mort |> filter(es_source == "two_by_two")
  if (nrow(es_2x2) >= 3)
    fits$fit_2x2 <- fit_one(
      es_2x2,
      formula = bf(yi | se(sei) ~ 1 + (1 | study_id)),
      priors  = c(.prior_intercept, .prior_tau),
      name = "fit_2x2", cache_dir = cache_dir)

  # 7. Reported only (back-calculated)
  es_rep <- es_mort |> filter(grepl("reported", es_source))
  if (nrow(es_rep) >= 3)
    fits$fit_reported <- fit_one(
      es_rep,
      formula = bf(yi | se(sei) ~ 1 + (1 | study_id)),
      priors  = c(.prior_intercept, .prior_tau),
      name = "fit_reported", cache_dir = cache_dir)

  # 8. Effect-size source meta-regression (two_by_two as reference)
  es_src <- es_mort |>
    filter(es_source %in% c("two_by_two", "reported_adjusted", "reported_unadjusted")) |>
    mutate(es_source_f = factor(es_source,
                                levels = c("two_by_two",
                                           "reported_adjusted",
                                           "reported_unadjusted")))
  if (length(unique(es_src$es_source_f)) >= 2 && nrow(es_src) >= 6)
    fits$fit_es_source <- fit_one(
      es_src,
      formula = bf(yi | se(sei) ~ 1 + es_source_f + (1 | study_id)),
      priors  = c(.prior_intercept, .prior_tau, .prior_mod),
      name    = "fit_es_source", cache_dir = cache_dir)

  # 9. Adjustment status
  es_adj <- es_mort |>
    filter(es_source %in% c("reported_adjusted", "reported_unadjusted")) |>
    mutate(adjusted_lbl = factor(if_else(es_source == "reported_adjusted",
                                         "adjusted", "unadjusted")))
  if (length(unique(es_adj$adjusted_lbl)) == 2 && nrow(es_adj) >= 4)
    fits$fit_adjustment <- fit_one(
      es_adj,
      formula = bf(yi | se(sei) ~ 1 + adjusted_lbl + (1 | study_id)),
      priors  = c(.prior_intercept, .prior_tau, .prior_mod),
      name = "fit_adjustment", cache_dir = cache_dir)

  # 9. Publication year trend
  es_yr <- es_mort |>
    filter(!is.na(pub_year)) |>
    mutate(year_c = pub_year - median(pub_year, na.rm = TRUE))
  fits$fit_year <- fit_one(
    es_yr,
    formula = bf(yi | se(sei) ~ 1 + year_c + (1 | study_id)),
    priors  = c(.prior_intercept, .prior_tau, .prior_mod),
    name = "fit_year", cache_dir = cache_dir)

  fits
}

# Helper: extract pooled OR posterior + I² from a brms fit (DerSimonian-style I²)
pooled_summary <- function(fit) {
  draws <- as_draws_df(fit)
  logOR <- draws$b_Intercept
  tau   <- draws$sd_study_id__Intercept
  # I² via Higgins method using a "typical" within-study variance
  vi_typical <- mean(fit$data$sei^2, na.rm = TRUE)
  i2 <- tau^2 / (tau^2 + vi_typical) * 100
  # 95% prediction interval: marginal predictive distribution for a new study
  pi_logOR <- rnorm(length(logOR), mean = logOR, sd = tau)
  list(
    logOR_draws = logOR,
    OR_draws    = exp(logOR),
    tau_draws   = tau,
    i2_draws    = i2,
    median_OR   = median(exp(logOR)),
    ci_OR       = quantile(exp(logOR), c(0.025, 0.975)),
    pi_OR       = quantile(exp(pi_logOR), c(0.025, 0.975)),
    p_harm      = mean(exp(logOR) > 1),
    median_tau  = median(tau),
    median_i2   = median(i2)
  )
}
