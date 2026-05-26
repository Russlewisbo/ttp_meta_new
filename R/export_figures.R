## export_figures.R
## Renders Figure 2 (forest plot) and Figure 3 (posterior density) as
## high-resolution TIFFs suitable for journal submission.
##
## Usage (from project root):
##   source("R/export_figures.R")
##
## Output: figures/fig2_forest.tiff, figures/fig3_posterior.tiff

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(metafor)
  library(posterior)
  library(ggdist)
})

source("R/01_prepare_data.R")
source("R/02_fit_models.R")

dir.create("figures", showWarnings = FALSE)

# ── Load data & models (uses brms file cache) ──────────────────────────────
es_all  <- prepare_effect_sizes("data/TTP_MetaAnalysis_Extraction_Complete_v4.xlsx",
                                include_inclusive_pool = FALSE,
                                drop_contaminated_ttp  = TRUE)
es_mort <- es_for_outcome(es_all, "mortality")
fits    <- fit_all_models(es_mort, cache_dir = "_brms_cache/primary")

base_post <- pooled_summary(fits$fit_base)

# ── Figure 2: Forest plot ───────────────────────────────────────────────────
forest_data <- es_mort |>
  arrange(yi) |>
  mutate(
    lbl = paste0(study_id,
                 ifelse(arm_id == "all", "", paste0(" [", arm_id, "]"))),
    lbl     = make.unique(lbl, sep = " #"),
    lbl     = factor(lbl, levels = lbl),
    or      = exp(yi),
    ci_low  = exp(yi - 1.96 * sei),
    ci_high = exp(yi + 1.96 * sei),
    type    = "Study"
  )

pooled_row <- tibble(
  lbl     = "POOLED (Bayesian)",
  or      = base_post$median_OR,
  ci_low  = base_post$ci_OR[1],
  ci_high = base_post$ci_OR[2],
  type    = "Pooled"
)
pi_row <- tibble(
  lbl     = "POOLED (Bayesian)",
  pi_low  = base_post$pi_OR[1],
  pi_high = base_post$pi_OR[2]
)

forest_data$lbl <- factor(forest_data$lbl,
                           levels = c("POOLED (Bayesian)", levels(forest_data$lbl)))
pooled_row$lbl  <- factor(pooled_row$lbl, levels = levels(forest_data$lbl))
pi_row$lbl      <- factor(pi_row$lbl,     levels = levels(forest_data$lbl))

fig2 <- ggplot(bind_rows(forest_data, pooled_row),
               aes(x = or, y = lbl, color = type)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_point(size = 2.5) +
  geom_errorbar(aes(xmin = pmax(ci_low, 0.05), xmax = pmin(ci_high, 50)),
                width = 0.25, linewidth = 0.7, orientation = "y") +
  geom_errorbar(data = pi_row,
                aes(y = lbl,
                    xmin = pmax(pi_low,  0.05),
                    xmax = pmin(pi_high, 50)),
                width = 0.15, linewidth = 0.7, linetype = "dashed",
                color = "darkred", inherit.aes = FALSE) +
  scale_color_manual(values = c("Study" = "steelblue", "Pooled" = "darkred")) +
  scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 8, 16, 32),
                limits  = c(0.05, 50)) +
  labs(title    = "Mortality forest plot — primary pool",
       subtitle = "Solid red bar = 95% CrI of pooled OR; dashed red bar = 95% prediction interval",
       x        = "OR, short vs long TTP (log scale)",
       y        = NULL,
       color    = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position    = "bottom",
        axis.text.y        = element_text(size = 7, lineheight = 0.85),
        panel.grid.major.y = element_line(color = "gray92"),
        panel.grid.minor   = element_blank())

fig2_h <- max(6, nrow(es_mort) * 0.55 + 2)

ggsave("figures/fig2_forest.tiff",
       plot     = fig2,
       width    = 12,
       height   = fig2_h,
       units    = "in",
       dpi      = 600,
       compression = "lzw")

message("Saved figures/fig2_forest.tiff  (", round(fig2_h, 1), " × 12 in, 600 dpi)")

# ── Figure 3: Posterior density ─────────────────────────────────────────────
fig3 <- tibble(OR = base_post$OR_draws) |>
  ggplot(aes(x = OR)) +
  stat_halfeye(fill = "#1f78b4", alpha = 0.6) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_x_continuous(limits = c(0.5, 6), breaks = seq(0.5, 6, 0.5)) +
  labs(x = "Pooled OR (short vs long TTP)", y = "Posterior density") +
  theme_minimal(base_size = 12)

ggsave("figures/fig3_posterior.tiff",
       plot     = fig3,
       width    = 8,
       height   = 5,
       units    = "in",
       dpi      = 600,
       compression = "lzw")

message("Saved figures/fig3_posterior.tiff  (5 × 8 in, 600 dpi)")
