# ==============================================================================
# task2.1_feature_set_comparison.R
# Compare demographic, cytokine, cytometry, and serology feature sets as
# predictors of Task 2.1 (Day-28 HAI antibody magnitude).
#
# The four feature sets do not share a common set of participants with the
# response (no participant has every feature set measured), so
# predictomics::compare_pipelines(option_type = "predictors") cannot be used
# directly - it requires every option's predictor matrix to share one fixed
# Y (i.e. the same n and row order across all options). Instead, each feature
# set is fit via its own predict_cv() call, restricted to the participants it
# shares with the response, and the resulting metrics are combined into a
# single comparison plot. The reference pipeline (z-score engineering, no
# feature selection, elastic net modelling) and a per-feature-set baseline
# (X = NULL, mean-only) are identical in spirit to what compare_pipelines()
# would have produced, just computed by hand across differing sample sizes.
# ==============================================================================

library(tidyverse)

if (!requireNamespace("predictomics", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("arthurhughes27/predictomics")
}
library(predictomics)

processed_data_path <- fs::path("data")
results_path        <- fs::path("analysis", "modelling", "results")
fs::dir_create(results_path)


# ==============================================================================
# Load processed feature sets and the task 2.1 response
# ==============================================================================

demographics_df <- readRDS(fs::path(processed_data_path, "demographics_train_df.rds"))
cytokine_df     <- readRDS(fs::path(processed_data_path, "cytokine_train_df.rds"))
cytometry_df    <- readRDS(fs::path(processed_data_path, "cytometry_train_df.rds"))
serology_df     <- readRDS(fs::path(processed_data_path, "serology_train_df.rds"))
response_df     <- readRDS(fs::path(processed_data_path, "task2.1_response_train.rds"))

feature_sets <- list(
  demographic = demographics_df,
  cytokines   = cytokine_df,
  cytometry   = cytometry_df,
  serology    = serology_df
)

reference_pipeline <- list(
  engineering_params = list(method = "engineer", col_transform = "z"),
  selection_params   = NULL,
  model_params       = list(method = "glmnet", impute = "mean")
)


# ==============================================================================
# Fit a baseline (X = NULL) and the reference pipeline for one feature set,
# restricted to the participants it shares with the response
# ==============================================================================

fit_feature_set <- function(feature_set_name, df, response_df,
                            reference_pipeline, folds = 10L, seed = 12345L) {

  common_ids <- intersect(response_df$participant_id, df$participant_id) %>%
    sort()

  response_aligned <- response_df %>%
    filter(participant_id %in% common_ids) %>%
    arrange(participant_id)

  Y <- response_aligned$value
  names(Y) <- response_aligned$participant_id

  X <- df %>%
    filter(participant_id %in% common_ids) %>%
    arrange(participant_id) %>%
    tibble::column_to_rownames("participant_id") %>%
    as.matrix()

  stopifnot(identical(rownames(X), names(Y)))

  n_folds <- min(folds, length(Y))

  baseline_fit <- suppressMessages(predict_cv(
    Y       = Y,
    X       = NULL,
    cv_type = "kfold",
    folds   = n_folds,
    seed    = seed,
    verbose = FALSE
  ))

  model_fit <- suppressMessages(predict_cv(
    Y                  = Y,
    X                  = X,
    cv_type            = "kfold",
    folds              = n_folds,
    seed               = seed,
    engineering_params = reference_pipeline$engineering_params,
    selection_params   = reference_pipeline$selection_params,
    model_params       = reference_pipeline$model_params,
    verbose            = FALSE
  ))

  bind_rows(
    tibble(feature_set = feature_set_name, role = "baseline", n = length(Y),
           !!!as.list(metrics(baseline_fit))),
    tibble(feature_set = feature_set_name, role = "model", n = length(Y),
           !!!as.list(metrics(model_fit)))
  )
}


# ==============================================================================
# Fit every feature set and combine metrics
# ==============================================================================

results <- purrr::imap(feature_sets, function(df, name) {
  message("[task2.1] Fitting feature set: ", name)
  fit_feature_set(name, df, response_df, reference_pipeline)
}) %>%
  bind_rows()

print(results)


# ==============================================================================
# Plot: metrics compared across feature sets, baseline vs. reference pipeline
# (z-score engineering, no selection, elastic net). Each feature set is fit
# on its own overlap with the response, so sample size (n) differs by set -
# this is annotated on the x-axis.
# ==============================================================================

plot_df <- results %>%
  mutate(feature_set_label = paste0(feature_set, "\n(n = ", n, ")")) %>%
  pivot_longer(
    cols      = c(RMSE, sRMSE, R2, SpearmanR),
    names_to  = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(metric, levels = c("RMSE", "sRMSE", "R2", "SpearmanR")),
    role   = factor(role, levels = c("baseline", "model"),
                    labels = c("Baseline (mean-only)",
                               "z-score + elastic net"))
  )

comparison_plot <- ggplot(plot_df, aes(x = feature_set_label, y = value, fill = role)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(
    aes(label = round(value, 3)),
    position = position_dodge(width = 0.7), vjust = -0.4, size = 3
  ) +
  facet_wrap(~metric, scales = "free_y") +
  labs(
    x        = "Feature set",
    y        = NULL,
    fill     = "Pipeline",
    title    = "Task 2.1 (antibody magnitude): feature set comparison",
    subtitle = paste0(
      "Each feature set fit on its own overlap of participants with the ",
      "response (n differs by set)"
    )
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(colour = "grey40", size = 10),
    axis.text.x      = element_text(size = 9),
    legend.position  = "bottom"
  )

comparison_plot


# ==============================================================================
# Save results
# ==============================================================================

saveRDS(results, fs::path(results_path, "task2.1_feature_set_comparison.rds"))

ggsave(
  filename = fs::path(results_path, "task2.1_feature_set_comparison.png"),
  plot     = comparison_plot,
  width    = 10, height = 7, dpi = 300
)

rm(list = ls())
