# ==============================================================================
# task2.1_feature_set_comparison.R
# Compare demographic, cytokine, cytometry, and serology feature sets as
# predictors of Task 2.1 (Day-28 HAI antibody magnitude), using
# predictomics::compare_pipelines(). Every pipeline (reference and each
# feature-set option) shares the same configuration: z-score engineering, no
# feature selection, elastic net (glmnet) modelling.
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


# ==============================================================================
# Restrict to participants with the response and every feature set available,
# and align row order across Y and every predictor matrix
# ==============================================================================

common_ids <- Reduce(
  intersect,
  c(list(response_df$participant_id),
    lapply(feature_sets, function(df) df$participant_id))
) %>%
  sort()

response_aligned <- response_df %>%
  filter(participant_id %in% common_ids) %>%
  arrange(participant_id)

Y <- response_aligned$value
names(Y) <- response_aligned$participant_id

to_aligned_matrix <- function(df) {
  df %>%
    filter(participant_id %in% common_ids) %>%
    arrange(participant_id) %>%
    tibble::column_to_rownames("participant_id") %>%
    as.matrix()
}

X_list <- lapply(feature_sets, to_aligned_matrix)

stopifnot(
  all(vapply(X_list, nrow, integer(1)) == length(Y)),
  all(vapply(X_list, function(m) identical(rownames(m), names(Y)), logical(1)))
)

# Combined predictor matrix (all feature sets together), used as the fixed
# "Reference" predictor set that each individual feature set is compared
# against; columns are prefixed by feature set to avoid name clashes.
X_combined <- do.call(cbind, Map(function(name, m) {
  colnames(m) <- paste0(name, "_", colnames(m))
  m
}, names(X_list), X_list))


# ==============================================================================
# Compare feature sets under a shared reference pipeline:
# z-score engineering, no feature selection, elastic net modelling.
# model_params$impute = "mean" guards against sporadic missing values within
# a participant's feature set (e.g. an analyte not measured for them).
# ==============================================================================

reference_pipeline <- list(
  engineering_params = list(method = "engineer", col_transform = "z"),
  selection_params   = NULL,
  model_params        = list(method = "glmnet", impute = "mean")
)

comparison <- compare_pipelines(
  Y                = Y,
  X                = X_combined,
  option_type      = "predictors",
  option_choices   = X_list,
  reference_params = reference_pipeline,
  cv_type          = "kfold",
  folds            = 10L,
  seed             = 12345L,
  metric           = "sRMSE",
  verbose          = TRUE
)

print(comparison)

comparison_plot <- plot(comparison, metric = "all")


# ==============================================================================
# Save results
# ==============================================================================

saveRDS(comparison, fs::path(results_path, "task2.1_feature_set_comparison.rds"))

ggplot2::ggsave(
  filename = fs::path(results_path, "task2.1_feature_set_comparison.png"),
  plot     = comparison_plot,
  width    = 10, height = 7, dpi = 300
)

rm(list = ls())
