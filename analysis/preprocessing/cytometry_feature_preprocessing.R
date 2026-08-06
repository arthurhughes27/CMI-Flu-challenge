library(tidyverse)

raw_data_path = fs::path("data-raw")
processed_data_path = fs::path("data")

cytometry_train <- read.delim(fs::path(raw_data_path, "publicData_ex_vivo_flow.tsv"))
cytometry_test <- read.delim(fs::path(raw_data_path, "2025LJI_ex_vivo_flow.tsv"))

cytometry_train_feature_df = cytometry_train %>% 
  dplyr::select(timepoint, name, unit) %>% 
  distinct()

cytometry_test_feature_df = cytometry_test %>% 
  dplyr::select(timepoint, name, unit) %>% 
  distinct()

cytometry_common_feature_df <- inner_join(
  cytometry_train_feature_df,
  cytometry_test_feature_df,
  by = c("timepoint", "name", "unit")
) %>% 
  filter(timepoint == "Pre-vacc")

cytometry_train_feature_df_filtered <- cytometry_train %>%
  semi_join(
    cytometry_common_feature_df,
    by = c("timepoint", "name", "unit")
  )

cytometry_train_df <- cytometry_train_feature_df_filtered %>%
  pivot_wider(
    id_cols = participant_id,
    names_from = c(timepoint, name, unit),
    values_from = value,
    names_sep = "_"
  ) %>% 
  janitor::clean_names() %>% 
  arrange(participant_id)


cytometry_test_feature_df_filtered <- cytometry_test %>%
  semi_join(
    cytometry_common_feature_df,
    by = c("timepoint", "name", "unit")
  )

cytometry_test_df <- cytometry_test_feature_df_filtered %>%
  pivot_wider(
    id_cols = participant_id,
    names_from = c(timepoint, name, unit),
    values_from = value,
    names_sep = "_"
  ) %>% 
  janitor::clean_names() %>% 
  arrange(participant_id)

saveRDS(cytometry_train_df, file = fs::path(processed_data_path, "cytometry_train_df.rds"))
saveRDS(cytometry_test_df, file = fs::path(processed_data_path, "cytometry_test_df.rds"))

rm(list = ls())
