library(tidyverse)

raw_data_path = fs::path("data-raw")
processed_data_path = fs::path("data")

serology_train <- read.delim(fs::path(raw_data_path, "publicData_serology.tsv"))
serology_test <- read.delim(fs::path(raw_data_path, "2025LJI_serology.tsv"))

serology_train_feature_df = serology_train %>% 
  dplyr::select(timepoint, virus_strain, assay) %>% 
  distinct()

serology_test_feature_df = serology_test %>% 
  dplyr::select(timepoint, virus_strain, assay) %>% 
  distinct()

serology_common_feature_df <- inner_join(
  serology_train_feature_df,
  serology_test_feature_df,
  by = c("timepoint", "virus_strain", "assay")
) %>% 
  filter(timepoint == "Pre-vacc")

serology_train_feature_df_filtered <- serology_train %>%
  semi_join(
    serology_common_feature_df,
    by = c("timepoint", "virus_strain", "assay")
  )

serology_train_df <- serology_train_feature_df_filtered %>%
  pivot_wider(
    id_cols = participant_id,
    names_from = c(timepoint, virus_strain, assay),
    values_from = value,
    names_sep = "_"
  ) %>% 
  janitor::clean_names() %>% 
  arrange(participant_id)


serology_test_feature_df_filtered <- serology_test %>%
  semi_join(
    serology_common_feature_df,
    by = c("timepoint", "virus_strain", "assay")
  )

serology_test_df <- serology_test_feature_df_filtered %>%
  pivot_wider(
    id_cols = participant_id,
    names_from = c(timepoint, virus_strain, assay),
    values_from = value,
    names_sep = "_"
  ) %>% 
  janitor::clean_names() %>% 
  arrange(participant_id)

saveRDS(serology_train_df, file = fs::path(processed_data_path, "serology_train_df.rds"))
saveRDS(serology_test_df, file = fs::path(processed_data_path, "serology_test_df.rds"))

rm(list = ls())
