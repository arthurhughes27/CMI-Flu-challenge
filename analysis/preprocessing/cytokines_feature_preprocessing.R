library(tidyverse)

raw_data_path = fs::path("data-raw")
processed_data_path = fs::path("data")

cytokine_train <- read.delim(fs::path(raw_data_path, "publicData_cytokine.tsv")) %>% 
  mutate(material = tolower(material))
cytokine_test <- read.delim(fs::path(raw_data_path, "2025LJI_cytokine.tsv"))

cytokine_train_feature_df = cytokine_train %>% 
  dplyr::select(timepoint, analyte, material) %>% 
  distinct()

cytokine_test_feature_df = cytokine_test %>% 
  dplyr::select(timepoint, analyte, material) %>% 
  distinct()

cytokine_common_feature_df <- inner_join(
  cytokine_train_feature_df,
  cytokine_test_feature_df,
  by = c("timepoint", "analyte", "material")
) %>% 
  filter(timepoint == "Pre-vacc")

cytokine_train_feature_df_filtered <- cytokine_train %>%
  semi_join(
    cytokine_common_feature_df,
    by = c("timepoint", "analyte", "material")
  )

cytokine_train_df <- cytokine_train_feature_df_filtered %>%
  pivot_wider(
    id_cols = participant_id,
    names_from = c(timepoint, analyte, material),
    values_from = value,
    names_sep = "_"
  ) %>% 
  janitor::clean_names() %>% 
  arrange(participant_id)

cytokine_test_feature_df_filtered <- cytokine_test %>%
  semi_join(
    cytokine_common_feature_df,
    by = c("timepoint", "analyte", "material")
  )

cytokine_test_df <- cytokine_test_feature_df_filtered %>%
  pivot_wider(
    id_cols = participant_id,
    names_from = c(timepoint, analyte, material),
    values_from = value,
    names_sep = "_"
  ) %>% 
  janitor::clean_names() %>% 
  arrange(participant_id)

saveRDS(cytokine_train_df, file = fs::path(processed_data_path, "cytokine_train_df.rds"))
saveRDS(cytokine_test_df, file = fs::path(processed_data_path, "cytokine_test_df.rds"))

rm(list = ls())
