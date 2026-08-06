library(tidyverse)

raw_data_path = fs::path("data-raw")
processed_data_path = fs::path("data")

cytokine_train <- read.delim(fs::path(raw_data_path, "publicData_cytokine.tsv"))

cytometry_train <- read.delim(fs::path(raw_data_path, "publicData_ex_vivo_flow.tsv"))

# aim_train <- read.delim("~/Desktop/Work/PhD/CMI-Flu-challenge/data-raw/2025LJI_aim.tsv")

serology_train <- read.delim(fs::path(raw_data_path, "publicData_serology.tsv"))
serology_test <- read.delim(fs::path(raw_data_path, "2025LJI_serology.tsv"))

# Task 1: early systems responses
# Task 1.1: IP10 (CXCL10) response
task1.1_response_train <- cytokine_train %>%
  filter(
    analyte == "IP10",
    timepoint %in% c("Pre-vacc", "1")
  ) %>%
  select(participant_id, timepoint, value) %>%
  tidyr::pivot_wider(
    names_from = timepoint,
    values_from = value
  ) %>%
  mutate(
    value = `1` / `Pre-vacc`
  ) %>%
  dplyr::select(participant_id, value) %>% 
  arrange(participant_id)
  

# Task 1.2: monocyte frequency
task1.2_response_train <- cytometry_train %>%
  filter(
    name == "Classical_monocytes",
    timepoint == 1
  ) %>%
  dplyr::select(participant_id, value) %>% 
  arrange(participant_id)

# Task 1.3: plasmablast frequency 
task1.3_response_train <- cytometry_train %>%
  filter(
    name == "Antibody-secreting_cells_(ASC)",
    timepoint == 7
  ) %>%
  dplyr::select(participant_id, value) %>% 
  arrange(participant_id)

# Task 1.4: frequency of antigen-reactive CD4 T cells


# Task 2: antibody responses
# Task 2.1: antibody (HAI) magnitude
strains_2025 = c("H1N1 A/Victoria/4897/2022",
                 "H3N2 A/District Of Columbia/27/2023_MDCK",
                 "Vic B/Austria/1359417/2021")

task2.1_response_train = serology_train %>%
  filter(virus_strain %in% strains_2025,
         timepoint == 28,
         assay == "hai") %>%
  group_by(participant_id) %>%
  summarise(
    value = exp(mean(log(value))),
    .groups = "drop"
  ) %>%
  arrange(participant_id)

# Task 2.2: antibody (HAI) breadth
strains_test = serology_test %>% 
  pull(virus_strain) %>% 
  unique()

task2.2_response_train = serology_train %>%
  filter(virus_strain %in% strains_test,
         timepoint == 28,
         assay == "hai") %>%
  group_by(participant_id) %>%
  summarise(
    value = exp(mean(log(value))),
    .groups = "drop"
  ) %>%
  arrange(participant_id)

# Task 2.3: antibody durability
task2.3_response_train = serology_train %>%
  filter(virus_strain %in% strains_test,
         timepoint == 365,
         assay == "hai") %>%
  group_by(participant_id) %>%
  summarise(
    value = exp(mean(log(value))),
    .groups = "drop"
  ) %>%
  arrange(participant_id)

saveRDS(task1.1_response_train, file = fs::path(processed_data_path, "task1.1_response_train.rds"))

saveRDS(task1.2_response_train, file = fs::path(processed_data_path, "task1.2_response_train.rds"))

saveRDS(task1.3_response_train, file = fs::path(processed_data_path, "task1.3_response_train.rds"))

# saveRDS(task1.4_response_train, file = fs::path(processed_data_path, "task1.4_response_train.rds"))

saveRDS(task2.1_response_train, file = fs::path(processed_data_path, "task2.1_response_train.rds"))

saveRDS(task2.2_response_train, file = fs::path(processed_data_path, "task2.2_response_train.rds"))

saveRDS(task2.3_response_train, file = fs::path(processed_data_path, "task2.3_response_train.rds"))

rm(list = ls())
