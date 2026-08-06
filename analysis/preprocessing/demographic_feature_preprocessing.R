# ==============================================================================
# Libraries
# ==============================================================================

library(tidyverse)
library(fastDummies)
library(janitor)


# ==============================================================================
# Paths
# ==============================================================================

raw_data_path <- fs::path("data-raw")
processed_data_path <- fs::path("data")


# ==============================================================================
# Load raw data
# ==============================================================================

demographics_all <- read.delim(
  fs::path(raw_data_path, "participants.tsv")
)

vaccine_history_all <- read.delim(
  fs::path(raw_data_path, "participant_vaccineHistory.tsv")
)

hla_all <- read.delim(
  fs::path(raw_data_path, "participant_hla.tsv")
) %>%
  mutate(
    participant_id = paste0(study_accession, ".", subject)
  )


# ==============================================================================
# Define train/test datasets
# ==============================================================================

demographics_train <- demographics_all %>%
  filter(study_accession != "2025LJI")

demographics_test <- demographics_all %>%
  filter(study_accession == "2025LJI")


train_ids <- unique(demographics_train$subject)
test_ids <- unique(demographics_test$subject)


vaccine_history_train <- vaccine_history_all %>%
  filter(subject %in% train_ids)

vaccine_history_test <- vaccine_history_all %>%
  filter(subject %in% test_ids)


hla_train <- hla_all %>%
  filter(subject %in% train_ids)

hla_test <- hla_all %>%
  filter(subject %in% test_ids)



# ==============================================================================
# Demographics processing
# ==============================================================================

process_demographics <- function(data) {
  
  data %>%
    select(
      participant_id,
      subject,
      study_accession,
      biological_sex,
      race,
      age,
      geolocation
    ) %>%
    fastDummies::dummy_cols(
      select_columns = c(
        "biological_sex",
        "race",
        "geolocation"
      ),
      remove_selected_columns = TRUE,
      ignore_na = TRUE
    )
}


demographics_train_encoded <- process_demographics(demographics_train)
demographics_test_encoded <- process_demographics(demographics_test)


demographic_features <- intersect(
  names(demographics_train_encoded),
  names(demographics_test_encoded)
)


demographics_train_df <- demographics_train_encoded %>%
  select(
    participant_id,
    age,
    all_of(setdiff(demographic_features,
                   c("participant_id",
                     "subject",
                     "study_accession")))
  )


demographics_test_df <- demographics_test_encoded %>%
  select(
    participant_id,
    age,
    all_of(setdiff(demographic_features,
                   c("participant_id",
                     "subject",
                     "study_accession")))
  )



# ==============================================================================
# Vaccine history processing
# ==============================================================================

process_vaccine_history <- function(data) {
  
  data %>%
    mutate(
      feature = paste0(
        "vacc_",
        vaccine_year,
        "_",
        received
      ),
      value = 1L
    ) %>%
    select(subject, feature, value) %>%
    distinct() %>%
    pivot_wider(
      id_cols = subject,
      names_from = feature,
      values_from = value,
      values_fill = 0
    )
}


vaccine_history_train_df <- process_vaccine_history(vaccine_history_train)

vaccine_history_test_df <- process_vaccine_history(vaccine_history_test)


vaccine_history_features <- intersect(
  names(vaccine_history_train_df),
  names(vaccine_history_test_df)
)


vaccine_history_train_df <- vaccine_history_train_df %>%
  select(subject, all_of(vaccine_history_features))

vaccine_history_test_df <- vaccine_history_test_df %>%
  select(subject, all_of(vaccine_history_features))


# ==============================================================================
# HLA processing
# ==============================================================================

process_hla <- function(data) {
  
  data %>%
    select(
      participant_id,
      locus_name,
      allele_1,
      allele_2
    ) %>%
    pivot_wider(
      id_cols = participant_id,
      names_from = locus_name,
      values_from = c(allele_1, allele_2),
      names_glue = "{locus_name}_{.value}"
    ) %>%
    fastDummies::dummy_cols(
      select_columns = setdiff(
        names(.),
        "participant_id"
      ),
      remove_selected_columns = TRUE,
      ignore_na = TRUE
    )
}


hla_train_df <- process_hla(hla_train)
hla_test_df <- process_hla(hla_test)


hla_features <- intersect(
  names(hla_train_df),
  names(hla_test_df)
)


hla_train_df <- hla_train_df %>%
  select(all_of(hla_features))


hla_test_df <- hla_test_df %>%
  select(all_of(hla_features))



# ==============================================================================
# Merge all feature blocks
# ==============================================================================

demographics_train_df <- demographics_train_df %>%
  left_join(
    hla_train_df,
    by = "participant_id"
  ) %>%
  clean_names()


demographics_test_df <- demographics_test_df %>%
  left_join(
    hla_test_df,
    by = "participant_id"
  ) %>%
  clean_names()



# ==============================================================================
# Final feature selection
# ==============================================================================

final_features <- intersect(
  names(demographics_train_df),
  names(demographics_test_df)
)


demographics_train_df <- demographics_train_df %>%
  select(all_of(final_features)) %>%
  mutate(
    across(
      where(is.numeric),
      ~replace_na(.x, 0)
    )
  )


demographics_test_df <- demographics_test_df %>%
  select(all_of(final_features)) %>%
  mutate(
    across(
      where(is.numeric),
      ~replace_na(.x, 0)
    )
  )



# ==============================================================================
# Save processed data
# ==============================================================================

saveRDS(
  demographics_train_df,
  fs::path(processed_data_path, "demographics_train_df.rds")
)

saveRDS(
  demographics_test_df,
  fs::path(processed_data_path, "demographics_test_df.rds")
)

rm(list = ls())