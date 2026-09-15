library(tidyverse)
library(immunarch)

raw_data_path = fs::path("data-raw")
processed_data_path = fs::path("data")

bcr_train <- read.delim(fs::path(raw_data_path, "publicData_bulkBCR.tsv"))
bcr_test <- read.delim(fs::path(raw_data_path, "2025LJI_bulkBCR.tsv"))

bcr_train = bcr_train %>% 
  dplyr::select(participant_id, subject, study_accession, 
                timepoint, ID, tissue,
                sequence_id, sequence, locus, 
                v_call, d_call, j_call, c_call, 
                cdr3, 
                v_identity, d_identity, j_identity, 
                clone_id)

bcr_test = bcr_test %>% 
  dplyr::select(participant_id, subject, study_accession, 
                timepoint, ID, tissue,
                sequence_id, sequence, locus, 
                v_call, d_call, j_call, c_call, 
                cdr3, 
                v_identity, d_identity, j_identity, 
                clone_id) %>% 
  mutate(tissue = toupper(tissue))
