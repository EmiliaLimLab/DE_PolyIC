#-------------------------------------------------------------------------------
# Deidentifying transcriptomic data 
#-------------------------------------------------------------------------------
metadata_GSC_2684 <- read.delim("data/List_for_sequencing_Carlsten_GSC-2684.csv", sep = ";", header = TRUE) %>%
  dplyr::select(External_ID, library_ID)
metadata_GSC_2684$External_ID <- metadata_GSC_2684$External_ID %>%
  gsub("_[0-9]{2}-[A-Za-z]+-[0-9]{4}$", "", .) %>%  # remove trailing portion
  gsub("_", ".", .) %>%                             # replace all underscores with periods
  gsub("-", "_", .)

# De-identify patients 
patient_deid_map <- read.csv("data/patient_ID_deidentified.csv") %>%
  mutate(deidentified = as.character(deidentified))
subject_prefix <- sub("\\..*$", "", metadata_GSC_2684$External_ID)
id_suffix <- sub("^[^.]*", "", metadata_GSC_2684$External_ID)
metadata_GSC_2684$External_ID <- paste0(patient_deid_map$deidentified[match(subject_prefix, patient_deid_map$patient_ID)], id_suffix)
write.csv(metadata_GSC_2684, "data/List_for_sequencing_Carlsten_GSC-2684_deidentified.csv", row.names = FALSE)
# 8 samples have NA as their deidentifier (corresponds to one patient's samples that had failed)
# 10 samples in total failed (excluding RNAseq_negative_control)
# 
#-------------------------------------------------------------------------------
# Deidentifying proteomic data 
#-------------------------------------------------------------------------------
cytokine_panel_samples_EveTechnologies <- read_tsv("/projects/lim_lab_scratch/taysia/atac-seq/proteomics/GSC2684_atac-seq_samples_EveTechnologies.tsv")

patient_deid_map <- read.csv("data/patient_ID_deidentified.csv") %>%
  mutate(deidentified = as.character(deidentified))
subject_prefix <- sub("\\..*$", "", cytokine_panel_samples_EveTechnologies$`Patient ID`)
id_suffix <- sub("^[^.]*", "", cytokine_panel_samples_EveTechnologies$`Patient ID`)
cytokine_panel_samples_EveTechnologies$deidentified_Patient_ID <- paste0(patient_deid_map$deidentified[match(subject_prefix, patient_deid_map$patient_ID)], id_suffix)

cytokine_panel_samples_EveTechnologies <- cytokine_panel_samples_EveTechnologies %>%
  mutate(`Library external Identifier` = paste0(
    deidentified_Patient_ID,
    sub("^[^_]+", "", `Library external Identifier`)
  ))
cytokine_panel_samples_EveTechnologies <- cytokine_panel_samples_EveTechnologies %>% select(-c(`Patient ID`, Notes))
cytokine_panel_samples_EveTechnologies %>% view()

write_tsv(cytokine_panel_samples_EveTechnologies, "data/cytokine_panel_labels.tsv")

#-------------------------------------------------------------------------------
# Deidentifying ATAC data 
#-------------------------------------------------------------------------------
DE_FA_M_1day <- read_csv("/projects/lim_lab_scratch/taysia/GSC-2684_updated/manuscript/data/ATAC-seq_DE_FA-M_1day_frip_filtered_samplesheet.csv")
patient_deid_map

DE_FA_M_1day <- DE_FA_M_1day %>%
  left_join(patient_deid_map, by = c("Factor" = "patient_ID")) %>%
  mutate(Factor = as.character(deidentified)) %>%
  dplyr::select(-deidentified)
write_csv(DE_FA_M_1day, "/projects/lim_lab_scratch/taysia/GSC-2684_updated/manuscript/data/ATAC-seq_DE_FA-M_1day_frip_filtered_samplesheet_deidentified.csv")

FA_P_FA_M_1day <- read_csv("/projects/lim_lab_scratch/taysia/GSC-2684_updated/manuscript/data/ATAC-seq_FA-P_FA-M_1day_frip_filtered_samplesheet.csv")
FA_P_FA_M_1day <- FA_P_FA_M_1day %>%
  left_join(patient_deid_map, by = c("Factor" = "patient_ID")) %>%
  mutate(Factor = as.character(deidentified)) %>%
  dplyr::select(-deidentified)
write_csv(FA_P_FA_M_1day, "/projects/lim_lab_scratch/taysia/GSC-2684_updated/manuscript/data/ATAC-seq_FA-P_FA-M_1day_frip_filtered_samplesheet_deidentified.csv")
