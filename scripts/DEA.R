out_dir <- "output/DEA"

#-------------------------------------------------------------------------------
# Load libraries
#-------------------------------------------------------------------------------
library(DESeq2)
library(tidyverse) 

#-------------------------------------------------------------------------------
# Read files
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

update_rna_seq_ids <- function(data, metadata){
  rename_map <- setNames(metadata$External_ID, metadata$library_ID)
  colnames(data) <- ifelse(colnames(data) %in% names(rename_map), 
                           rename_map[colnames(data)],
                           colnames(data))
  return (data)
}

gene_counts_length_scaled_GSC_2684 <- read.delim("data/salmon.merged.gene_counts_length_scaled.tsv",
                                                 check.names=FALSE) %>% column_to_rownames("gene_id")
gene_counts_length_scaled_GSC_2684 <- update_rna_seq_ids(gene_counts_length_scaled_GSC_2684, metadata_GSC_2684)

cibersort_results <- read.csv("data/CIBERSORTx_cell_proportions.csv")

cibersort_meta <- cibersort_results %>%
  left_join(metadata_GSC_2684, by = c("Mixture" = "library_ID")) %>%
  dplyr::select(External_ID, basal_prop = Basal, multiciliated_prop = Multiciliated.lineage, secretory_prop = Secretory) %>%
  distinct()

#-------------------------------------------------------------------------------
# Data wrangling
#-------------------------------------------------------------------------------
# removing 2.DE_P.1week and 1.DE_M.1week (not used in downstream analysis); IC samples kept
DE_FA_with_IC <- gene_counts_length_scaled_GSC_2684 %>% dplyr::select(!matches("2\\.DE_P\\.1week|1\\.DE_M\\.1week"))

extract_sample_metadata <- function(df, prop_df) {
  sample_df <- df %>%
    dplyr::select(-dplyr::any_of("gene_name")) %>%
    colnames() %>%
    as_tibble_col(column_name = "External_ID") %>%
    tidyr::separate_wider_delim(
      cols = External_ID,
      delim = ".",
      names = c("subject", "condition", "time"),
      cols_remove = FALSE
    ) %>%
    dplyr::left_join(
      prop_df,            
      by = "External_ID"
    )
  
  sample_df %>%
    column_to_rownames("External_ID")
}

#create combined condition variable: 
DE_FA_with_IC_factors <- extract_sample_metadata(DE_FA_with_IC, cibersort_meta)
DE_FA_with_IC_factors$group <- paste(DE_FA_with_IC_factors$condition, DE_FA_with_IC_factors$time, sep = "_")

#-------------------------------------------------------------------------------
# Differential expression analysis
#-------------------------------------------------------------------------------
# Function definition: create DESEQDataSet (DDS) object
# * Note: Length-scaled gene counts are not integers, hence need to round
create_dds <- function(counts_matrix, metadata) {
  dds <- DESeqDataSetFromMatrix(countData=round(counts_matrix), 
                                colData=metadata, 
                                ~ basal_prop + multiciliated_prop + secretory_prop + subject + group)
  
  # Pre-filtering low count genes
  # * Recommended starting values from developer: X=10, Y=smallest_group_sample_size
  # * May want to adjust based on data
  X <- 10
  Y <- metadata %>% 
    group_by(condition) %>% 
    summarise(n = n()) %>%
    pull(n) %>% 
    min()
  keep <- rowSums(counts(dds) >= X) >= Y
  dds_filtered_lc <- dds[keep, ]
  
  return(dds_filtered_lc)
}

deseq_DE_FA_with_IC <- DESeq(create_dds(dplyr::select(DE_FA_with_IC, -gene_name), DE_FA_with_IC_factors))

deseq_DE_FA_1day_M_with_IC <- results(deseq_DE_FA_with_IC, contrast=c("group", "DE_M_1day", "FA_M_1day")) %>% 
  data.frame() %>% 
  arrange(padj)

deseq_DE_FA_1day_P_with_IC <- results(deseq_DE_FA_with_IC, contrast=c("group", "DE_P_1day", "FA_P_1day")) %>% 
  data.frame() %>% 
  arrange(padj)

deseq_DE_FA_1week_P_with_IC <- results(deseq_DE_FA_with_IC, contrast=c("group", "DE_P_1week", "FA_P_1week")) %>% 
  data.frame() %>% 
  arrange(padj)

deseq_M_P_IC <- results(deseq_DE_FA_with_IC, contrast=c("group", "IC_P_1day", "IC_M_1day")) %>% 
  data.frame() %>% 
  arrange(padj)

deseq_FA_1day_M_P_with_IC <- results(deseq_DE_FA_with_IC, contrast=c("group", "FA_P_1day", "FA_M_1day")) %>% 
  data.frame() %>% 
  arrange(padj) 

deseq_DE_1day_M_P_with_IC <- results(deseq_DE_FA_with_IC, contrast=c("group", "DE_P_1day", "DE_M_1day")) %>% 
  data.frame() %>% 
  arrange(padj)

deseq_DE_1week_1day_with_IC <- results(deseq_DE_FA_with_IC, contrast = c("group", "DE_P_1week", "DE_P_1day")) %>%
  data.frame() %>%
  arrange(padj)

#-------------------------------------------------------------------------------
# Write DESeq2 results to file
#-------------------------------------------------------------------------------
# Create ouptut directory if it doesn't exist
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

write_deseq_results <- function(deseq_results, out_fn, count_data = NULL) {
  if (!"gene_id" %in% colnames(deseq_results)) {
    deseq_results <- deseq_results %>% rownames_to_column(var = "gene_id")
  }
  if (!is.null(count_data) && !"gene_name" %in% colnames(deseq_results)) {
    if (!"gene_id" %in% colnames(count_data)) {
      count_data <- count_data %>% rownames_to_column(var = "gene_id")
    }
    deseq_results <- deseq_results %>%
      left_join(dplyr::select(count_data, gene_id, gene_name), by = "gene_id") %>%
      relocate(gene_name, .after = gene_id)
  }
  write.table(deseq_results, out_fn, row.names = FALSE, col.names = TRUE,
              sep = "\t", quote = FALSE)
}

# Write tables
write_deseq_results(deseq_DE_FA_1day_M_with_IC, paste0(out_dir, "/GSC-2684_deseq_DE_FA_1day_M_IC_celltype_props_in_model.tsv"), DE_FA_with_IC)
write_deseq_results(deseq_DE_FA_1day_P_with_IC, paste0(out_dir, "/GSC-2684_deseq_DE_FA_1day_P_IC_celltype_props_in_model.tsv"), DE_FA_with_IC)
write_deseq_results(deseq_DE_FA_1week_P_with_IC, paste0(out_dir, "/GSC-2684_deseq_DE_FA_1week_P_IC_celltype_props_in_model.tsv"), DE_FA_with_IC)

write_deseq_results(deseq_M_P_IC, paste0(out_dir, "/GSC-2684_deseq_M_P_IC_celltype_props_in_model.tsv"), DE_FA_with_IC)
write_deseq_results(deseq_DE_1day_M_P_with_IC, paste0(out_dir, "/GSC-2684_deseq_DE_1day_M_P_IC_celltype_props_in_model.tsv"), DE_FA_with_IC)
write_deseq_results(deseq_FA_1day_M_P_with_IC, paste0(out_dir, "/GSC-2684_deseq_FA_1day_M_P_IC_celltype_props_in_model.tsv"), DE_FA_with_IC)

write_deseq_results(deseq_DE_1week_1day_with_IC, paste0(out_dir, "/GSC-2684_deseq_DE_P_1week_vs_1day_celltype_props_in_model.tsv"), DE_FA_with_IC)
#-------------------------------------------------------------------------------
# Updating gene names
#-------------------------------------------------------------------------------
# Updating outdatted gene IDs that now have been given gene names 
updated_ensembl_gene_ids <- read.delim("data/updated_ensembl_id_gene_names.tsv")

rename_genes <- function (data, updated_gene_names){
  gene_counts_updated <- data %>%
    left_join(updated_gene_names, by = c("gene_id" = "Gene.stable.ID")) %>%
    mutate(gene_name = ifelse(is.na(Gene.name) | Gene.name == "", gene_name, Gene.name)) %>%
    dplyr::select(-Gene.description, -Gene.name)
  return(gene_counts_updated)
}

deseq_DE_FA_M_1day <- rename_genes(read.delim(file.path(out_dir, "GSC-2684_deseq_DE_FA_1day_M_IC_celltype_props_in_model.tsv")), updated_ensembl_gene_ids)
deseq_DE_FA_P_1day <- rename_genes(read.delim(file.path(out_dir, "GSC-2684_deseq_DE_FA_1day_P_IC_celltype_props_in_model.tsv")), updated_ensembl_gene_ids)
deseq_DE_FA_P_1week <- rename_genes(read.delim(file.path(out_dir, "GSC-2684_deseq_DE_FA_1week_P_IC_celltype_props_in_model.tsv")), updated_ensembl_gene_ids)
deseq_M_P_IC <- rename_genes(read.delim(file.path(out_dir, "GSC-2684_deseq_M_P_IC_celltype_props_in_model.tsv")), updated_ensembl_gene_ids)
deseq_M_P_DE <- rename_genes(read.delim(file.path(out_dir, "GSC-2684_deseq_DE_1day_M_P_IC_celltype_props_in_model.tsv")), updated_ensembl_gene_ids)
deseq_M_P_FA <- rename_genes(read.delim(file.path(out_dir, "GSC-2684_deseq_FA_1day_M_P_IC_celltype_props_in_model.tsv")), updated_ensembl_gene_ids)
deseq_DE_P_1week_1day <- rename_genes(read.delim(file.path(out_dir, "GSC-2684_deseq_DE_P_1week_vs_1day_celltype_props_in_model.tsv")), updated_ensembl_gene_ids)

#-------------------------------------------------------------------------------
# Filtering out high lfcSE genes
#-------------------------------------------------------------------------------
deseq_DE_FA_M_1day_lfcSE_cutoff <- deseq_DE_FA_M_1day %>% filter(lfcSE < 1 )
deseq_DE_FA_P_1day_lfcSE_cutoff <- deseq_DE_FA_P_1day %>% filter(lfcSE < 1 )
deseq_DE_FA_P_1week_lfcSE_cutoff <- deseq_DE_FA_P_1week %>% filter(lfcSE < 1 )

deseq_M_P_IC_lfcSE_cutoff <- deseq_M_P_IC %>% filter(lfcSE < 1 )
deseq_M_P_DE_lfcSE_cutoff <- deseq_M_P_DE %>% filter(lfcSE < 1 )
deseq_M_P_FA_lfcSE_cutoff <- deseq_M_P_FA %>% filter(lfcSE < 1 )

deseq_DE_P_1week_1day_lfcSE_cutoff <- deseq_DE_P_1week_1day %>% filter(lfcSE < 1 ) %>% tibble()

#-------------------------------------------------------------------------------
# Saving filtered files
#-------------------------------------------------------------------------------
write_deseq_results(deseq_DE_FA_M_1day_lfcSE_cutoff, paste0(out_dir, "/deseq_DE_FA_M_1day_lfcSE_cutoff.tsv"))
write_deseq_results(deseq_DE_FA_P_1day_lfcSE_cutoff, paste0(out_dir, "/deseq_DE_FA_P_1day_lfcSE_cutoff.tsv"))
write_deseq_results(deseq_DE_FA_P_1week_lfcSE_cutoff, paste0(out_dir, "/deseq_DE_FA_P_1week_lfcSE_cutoff.tsv"))
write_deseq_results(deseq_M_P_IC_lfcSE_cutoff, paste0(out_dir, "/deseq_M_P_IC_lfcSE_cutoff.tsv"))
write_deseq_results(deseq_M_P_FA_lfcSE_cutoff, paste0(out_dir, "/deseq_M_P_FA_lfcSE_cutoff.tsv"))
write_deseq_results(deseq_M_P_DE_lfcSE_cutoff, paste0(out_dir, "/deseq_M_P_DE_lfcSE_cutoff.tsv"))
write_deseq_results(deseq_DE_P_1week_1day_lfcSE_cutoff, paste0(out_dir, "/deseq_DE_P_1week_1day_lfcSE_cutoff.tsv"))
