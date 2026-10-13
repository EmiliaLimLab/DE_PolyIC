out_dir <- "output/ATAC/DEA_GSEA_patients_with_ATAC_samples"

#-------------------------------------------------------------------------------
# Load libraries
#-------------------------------------------------------------------------------
library(tidyverse)
library(DESeq2) 
library(ggtext)
library(patchwork)
library(ggplot2)
library(biomaRt)
library(ggbeeswarm)
library(scales)
library(ggrepel)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(stringr)
library(colorspace)
library(tidytext)
#-------------------------------------------------------------------------------
# Read files
#-------------------------------------------------------------------------------
metadata_GSC_2684 <- read_csv("data/List_for_sequencing_Carlsten_GSC-2684_deidentified.csv")

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

columns_to_keep <- c(
  "gene_name",  # Always keep gene_name
  "4.FA_M.1day", 
  "4.FA_P.1day", 
  "2.FA_M.1day", 
  "2.FA_P.1day",  
  "4.DE_M.1day", 
  "2.DE_M.1day"
)

patient_2_and_4 <- gene_counts_length_scaled_GSC_2684 %>%
  dplyr::select(all_of(columns_to_keep))

#-------------------------------------------------------------------------------
# Data wrangling
#-------------------------------------------------------------------------------
extract_sample_metadata <- function(df) {
  df %>%
    dplyr::select(-gene_name) %>%
    colnames() %>%
    as_tibble_col(column_name = "sample") %>%
    separate_wider_delim(
      cols = sample,
      delim = ".",
      names = c("subject", "condition", "time"),
      cols_remove = FALSE
    ) %>%
    column_to_rownames("sample")
}

#create combined condition variable: 
patient_2_and_4_factors <- extract_sample_metadata(patient_2_and_4)
patient_2_and_4_factors$group <- paste(patient_2_and_4_factors$condition, patient_2_and_4_factors$time, sep = "_")

#-------------------------------------------------------------------------------
# Differential expression analysis
#-------------------------------------------------------------------------------
# Function definition: create DESEQDataSet (DDS) object
# * Note: Length-scaled gene counts are not integers, hence need to round

# Cell-type proportions are NOT included as covariates in the DESeq2 design.
# This subset has only 6 samples (2 subjects x 3 groups). 
# The design ~ subject + group already uses 4 parameters (intercept, subject, 2 group terms), 
# leaving 2 residual degrees of freedom for dispersion estimation.
# Adding basal, multiciliated and secretory proportions would require 7 parameters for 6 samples, 
# giving a rank-deficient model matrix that DESeq2 cannot fit.

create_dds <- function(counts_matrix, metadata) {
  dds <- DESeqDataSetFromMatrix(countData=round(counts_matrix), 
                                colData=metadata, 
                                ~ subject + group)
  
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

# Run DESeq: 
deseq_patient_2_and_4 <- DESeq(create_dds(dplyr::select(patient_2_and_4, -gene_name), patient_2_and_4_factors))

deseq_FA_PvsFA_M_patient_2_and_4 <- results(deseq_patient_2_and_4, contrast=c("group", "FA_P_1day", "FA_M_1day")) %>% 
  data.frame() %>% 
  arrange(padj) 

deseq_DE_MvsFA_M_patient_2_and_4 <- results(deseq_patient_2_and_4, contrast=c("group", "DE_M_1day", "FA_M_1day")) %>% 
  data.frame() %>% 
  arrange(padj)

#-------------------------------------------------------------------------------
# Updating gene names
#-------------------------------------------------------------------------------
# Updating outdatted gene IDs that now have been given gene names 
updated_ensembl_gene_ids <- read.delim("data/updated_ensembl_id_gene_names.tsv")

format_deseq_results <- function(deseq_results, count_data) {
  # add gene_id only if it doesn't already exist
  if (!"gene_id" %in% colnames(deseq_results)) {
    deseq_results <- deseq_results %>% rownames_to_column(var = "gene_id")
  }
  if (!"gene_id" %in% colnames(count_data)) {
    count_data <- count_data %>% rownames_to_column(var = "gene_id")
  }
  
  deseq_results %>% 
    left_join(
      count_data %>% dplyr::select(gene_id, gene_name),
      by = "gene_id"
    ) %>% 
    relocate(gene_name, .after = gene_id)
}

rename_genes <- function (data, updated_gene_names){
  gene_counts_updated <- data %>%
    left_join(updated_gene_names, by = c("gene_id" = "Gene.stable.ID")) %>%
    mutate(gene_name = ifelse(is.na(Gene.name) | Gene.name == "", gene_name, Gene.name)) %>%
    dplyr::select(-Gene.description, -Gene.name)
  return(gene_counts_updated)
}

deseq_FA_PvsFA_M_patient_2_and_4 <- deseq_FA_PvsFA_M_patient_2_and_4 %>%
  format_deseq_results(patient_2_and_4) %>%
  rename_genes(updated_ensembl_gene_ids)

deseq_DE_MvsFA_M_patient_2_and_4 <- deseq_DE_MvsFA_M_patient_2_and_4 %>%
  format_deseq_results(patient_2_and_4) %>%
  rename_genes(updated_ensembl_gene_ids)
#-------------------------------------------------------------------------------
# Filtering out high lfcSE genes
#-------------------------------------------------------------------------------
deseq_FA_PvsFA_M_patient_2_and_4_lfcSE_cutoff <- deseq_FA_PvsFA_M_patient_2_and_4 %>% filter(lfcSE < 1 )
deseq_DE_MvsFA_M_patient_2_and_4_lfcSE_cutoff <- deseq_DE_MvsFA_M_patient_2_and_4 %>% filter(lfcSE < 1 )

#-------------------------------------------------------------------------------
# Saving DESeq filtered files
#-------------------------------------------------------------------------------
write_deseq_results <- function(deseq_results, out_fn) {
  deseq_results %>% 
    write.table(out_fn, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}

write_deseq_results(deseq_FA_PvsFA_M_patient_2_and_4_lfcSE_cutoff, paste0(out_dir, "/deseq_FA_PvsFA_M_patient_2_and_4_lfcSE_cutoff.tsv"))
write_deseq_results(deseq_DE_MvsFA_M_patient_2_and_4_lfcSE_cutoff, paste0(out_dir, "/deseq_DE_MvsFA_M_patient_2_and_4_lfcSE_cutoff.tsv"))

#-------------------------------------------------------------------------------
# Gene set enrichment analysis (GSEA) using gene lists ranked by DESeq2 log2FC
#-------------------------------------------------------------------------------
# Create gene lists inputs and sort in decreasing order
create_sorted_gene_list <- function(df) {
  gene_list <- df$log2FoldChange
  names(gene_list) <- df$gene_id
  gene_list <- na.omit(gene_list)
  gene_list_sorted <- sort(gene_list, decreasing = TRUE)
  return(gene_list_sorted)
}

FA_PvsFA_M_patient_2_and_4_gene_list_sorted <- create_sorted_gene_list(deseq_FA_PvsFA_M_patient_2_and_4_lfcSE_cutoff)
DE_MvsFA_M_patient_2_and_4_gene_list_sorted <- create_sorted_gene_list(deseq_DE_MvsFA_M_patient_2_and_4_lfcSE_cutoff)

# remove NA Entrez IDs
get_entrez_ids <- function(ensembl_ids) {
  deseq_DE_gene_list_entrez <- mapIds(org.Hs.eg.db,
                                      keys = names(ensembl_ids),
                                      column = "ENTREZID",
                                      keytype = "ENSEMBL",
                                      multiVals = "first") %>%
    enframe() %>%
    dplyr::filter(!is.na(value)) %>% 
    left_join(enframe(ensembl_ids), by = "name") %>%
    dplyr::select(-name) %>%
    deframe()
  
  return(deseq_DE_gene_list_entrez)
}

FA_PvsFA_M_patient_2_and_4_gene_list_entrez_sorted <- get_entrez_ids(FA_PvsFA_M_patient_2_and_4_gene_list_sorted) %>% sort(decreasing = TRUE)
DE_MvsFA_M_patient_2_and_4_gene_list_entrez_sorted <- get_entrez_ids(DE_MvsFA_M_patient_2_and_4_gene_list_sorted) %>% sort(decreasing = TRUE)

# Non unique Entrez IDs must be removed prior to GSEA
FA_PvsFA_M_patient_2_and_4_gene_list_entrez_sorted <- FA_PvsFA_M_patient_2_and_4_gene_list_entrez_sorted[!duplicated(names(FA_PvsFA_M_patient_2_and_4_gene_list_entrez_sorted))]
DE_MvsFA_M_patient_2_and_4_gene_list_entrez_sorted <- DE_MvsFA_M_patient_2_and_4_gene_list_entrez_sorted[!duplicated(names(DE_MvsFA_M_patient_2_and_4_gene_list_entrez_sorted))]

go_enrichment <- function(sorted_entrez_gene_list, 
                          ont = "BP", 
                          pvalueCutoff = 1,
                          minGSSize = 3,
                          maxGSSize = 800) {
  
  # Run GSEA for GO
  set.seed(1234)
  res <- gseGO(geneList = sorted_entrez_gene_list,
               OrgDb = org.Hs.eg.db,
               ont = ont,
               keyType = "ENTREZID",
               pvalueCutoff = pvalueCutoff,
               minGSSize = minGSSize,
               maxGSSize = maxGSSize,
               verbose = TRUE)
  
  # Convert to data frame
  res_df <- as.data.frame(res@result)
  
  # If there are results, add gene_ratio and direction
  if (nrow(res_df) > 0) {
    res_df <- res_df %>%
      mutate(gene_ratio = (str_count(core_enrichment, "/") + 1) / setSize,
             direction = ifelse(NES < 0, "Suppressed", "Activated"))
  }
  
  return(res_df)
}

FA_PvsFA_M_patient_2_and_4_gseGO <- go_enrichment(FA_PvsFA_M_patient_2_and_4_gene_list_entrez_sorted, ont = "BP")
DE_MvsFA_M_patient_2_and_4_gseGO <- go_enrichment(DE_MvsFA_M_patient_2_and_4_gene_list_entrez_sorted, ont = "BP")

#-------------------------------------------------------------------------------
# Saving files
#-------------------------------------------------------------------------------
write_GO_results <- function(deseq_results, out_fn) {
  deseq_results %>% 
    write.table(out_fn, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}

write_GO_results(FA_PvsFA_M_patient_2_and_4_gseGO, paste0(out_dir, "/FA_PvsFA_M_patient_2_and_4_gseGO.tsv"))
write_GO_results(DE_MvsFA_M_patient_2_and_4_gseGO, paste0(out_dir, "/DE_MvsFA_M_patient_2_and_4_gseGO.tsv"))
