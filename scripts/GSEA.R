out_dir <- "output/GSEA"

#-------------------------------------------------------------------------------
# Load libraries
#-------------------------------------------------------------------------------
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyverse)
library(clusterProfiler)
library(stringr)
library(colorspace)
library(tidytext) 
#-------------------------------------------------------------------------------
# Read files
#-------------------------------------------------------------------------------
deseq_DE_FA_M_1day_lfcSE_cutoff <- read_tsv("output/DEA/deseq_DE_FA_M_1day_lfcSE_cutoff.tsv")
deseq_DE_FA_P_1day_lfcSE_cutoff <- read_tsv("output/DEA/deseq_DE_FA_P_1day_lfcSE_cutoff.tsv")
deseq_DE_FA_P_1week_lfcSE_cutoff <- read_tsv("output/DEA/deseq_DE_FA_P_1week_lfcSE_cutoff.tsv")
deseq_M_P_IC_lfcSE_cutoff <- read_tsv("output/DEA/deseq_M_P_IC_lfcSE_cutoff.tsv")
deseq_M_P_FA_lfcSE_cutoff <- read_tsv("output/DEA/deseq_M_P_FA_lfcSE_cutoff.tsv")
deseq_M_P_DE_lfcSE_cutoff <- read_tsv("output/DEA/deseq_M_P_DE_lfcSE_cutoff.tsv")
deseq_DE_P_1week_1day_lfcSE_cutoff <- read_tsv("output/DEA/deseq_DE_P_1week_1day_lfcSE_cutoff.tsv")

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

DE_FA_M_1day_gene_list_sorted <- create_sorted_gene_list(deseq_DE_FA_M_1day_lfcSE_cutoff)
DE_FA_P_1day_gene_list_sorted <- create_sorted_gene_list(deseq_DE_FA_P_1day_lfcSE_cutoff)
DE_FA_P_1week_gene_list_sorted <- create_sorted_gene_list(deseq_DE_FA_P_1week_lfcSE_cutoff)

M_P_IC_gene_list_sorted <- create_sorted_gene_list(deseq_M_P_IC_lfcSE_cutoff)
M_P_DE_gene_list_sorted <- create_sorted_gene_list(deseq_M_P_DE_lfcSE_cutoff)
M_P_FA_gene_list_sorted <- create_sorted_gene_list(deseq_M_P_FA_lfcSE_cutoff)

DE_P_1week_1day_gene_list_sorted <- create_sorted_gene_list(deseq_DE_P_1week_1day_lfcSE_cutoff)

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

# DE vs FA
DE_FA_M_1day_gene_list_entrez_sorted <- get_entrez_ids(DE_FA_M_1day_gene_list_sorted) %>% sort(decreasing = TRUE)
DE_FA_P_1day_gene_list_entrez_sorted <- get_entrez_ids(DE_FA_P_1day_gene_list_sorted) %>% sort(decreasing = TRUE)
DE_FA_P_1week_gene_list_entrez_sorted <- get_entrez_ids(DE_FA_P_1week_gene_list_sorted) %>% sort(decreasing = TRUE)
# Non unique Entrez IDs must be removed prior to GSEA
DE_FA_M_1day_gene_list_entrez_sorted <- DE_FA_M_1day_gene_list_entrez_sorted[!duplicated(names(DE_FA_M_1day_gene_list_entrez_sorted))]
DE_FA_P_1day_gene_list_entrez_sorted <- DE_FA_P_1day_gene_list_entrez_sorted[!duplicated(names(DE_FA_P_1day_gene_list_entrez_sorted))]
DE_FA_P_1week_gene_list_entrez_sorted <- DE_FA_P_1week_gene_list_entrez_sorted[!duplicated(names(DE_FA_P_1week_gene_list_entrez_sorted))]

# P vs M
M_P_IC_gene_list_entrez_sorted <- get_entrez_ids(M_P_IC_gene_list_sorted) %>% sort(decreasing = TRUE)
M_P_DE_gene_list_entrez_sorted <- get_entrez_ids(M_P_DE_gene_list_sorted) %>% sort(decreasing = TRUE)
M_P_FA_gene_list_entrez_sorted <- get_entrez_ids(M_P_FA_gene_list_sorted) %>% sort(decreasing = TRUE)
# Non unique Entrez IDs must be removed prior to GSEA  
M_P_IC_gene_list_entrez_sorted <- M_P_IC_gene_list_entrez_sorted[!duplicated(names(M_P_IC_gene_list_entrez_sorted))]
M_P_DE_gene_list_entrez_sorted <- M_P_DE_gene_list_entrez_sorted[!duplicated(names(M_P_DE_gene_list_entrez_sorted))]
M_P_FA_gene_list_entrez_sorted <- M_P_FA_gene_list_entrez_sorted[!duplicated(names(M_P_FA_gene_list_entrez_sorted))]

# 1 week vs 1 day
DE_P_1week_1day_gene_list_entrez_sorted <- get_entrez_ids(DE_P_1week_1day_gene_list_sorted) %>% sort(decreasing = TRUE)
DE_P_1week_1day_gene_list_entrez_sorted <- DE_P_1week_1day_gene_list_entrez_sorted[!duplicated(names(DE_P_1week_1day_gene_list_entrez_sorted))]

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

DE_FA_M_1day_gseGO <- go_enrichment(DE_FA_M_1day_gene_list_entrez_sorted, ont = "BP")
DE_FA_P_1day_gseGO <- go_enrichment(DE_FA_P_1day_gene_list_entrez_sorted, ont = "BP")
DE_FA_P_1week_gseGO <- go_enrichment(DE_FA_P_1week_gene_list_entrez_sorted, ont = "BP")

M_P_IC_gseGO <- go_enrichment(M_P_IC_gene_list_entrez_sorted, ont = "BP")
M_P_DE_gseGO <- go_enrichment(M_P_DE_gene_list_entrez_sorted, ont = "BP")
M_P_FA_gseGO <- go_enrichment(M_P_FA_gene_list_entrez_sorted, ont = "BP")

DE_P_1week_1day_gseGO <- go_enrichment(DE_P_1week_1day_gene_list_entrez_sorted, ont = "BP")
#-------------------------------------------------------------------------------
# Saving files
#-------------------------------------------------------------------------------
write_GO_results <- function(deseq_results, out_fn) {
  deseq_results %>% 
    write.table(out_fn, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}

write_GO_results(DE_FA_M_1day_gseGO, paste0(out_dir, "/DE_FA_M_1day_gseGO.tsv"))
write_GO_results(DE_FA_P_1day_gseGO, paste0(out_dir, "/DE_FA_P_1day_gseGO.tsv"))
write_GO_results(DE_FA_P_1week_gseGO, paste0(out_dir, "/DE_FA_P_1week_gseGO.tsv"))
write_GO_results(M_P_IC_gseGO, paste0(out_dir, "/M_P_IC_gseGO.tsv"))
write_GO_results(M_P_DE_gseGO, paste0(out_dir, "/M_P_DE_gseGO.tsv"))
write_GO_results(M_P_FA_gseGO, paste0(out_dir, "/M_P_FA_gseGO.tsv"))
write_GO_results(DE_P_1week_1day_gseGO, paste0(out_dir, "/DE_P_1week_1day_gseGO.tsv"))

DE_FA_M_1day_gseGO <- read.delim(paste0(out_dir, "/DE_FA_M_1day_gseGO.tsv"))
DE_FA_P_1day_gseGO <- read.delim(paste0(out_dir, "/DE_FA_P_1day_gseGO.tsv"))
DE_FA_P_1week_gseGO <- read.delim(paste0(out_dir, "/DE_FA_P_1week_gseGO.tsv"))
M_P_FA_gseGO <- read.delim(paste0(out_dir, "/M_P_FA_gseGO.tsv"))
M_P_DE_gseGO <- read.delim(paste0(out_dir, "/M_P_DE_gseGO.tsv"))
M_P_IC_gseGO <- read.delim(paste0(out_dir, "/M_P_IC_gseGO.tsv"))
#-------------------------------------------------------------------------------
# GSEA GO Plots
#-------------------------------------------------------------------------------
df_list_DEvsFA <- list(DE_FA_M_1day_gseGO, 
                       DE_FA_P_1day_gseGO,      
                       DE_FA_P_1week_gseGO)  

# Create titles for each facet
titles_DEvsFA <- c("Media Control, 1 Day\nPost-Exposure",
                   "Poly(I:C), 1 Day\nPost-Exposure", 
                   "Poly(I:C), 1 Week\nPost-Exposure")

plot_go_enrichment_faceted_top <- function(df_list, titles,
                                           title_text = "GO Enrichment Analysis",
                                           caption_text = NULL,
                                           n_terms = 8,
                                           padj_cutoff = 0.05) {
  
  combined_df <- bind_rows(
    df_list[[1]] %>% mutate(dataset = titles[1]),
    df_list[[2]] %>% mutate(dataset = titles[2]),
    df_list[[3]] %>% mutate(dataset = titles[3])
  ) %>%
    mutate(dataset = factor(dataset, levels = titles))
  
  top_pathways <- combined_df %>%
    filter(!is.na(p.adjust), p.adjust <= padj_cutoff) %>%
    group_by(dataset) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = n_terms) %>%
    ungroup() %>%
    mutate(
      Description = str_wrap(Description, width = 65),
      Description_panel = reorder_within(
        Description,
        gene_ratio,
        interaction(dataset, direction, sep = "__")
      )
    )
  
  ggplot(
    top_pathways,
    aes(x = gene_ratio, y = Description_panel, fill = p.adjust, size = setSize)
  ) +
    geom_point(pch = 21, color = "black") +
    facet_grid(dataset ~ direction, scales = "free_y", space = "free_y") +
    scale_y_reordered() +
    scale_x_continuous(
      limits = c(0, 0.80),
      breaks = seq(0, 0.8, by = 0.1)
    ) +
    scale_fill_continuous_divergingx(
      palette = "RdBu",
      mid = 0.05,
      limits = c(0, 0.05),
      oob = scales::squish
    ) +
    scale_size_continuous(
      range = c(4, 14),
      breaks = c(50, 100, 150, 200)
    ) +
    labs(
      x = "Gene ratio",
      y = "Description",
      fill = "*P*<sub>adj</sub>",
      size = "Set size",
      title = title_text,
      subtitle = caption_text
    ) +
    theme_bw() +
    theme(
      text = element_text(size = 20),
      plot.title = element_text(size = 25, hjust = 0.5, face = "bold",
                                margin = margin(b = 15)),
      axis.text.y = ggtext::element_markdown(size = 24, lineheight = 0.5,
                                             margin = margin(l = 10)),
      axis.text.x = element_text(size = 24),
      axis.title.x = element_text(size = 22),
      axis.title.y = element_text(size = 24),
      legend.title = ggtext::element_markdown(size = 20),
      legend.text = element_text(size = 20),
      strip.text.x = element_text(face = "bold", size = 24),
      strip.text.y = element_text(face = "bold", size = 18, angle = 270),
      panel.spacing.x = unit(0.4, "cm"),
      panel.spacing.y = unit(0.2, "cm"),
      plot.subtitle = element_text(size = 16, hjust = 0.5,
                                   margin = margin(b = 5))
    )
}

go_enrichment_DEvsFA <- plot_go_enrichment_faceted_top(df_list_DEvsFA, titles_DEvsFA, 
                                                       title_text = "GO Enrichment Analysis:\nDE vs FA Samples Across Conditions",
                                                       caption_text = "Top activated and suppressed GO terms by adjusted p-value")

ggsave(
  filename = paste0(out_dir, "/go_enrichment_DEvsFA_plot.png"),
  plot = go_enrichment_DEvsFA,
  width = 28,
  height = 12,
  dpi = 300
)


df_list_MP <- list(M_P_FA_gseGO,
                   M_P_DE_gseGO,
                   M_P_IC_gseGO) 

# Create titles for each facet
titles_MP <- c("FA-Exposed Samples,\n1 Day Post-Exposure",
               "DE-Exposed Samples,\n1 Day Post-Exposure",
               "Incubator Control")

# Generate the faceted plot
go_enrichment_PvsM <- plot_go_enrichment_faceted_top(df_list_MP, titles_MP, 
                                                     title_text = "GO Enrichment Analysis:\nPolyI:C vs Control Media Samples Across Conditions", 
                                                     caption_text = "Top activated and suppressed GO terms by adjusted p-value")


ggsave(
  filename = paste0(out_dir, "/go_enrichment_PvsM_plot.png"),
  plot = go_enrichment_PvsM,
  width = 28,
  height = 12,
  dpi = 300
)

