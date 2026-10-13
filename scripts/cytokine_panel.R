out_dir <- "output/cytokine_panel/"

#-------------------------------------------------------------------------------
# library
#-------------------------------------------------------------------------------
library(readxl)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(colorspace) 
library(ggplot2)

#-------------------------------------------------------------------------------
# Data wrangling
#-------------------------------------------------------------------------------
cytokine_panel <- read_xlsx("data/Aug-18-2025-Human-Cytokine-96-Plex-Discovery-Assay-Jha-67551_reformatted.xlsx", skip=1)
cytokine_panel_title <- read_xlsx("data/Aug-18-2025-Human-Cytokine-96-Plex-Discovery-Assay-Jha-67551_reformatted.xlsx")
cytokine_panel_title <- colnames(cytokine_panel_title)
cytokine_names <- cytokine_panel_title[!grepl("^\\.\\.\\.", cytokine_panel_title)]

obs_cols <- grep("Obs Conc", names(cytokine_panel), value = TRUE)
stopifnot(length(obs_cols) == length(cytokine_names)) 

cytokine_panel <- cytokine_panel %>%
  dplyr::select(1, all_of(obs_cols))

eve_tech_sample_names <- read_tsv("data/cytokine_panel_labels.tsv")
#-------------------------------------------------------------------------------
# Dataframe Cleaning
#-------------------------------------------------------------------------------
names(cytokine_panel)[2:(length(cytokine_names)+1)] <- cytokine_names

# replaces OOR (out of range) values, following Eve Technologies guidance
# if OOR >: replaces with the highest value obtained for that analyte
# if OOR <: replaces with half the lowest value obtained for that analyte
# "*" (Eve Technologies extrapolated-value flag) is stripped so those values are kept
# analytes with no numeric values at all become NA (no basis for imputation)
replace_OOR <- function(df) {
  df[, -1] <- lapply(df[, -1], function(col) {
    col_char <- str_remove(as.character(col), "^\\*")
    numeric_col <- suppressWarnings(as.numeric(col_char))
    
    if (all(is.na(numeric_col))) return(numeric_col)
    
    col_char[grep("OOR >", col_char)] <- max(numeric_col, na.rm = TRUE)
    col_char[grep("OOR <", col_char)] <- min(numeric_col, na.rm = TRUE) / 2
    
    suppressWarnings(as.numeric(col_char))
  })
  df
}

cytokine_panel_formatted <- cytokine_panel %>%
  # Join metadata to map labels
  left_join(
    eve_tech_sample_names %>%
      dplyr::select(`Eve Technologies Label ID`, `Library external Identifier`) %>%
      mutate(
        # Remove the date from the end (everything after last underscore + date)
        `Library external Identifier` = str_replace(`Library external Identifier`, "_\\d{2}-[A-Za-z]+-\\d{4}$", "")
      ),
    by = c("Assay Buffer...1" = "Eve Technologies Label ID")
  ) %>%
  # Replace Assay Buffer...1 with cleaned Library external Identifier
  mutate(`Assay Buffer...1` = `Library external Identifier`) %>%
  # Drop the extra column from join
  dplyr::select(-`Library external Identifier`) %>%
  filter(!is.na(`Assay Buffer...1`))                                 

cytokine_panel_formatted <- replace_OOR(cytokine_panel_formatted)    

colnames(cytokine_panel_formatted)[1] <- "Sample_ID"
cytokine_panel_formatted <- cytokine_panel_formatted %>%
  mutate(Group = case_when(
    str_detect(Sample_ID, "FA") ~ "FA",
    str_detect(Sample_ID, "DE") ~ "DE",
    str_detect(Sample_ID, "IC") ~ "IC",
    TRUE ~ NA_character_
  )) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "_1Day$", "_1day")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "_1Week$", "_1week"))

cytokine_panel_log <- cytokine_panel_formatted %>%
  mutate(across(where(is.numeric), ~ log2(. + 1))) # log2 transformation + 1

keep_remove_cols <- function(panel_log, keep, remove){
  panel_log <- panel_log %>% filter(
    str_detect(Sample_ID, keep),  # keep samples
    !str_detect(Sample_ID, remove)        # remove samples
  )
  return (panel_log)
}

samples_IC_1day <- keep_remove_cols(cytokine_panel_log, "-M_1day|-P_1day", "DE|FA")
samples_FA_1day <- keep_remove_cols(cytokine_panel_log, "FA-M_1day|FA-P_1day", "IC")
samples_DE_1day <- keep_remove_cols(cytokine_panel_log, "DE-M_1day|DE-P_1day", "IC")

samples_DE_1day_1week <- cytokine_panel_log %>% # Comparing DE exposed P at day 1 vs week 1
  filter(
    str_detect(Sample_ID, "DE-P_1day|DE-P_1week"),  # keep 1day samples
    !str_detect(Sample_ID, c("IC"))        # remove IC samples
  )

pivot_cytokine <- function(df, id_cols = c("Sample_ID", "Group")) {
  df %>%
    pivot_longer(
      cols = -all_of(id_cols),
      names_to = "Cytokine",
      values_to = "Value"
    )
}

samples_long_IC_1day <- pivot_cytokine(samples_IC_1day)
samples_long_FA_1day <- pivot_cytokine(samples_FA_1day)
samples_long_DE_1day <- pivot_cytokine(samples_DE_1day)
samples_long_DE_1day_1week <- pivot_cytokine(samples_DE_1day_1week)

samples_long_IC_1day <- samples_long_IC_1day %>% mutate(Treatment = case_when(
  str_detect(Sample_ID, "-M_") ~ "Media",
  str_detect(Sample_ID, "-P_") ~ "PolyIC"))
samples_long_IC_1day$Subject <- str_extract(samples_long_IC_1day$Sample_ID, "^[0-9]+")

samples_long_FA_1day <- samples_long_FA_1day %>% mutate(Treatment = case_when(
  str_detect(Sample_ID, "-M_") ~ "Media",
  str_detect(Sample_ID, "-P_") ~ "PolyIC"))
samples_long_FA_1day$Subject <- str_extract(samples_long_FA_1day$Sample_ID,  "^[0-9]+")

samples_long_DE_1day <- samples_long_DE_1day %>% mutate(Treatment = case_when(
  str_detect(Sample_ID, "-M_") ~ "Media",
  str_detect(Sample_ID, "-P_") ~ "PolyIC"))
samples_long_DE_1day$Subject <- str_extract(samples_long_DE_1day$Sample_ID, "^[0-9]+")

samples_long_DE_1day_1week <- samples_long_DE_1day_1week %>% mutate(Treatment = case_when(
  str_detect(Sample_ID, "_1day") ~ "1day",
  str_detect(Sample_ID, "_1week") ~ "1week"))
samples_long_DE_1day_1week$Subject <- str_extract(samples_long_DE_1day_1week$Sample_ID, "^[0-9]+")
#-------------------------------------------------------------------------------
# Media Normalized (to view how PolyIC induces viral response)
#-------------------------------------------------------------------------------
media_normalization <- function(samples_long, log_name){
  stopifnot(!any(is.na(samples_long$Subject)))                       
  
  out <- samples_long %>% filter(Treatment %in% c("Media", "PolyIC")) %>%
    dplyr::select(Subject, Group, Cytokine, Treatment, Value) %>%
    pivot_wider(
      names_from = Treatment,
      values_from = Value
    )
  
  # duplicate wells per subject would create list-columns here
  stopifnot(is.numeric(out$Media), is.numeric(out$PolyIC))
  
  out %>%
    filter(!is.na(Media), !is.na(PolyIC)) %>%
    mutate(!!sym(log_name) := PolyIC - Media) # mutate(log2FC_M_P_DE = PolyIC - Media)
}

#Samples log 2 transformed and then normalized to media 
M_P_IC_media_normalized <- media_normalization(samples_long_IC_1day, log_name = "log2FC_M_P_IC")
M_P_FA_media_normalized <- media_normalization(samples_long_FA_1day, log_name = "log2FC_M_P_FA")
M_P_DE_media_normalized <- media_normalization(samples_long_DE_1day, log_name = "log2FC_M_P_DE")

all_M_P_normalized <- bind_rows(M_P_IC_media_normalized, M_P_FA_media_normalized, M_P_DE_media_normalized)
all_M_P_normalized_long <- all_M_P_normalized %>%
  pivot_longer(
    cols = starts_with("log2FC_"),    # all log2FC columns
    names_to = "Condition",
    values_to = "log2FC"
  ) %>%
  filter(!is.na(log2FC)) 

# Taking the average cytokine level (grouping patients per condition)
all_M_P_summary_normalized_df <- all_M_P_normalized_long %>%
  group_by(Condition, Cytokine) %>%
  summarise(
    mean_log2FC = mean(log2FC, na.rm = TRUE),
    .groups = "drop"
  )

#-------------------------------------------------------------------------------
# Paired Wilcox signed-rank Test
#------------------------------------------------------------------------------
paired_wilcox <- function(df, x, y) {
  df %>%
    group_by(Cytokine) %>%
    summarise(
      n           = sum(!is.na(.data[[x]]) & !is.na(.data[[y]])),
      median_diff = median(.data[[y]] - .data[[x]], na.rm = TRUE),  # y - x, log2 scale
      p_value     = tryCatch(
        suppressWarnings(wilcox.test(.data[[y]], .data[[x]], paired = TRUE)$p.value),
        error = function(e) NA_real_   # e.g. all differences are zero
      ),
      .groups = "drop"
    ) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    arrange(p_value)
}

# PolyIC vs Media within each exposure
wilcox_IC <- paired_wilcox(M_P_IC_media_normalized, "Media", "PolyIC")
wilcox_FA <- paired_wilcox(M_P_FA_media_normalized, "Media", "PolyIC")
wilcox_DE <- paired_wilcox(M_P_DE_media_normalized, "Media", "PolyIC")

# Does DE change the PolyIC response vs FA
de_vs_fa <- inner_join(
  M_P_FA_media_normalized %>% dplyr::select(Subject, Cytokine, FA = log2FC_M_P_FA),
  M_P_DE_media_normalized %>% dplyr::select(Subject, Cytokine, DE = log2FC_M_P_DE),
  by = c("Subject", "Cytokine")
)
wilcox_DE_vs_FA <- paired_wilcox(de_vs_fa, "FA", "DE")

# DE PolyIC: 1 day vs 1 week
day_week <- samples_long_DE_1day_1week %>%
  dplyr::select(Subject, Cytokine, Treatment, Value) %>%
  pivot_wider(names_from = Treatment, values_from = Value)
wilcox_day_week <- paired_wilcox(day_week, "1day", "1week")

#-------------------------------------------------------------------------------
# Heatmap of Poly(I:C) Induction at 1 day Post-exposure (media-normalized)
#-------------------------------------------------------------------------------
heatmap_df <- all_M_P_summary_normalized_df %>%
  pivot_wider(
    names_from = Condition,
    values_from = mean_log2FC
  )
new_order_media_normalized <- c("Cytokine", "log2FC_M_P_IC", "log2FC_M_P_FA", "log2FC_M_P_DE")
heatmap_df <- heatmap_df[, new_order_media_normalized]

colnames(heatmap_df) <- c(
  "Cytokine",
  "IC + PolyI:C",
  "FA + PolyI:C",
  "DE + PolyI:C"
)

heatmap_df_rename <- heatmap_df

plot_cytokine_heatmap_media_normalized <- function(df,
                                                   title = NULL,
                                                   cytokine_order = NULL,
                                                   cluster_cols = FALSE,
                                                   zscore_rows = TRUE,
                                                   show_legend_bool = TRUE,
                                                   colors = viridis(3),
                                                   clustering_distance = "euclidean", 
                                                   clustering_method = "ward.D2") {
  # Ensure numeric columns
  df_numeric <- df %>%
    mutate(across(-Gene, as.numeric))
  
  # Only keep cytokines that exist in df
  if (!is.null(cytokine_order)) {
    cytokine_order <- intersect(cytokine_order, df_numeric$Gene)
  }
  
  # Convert to matrix with Gene as rownames
  expr_matrix <- df_numeric %>%
    column_to_rownames("Gene")
  
  if (!is.null(cytokine_order)) {
    expr_matrix <- expr_matrix[cytokine_order, , drop = FALSE]  # reorder rows
  }
  
  expr_matrix <- as.matrix(expr_matrix)
  
  # Z-score rows if requested
  if (zscore_rows) {
    missing_mask <- is.na(expr_matrix)  
    expr_matrix <- t(scale(t(expr_matrix)))
    expr_matrix[is.na(expr_matrix) | is.infinite(expr_matrix)] <- 0    # zero-variance rows -> 0
    expr_matrix[missing_mask] <- NA     
  }
  
  # Transpose so cytokines are on x-axis
  expr_matrix <- t(expr_matrix)

  col_fun <- colorRamp2(c(-2, 0, 2), colors)
  
  # Plot heatmap
  Heatmap(expr_matrix,
          name = ifelse(zscore_rows, "Z-score", "log2FC"),
          col = col_fun,
          na_col = "grey80", 
          cluster_columns = cluster_cols,
          clustering_distance_columns = clustering_distance,
          clustering_method_columns = clustering_method,  
          column_dend_height = unit(3, "cm"),
          column_dend_gp = gpar(lwd = 2),  
          cluster_rows = FALSE,
          show_column_names = TRUE,
          row_names_gp = gpar(fontsize = 30, fontface = "bold"),
          column_title = title,
          column_names_gp = gpar(fontsize = 18),   
          column_title_gp = gpar(fontsize = 25, fontface = "bold"),
          column_names_rot = -45,
          show_heatmap_legend = show_legend_bool
  )
} 

# renaming proteins to their gene names 
cytokine_96_panel_proteins <- read_tsv("/projects/lim_lab_scratch/taysia/atac-seq/proteomics/cytokine_96_panel_proteins.tsv") %>%
  dplyr::select(-HGNC_symbol) %>% distinct()      
heatmap_df_renamed_UNIPROT <- heatmap_df_rename %>%
  left_join(cytokine_96_panel_proteins, by = "Cytokine") %>%
  dplyr::select(-Cytokine) %>% dplyr::rename(Gene = Protein)  

p1 <- plot_cytokine_heatmap_media_normalized(
  heatmap_df_renamed_UNIPROT, 
  cluster_cols = TRUE,
  zscore_rows = TRUE,
  clustering_distance = "euclidean",
  clustering_method = "ward.D2"
)
draw(p1) 

pdf(paste0(out_dir, "cytokine_heatmap.pdf"), width = 30, height = 8)
draw(p1)
dev.off()

#-------------------------------------------------------------------------------
# Cluster 1 and 2 GO Plots
#-------------------------------------------------------------------------------
cluster1_enrichment <- read_tsv("data/cluster1_enrichment.Process.tsv")
cluster2_enrichment <- read_tsv("data/cluster2_enrichment.Process.tsv")

plot_go_enrichment <- function(df, title, upper_lim) {
  df_clean <- df %>%
    mutate(
      signal_value = signal,
      p.adjust = `false discovery rate`,
      Description = `term description`,
      rank = row_number(),
      is_significant = p.adjust < 0.05
    ) %>%
    arrange(p.adjust) %>%
    slice_head(n = 10) %>%
    mutate(
      Description = factor(Description, levels = rev(unique(Description))),
      gene_count_factor = as.factor(`observed gene count`)
    )
  
  ggplot(df_clean, aes(x = signal_value, y = Description, size = gene_count_factor, fill = p.adjust)) +
    geom_point(shape = 21, color = "black") +
    scale_x_continuous(limits = c(0, upper_lim), breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8)) +
    scale_size_manual(  # Changed to manual scale
      values = setNames(
        seq(5, 20, length.out = length(unique(df_clean$gene_count_factor))), 
        levels(df_clean$gene_count_factor)
      ),
      name = "Gene Count"
    ) +
    scale_fill_continuous_divergingx(
      palette = "RdBu",
      mid = 0.05,
      limits = c(0, 0.1),
      oob = scales::squish,
      name = expression(P[adj])
    ) +
    labs(
      x = "Signal",
      y = "GO Term",
      title = title
    ) +
    theme_bw() +
    theme(
      text = element_text(size = 20),
      plot.title = element_text(size = 25, hjust = 0.5, face = "bold", margin = margin(b = 15)),
      axis.text.y = ggtext::element_markdown(size = 28, lineheight = 0.5, margin = margin(l = 10)),
      axis.text.x = element_text(size = 28),
      axis.title.x = element_text(size = 28),
      axis.title.y = element_text(size = 22),
      legend.title = element_text(size = 24),
      legend.text = element_text(size = 20),
      legend.key = element_blank()
    ) 
}

clust1_GO <- plot_go_enrichment(cluster1_enrichment, "Enriched GO Terms in Cluster 1", upper_lim = 3)
clust2_GO <- plot_go_enrichment(cluster2_enrichment, "Enriched GO Terms in Cluster 2", upper_lim = 5)

ggsave(
  filename = (paste0(out_dir,"clust1_GO.pdf")),
  plot = clust1_GO,
  device = cairo_pdf,
  width = 14,
  height = 10,
  units = "in"
)

ggsave(
  filename = (paste0(out_dir,"clust2_GO.pdf")),
  plot = clust2_GO,
  device = cairo_pdf,
  width = 14,
  height = 10,
  units = "in"
)

week1_down <- read_tsv("data/Week1_down_enrichment.Process.tsv")
week1_up <- read_tsv("data/Week1_up_enrichment.Process.tsv")

week1_down_GO <- plot_go_enrichment(week1_down, "Enriched GO Terms in Cytokines with\nlog2FC < -1 at 1 Week Post-Exposure", upper_lim = 5)
week1_up_GO <- plot_go_enrichment(week1_up, "Enriched GO Terms in Cytokines with\nlog2FC > 1 at 1 Week Post-Exposure", upper_lim = 2)

ggsave(
  filename = (paste0(out_dir,"week1_down_GO.pdf")),
  plot = week1_down_GO,
  device = cairo_pdf,
  width = 14,
  height = 10,
  units = "in"
)

ggsave(
  filename = (paste0(out_dir,"week1_up_GO.pdf")),
  plot = week1_up_GO,
  device = cairo_pdf,
  width = 14,
  height = 10,
  units = "in"
)

#-------------------------------------------------------------------------------
# Effect of time on polyI:C response - 1 day vs 1 week Barplot
#-------------------------------------------------------------------------------
#summarize as a patient average:

# Filter and compute DE − FA per patient (logFC values)
DE_1day_1week_df_summary <- samples_long_DE_1day_1week %>%
  filter(Treatment %in% c("1day", "1week")) %>%
  group_by(Subject, Cytokine) %>%
  filter(n() == 2) %>%  # ensure patient-paired
  reframe(
    DE_minus_FA = Value[Treatment == "1week"] - Value[Treatment == "1day"])

# Summarize across patients: mean difference per cytokine
DE_1day_1week_df_summary_avg <- DE_1day_1week_df_summary %>%
  group_by(Cytokine) %>%
  summarise(
    mean_DE_minus_FA = mean(DE_minus_FA, na.rm = TRUE),
    .groups = "drop"
  )

DE_1day_1week_df_summary_avg_renamed <- DE_1day_1week_df_summary_avg %>%
  left_join(cytokine_96_panel_proteins %>% dplyr::distinct(Cytokine, Protein), by = "Cytokine") %>%
  mutate(Cytokine = Protein) # Protein is the gene name 

DE_1day_1week_df_summary_avg_renamed_abs1 <- DE_1day_1week_df_summary_avg_renamed %>% filter(abs(mean_DE_minus_FA) > 1)
cluster1 <- c("IL15", "TGFA", "PDGFB", "KITLG", "IFNB1", "TUBGCP2", "IL31", "CXCL12", "IL6", "IL12A/EBI3")
cluster2 <- c("VEGFA", "CSF3", "IL27", "LTA", "CCL4", "IL25", "IL13", "CCL8", "IL17F", "TSLP", "IL4", "CCL24", "FASLG", "IL1B", "TNFSF10", "CCL3", "FLT3LG", "CXCL8", "CXCL2", "CCL2", "IFNL2", "IL22", "TNF", "CCL19")

DE_1day_1week_df_summary_avg_renamed_abs1 <- DE_1day_1week_df_summary_avg_renamed_abs1 %>%
  mutate(
    `Cluster` = case_when(
      Protein %in% cluster1 ~ "Cluster 1",
      Protein %in% cluster2 ~ "Cluster 2",
      TRUE ~ "Cluster 3-5"
    ))

DE_1day_1week_waterfall <- ggplot(DE_1day_1week_df_summary_avg_renamed_abs1, aes(x = reorder(Cytokine, mean_DE_minus_FA), y = mean_DE_minus_FA,
                                                      fill = `Cluster`)) +
  geom_col() +
  scale_fill_manual(values = c("Cluster 1" = "#6F8596", "Cluster 2" = "#DE7576", "Cluster 3-5" = "#D2D2D2")) +
  labs(
    x = "Cytokine",
    y = "Average Patient-Paired log2FC\nin Cytokine Levels",
    title = "Effect of Time on Poly(I:C)-Induced Cytokine Responses\nin DE-Exposed Samples",
    subtitle = "Only cytokines with |log2FC| > 1 shown") +
  theme_bw() + 
  theme(
    plot.title = element_text(face = "bold", size = 25), 
    plot.subtitle = element_text(size = 20),
    axis.title = element_text(face = "bold", size = 20),
    axis.text.x = element_text(angle = -45, hjust = 0, size = 25),
    axis.text.y = element_text(size = 20),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 20)
  )

ggsave(
  filename = (paste0(out_dir,"DE_1day_1week_waterfall.pdf")),
  plot = DE_1day_1week_waterfall,
  device = cairo_pdf,
  width = 12,
  height = 8,
  units = "in"
)
