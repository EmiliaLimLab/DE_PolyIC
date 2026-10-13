#-------------------------------------------------------------------------------
# library
#-------------------------------------------------------------------------------
library(GenomicRanges)
library(IRanges)
library(rtracklayer)
library(EnrichedHeatmap)
library(dplyr)
library(tidyr)
library(ggplot2)
library(circlize)       
library(viridis)  
library(clusterProfiler)   
library(org.Hs.eg.db)     
library(readr)  
library(stringr)
library(ggrepel)
library(ggtext)  

#-------------------------------------------------------------------------------
# Average signal plot across all 3 conditions
#-------------------------------------------------------------------------------
bw_dir  <- "/projects/lim_lab_scratch/tlo/projects/IPPL/nfcore-atacseq/GSC-2684/bwa/merged_library/bigwig"
rds_dir <- "output/ATAC"
out_dir <- "output/ATAC"

# Build a GRanges of all DARs from an annotated DAR table
make_dar_gr <- function(dar_df) {
  gr <- GRanges(
    seqnames = dar_df$seqnames,
    ranges   = IRanges(start = dar_df$start, end = dar_df$end),
    strand   = dar_df$strand
  )
  mcols(gr) <- dar_df %>%
    dplyr::select(logFC, FDR, SYMBOL, annotation, distanceToTSS)
  gr
}

# Return 1 bp centers only; normalizeToMatrix's `extend` adds the +/- 2500 bp
center_extend <- function(gr) {
  resize(gr, width = 1, fix = "center")
}

create_profile_matrix <- function(bigwig_path, regions, upstream = 2500, downstream = 2500, bin_size = 50) {
  
  bw <- import(bigwig_path, format = "BigWig")
  
  normalizeToMatrix(
    bw,
    regions,
    extend = c(upstream, downstream),
    w = bin_size,
    mean_mode = "w0",
    value_column = "score",
    background = 0,
    smooth = TRUE,
    keep = c(0, 0.99),
    limit = NA,
    verbose = FALSE
  )
}

# Average profile per sample -> long format -> average per condition
make_condition_avg <- function(mat_list, sample_metadata) {
  avg_mat <- sapply(mat_list, colMeans, na.rm = TRUE)
  
  df <- as.data.frame(avg_mat)
  df$Bin <- seq_len(nrow(df))
  df <- pivot_longer(df, -Bin, names_to = "Sample", values_to = "Signal")
  df$Condition <- sample_metadata$Condition[match(df$Sample, sample_metadata$Sample)]
  df$Distance  <- seq(-2475, 2475, by = 50)[df$Bin]   # midpoint of each 50 bp bin
  
  df %>%
    group_by(Condition, Bin, Distance) %>%
    summarise(Signal = mean(Signal, na.rm = TRUE), .groups = "drop")
}

# Profile all 3 conditions over a shared region set
dar_PvsM   <- readRDS(file.path(rds_dir, "FA_P_FA_M_1day_DARs_annotated.rds"))
dar_DEvsFA <- readRDS(file.path(rds_dir, "DE_FA_M_1day_DARs_annotated.rds"))

common_regions <- center_extend(unique(c(granges(make_dar_gr(dar_PvsM)),
                                         granges(make_dar_gr(dar_DEvsFA)))))

bw_files_all <- c(
  Patient_4_DE_M = file.path(bw_dir, "F147887_REP1_T1_REP1.mLb.clN.bigWig"),
  Patient_2_DE_M = file.path(bw_dir, "F147886_REP1_T1_REP1.mLb.clN.bigWig"),
  Patient_4_FA_M = file.path(bw_dir, "F147881_REP1_T1_REP1.mLb.clN.bigWig"),
  Patient_2_FA_M = file.path(bw_dir, "F147880_REP1_T1_REP1.mLb.clN.bigWig"),
  Patient_4_FA_P = file.path(bw_dir, "F147884_REP1_T1_REP1.mLb.clN.bigWig"),
  Patient_2_FA_P = file.path(bw_dir, "F147883_REP1_T1_REP1.mLb.clN.bigWig")
)
mat_all <- lapply(bw_files_all, create_profile_matrix, regions = common_regions)
saveRDS(mat_all, file.path(out_dir, "DARs_3cond_mat_list_common.rds"))

sample_metadata <- data.frame(Sample = names(bw_files_all),
                              Condition = rep(c("DE_M", "FA_M", "FA_P"), each = 2))
plot_avg_3cond <- make_condition_avg(mat_all, sample_metadata)

condition_colors <- c("DE_M" = "darkred", "FA_M" = "darkblue", "FA_P" = "darkgreen")

signal_plot_3cond <- ggplot(plot_avg_3cond, aes(x = Distance, y = Signal, color = Condition)) +
  geom_line(linewidth = 1.2) +
  guides(color = guide_legend(override.aes = list(linewidth = 3))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  ylim(0.01, 0.28) +   
  scale_x_continuous(
    breaks = c(-2500, -1250, 0, 1250, 2500),
    labels = c("-2.5kb", "-1.25kb", "Center", "+1.25kb", "+2.5kb")
  ) +
  scale_color_manual(
    values = condition_colors,
    labels = c("DE_M" = "DE + Media Control", "FA_M" = "FA + Media Control", "FA_P" = "FA + Poly(I:C)")
  ) +
  labs(
    x = "Distance from peak center (bp)",
    y = "Mean ATAC-seq signal (RPM)",
    title = "ATAC-seq Profile Across DARs in 3 Conditions",
    color = "Condition"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 28),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15)
  )

signal_plot_3cond <- signal_plot_3cond +
  theme(plot.margin = margin(t = 5.5, r = 30, b = 5.5, l = 5.5, unit = "pt"))

ggsave(
  filename = (paste0(out_dir,"/signal_plot_3cond.pdf")),
  plot = signal_plot_3cond,
  device = cairo_pdf,
  width = 9,
  height = 10,
  units = "in"
)


#-------------------------------------------------------------------------------
# partitioned heatmap setup and reading input files
#-------------------------------------------------------------------------------
bw_dir   <- "/projects/lim_lab_scratch/tlo/projects/IPPL/nfcore-atacseq/GSC-2684/bwa/merged_library/bigwig"
data_dir <- "output/ATAC"

bw_files <- c(
  Patient_2_DE_M = file.path(bw_dir, "F147886_REP1_T1_REP1.mLb.clN.bigWig"),
  Patient_4_DE_M = file.path(bw_dir, "F147887_REP1_T1_REP1.mLb.clN.bigWig"),
  Patient_2_FA_M = file.path(bw_dir, "F147880_REP1_T1_REP1.mLb.clN.bigWig"),
  Patient_4_FA_M = file.path(bw_dir, "F147881_REP1_T1_REP1.mLb.clN.bigWig"),
  Patient_2_FA_P = file.path(bw_dir, "F147883_REP1_T1_REP1.mLb.clN.bigWig"),
  Patient_4_FA_P = file.path(bw_dir, "F147884_REP1_T1_REP1.mLb.clN.bigWig")
)

DE_FA_M_1day_DARs_merged   <- readRDS(file.path(data_dir, "DE_FA_M_1day_DARs_merged.rds"))
DE_FA_M_1day_DARs_merged <- DE_FA_M_1day_DARs_merged$merged_all

FA_P_FA_M_1day_DARs_merged <- readRDS(file.path(data_dir, "FA_P_FA_M_1day_DARs_merged.rds"))
FA_P_FA_M_1day_DARs_merged <- FA_P_FA_M_1day_DARs_merged$merged_all

# Significant DARs, kept as GRanges 
sig_DARs    <- DE_FA_M_1day_DARs_merged[DE_FA_M_1day_DARs_merged$FDR < 0.05]
sig_DARs_PM <- FA_P_FA_M_1day_DARs_merged[FA_P_FA_M_1day_DARs_merged$FDR < 0.05]

# Shared function: signal matrix for one bigWig over a set of DARs
# DAR body is scaled; +/- 1500 bp flanks in 10 bp bins
make_mat <- function(bw_file, regions) {
  bw <- import(bw_file, format = "BigWig")
  bw <- subsetByOverlaps(bw, regions + 5000)
  normalizeToMatrix(bw, regions, value_column = "score",
                    extend = 1500, mean_mode = "w0", w = 10, keep = c(0, 0.99))
}

title_size      <- 25   # column titles ("PolyI:C (avg)", "Media (avg)")
row_title_size  <- 20   # row split titles ("More open in ...")
axis_text_size  <- 35  # axis labels on line plots and heatmap bottom
legend_title_sz <- 15
legend_label_sz <- 12
line_width      <- 2.5    # thickness of lines in the line plots

heatmap_legend_style <- list(
  title_gp     = gpar(fontsize = legend_title_sz, fontface = "bold"),
  labels_gp    = gpar(fontsize = legend_label_sz),
  legend_height = unit(4, "cm"),
  grid_width   = unit(6, "mm")
)

#-------------------------------------------------------------------------------
# DE vs FA (Media) - partitioned heatmap
#-------------------------------------------------------------------------------
mat_DE_P2 <- make_mat(bw_files["Patient_2_DE_M"], sig_DARs)
mat_DE_P4 <- make_mat(bw_files["Patient_4_DE_M"], sig_DARs)
mat_FA_P2 <- make_mat(bw_files["Patient_2_FA_M"], sig_DARs)
mat_FA_P4 <- make_mat(bw_files["Patient_4_FA_M"], sig_DARs)

# Average patients within each condition
mat_DE_avg <- (mat_DE_P2 + mat_DE_P4) / 2
mat_FA_avg <- (mat_FA_P2 + mat_FA_P4) / 2

# # Shared color scale: 99th percentile across the individual matrices
# global_max_bw_DE_FA_M <- quantile(c(mat_DE_P2, mat_DE_P4, mat_FA_P2, mat_FA_P4), 0.99, na.rm = TRUE)
# col_DE_FA <- colorRamp2(
#   breaks = c(0, global_max_bw_DE_FA_M / 2, global_max_bw_DE_FA_M),
#   colors = viridis(3, option = "plasma")
# )

col_DE_FA <- colorRamp2(c(0, 0.4), c("white", "red"))

# Split rows by direction of change
group <- factor(ifelse(mcols(sig_DARs)$logFC > 0, "More open in DE", "More open in FA"),
                levels = c("More open in DE", "More open in FA"))
group_colors <- c("More open in DE" = "#76150D", "More open in FA" = "#0700B2")

lgd_lines <- Legend(
  at = levels(group),
  title = "Chromatin\nAccessibility",
  type = "lines",
  legend_gp = gpar(col = group_colors, lwd = 5),
  title_gp = gpar(fontsize = legend_title_sz, fontface = "bold"),
  labels_gp = gpar(fontsize = legend_label_sz),
  grid_height = unit(8, "mm"),
  grid_width = unit(12, "mm")
)

ht_DE_bw <- EnrichedHeatmap(mat_DE_avg,
                            name = "Signal",
                            col = col_DE_FA,
                            row_split = group,
                            row_title = NULL,
                            row_title_gp = gpar(fontsize = row_title_size),
                            column_title_gp = gpar(fontsize = title_size, fontface = "bold"),
                            axis_name_gp = gpar(fontsize = axis_text_size),
                            heatmap_legend_param = heatmap_legend_style,
                            top_annotation = HeatmapAnnotation(
                              lines = anno_enriched(
                                gp = gpar(col = group_colors, lwd = line_width),
                                ylim = c(0, 0.4),
                                height = unit(4, "cm"),
                                axis_param = list(side = "left", facing = "outside",
                                                  gp = gpar(fontsize = axis_text_size))
                              )
                            ),
                            column_title = "DE (avg)")

ht_FA_bw <- EnrichedHeatmap(mat_FA_avg,
                            name = "FA",
                            col = col_DE_FA,
                            row_split = group,
                            row_title = NULL,
                            row_title_gp = gpar(fontsize = row_title_size),
                            column_title_gp = gpar(fontsize = title_size, fontface = "bold"),
                            axis_name_gp = gpar(fontsize = axis_text_size),
                            show_heatmap_legend = FALSE,
                            top_annotation = HeatmapAnnotation(
                              lines = anno_enriched(
                                gp = gpar(col = group_colors, lwd = line_width),
                                ylim = c(0, 0.4),
                                height = unit(4, "cm"),
                                axis = FALSE
                              )
                            ),
                            column_title = "FA (avg)")

row_anno_DE_FA <- rowAnnotation(
  foo = anno_simple(as.character(group), col = group_colors, height = unit(2, "mm")),
  show_annotation_name = FALSE,
  width = unit(6, "mm")
)

pdf(file.path(out_dir, "DARs_signal_enrichment_DEvsFA.pdf"), width = 16, height = 14)
draw(
  ht_DE_bw + ht_FA_bw + row_anno_DE_FA,
  heatmap_legend_side = "right",
  annotation_legend_list = list(lgd_lines),
  merge_legend = TRUE,
  legend_gap = unit(6, "mm"),
  column_title = paste("DARs Signal Enrichment (n =", nrow(mat_DE_avg), "regions)"),
  column_title_gp = gpar(fontsize = 20, fontface = "bold"),
  padding = unit(c(10, 2, 2, 2), "mm")
)
dev.off()
#-------------------------------------------------------------------------------
# Poly(I:C) vs Media (FA) - partitioned heatmap
#-------------------------------------------------------------------------------
mat_FA_M_P2 <- make_mat(bw_files["Patient_2_FA_M"], sig_DARs_PM)
mat_FA_M_P4 <- make_mat(bw_files["Patient_4_FA_M"], sig_DARs_PM)
mat_FA_P_P2 <- make_mat(bw_files["Patient_2_FA_P"], sig_DARs_PM)
mat_FA_P_P4 <- make_mat(bw_files["Patient_4_FA_P"], sig_DARs_PM)

# Average patients within each condition
mat_FA_M_avg <- (mat_FA_M_P2 + mat_FA_M_P4) / 2
mat_FA_P_avg <- (mat_FA_P_P2 + mat_FA_P_P4) / 2

# # Shared color scale: 99th percentile across the individual matrices
# global_max_bw_FA_M_P <- quantile(c(mat_FA_M_P2, mat_FA_M_P4, mat_FA_P_P2, mat_FA_P_P4), 0.99, na.rm = TRUE)
# col_PvsM <- colorRamp2(
#   breaks = c(0, global_max_bw_FA_M_P / 2, global_max_bw_FA_M_P),
#   colors = viridis(3, option = "plasma")
# )

col_PvsM <- colorRamp2(c(0, 0.4), c("white", "red"))

# Split rows by direction of change
group_PM <- factor(ifelse(mcols(sig_DARs_PM)$logFC > 0, "More open in PolyI:C", "More open in media control"),
                   levels = c("More open in PolyI:C", "More open in media control"))
group_colors_PM <- c("More open in PolyI:C" = "#2A6218", "More open in media control" = "#0700B2")

lgd_lines_PM <- Legend(
  at = levels(group_PM),
  title = "Chromatin\nAccessibility",
  type = "lines",
  legend_gp = gpar(col = group_colors_PM, lwd = 5),
  title_gp = gpar(fontsize = legend_title_sz, fontface = "bold"),
  labels_gp = gpar(fontsize = legend_label_sz),
  grid_height = unit(8, "mm"),
  grid_width = unit(12, "mm")
)

ht_P_bw <- EnrichedHeatmap(mat_FA_P_avg,
                           name = "Signal",
                           col = col_PvsM,
                           row_split = group_PM,
                           row_title = NULL,
                           row_title_gp = gpar(fontsize = row_title_size),
                           column_title_gp = gpar(fontsize = title_size, fontface = "bold"),
                           axis_name_gp = gpar(fontsize = axis_text_size),
                           heatmap_legend_param = heatmap_legend_style,
                           top_annotation = HeatmapAnnotation(
                             lines = anno_enriched(
                               gp = gpar(col = group_colors_PM, lwd = line_width),
                               ylim = c(0, 0.4),
                               height = unit(4, "cm"),
                               axis_param = list(side = "left", facing = "outside",
                                                 gp = gpar(fontsize = axis_text_size))
                             )
                           ),
                           column_title = "PolyI:C (avg)")

ht_M_bw <- EnrichedHeatmap(mat_FA_M_avg,
                           name = "Media",
                           col = col_PvsM,
                           row_split = group_PM,
                           row_title = NULL,
                           row_title_gp = gpar(fontsize = row_title_size),
                           column_title_gp = gpar(fontsize = title_size, fontface = "bold"),
                           axis_name_gp = gpar(fontsize = axis_text_size),
                           show_heatmap_legend = FALSE,
                           top_annotation = HeatmapAnnotation(
                             lines = anno_enriched(
                               gp = gpar(col = group_colors_PM, lwd = line_width),
                               ylim = c(0, 0.4),
                               height = unit(4, "cm"),
                               axis = FALSE
                             )
                           ),
                           column_title = "Media (avg)")
row_anno_PM <- rowAnnotation(
  foo = anno_simple(as.character(group_PM), col = group_colors_PM, height = unit(2, "mm")),
  show_annotation_name = FALSE,
  width = unit(6, "mm")
)

pdf(file.path(out_dir, "DARs_signal_enrichment_PvsM.pdf"), width = 16, height = 14)
draw(
  ht_P_bw + ht_M_bw + row_anno_PM,
  heatmap_legend_side = "right",
  annotation_legend_list = list(lgd_lines_PM),
  merge_legend = TRUE,
  legend_gap = unit(6, "mm"),
  column_title = paste("DARs Signal Enrichment (n =", nrow(mat_FA_P_avg), "regions)"),
  column_title_gp = gpar(fontsize = 20, fontface = "bold"),
  padding = unit(c(10, 2, 2, 2), "mm")
)
dev.off()


#-------------------------------------------------------------------------------
# volcano plot setup and reading input files
#-------------------------------------------------------------------------------
DE_FA_M_1day_all_regions_annotated   <- readRDS("output/ATAC/DE_FA_M_1day_all_regions_annotated.rds")
FA_P_FA_M_1day_all_regions_annotated <- readRDS("output/ATAC/FA_P_FA_M_1day_all_regions_annotated.rds")

DE_MvsFA_M_patient_2_and_4_gseGO <- read.delim("output/ATAC/DEA_GSEA_patients_with_ATAC_samples/DE_MvsFA_M_patient_2_and_4_gseGO.tsv")
FA_PvsFA_M_patient_2_and_4_gseGO <- read.delim("output/ATAC/DEA_GSEA_patients_with_ATAC_samples/FA_PvsFA_M_patient_2_and_4_gseGO.tsv")

deseq_DE_MvsFA_M_patient_2_and_4_lfcSE_cutoff <- read_tsv("output/ATAC/DEA_GSEA_patients_with_ATAC_samples/deseq_DE_MvsFA_M_patient_2_and_4_lfcSE_cutoff.tsv")
deseq_FA_PvsFA_M_patient_2_and_4_lfcSE_cutoff <- read_tsv("output/ATAC/DEA_GSEA_patients_with_ATAC_samples/deseq_FA_PvsFA_M_patient_2_and_4_lfcSE_cutoff.tsv")

deseq_DE_MvsFA_M_patient_2_and_4_sig <- deseq_DE_MvsFA_M_patient_2_and_4_lfcSE_cutoff %>% filter(padj < 0.05)
deseq_FA_PvsFA_M_patient_2_and_4_sig <- deseq_FA_PvsFA_M_patient_2_and_4_lfcSE_cutoff %>% filter(padj < 0.05)

# Detecting Leading Edge Genes from Defense Response to Virus Pathway 
get_leading_edge <- function(gsea_df, go_id) {
  ids <- gsea_df %>%
    dplyr::filter(ID == go_id) %>%
    dplyr::pull(core_enrichment) %>%
    strsplit("/") %>%
    unlist()
  if (length(ids) == 0) {
    warning(go_id, " not found in GSEA results")
    return(character(0))
  }
  bitr(ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)$SYMBOL
}

DE_FA_1day_M_GO_viral_response_genes <- get_leading_edge(DE_MvsFA_M_patient_2_and_4_gseGO, "GO:0051607")
P_M_FA_GO_viral_response_genes       <- get_leading_edge(FA_PvsFA_M_patient_2_and_4_gseGO, "GO:0051607")

#-------------------------------------------------------------------------------
# Plotting volcano plot
#-------------------------------------------------------------------------------
# Rename ATAC columns to match the plotting function
DE_FA_M_1day_DARs_annotated_relabel <- DE_FA_M_1day_all_regions_annotated %>%
  dplyr::rename(padj = FDR, gene_name = SYMBOL, log2FoldChange = logFC)
FA_P_FA_M_1day_DARs_annotated_relabel <- FA_P_FA_M_1day_all_regions_annotated %>%
  dplyr::rename(padj = FDR, gene_name = SYMBOL, log2FoldChange = logFC)

# Viral leading-edge gene DARs highlighted; optional RNA-seq DEG labels in bold blue
volcano_plot_atac_viral <- function(volcano_data, plot_title, viral_genes,
                                    deg_genes = NULL, caption_text = NULL,
                                    genes_to_label_num = 30, ylim = 5,
                                    viral_color = "#a7c957") {
  
  volcano_data <- volcano_data %>%
    mutate(
      log10_padj = -log10(padj),
      viral = gene_name %in% viral_genes,
      plot_group = case_when(
        viral ~ "Viral response",
        padj < 0.05 & log2FoldChange > 0 ~ "⬆Upregulated",
        padj < 0.05 & log2FoldChange < 0 ~ "⬇Downregulated",
        TRUE ~ "Not significant"
      )
    ) %>%
    filter(is.finite(log2FoldChange), is.finite(log10_padj))
  
  # Top N most significant DARs get labels; RNA-seq DEGs are bold blue
  genes_to_label <- volcano_data %>%
    filter(padj < 0.05) %>%
    arrange(padj) %>%
    slice_head(n = genes_to_label_num) %>%
    mutate(
      is_deg_gene = if (is.null(deg_genes)) FALSE else gene_name %in% deg_genes,
      label_fontface = ifelse(is_deg_gene, "bold", "plain")
    )
  
  message("n_sig_total = ", sum(volcano_data$padj < 0.05, na.rm = TRUE))
  message("n_genes_to_label = ", nrow(genes_to_label))
  message("n_viral_points = ", sum(volcano_data$viral))
  
  ggplot(volcano_data, aes(x = log2FoldChange, y = log10_padj)) +
    
    # non-viral points (bottom layer)
    geom_point(
      data = subset(volcano_data, !viral),
      aes(fill = plot_group),
      shape = 21, size = 6, alpha = 0.8, color = "black", stroke = 0.5
    ) +
    
    # viral leading-edge points (top layer)
    geom_point(
      data = subset(volcano_data, viral),
      fill = viral_color,
      shape = 21, size = 6, alpha = 0.95, color = "black", stroke = 0.5
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    scale_fill_manual(
      values = c(
        "⬇Downregulated"  = "#3D5F75",
        "⬆Upregulated"    = "#C7828C",
        "Not significant" = "#A5A5A5"
      ),
      breaks = c("⬇Downregulated", "⬆Upregulated"),
      labels = c("Downregulated", "Upregulated")
    ) +
    geom_label_repel(
      data = genes_to_label,
      aes(x = log2FoldChange, y = log10_padj,
          label = gene_name, fontface = label_fontface),
      size = 10,
      box.padding = 0.25,
      point.padding = 0.2,
      max.overlaps = Inf,
      min.segment.length = 0,
      force = 85, 
      seed = 42,
      segment.alpha = 0.5,
      segment.size = 0.8,
      linewidth = 0,
      segment.color = "black",
      fill = alpha("white", 0.6),
      color = ifelse(genes_to_label$is_deg_gene, "#2980B9", "black"),
      label.r = 0,
      label.size = 0
    ) +
    labs(
      x = "log<sub>2</sub>FC",
      y = "-log<sub>10</sub>*P*<sub>adj</sub>",
      title = plot_title,
      caption = caption_text
    ) +
    xlim(-6, 6) +
    ylim(0, ylim) +
    theme_bw() +
    theme(
      text = element_text(size = 40),
      axis.text = element_text(size = 40),
      axis.title = element_markdown(size = 40),
      axis.title.x = element_markdown(),
      axis.title.y = element_markdown(),
      legend.title = element_blank(),
      legend.position = "top",
      legend.text = element_markdown(),
      legend.spacing.x = unit(0.7, "cm"),
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
      plot.caption = element_text(hjust = 0, size = 14)
    )
}

DAR_DE_FA_M_volcano <- volcano_plot_atac_viral(
  DE_FA_M_1day_DARs_annotated_relabel,
  "DARs in DE vs FA in Media Control Samples, 1 Day Post-Exposure",
  viral_genes = DE_FA_1day_M_GO_viral_response_genes,
  deg_genes = deseq_DE_MvsFA_M_patient_2_and_4_sig$gene_name, ylim = 4)

DAR_P_M_FA_volcano <- volcano_plot_atac_viral(
  FA_P_FA_M_1day_DARs_annotated_relabel,
  "DARs in PolyI:C vs Media in FA-Exposed Samples,\n1 Day Post-Exposure",
  viral_genes = P_M_FA_GO_viral_response_genes,
  deg_genes = deseq_FA_PvsFA_M_patient_2_and_4_sig$gene_name)

ggsave(
  filename = "output/ATAC/DAR_DE_FA_M_volcano.pdf",
  plot = DAR_DE_FA_M_volcano,
  device = cairo_pdf, 
  width = 10,
  height = 12,
  dpi = 300
)

ggsave(
  filename = "output/ATAC/DAR_P_M_FA_volcano.pdf",
  plot = DAR_P_M_FA_volcano,
  device = cairo_pdf, 
  width = 10,
  height = 12,
  dpi = 300
)


#-------------------------------------------------------------------------------
# Waterfall plot of DEGs Associated with significant DARs
#-------------------------------------------------------------------------------
DE_FA_M_1day_DARs_annotated_sig <- readRDS("output/ATAC/DE_FA_M_1day_DARs_annotated.rds") %>%
  filter(FDR < 0.05, !is.na(SYMBOL))

deseq_DE_MvsFA_M_patient_2_and_4_lfcSE_cutoff <- read_tsv("output/ATAC/DEA_GSEA_patients_with_ATAC_samples/deseq_DE_MvsFA_M_patient_2_and_4_lfcSE_cutoff.tsv")

# RNA-seq results for genes with a significant DAR
deseq_DE_MvsFA_M_RNA <- deseq_DE_MvsFA_M_patient_2_and_4_lfcSE_cutoff %>%
  filter(gene_name %in% DE_FA_M_1day_DARs_annotated_sig$SYMBOL) %>%
  arrange(padj) %>%
  distinct(gene_name, .keep_all = TRUE) %>%   
  mutate(direction = ifelse(log2FoldChange > 0, "up", "down")) %>%
  dplyr::rename(rna_logFC = log2FoldChange)

# One representative DAR per gene (largest |logFC|)
DE_FA_M_1day_DARs_annotated_sig_RNA <- DE_FA_M_1day_DARs_annotated_sig %>%
  filter(SYMBOL %in% deseq_DE_MvsFA_M_RNA$gene_name) %>%
  group_by(SYMBOL) %>%
  slice_max(order_by = abs(logFC), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    SYMBOL,
    atac_logFC = logFC,
    atac_FDR   = FDR,
    direction  = ifelse(atac_logFC > 0, "up", "down")
  )

# Join RNA + ATAC, add significance stars
DE_MvsFA_M_DARs_RNA_logFC_plot_df <- deseq_DE_MvsFA_M_RNA %>%
  inner_join(
    DE_FA_M_1day_DARs_annotated_sig_RNA,
    by = c("gene_name" = "SYMBOL"),
    suffix = c("_rna", "_atac")
  ) %>%
  mutate(
    concordance = ifelse(direction_rna == direction_atac, "same", "opposite"),
    sig_star = ifelse(!is.na(padj) & padj < 0.05, "*", NA),
    star_y = case_when(
      sig_star == "*" & rna_logFC > 0 ~ rna_logFC + 0.1,
      sig_star == "*" & rna_logFC < 0 ~ rna_logFC - 0.2,
      TRUE ~ NA_real_
    )
  )

# Plot
DE_MvsFA_M_DARs_RNA_logFC_plot <- ggplot(
  DE_MvsFA_M_DARs_RNA_logFC_plot_df,
  aes(x = reorder(gene_name, rna_logFC), y = rna_logFC, fill = direction_atac)) +
  geom_col(width = 1) +
  geom_text(
    aes(y = star_y, label = sig_star),
    size = 20,
    na.rm = TRUE,
    alpha = 0.7
  ) +
  scale_fill_manual(
    values = c("up" = "#C7828C", "down" = "#3D5F75"),
    name = "DAR Chromatin Activity",
    labels = c("up" = "More accessible in DE", "down" = "Less accessible in DE")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 30),
    axis.text.x = element_text(angle = -45, hjust = 0, size = 25),
    axis.text = element_text(colour = "black", size = 30),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    plot.subtitle = element_text(size = 14)
  ) +
  labs(
    x = "Gene",
    y = "LogFC",
    title = "LogFC of Genes (RNA-Seq) Associated with Significant DARs (ATAC-Seq)",
    subtitle = "DE vs FA, 1 Day Post-Exposure with Media,\nResults based on paired analysis of n=2 patients"
  )

ggsave(
  filename = "output/ATAC/DE_MvsFA_M_DARs_RNA_logFC_waterfall_plot.pdf",
  plot = DE_MvsFA_M_DARs_RNA_logFC_plot,
  device = cairo_pdf, 
  width = 18,
  height = 8,
  dpi = 300
)
