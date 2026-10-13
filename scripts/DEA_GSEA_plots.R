#-------------------------------------------------------------------------------
# Load libraries
#-------------------------------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(readr)
library(stringr)
library(tidytext)
library(colorspace)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(ggtext)

#-------------------------------------------------------------------------------
# GSEA GO Plots 
#-------------------------------------------------------------------------------
DE_FA_M_1day_gseGO <- read.delim("output/GSEA/DE_FA_M_1day_gseGO.tsv")
DE_FA_P_1day_gseGO <- read.delim("output/GSEA/DE_FA_P_1day_gseGO.tsv")
DE_FA_P_1week_gseGO <- read.delim("output/GSEA/DE_FA_P_1week_gseGO.tsv")
M_P_FA_gseGO <- read.delim("output/GSEA/M_P_FA_gseGO.tsv")
M_P_DE_gseGO <- read.delim("output/GSEA/M_P_DE_gseGO.tsv")
M_P_IC_gseGO <- read.delim("output/GSEA/M_P_IC_gseGO.tsv")

DE_P_1week_1day_gseGO <- read.delim("output/GSEA/DE_P_1week_1day_gseGO.tsv")

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
    guides(
      fill = guide_colourbar(
        theme = theme(
          legend.key.height = unit(12, "lines"),
          legend.key.width  = unit(1.5, "lines")
        )
      )
    ) +
    scale_size_continuous(
      range = c(5, 18),
      breaks = c(50, 100, 150, 200, 250, 300)
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
      axis.text.y = ggtext::element_markdown(size = 35, lineheight = 0.5,
                                             margin = margin(l = 10)),
      axis.text.x = element_text(size = 30),
      axis.title.x = element_text(size = 35),
      axis.title.y = element_text(size = 35),
      legend.title = ggtext::element_markdown(size = 30),
      legend.text = element_text(size = 30),
      strip.text.x = element_text(face = "bold", size = 24),
      strip.text.y = element_text(face = "bold", size = 18, angle = 270),
      panel.spacing.x = unit(0.4, "cm"),
      panel.spacing.y = unit(0.2, "cm"),
      plot.subtitle = element_text(size = 20, hjust = 0.5,
                                   margin = margin(b = 5))
    )
}

go_enrichment_DEvsFA <- plot_go_enrichment_faceted_top(df_list_DEvsFA, titles_DEvsFA, 
                                                       title_text = "GO Enrichment Analysis:\nDE vs FA Samples Across Conditions",
                                                       caption_text = "Top activated and suppressed GO terms by adjusted p-value")

ggsave(
  filename = "output/GSEA/go_enrichment_DEvsFA_plot.pdf",
  plot = go_enrichment_DEvsFA,
  device = cairo_pdf, 
  width = 28,
  height = 15,
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
  filename = "output/GSEA/go_enrichment_PvsM_plot.pdf",
  plot = go_enrichment_PvsM,
  device = cairo_pdf, 
  width = 28,
  height = 15,
  dpi = 300
)
#-------------------------------------------------------------------------------
# Quadrant Plots (Leading Edge Genes from Viral Response GO Term GO:0051607) 
#-------------------------------------------------------------------------------
deseq_DE_FA_M_1day_lfcSE_cutoff <- read_tsv("output/DEA/deseq_DE_FA_M_1day_lfcSE_cutoff.tsv")
deseq_DE_FA_P_1day_lfcSE_cutoff <- read_tsv("output/DEA/deseq_DE_FA_P_1day_lfcSE_cutoff.tsv")
deseq_DE_FA_P_1week_lfcSE_cutoff <- read_tsv("output/DEA/deseq_DE_FA_P_1week_lfcSE_cutoff.tsv")
deseq_M_P_IC_lfcSE_cutoff <- read_tsv("output/DEA/deseq_M_P_IC_lfcSE_cutoff.tsv")
deseq_M_P_FA_lfcSE_cutoff <- read_tsv("output/DEA/deseq_M_P_FA_lfcSE_cutoff.tsv")
deseq_M_P_DE_lfcSE_cutoff <- read_tsv("output/DEA/deseq_M_P_DE_lfcSE_cutoff.tsv")
deseq_DE_P_1week_1day_lfcSE_cutoff <- read_tsv("output/DEA/deseq_DE_P_1week_1day_lfcSE_cutoff.tsv")

get_core_symbols <- function(go_res, go_id = "GO:0051607") {  # defense response to virus 
  go_res %>%
    as.data.frame() %>%  
    filter(ID == go_id) %>%
    pull(core_enrichment) %>%
    strsplit("/") %>%
    unlist() %>%
    bitr(fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db) %>%
    pull(SYMBOL)
}

M_P_FA_gseGO_viral_response_genes <- get_core_symbols(M_P_FA_gseGO)
M_P_DE_GO_viral_response_genes <- get_core_symbols(M_P_DE_gseGO)
M_P_IC_GO_viral_response_genes <- get_core_symbols(M_P_IC_gseGO)
DE_FA_M_1day_gseGO_viral_response_genes <- get_core_symbols(DE_FA_M_1day_gseGO)
DE_FA_P_1day_gseGO_viral_response_genes <- get_core_symbols(DE_FA_P_1day_gseGO)
DE_FA_P_1week_gseGO_viral_response_genes <- get_core_symbols(DE_FA_P_1week_gseGO)
DE_P_1week_1day_gseGO_viral_response_genes <- get_core_symbols(DE_P_1week_1day_gseGO) 

PvsM_inFA_vs_DEvsFA_1day_P <- union(M_P_FA_gseGO_viral_response_genes, DE_FA_P_1day_gseGO_viral_response_genes)
PvsM_inFA_vs_DEvsFA_1week_P <- union(M_P_FA_gseGO_viral_response_genes, DE_FA_P_1week_gseGO_viral_response_genes)
PvsM_inFA_vs_DE_P_1week_1day <- union(M_P_FA_gseGO_viral_response_genes, DE_P_1week_1day_gseGO_viral_response_genes)

viral_label <- "Leading Edge Gene from Defense Response to Virus (GO:0051607) GO Term"
x_lab_PvsM  <- expression(log[2] * "FC in PolyI:C vs Media in FA-Exposed Samples")

quadrant_theme <- theme_bw() +
  theme(
    text = element_text(size = 20),
    legend.position = "top",
    legend.justification = "center",
    legend.title = element_blank(),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 35),
    panel.spacing.x = unit(0.4, "cm"),
    panel.spacing.y = unit(0.2, "cm"),
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

# function builds the plot, linear model and correlation for one comparison
plot_fc_quadrant <- function(y_res, genes, y_lab,
                             x_res = deseq_M_P_FA_lfcSE_cutoff,
                             x_lab = x_lab_PvsM,
                             extra_labels = NULL,
                             n_labels = 30) {
  
  # Foreground: leading edge genes
  fc_combined <- x_res %>%
    dplyr::select(gene_name, x = log2FoldChange) %>%
    filter(gene_name %in% genes) %>%
    left_join(dplyr::select(y_res, gene_name, y = log2FoldChange), by = "gene_name") %>%
    mutate(viral = viral_label)
  
  # Background: all genes present in both comparisons
  background <- distinct(x_res, gene_name, .keep_all = TRUE) %>%
    inner_join(distinct(y_res, gene_name, .keep_all = TRUE), by = "gene_name") %>%
    transmute(x = log2FoldChange.x, y = log2FoldChange.y)
  
  # Labels: top n by combined |FC|, plus any extra genes requested
  label_genes <- fc_combined %>%
    arrange(desc(abs(x) + abs(y))) %>%
    slice_head(n = n_labels) %>%
    bind_rows(filter(fc_combined, gene_name %in% extra_labels)) %>%
    distinct(gene_name, .keep_all = TRUE)
  
  p <- ggplot(mapping = aes(x = x, y = y)) +
    geom_point(data = background, shape = 21, size = 5, alpha = 0.7,
               color = "black", stroke = 0.5, fill = "#A5A5A5") +
    geom_point(data = fc_combined, aes(fill = viral), shape = 21, size = 5,
               alpha = 0.9, color = "black", stroke = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_label_repel(
      data = label_genes, aes(label = gene_name),
      box.padding = 0.5, point.padding = 0.5, max.overlaps = Inf,
      min.segment.length = 0.3, force = 45, size = 9, seed = 42,
      segment.alpha = 0.5, segment.size = 0.5, segment.color = "black",
      fill = alpha("white", 0.6), label.r = 0, label.size = 0, linewidth = 0
    ) +
    coord_cartesian(xlim = c(-3, 4), ylim = c(-4, 4)) +
    scale_fill_manual(values = setNames("#a7c957", viral_label)) +
    labs(x = x_lab, y = y_lab) +
    quadrant_theme
  
  list(
    plot  = p,
    data  = fc_combined,
    model = lm(y ~ x, data = fc_combined),
    cor   = cor.test(fc_combined$x, fc_combined$y)
  )
}

res_1day <- plot_fc_quadrant(
  y_res = deseq_DE_FA_P_1day_lfcSE_cutoff,
  genes = PvsM_inFA_vs_DEvsFA_1day_P,
  y_lab = expression(log[2] * "FC in DE vs FA in PolyI:C-Treated Samples, 1 Day Post-Exposure")
)
res_1day$cor

res_1week <- plot_fc_quadrant(
  y_res = deseq_DE_FA_P_1week_lfcSE_cutoff,
  genes = PvsM_inFA_vs_DEvsFA_1week_P, 
  y_lab = expression(log[2] * "FC in DE vs FA in PolyI:C-Treated Samples, 1 Week Post-Exposure"),
  extra_labels = "NLRP3"
)
res_1week$cor

res_1week_1day <- plot_fc_quadrant(
  y_res = deseq_DE_P_1week_1day_lfcSE_cutoff,
  genes = PvsM_inFA_vs_DE_P_1week_1day,
  y_lab = expression(log[2] * "FC in DE Samples, 1 Week vs 1 Day Post-Exposure"),
  extra_labels = "NLRP3"
)
res_1week_1day$cor

# Quadrant Plots Together
DEA_antiviral_responses_quadrant_plot <- 
  res_1day$plot + 
  plot_spacer() +  # Adds empty space
  res_1week$plot +
  plot_spacer() +
  res_1week_1day$plot +
  
  plot_layout(
    ncol = 5, 
    guides = "collect",
    widths = c(5, 0.5, 5, 0.5, 5) 
  ) &
  plot_annotation(
    title = "Comparing Fold Changes of DEGs adjusted for basal, secretory, and multiciliated epithelial cell proportions",
    theme = theme(
      legend.position = "top",
      legend.justification = "center",
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
    )
  )

ggsave(
  filename = "output/DEA/DEA_antiviral_responses_quadrant_plot.pdf",
  plot = DEA_antiviral_responses_quadrant_plot,
  device = cairo_pdf, 
  width = 30,
  height = 10,
  dpi = 300
)


#-------------------------------------------------------------------------------
# DEA volcanos with Viral-Associated dots coloured
#-------------------------------------------------------------------------------
viral_LEG <- function(df, viral_LEGs){
  df <- df %>% mutate(
    viral = gene_name %in% viral_LEGs
  )
  return(df)
}

deseq_DE_FA_M_1day_viral_LEG <- viral_LEG(deseq_DE_FA_M_1day_lfcSE_cutoff, DE_FA_M_1day_gseGO_viral_response_genes)
deseq_DE_FA_P_1day_viral_LEG <- viral_LEG(deseq_DE_FA_P_1day_lfcSE_cutoff, DE_FA_P_1day_gseGO_viral_response_genes)
deseq_DE_FA_P_1week_viral_LEG <- viral_LEG(deseq_DE_FA_P_1week_lfcSE_cutoff, DE_FA_P_1week_gseGO_viral_response_genes)
deseq_P_M_FA_viral_LEG <- viral_LEG(deseq_M_P_FA_lfcSE_cutoff, M_P_FA_gseGO_viral_response_genes)
deseq_P_M_DE_viral_LEG <- viral_LEG(deseq_M_P_DE_lfcSE_cutoff, M_P_DE_GO_viral_response_genes)
deseq_P_M_IC_viral_LEG <- viral_LEG(deseq_M_P_IC_lfcSE_cutoff, M_P_IC_GO_viral_response_genes)

volcano_plot_viral <- function(volcano_data, plot_title, caption_text = NULL, genes_to_label_num = 30) {
  
  volcano_data <- volcano_data %>%
    mutate(
      log10_padj = -log10(padj),
      plot_group = case_when(
        viral ~ "Viral response",
        padj < 0.05 & log2FoldChange > 0 ~ "⬆Upregulated",
        padj < 0.05 & log2FoldChange < 0 ~ "⬇Downregulated",
        TRUE ~ "Not significant"
      )
    ) %>%
    filter(
      is.finite(log2FoldChange),
      is.finite(log10_padj)
    )
  genes_to_label <- volcano_data %>%
    filter(padj < 0.05) %>%
    arrange(padj) %>%  
    slice_head(n = genes_to_label_num)
  
  
  message("n_sig_total = ", sum(volcano_data$padj < 0.05, na.rm = TRUE))
  message("n_genes_to_label = ", nrow(genes_to_label))
  
  # Filter out non-finite values 
  ggplot(volcano_data,
         aes(x = log2FoldChange, y = log10_padj)) +
    
    # non-viral points (bottom layer)
    geom_point(
      data = subset(volcano_data, !viral),
      aes(fill = plot_group),
      shape = 21,
      size = 4,
      alpha = 0.8,
      color = "black",
      stroke = 0.5
    ) +
    
    # viral points (top layer)
    geom_point(
      data = subset(volcano_data, viral),
      fill = "#a7c957",
      shape = 21,
      size = 4,
      alpha = 0.95,
      color = "black",
      stroke = 0.5
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    scale_fill_manual(
      values = c(
        "⬇Downregulated" = "#3D5F75",
        "⬆Upregulated"   = "#C7828C",
        "Not significant" = "#A5A5A5"
      ),
      breaks = c("⬇Downregulated", "⬆Upregulated"),
      labels = c("Downregulated", "Upregulated")
    )+
    geom_label_repel(
      data = genes_to_label,
      aes(
        x = log2FoldChange,
        y = log10_padj,
        label = gene_name
      ),
      size = 9, 
      box.padding = 0.25,
      point.padding = 0.2,
      max.overlaps = Inf,
      min.segment.length = 0,
      force = 50,
      seed = 42,
      segment.alpha = 0.5,
      segment.size = 0.4,
      linewidth = 0,
      segment.color = "black",
      fill = alpha("white", 0.6),
      label.r = 0,
      label.size = 0
    ) +
    labs(
      x = "log<sub>2</sub>FC",
      y = "-log<sub>10</sub>*P*<sub>adj</sub>",
      title = plot_title,
      caption = caption_text
    ) +
    coord_cartesian(xlim = c(-5, 5), ylim = c(0, 5)) +
    theme_bw() +
    theme(
      text = element_text(size = 27),
      axis.text = element_text(size = 30),
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

deseq_P_M_FA_viral_LEG_volcano <- volcano_plot_viral(deseq_P_M_FA_viral_LEG, "PolyI:C vs Media in FA-Exposed Samples", genes_to_label_num = 30)

ggsave(
  filename = "output/DEA/deseq_P_M_FA_viral_LEG_volcano.pdf",
  plot = deseq_P_M_FA_viral_LEG_volcano,
  device = cairo_pdf, 
  width = 10,
  height = 12,
  dpi = 300
)

deseq_DE_FA_M_1day_viral_LEG_volcano <- volcano_plot_viral(deseq_DE_FA_M_1day_viral_LEG, "DE vs FA in Media Control Samples 1 Day PE")
deseq_DE_FA_P_1day_viral_LEG_volcano <- volcano_plot_viral(deseq_DE_FA_P_1day_viral_LEG, "DE vs FA in PolyI:C Samples 1 Day PE")
deseq_DE_FA_P_1week_viral_LEG_volcano <- volcano_plot_viral(deseq_DE_FA_P_1week_viral_LEG, "DE vs FA in PolyI:C Samples 1 Week PE", genes_to_label_num = 15)

deseq_P_M_DE_viral_LEG_volcano <- volcano_plot_viral(deseq_P_M_DE_viral_LEG, "PolyI:C vs Media in DE-Exposed Samples", genes_to_label_num = 15)
deseq_P_M_IC_viral_LEG_volcano <- volcano_plot_viral(deseq_P_M_IC_viral_LEG, "PolyI:C vs Media in Incubator Control Samples", genes_to_label_num = 30)

suppl_volcanos <- deseq_DE_FA_M_1day_viral_LEG_volcano + deseq_DE_FA_P_1day_viral_LEG_volcano + deseq_DE_FA_P_1week_viral_LEG_volcano + 
  deseq_P_M_IC_viral_LEG_volcano + deseq_P_M_DE_viral_LEG_volcano + plot_layout(ncol = 3, nrow = 2)

ggsave(
  filename = "output/DEA/supplementary_volcanos.pdf",
  plot = suppl_volcanos,
  device = cairo_pdf,
  width = 18,
  height = 17,
  dpi = 300
)
