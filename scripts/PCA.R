out_dir <- "output/PCA"
ntop    <- 500  

#-------------------------------------------------------------------------------
# Load libraries
#-------------------------------------------------------------------------------
library(DESeq2)
library(tidyverse)

#-------------------------------------------------------------------------------
# Read files
#-------------------------------------------------------------------------------
# De-identified IDs are already formatted as Participant.Exposure.Time
metadata_names <- read.csv("data/List_for_sequencing_Carlsten_GSC-2684_deidentified.csv", header = TRUE) %>%
  dplyr::select(External_ID, library_ID)

update_ids <- function(data, metadata){
  rename_map <- setNames(metadata$External_ID, metadata$library_ID)
  colnames(data) <- ifelse(colnames(data) %in% names(rename_map),
                           rename_map[colnames(data)],
                           colnames(data))
  return(data)
}

length_scaled_gene_counts <- read.delim("data/salmon.merged.gene_counts_length_scaled.tsv",
                                        check.names = FALSE) %>%
  column_to_rownames(var = "gene_id")

length_scaled_gene_counts <- update_ids(length_scaled_gene_counts, metadata_names)

#-------------------------------------------------------------------------------
# Select samples
#-------------------------------------------------------------------------------
# Keep only Participant.Exposure.Time IDs; drops the negative control, any
# un-renamed library IDs, and the failed patient (NA de-identifier)
all_samples  <- setdiff(colnames(length_scaled_gene_counts), "gene_name")
keep_samples <- all_samples[grepl("^[^.]+\\.[^.]+\\.[^.]+$", all_samples) &
                              !startsWith(all_samples, "NA.")]

message("Samples kept: ", length(keep_samples), " of ", length(all_samples))
message("Samples dropped: ", paste(setdiff(all_samples, keep_samples), collapse = ", "))

counts_pca <- length_scaled_gene_counts[, keep_samples]

#-------------------------------------------------------------------------------
# PCA on rlog-transformed length-scaled gene counts
#-------------------------------------------------------------------------------
# Approach follows nf-core/rnaseq:
# https://github.com/nf-core/rnaseq/blob/master/bin/deseq2_qc.r
# * Length-scaled gene counts are not integers, hence need to round
# * design = ~1 (intercept only) makes the transformation blind to sample groups
coldata <- data.frame(sample = keep_samples, row.names = keep_samples)

dds_pca <- DESeqDataSetFromMatrix(countData = round(counts_pca),
                                  colData   = coldata,
                                  design    = ~1)
dds_pca <- estimateSizeFactors(dds_pca)
rld     <- rlog(dds_pca, blind = TRUE)

run_pca <- function(transformed, ntop = 500) {
  mat        <- assay(transformed)
  rv         <- rowVars(mat)
  select     <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca        <- prcomp(t(mat[select, ]), center = TRUE, scale. = FALSE)
  percentVar <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  list(scores = as.data.frame(pca$x), percentVar = percentVar)
}

pca_res <- run_pca(rld, ntop = ntop)

#-------------------------------------------------------------------------------
# Write results
#-------------------------------------------------------------------------------
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# PC1 and PC2 scores per sample
pca_res$scores %>%
  dplyr::select(PC1, PC2) %>%
  rownames_to_column(var = "sample") %>%
  separate_wider_delim(sample, delim = ".",
                       names = c("Participant", "Exposure", "Time")) %>%
  write.table(file.path(out_dir, "GSC-2684_PCA.tsv"),
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

# Percent variance explained per PC
tibble(PC = paste0("PC", seq_along(pca_res$percentVar)),
       percent_variance = round(pca_res$percentVar, 2),
       ntop_genes = ntop) %>%
  write.table(file.path(out_dir, "GSC-2684_PCA_percent_variance.tsv"),
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

message("PC1: ", round(pca_res$percentVar[1]), "% | PC2: ",
        round(pca_res$percentVar[2]), "% variance (top ", ntop, " genes)")
#-------------------------------------------------------------------------------
# RNA PCA Plot
#-------------------------------------------------------------------------------
PCA <- read_tsv("output/PCA/GSC-2684_PCA.tsv")
PCA <- PCA %>% mutate(Participant = factor(Participant)) %>% mutate(Exposure = str_replace_all(Exposure, "_", "-"))

PVA_percent_variance <- read_tsv("output/PCA/GSC-2684_PCA_percent_variance.tsv")
PCA1 <- PVA_percent_variance %>% filter(PC == "PC1") %>% pull(percent_variance) 
PCA2 <- PVA_percent_variance %>% filter(PC == "PC2") %>% pull(percent_variance) 

shape_vals <- c(
  "IC-M" = 3,   
  "IC-P" = 4,
  "DE-M" = 24,
  "DE-P" = 23,
  "FA-M" = 25,
  "FA-P" = 22
)

patient_cols <- c(
  "1" = "#f2cc8f",
  "2" = "#81b29a",
  "3" = "#3d405b",
  "4" = "#219ebc",
  "5" = "#e07a5f"
)

PCA_RNA_plot <- ggplot(
  PCA,
  aes(x = PC1,
      y = PC2,
      colour = Participant,
      shape  = Exposure,
      fill   = Time)) +
  geom_point(size = 10, stroke = 3, alpha = 0.8) +
  scale_shape_manual(values = shape_vals, name = "Exposure") +
  scale_colour_manual(
    values = patient_cols,
    name   = "Participant"
  ) +
  scale_fill_manual(
    values = c("1day" = "white", "1week" = "black"),
    name   = "Time",
    labels = c("1day" = "1 Day", "1week" = "1 Week")
  ) +
  guides(
    fill  = guide_legend(override.aes = list(shape = 21, colour = "black", size = 12, stroke = 1)),
    shape = guide_legend(override.aes = list(fill = "white", size = 12, stroke = 1)),
    colour = guide_legend(override.aes = list(size = 12, stroke = 1))
  ) +
  labs(
    x = paste0("PC1: ", round(PCA1, 1), "% variance"),
    y = paste0("PC2: ", round(PCA2, 1), "% variance")
  ) +
  theme_bw() +
  theme(axis.title = element_text(size = 45),
        legend.text = element_text(size = 35),
        legend.title = element_text(size = 40),
        axis.text = element_text(size = 45),
        legend.key.size = unit(1.5, "cm")
  ) 

ggsave(
  filename = "output/PCA/PCA_RNA_plot.pdf",
  plot = PCA_RNA_plot,
  device = cairo_pdf, 
  width = 16,
  height = 12,
  dpi = 300
)
