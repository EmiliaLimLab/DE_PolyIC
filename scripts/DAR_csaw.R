# Pipeline overview
# 
# Step 1  count_and_test_windows():
#   - Counts ATAC-seq fragments in 500 bp sliding windows (250 bp spacing)
#     (MAPQ >= 30, duplicates skipped, paired-end fragments)
#   - Filters low-abundance windows (average logCPM > -1, ~0.5 CPM)
#   - Library sizes = reads in the kept windows (keep.lib.sizes = FALSE)
#   - Paired design (~ patient + condition)
#   - Saves window-level results + counts
#
# Step 2  merge_csaw_results():
#   - Merges all filtered windows into regions (merge first, filter last)
#   - Combines window p-values per region (Simes' method), region-level BH FDR
#   - Returns all regions and significant DARs (FDR < 0.05)

out_dir <- "output/ATAC"
#-------------------------------------------------------------------------------
# library
#-------------------------------------------------------------------------------
library(csaw)
library(BiocParallel)
library(readr)
library(dplyr)
library(edgeR)
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
#-------------------------------------------------------------------------------
# Read files
#-------------------------------------------------------------------------------
DE_FA_M_1day_samplesheet_filtered <- read_csv("data/ATAC-seq_DE_FA-M_1day_frip_filtered_samplesheet_deidentified.csv")
FA_P_FA_M_1day_samplesheet_filtered <- read_csv("data/ATAC-seq_FA-P_FA-M_1day_frip_filtered_samplesheet_deidentified.csv")

#-------------------------------------------------------------------------------
# Step 1: count fragments in windows and test each window
#-------------------------------------------------------------------------------
count_and_test_windows <- function(samplesheet, comparison,
                                   condition_col = "condition",
                                   patient_col = "patient",
                                   bam_col = "bam",
                                   colname_coef,
                                   width = 500, spacing = 250,
                                   min_logCPM = -1,
                                   baseline,
                                   exposure,
                                   out_dir,
                                   workers = 8) {
  
  # Prepare metadata
  samples <- samplesheet %>%
    transmute(
      sample_id = SampleID,
      bam       = !!sym(bam_col),
      condition = !!sym(condition_col),
      patient   = factor(!!sym(patient_col))
    )
  samples$condition <- factor(samples$condition, levels = c(baseline, exposure))
  
  # Read filtering: MAPQ >= 30, skip duplicate-flagged reads, paired-end fragments
  param <- readParam(minq = 30, dedup = TRUE, pe = "both")
  
  # Count fragments in sliding windows
  win_counts <- windowCounts(
    bam.files = samples$bam,
    width     = width,
    spacing   = spacing,
    param     = param,
    BPPARAM   = MulticoreParam(workers = workers)
  )
  
  # Convert to DGEList
  y <- asDGEList(win_counts)
  
  # Pre-filter low-abundance windows; library sizes recomputed from kept windows
  keep <- aveLogCPM(y) > min_logCPM
  y <- y[keep, , keep.lib.sizes = FALSE]
  message(comparison, ": kept ", sum(keep), " of ", length(keep), " windows")
  
  # Design matrix (paired by patient)
  design <- model.matrix(~ patient + condition, data = samples)
  message(comparison, ": residual df = ", nrow(design) - ncol(design))
  
  # Fit GLM and test
  y   <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  res <- glmQLFTest(fit, coef = colname_coef)
  
  # Window-level FDR
  res$table$FDR <- p.adjust(res$table$PValue, method = "BH")
  
  # Filtered window coordinates with stats
  win_ranges_filt <- rowRanges(win_counts)[keep]
  win_ranges_filt$logFC  <- res$table$logFC
  win_ranges_filt$PValue <- res$table$PValue
  win_ranges_filt$FDR    <- res$table$FDR
  
  # Raw count matrix for kept windows
  count_matrix <- assay(win_counts)[keep, ]
  colnames(count_matrix) <- samples$sample_id
  
  results <- list(
    ranges           = win_ranges_filt,
    counts           = count_matrix,
    samples          = samples,
    design           = design,
    filtered_dgelist = y,
    test_results     = res$table
  )
  
  saveRDS(results, file = file.path(out_dir, paste0(comparison, ".rds")))
  
  return(results)
}

#-------------------------------------------------------------------------------
# Step 2: merge windows into regions and call DARs
#-------------------------------------------------------------------------------
merge_csaw_results <- function(csaw_list, tol = 500, max_width = 5000, fdr_cutoff = 0.05) {
  
  win_ranges_filt <- csaw_list$ranges
  test_results    <- csaw_list$test_results
  
  # Approximate lfcSE from the F statistic
  test_results$lfcSE <- abs(test_results$logFC) / sqrt(pmax(test_results$F, 0))
  test_results$lfcSE[!is.finite(test_results$lfcSE)] <- NA
  
  # Merge ALL filtered windows (merge first, filter last)
  merged <- mergeWindows(win_ranges_filt, tol = tol, max.width = max_width)
  
  # Combine window stats into region stats (Simes' method; BH applied)
  combined <- combineTests(ids = merged$ids, tab = test_results,
                           pval.col = "PValue", fc.col = "logFC")
  
  # Region coordinates with combined stats
  dar_gr <- merged$region
  mcols(dar_gr) <- combined
  # Columns: num.tests, num.up.logFC, num.down.logFC, PValue, FDR,
  #          direction, rep.test, rep.logFC
  
  # Stats from the representative (best) window of each region
  dar_gr$best_lfcSE <- test_results$lfcSE[combined$rep.test]
  dar_gr$logFC      <- combined$rep.logFC
  
  # Significant DARs
  dar_gr_sig <- dar_gr[dar_gr$FDR < fdr_cutoff]
  
  csaw_list$merged_all   <- dar_gr
  csaw_list$merged_sig   <- dar_gr_sig
  csaw_list$test_results <- test_results
  
  return(csaw_list)
}

#-------------------------------------------------------------------------------
# Run: P vs M (FA, 1 day)
#-------------------------------------------------------------------------------
FA_P_FA_M_1day_DAR_csaw_counts <- count_and_test_windows(
  samplesheet   = FA_P_FA_M_1day_samplesheet_filtered,
  comparison    = "FA_P_FA_M_1day_DAR_csaw_counts",
  condition_col = "Treatment",
  patient_col   = "Factor",
  bam_col       = "bamReads",
  baseline      = "M",
  exposure      = "P",
  colname_coef  = "conditionP",
  out_dir       = out_dir,
  workers       = 8
)

FA_P_FA_M_1day_DARs_merged <- merge_csaw_results(FA_P_FA_M_1day_DAR_csaw_counts)
saveRDS(FA_P_FA_M_1day_DARs_merged, file.path(out_dir, "FA_P_FA_M_1day_DARs_merged.rds"))

#-------------------------------------------------------------------------------
# Run: DE vs FA (M, 1 day)
#-------------------------------------------------------------------------------
DE_FA_M_1day_DAR_csaw_counts <- count_and_test_windows(
  samplesheet   = DE_FA_M_1day_samplesheet_filtered,
  comparison    = "DE_FA_M_1day_DAR_csaw_counts",
  condition_col = "Condition",
  patient_col   = "Factor",
  bam_col       = "bamReads",
  baseline      = "FA",
  exposure      = "DE",
  colname_coef  = "conditionDE",
  out_dir       = out_dir,
  workers       = 8
)
DE_FA_M_1day_DARs_merged <- merge_csaw_results(DE_FA_M_1day_DAR_csaw_counts)
saveRDS(DE_FA_M_1day_DARs_merged, file.path(out_dir, "DE_FA_M_1day_DARs_merged.rds"))

#-------------------------------------------------------------------------------
# Annotation 
#-------------------------------------------------------------------------------
annotate_DARs <- function(dar_gr, output_file,
                          txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                          tss_region = c(-3000, 3000),
                          verbose = FALSE) {
  
  if (!inherits(dar_gr, "GRanges")) stop("dar_gr must be a GRanges object!")
  
  # Annotate regions: nearest gene, genomic feature (promoter, intron, etc.),
  # distance to TSS, gene symbols via org.Hs.eg.db
  peak_anno <- annotatePeak(
    dar_gr,
    TxDb      = txdb,
    tssRegion = tss_region,
    annoDb    = "org.Hs.eg.db",
    verbose   = verbose
  )

  # alongside the annotation columns
  anno_df <- as.data.frame(peak_anno)
  
  saveRDS(anno_df, file = output_file)
  write.csv(anno_df, sub("\\.rds$", ".csv", output_file), row.names = FALSE)
  
  message("Annotation complete! Saved to: ", output_file, " (+ .csv)")
  
  return(anno_df)
}

# Run for each comparison
comparisons <- c("DE_FA_M_1day", "FA_P_FA_M_1day")

DE_FA_M_1day_DARs_merged <- readRDS(file.path(out_dir, "DE_FA_M_1day_DARs_merged.rds"))
FA_P_FA_M_1day_DARs_merged <- readRDS(file.path(out_dir, "FA_P_FA_M_1day_DARs_merged.rds"))

for (cmp in comparisons) {
  
  res <- readRDS(file.path(out_dir, paste0(cmp, "_DARs_merged.rds")))
  
  # Significant DARs (region-level FDR < 0.05)
  dar_gr <- res$merged_sig
  if (length(dar_gr) == 0) {
    message(cmp, ": no DARs to annotate, skipping")
  } else {
    annotate_DARs(dar_gr, file.path(out_dir, paste0(cmp, "_DARs_annotated.rds")))
  }
  
  # All tested regions: background set for enrichment comparisons
  annotate_DARs(res$merged_all,
                file.path(out_dir, paste0(cmp, "_all_regions_annotated.rds")))
}
