# DE_PolyIC (GSC-2684)

Manuscript analysis for the diesel exhaust exposure + poly(I:C) challenge study, integrating
RNA-seq, ATAC-seq, and cytokine (Eve Technologies 96-plex) data.

## Directory structure

```
data/      input data (sample sheets, counts, deidentification maps, cytokine panel)
scripts/   analysis scripts (R)
output/    generated results, organized by analysis (DEA/, GSEA/, ATAC/, PCA/, cytokine_panel/)
```

## Data

- `List_for_sequencing_Carlsten_GSC-2684_deidentified.csv` — sample
  metadata linking library IDs to subject/condition/timepoint.
- `salmon.merged.gene_counts_length_scaled.tsv` — RNA-seq gene-level counts (nf-core/rnaseq output).
- `ATAC-seq_*_frip_filtered_samplesheet_deidentified.csv` — ATAC-seq sample
  sheets (FRiP-filtered) for the DE+M vs FA+M and FA+P vs FA+M comparisons.
- `Aug-18-2025-Human-Cytokine-96-Plex-Discovery-Assay-Jha-67551_reformatted.xlsx`,
  `cytokine_panel_labels.tsv` — sample metadata linking cytokine panel IDs to subject/condition/timepoint.
- `CIBERSORTx_cell_proportions.csv` — estimated cell-type proportions used as
  covariates in some DESeq2 models.
- `cluster1/2_enrichment.Process.tsv`, `Week1_up/down_enrichment.Process.tsv` —
  precomputed GO enrichment tables (from StringDB) used by `cytokine_panel.R`.


## Scripts (run order)

1. **`DEA.R`** — deidentifies sample metadata, runs DESeq2 differential
   expression across conditions/timepoints, writes `output/DEA/*.tsv`.
2. **`PCA.R`** — PCA on (variance-stabilized) RNA-seq counts using
   deidentified sample IDs; writes `output/PCA/`.
3. **`GSEA.R`** — gene set enrichment (GO) on DESeq2 log2FC-ranked gene lists
   from `output/DEA/`; writes `output/GSEA/*.tsv`.
4. **`DEA_GSEA_plots.R`** — figures from the DEA/GSEA results.
5. **`DAR_csaw.R`** — ATAC-seq differential accessible region (DAR) calling
   with `csaw`: counts fragments in sliding windows, filters, tests with a
   paired (patient + condition) design, merges windows into regions (Simes'
   method, BH FDR < 0.05). Writes `output/ATAC/*DARs*`.
6. **`ATAC_plots.R`** — ATAC signal/enrichment and volcano plots from the
   `DAR_csaw.R` outputs; reads bigwigs from the nf-core/atacseq run.
7. **`DEA_GSEA_of_Patients_with_ATAC_samples.R`** — DEA/GSEA restricted to the
   subset of patients with matched ATAC-seq samples, for cross-omic
   comparison; writes `output/ATAC/DEA_GSEA_patients_with_ATAC_samples/`.
8. **`cytokine_panel.R`** — wrangles the 96-plex cytokine panel, heatmaps,
   waterfall plots, and GO enrichment on clustered cytokine groups; writes
   `output/cytokine_panel/`.
