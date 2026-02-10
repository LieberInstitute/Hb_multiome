#   Annotate DAR peaks with nearest genes and do GO by cell type on those genes.
#   Not technically TF-related

library(tidyverse)
library(TFBSTools)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # or species-appropriate
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg38)
library(clusterProfiler)
library(here)
library(sessioninfo)
library(GenomicRanges)
library(Signac)

peak_path = here(
    'processed-data', '06_peak_calling', '21_overlaping_FDRscores_TopHeatmap',
    'all_peaks_categorized.csv.gz'
)
seur_path = here(
    "processed-data", "06_peak_calling", "12_pseudobulk_MACS2",
    "Mid_pseudobulk.spearman.5e5_merged_peaks.rds"
)
promoter_window = 2000  # +/- around TSS

#   Load in DARs
peak_df = read_csv(peak_path, show_col_types = FALSE) |>
    filter(grepl('DAR', category)) |>
    distinct(peak_id, cell_type, .keep_all = TRUE)

peak_gr = peak_df |>
    pull(peak_id) |>
    unique() |>
    StringToGRanges(sep = c("-", "-"))

#   Add nearest gene + basic annotation
ann = annotatePeak(
        peak_gr, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
        tssRegion = c(-1 * promoter_window, promoter_window),
        annoDb = "org.Hs.eg.db"
    ) |>
    as.data.frame() |>
    as_tibble() |>
    dplyr::rename(
        nearest_gene_id = ENSEMBL, nearest_gene_name = SYMBOL
    ) |>
    mutate(peak_id = paste(seqnames, start, end, sep = "-")) |>
    select(peak_id, annotation, nearest_gene_id, nearest_gene_name)

peak_df = peak_df |>
    left_join(ann, by = "peak_id") |>
    #   For GO, an empirically linked gene is stronger evidence than using the
    #   nearest gene. Use whichever is available though
    mutate(gene_for_go = coalesce(link_gene_name, nearest_gene_name))

#   This exact set of genes was tested for linkage, and is also the set from
#   which nearest genes are drawn, so should be the appropriate universe for GO
seur = readRDS(seur_path)
back_universe = rownames(seur[['RNA']])

## Functional enrichment (GO/Pathways) for linked genes

for (this_cell_type in unique(peak_df$cell_type)) {
    gene_set = peak_df |>
        filter(cell_type == this_cell_type) |>
        pull(gene_for_go) |>
        unique()

    ego = enrichGO(
        gene = gene_set, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
        universe = back_universe, pAdjustMethod= "BH", pvalueCutoff = 1,
        qvalueCutoff = 0.05
    )

    dotplot(ego, showCategory = 5) + ggtitle(this_cell_type)
}

# Optional: rrvgo to reduce terms

##########  QA checks:
# Replicate structure (≥2 per group) for DAR calling
# GC/length matching for motif background
# Blacklist and low-mappability filters
# Sensitivity analysis on FDR/logFC thresholds
# Window size for LinkPeaks (±100–500 kb) and per-cluster stability.

## Desired Deliverables: 
# DARs_<celltype>_vs_rest.csv (full table)
# DARs_<celltype>_topN.bed (browser tracks)
# motif_enrichment_<celltype>.csv
# GO_<celltype>_linked_genes.csv
# QC plots: MA/volcano, region annotations pie/bar, distance-to-TSS, motif volcano, GO dotplot, link distance distribution.
