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
library(rtracklayer)

peak_path = here(
    'processed-data', '06_peak_calling', '21_overlaping_FDRscores_TopHeatmap',
    'all_peaks_categorized.csv.gz'
)
gtf_path = '/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz'
promoter_window = 2000  # +/- around TSS

#   Load in DARs
peak_df = read_csv(peak_path, show_col_types = FALSE) |>
    dplyr::filter(grepl('DAR', category)) |>
    distinct(peak_id, cell_type, .keep_all = TRUE)

peak_gr = peak_df |>
    pull(peak_id) |>
    unique() |>
    StringToGRanges(sep = c("-", "-"))

gtf = import(gtf_path) |>
    as.data.frame() |>
    as_tibble() |>
    filter(type == "gene") |>
    select(gene_id, gene_name)

#   Add nearest gene + basic annotation
ann = annotatePeak(
        peak_gr, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
        tssRegion = c(-1 * promoter_window, promoter_window),
        annoDb = "org.Hs.eg.db"
    ) |>
    as.data.frame() |>
    as_tibble() |>
    dplyr::rename(nearest_gene_name = SYMBOL) |>
    mutate(
        peak_id = paste(seqnames, start, end, sep = "-"),
        nearest_gene_id = gtf$gene_id[match(nearest_gene_name, gtf$gene_name)]
    ) |>
    select(peak_id, annotation, nearest_gene_id, nearest_gene_name)

peak_df = peak_df |>
    left_join(ann, by = "peak_id") |>
    #   For GO, an empirically linked gene is stronger evidence than using the
    #   nearest gene. Use whichever is available though
    mutate(gene_for_go = coalesce(link_gene_id, nearest_gene_id)) |>
    #   It's still a bit unclear why a considerable fraction of nearest genes
    #   don't have Ensembl IDs or symbols in the GTF. We'll only consider
    #   genes in the GTF for GO
    filter(!is.na(gene_for_go))

## Functional enrichment (GO/Pathways) for linked genes

for (this_cell_type in unique(peak_df$cell_type)) {
    gene_set = peak_df |>
        filter(cell_type == this_cell_type) |>
        pull(gene_for_go) |>
        unique()

    #   Note the universe here-- we're constraining nearest genes to those in
    #   the GTF
    ego = enrichGO(
        gene = gene_set, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
        universe = gtf$gene_id, pAdjustMethod= "BH", pvalueCutoff = 1,
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
