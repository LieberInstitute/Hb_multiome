## Scan DAR sequences for TF motifs, test enrichment vs background, and optionally footprint

library(tidyverse)
library(JASPAR2024)
library(TFBSTools)
library(motifmatchr)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg38)
library(clusterProfiler)
library(here)
library(sessioninfo)
library(GenomicRanges)
library(Signac)
library(Seurat)
library(ggrepel)

peak_path = here(
    'processed-data', '06_peak_calling', '21_overlaping_FDRscores_TopHeatmap',
    'all_peaks_categorized.csv.gz'
)
seur_path = here('processed-data', '11_TF', 'minimal_seur.rds')
plot_dir = here('plots', '11_TF', '02_motif_enrichment')
promoter_window = 2000  # +/- around TSS
volcano_fold_cutoff = 1.5
volcano_fdr_cutoff = 0.05

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

## Inputs & goal

## Confirm DARs are pseudobulk data (confirm)

#   Consider filtering DARs by logFC

#   Load in and keep standard chromosomes only to match BSgenome
seur = readRDS(seur_path)
seur = seur[as.vector(grepl('^chr', seqnames(granges(seur)))), ]

#   Grab unique DARs on the standard chromosomes
peak_df = read_csv(peak_path, show_col_types = FALSE) |>
    filter(grepl('DAR', category), peak_id %in% rownames(seur)) |>
    select(peak_id, cell_type) |>
    distinct()
stopifnot(nrow(peak_df) > 0)

#   RegionStats was already computed. Now attach TF motifs and test for
#   enrichment by cell type

jaspar_db = JASPAR2024()
pfm = getMatrixSet(
    jaspar_db@db, opts = list(species = 9606, collection = "CORE")
)
seur = AddMotifs(
    object = seur, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm
)

motif_df_list = list()
for (this_cell_type in unique(peak_df$cell_type)) {
    dar_peaks = peak_df |>
        filter(cell_type == this_cell_type) |>
        pull(peak_id)

    motif_df_list[[this_cell_type]] = FindMotifs(seur, features = dar_peaks) |>
        mutate(cell_type = this_cell_type) |>
        as_tibble()
}
motif_df = bind_rows(motif_df_list)

p = motif_df |>
    mutate(
        neg_log_p = -log10(p.adjust),
        log_fold = log2(fold.enrichment),
        significant = (p.adjust < volcano_fdr_cutoff) &
            (fold.enrichment > volcano_fold_cutoff)
    ) |>
    ggplot(
        aes(
            x = log_fold, y = neg_log_p, color = significant, label = motif.name
        )
    ) +
    geom_point(alpha = 0.6, size = 2) +
    geom_text_repel(
        data = \(x) x |>
            filter(significant) |>
            group_by(cell_type) |>
            slice_max(fold.enrichment, n = 5),
        size = 5
    ) +
    scale_color_manual(values = c("grey50", "red")) +
    geom_hline(
        yintercept = -log10(volcano_fdr_cutoff), linetype = "dashed",
        color = "grey30"
    ) +
    geom_vline(
        xintercept = log2(volcano_fold_cutoff), linetype = "dashed",
        color = "grey30"
    ) +
    facet_wrap(~cell_type, ncol = 4, scales = "free") +
    theme_bw(base_size = 20) +
    labs(
        x = "Log2 Fold Enrichment", y = "-Log10 Adjusted P-Value",
        color = "Is Significant"
    )
pdf(file.path(plot_dir, "motif_volcano.pdf"), width = 16, height = 16)
print(p)
dev.off()

## Maybe consider external data
# ENCODE cCREs, FANTOM5 enhancers, Vista, DHS, blacklist → annotate class (promoter/enhancer), confidence tiers.
# ChromHMM states if you have matched samples/tissues.

## get DAR and background (size/GC-matched) peak sets



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
