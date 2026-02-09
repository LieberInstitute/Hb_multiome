
## Scan DAR sequences for TF motifs, test enrichment vs background, and optionally footprint

library(tidyverse)
library(JASPAR2024)     # or latest available JASPAR set
library(TFBSTools)
library(motifmatchr)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # or species-appropriate
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg38)
library(clusterProfiler)
library(here)
library(sessioninfo)
library(GenomicRanges)
library(Signac)
library(Seurat)

peak_path = here(
    'processed-data', '06_peak_calling', '21_overlaping_FDRscores_TopHeatmap',
    'all_peaks_categorized.csv.gz'
)
seur_path = here(
    'processed-data', '06_peak_calling', '11_peaks_merge_MACS2',
    'Seurat_peaks_merged_cell_level_Mid_resolution.rds'
)
promoter_window = 2000  # +/- around TSS

## Inputs & goal

# Peak set: use merged MACS2 peaks (across samples) or Signac “peaks” assay. (all done)


## Confirm DARs are pseudobulk data (confirm) 

#   Consider filtering DARs by logFC


## Genomic annotation of DARs
# Add basic region context and nearest/overlapping genes

seur = readRDS(seur_path)

#   Grab unique DARs
peak_df = read_csv(peak_path, show_col_types = FALSE) |>
    filter(grepl('DAR', category)) |>
    select(peak_id, cell_type) |>
    distinct()
stopifnot(all(peak_df$peak_id %in% rownames(seur)))

#   RegionStats was already computed. Now attach TF motifs and test for
#   enrichment by cell type

jaspar_db = JASPAR2024()
pfm = getMatrixSet(
    jaspar_db@db, opts = list(species = 9606, collection = "CORE")
)
seur = AddMotifs(
    object = seur, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm
)

for (this_cell_type in unique(peak_df$cell_type)) {
    dar_peaks = peak_df |>
        filter(cell_type == this_cell_type) |>
        pull(peak_id)
}

## Maybe consider external data
# ENCODE cCREs, FANTOM5 enhancers, Vista, DHS, blacklist → annotate class (promoter/enhancer), confidence tiers.
# ChromHMM states if you have matched samples/tissues.

## get DAR and background (size/GC-matched) peak sets

dar_ids <- DARs %>% rownames()
bg_ids  <- setdiff(rownames(dar), dar_ids)

dar_gr <- peak_gr[names(peak_gr) %in% dar_ids]
bg_gr  <- peak_gr[names(peak_gr) %in% bg_ids]

pwm_list <- getMatrixSet(JASPAR2024, opts = list(species=9606, collection="CORE"))
matches_dar <- matchMotifs(pwm_list, dar_gr, genome=BSgenome.Hsapiens.UCSC.hg38)
matches_bg  <- matchMotifs(pwm_list, bg_gr,  genome=BSgenome.Hsapiens.UCSC.hg38)

# simple enrichment (Fisher) per motif
enrich <- lapply(seq_len(length(pwm_list)), function(i){
    a <- sum(assay(matches_dar))[i]; b <- length(dar_gr) - a
    c <- sum(assay(matches_bg))[i];  d <- length(bg_gr) - c
    ft <- fisher.test(matrix(c(a,b,c,d), nrow=2))
    data.frame(motif=names(pwm_list)[i], OR=ft$estimate, p=ft$p.value)
}) |> dplyr::bind_rows() |> mutate(FDR=p.adjust(p,"BH"))


## Functional enrichment (GO/Pathways) for linked genes

genes_for_go <- unique(na.omit(DARs_annot$link_gene))  # or SYMBOL from nearest/overlap
ego <- enrichGO(gene         = genes_for_go,
                OrgDb        = org.Hs.eg.db,
                keyType      = "SYMBOL",
                ont          = "BP",
                pAdjustMethod= "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2)
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
