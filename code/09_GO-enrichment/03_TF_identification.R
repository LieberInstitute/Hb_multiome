
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


## Inputs & goal

## Pre-processing steps completed:
# Filter cells (both assays): min fragments in peaks, TSS enrichment, blacklist removal, nucleosome signal; RNA nFeature/nCount/mito%. (all done)
# Normalization: RNA (SCTransform/LogNormalize); ATAC (TF-IDF + SVD/LSI). (all done)
# WNN integration & clustering: derive final clusters/cell-types to compare. (all done)
# Peak set: use merged MACS2 peaks (across samples) or Signac “peaks” assay. (all done)


## Confirm DARs are pseudobulk data (confirm) 

## Keep significant DARs (tune thresholds)
# DARs <- dar %>% filter(FDR <= 0.05, abs(logFC) >= 0.25)
# (ArchR: getMarkerFeatures with useGroups and bgdGroups, then markerTest.)


## Genomic annotation of DARs
# Add basic region context and nearest/overlapping genes

peak_gr <- StringToGRanges(rownames(DARs), sep=c("-", "-"))  # "chr-start-end" to GRanges

ann <- annotatePeak(
    peak_gr,
    TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,
    tssRegion=c(-2000, 2000),  # promoter window (adjust as needed)
    annoDb="org.Hs.eg.db"
) %>% as.data.frame()


DARs_annot <- DARs %>%
    tibble::rownames_to_column("peak_id") %>%
    left_join(ann %>% 
                  transmute(peak_id=paste(seqnames, start, end, sep="-"),
                            annotation, geneId=geneId, SYMBOL=SYMBOL, 
                            distanceToTSS=distanceToTSS))

## Maybe consider external data
# ENCODE cCREs, FANTOM5 enhancers, Vista, DHS, blacklist → annotate class (promoter/enhancer), confidence tiers.
# ChromHMM states if you have matched samples/tissues.

## We could merge with Peak-to-gene linkage, but its not considered by now


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

## Minimal skeleton
annotate_DARs_multiome <- function(obj, contrast_celltype, fdr_atac=0.05, lfc_atac=0.25) {
    # assume we have merged peaks
    pb_list <- make_pseudobulk(obj, assay="peaks", group=c("sample_id","cell_type"))
    dar_tbl <- run_edgeR(pb_list, contrast_celltype, fdr=fdr_atac, lfc=lfc_atac)
    dar_tbl <- add_genomic_annotations(dar_tbl, genome="hg38")
    links   <- link_peaks_to_genes(obj, distance=5e5)
    dar_tbl <- dar_tbl |> left_join(links, by="peak_id")
    degs    <- get_deg_table(obj, contrast_celltype)
    dar_tbl <- dar_tbl |> left_join(degs, by=c("SYMBOL"="gene"))
    dar_tbl <- add_concordance(dar_tbl)
    dar_tbl <- add_motif_enrichment(dar_tbl, genome="hg38", db="JASPAR2024")
    dar_tbl
}


## Desired Deliverables: 
# DARs_<celltype>_vs_rest.csv (full table)
# DARs_<celltype>_topN.bed (browser tracks)
# motif_enrichment_<celltype>.csv
# GO_<celltype>_linked_genes.csv
# QC plots: MA/volcano, region annotations pie/bar, distance-to-TSS, motif volcano, GO dotplot, link distance distribution.
