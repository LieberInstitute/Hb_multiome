########################################################################
## EDA: Compute LinkPeaks() and filter High-Confident Peaks 
## INPUT: Peaks generated with CallPeaks() - MACS2
##
## CVS tables with links peaks "global" and "local" with
## - Spearman at 5e5 open-windows sized (check below details) 
##
## Authors. CSC
## Date. August 18, 2025
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
## Note. Seurat objects were created with module load conda_R/4.3.x
########################################################################

library("Seurat")
library("Signac")
## I use BSgenome.Hsapiens.UCSC.hg38 for extracting DNA motifs, k-mers, sequence-based features and compute Tn5 bias correction
library("BSgenome.Hsapiens.UCSC.hg38")  # full reference genome sequence / actual DNA bases (A/T/C/G) for each chromosome
library("tidyverse")
library("tidyr")
library("stringr")
library("here")

#===============================================================================
# resolution_level = "Fine"     # 42 clusters (small clusters - not run)
# resolution_level = "Broad"    # 8 cell-types
# resolution_level = "Mid"      # 18 cell-types
#===============================================================================

p_met = "spearman"
w_size = "5e5"

## read input arguments
args = commandArgs(trailingOnly = TRUE)
resolution_level <- args[2]

## for testing:
# resolution_level = "Mid" 

if (length(resolution_level)) {
    message(
        "Processing job for resolution_level:\n",
        resolution_level
    )
    f_sufix <- paste0(".", p_met, ".", w_size, ".cells_filtered_2perc")
} else {
    message("Input argument missed")
    stop()
}


# Check/create directories

## clusters renamed for Spatial-Registration on Visium project
inputRDS_Dir <- here(
  "processed-data",
  "05_Clustering_ARCr",
  "22_add_mid_level_clustering" # Seurat multiome final version with annotations
)
cvsDir <- here(
    "processed-data",
    "06_peak_calling",
    "02_link_peaks_MACS2" # local peaks redo with Signac::CallPeaks()
)

if (!dir.exists(cvsDir)) {
    dir.create(cvsDir)
}
   

##==============================================================================
## Set desired meta-data as current level: Mid or Broad

# Set Seurat identities from a metadata column
set_idents_from_meta <- function(seurat_obj, meta_col, level_order = NULL, na_fill = "Unknown") {
    ## double check level exist on meta-data
    if (!meta_col %in% colnames(seurat_obj@meta.data)) {
        stop("Meta column '", meta_col, "' not found in SeuratOBJ@meta.data")
    }
    # extract target vector
    target_vec <- as.character(seurat_obj[[meta_col]][, 1])
    
    # decide levels and keep appearance order
    if (is.null(level_order)) {
        level_order <- sort(unique(target_vec))
    }
    
    # only update if different from current Idents
    current_idents <- as.character(Idents(seurat_obj))
    if (!identical(current_idents, target_vec)) {
        seurat_obj <- SetIdent(seurat_obj, value = factor(target_vec, levels = level_order))
        message("Idents set from meta column '", meta_col, "'.")
    } else {
        message("Idents already match '", meta_col, "', nothing to do.")
    }
    
    return(seurat_obj)
}


##==============================================================================

## Load Seurat and make verification

# Use Seurat with clusters renamed for Spatial-Registration on Visium project
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_HD.rds"
seurat_name <- here(inputRDS_Dir, Seurat_base_name)

SeuratOBJ <- readRDS(here(inputRDS_Dir, Seurat_base_name))
SeuratOBJ
DefaultAssay(SeuratOBJ) <- "RNA"

# colnames(SeuratOBJ@meta.data)
msg <- switch(resolution_level,
    "Fine" = paste0("Fine-level Cell-Types:\n ", paste(sort(unique(SeuratOBJ$cluster_ann)), collapse = "\n")),
    "Broad" = paste0("Broad-level Cell-Types:\n", paste(sort(unique(SeuratOBJ$merged_cluster)), collapse = "\n")),
    "Mid" = paste0("Mid-level Cell-Types:\n", paste(sort(unique(SeuratOBJ$mid_cluster)), collapse = "\n"))
)

message(msg)
message("Processing ", length(Cells(SeuratOBJ)), " cells")

## set resolution_level
meta_col <- case_when(
    resolution_level=="Fine" ~ "cluster_ann",
    resolution_level=="Broad" ~ "merged_cluster",
    resolution_level=="Mid" ~ "mid_cluster"
)

## set desired idents as current level
SeuratOBJ <- set_idents_from_meta(SeuratOBJ, meta_col = meta_col)
levels(SeuratOBJ)

message("Seurat loaded and ready!")


##==============================================================================
## filter genes to those expressed in 2% of cells

rna_counts <- GetAssayData(SeuratOBJ, assay="RNA", layer="data")
length(rownames(rna_counts)) # [1] 36601
#length(rownames(rna_counts)[Matrix::rowSums(rna_counts > 0)]) # 34738
#length(rownames(rna_counts)[Matrix::rowSums(rna_counts > 0) > 0.02]) # 34738

# filter rna count expressed in at least 2% of the cells
keep_genes <- rownames(rna_counts)[Matrix::rowSums(rna_counts > 0) > 0.02 * ncol(rna_counts)]
length(keep_genes) # in count: [1] 14526
# Note. By default behavior on Seurat v5, after remove cells, ATAC assay was removed, which will be attached later
SeuratOBJ_subset <- SeuratOBJ
SeuratOBJ_subset <- SeuratOBJ_subset[keep_genes, ]  
Assays(SeuratOBJ_subset)

##==============================================================================

## Pre-processing to identifies cis-regulatory elements by linking chromatin-accessible peaks

## Set ATAC assay 
DefaultAssay(SeuratOBJ) <- "ATAC"
class(SeuratOBJ[["ATAC"]])

message("Chromatin loaded!")

## pre-processed peaks / QC
message("Starting GC content correction ... ")

## GC content correction
genome <- BSgenome.Hsapiens.UCSC.hg38

SeuratOBJ <- RegionStats(
    object = SeuratOBJ,
    genome = genome,
    assay = "ATAC"
)
SeuratOBJ

message("GC content correction done!")

# ## For documentation purposes, I fix this chunk in case we need to filter peaks on Seurat v5
# ## In this case, we have only Std Chromosomes - just skip it 
# 
# table(seqnames(granges(SeuratOBJ)))
# 
# # pull assay
# atac <- SeuratOBJ[["ATAC"]]     # Get ChromatinAssay 
# gr_all <- granges(SeuratOBJ)        # GRanges of ATAC features (peaks)
# frag_list  <- Fragments(SeuratOBJ)              # carry fragment(s)
# annot <- tryCatch(Annotation(SeuratOBJ), error = function(e) NULL)
# 
# # normalize
# seqlevelsStyle(gr_all) <- "UCSC"    # normalize naming style
# gr_std <-  keepStandardChromosomes(gr_all, pruning.mode = "coarse")
# 
# # give GRanges canonical IDs that match the assay rownames - features Seurat uses
# ids_from_assay <- rownames(atac)
# ids_from_gr <- Signac::GRangesToString(gr_std)
# # build the keep list by intersecting with the assay’s rownames
# peaks.keep <- intersect(ids_from_gr, ids_from_assay)
# # check returned features > 0
# length(peaks.keep) # 262951
# stopifnot(length(peaks.keep) > 0) 
# 
# # ensure the GRanges names match the assay IDs (prefer the assay’s own IDs)
# if (length(peaks.keep) == length(ids_from_assay)) {
#     message("All peaks are already on standard chromosomes; no subsetting needed.")
# } else {
#     idx <- match(peaks.keep, ids_from_assay)                # integer indices
#     counts_mat <- GetAssayData(SeuratOBJ, assay = "ATAC", layer = "counts")[idx, , drop = FALSE]
# }
# 
# head(SeuratOBJ[["ATAC"]]@meta.features, n=3)
# # count percentile AA AC AG AT CA  CC  CG  CT GA  GC  GG GT TA
# # chr1-180813-181799   893  0.6523116 48 67 76  6 59 115 107  67 49 142 113 28 41
# # chr1-182478-183337   137  0.0943788 24 32 86 28 66  77  19  71 63  73  90 53 17
# # chr1-183785-184772   592  0.5414279 37 50 62 34 73 131  16 107 51  72  61 44 22
# # TC TG TT GC.percent sequence.length
# # chr1-180813-181799 24 36  8   68.99696             987
# # chr1-182478-183337 51 84 25   59.65116             860
# # chr1-183785-184772 74 89 64   56.17409             988
# 
# # SeuratOBJ[["ATAC"]] <- SeuratOBJ[["ATAC"]][peaks.keep, ] --> easy step fails
# 
# # Subset counts *matrix* from the assay
# counts_mat <- GetAssayData(SeuratOBJ, assay = "ATAC", layer = "counts")[peaks.keep, , drop = FALSE]
# 
# # Align ranges to the same order as counts
# names(gr_all) <- Signac::GRangesToString(gr_all)
# new_ranges <- gr_all[peaks.keep]
# 
# ## Build a new ChromatinAssay with counts + ranges (+ fragments/annotation)
# subsetted_atac_assay <- CreateChromatinAssay(
#     counts     = counts_mat,
#     ranges     = new_ranges,
#     fragments  = if (length(frag_list) > 0) frag_list else NULL,
#     annotation = annot
# )
# 
# # carry over per-peak meta.features for kept peaks (same row order!)
# mf_old <- tryCatch(atac@meta.features, error = function(e) NULL)
# if (!is.null(mf_old)) {
#     mf_new <- mf_old[peaks.keep, , drop = FALSE]
#     subsetted_atac_assay@meta.features <- mf_new
# }
# 
# # Replace the original "ATAC" assay with the subsetted one
# SeuratOBJ[["ATAC"]] <- subsetted_atac_assay
# DefaultAssay(SeuratOBJ) <- "ATAC"

## Re-atach chromatin 
SeuratOBJ_subset[["ATAC"]]  <- SeuratOBJ[["ATAC"]]
SeuratOBJ <- SeuratOBJ_subset 

## set a subset of genes to test
#hb_cannonical_genes <- c("GPR151",  "POU4F1", "TAC3")

message("Pre-processing ready ...")


##==============================================================================
## Process "global" Link peak-genes

message("Computing global link-peaks correlations ...")

## find peaks that are correlated with the expression of nearby genes 
atac <- LinkPeaks(
    object = SeuratOBJ,
    peak.assay = "ATAC",
    expression.assay = "RNA",
    genes.use = keep_genes,
    method = p_met,
    distance = as.numeric(w_size)             # Only consider peaks within ±500 kb of gene TSS (cis-window)
)

message("Global link-peaks correlations completed!")

## inspect data
head(Links(atac), n=3)
# GRanges object with 5 ranges and 5 metadata columns:
#     seqnames              ranges strand |     score        gene
#        <Rle>           <IRanges>  <Rle> | <numeric> <character>

link_df <- as.data.frame(Links(atac))
summary(link_df$score)

write.csv(
    link_df,
    file = here(cvsDir, paste0(resolution_level, "_global_link_peak_genes", f_sufix, ".csv")),
    row.names = FALSE
)

message("Global peaks saved!")



##==============================================================================
## Process "local" Link peak-genes (by cluster)

## find peaks by cluster correlated with the expression of nearby genes 

clusters <- levels(SeuratOBJ)

# make a list to store all subsets: cluster level
seurat_subsets <- list()

for (clust in clusters) {
    message("Subsetting cluster: ", clust)
    
    seurat_subsets[[clust]] <- subset(
        SeuratOBJ,
        idents = clust
    )
}

# now access each subset by name

for (seurat_cluster in names(seurat_subsets)) {
    # seurat_cluster = "Endo"
    
    message("Processing seurat cluster: ", seurat_cluster)
    
    seurat_subset <- seurat_subsets[[seurat_cluster]]
    seurat_subset
    
    message("Total cells in cluster ", seurat_cluster, ": ", length(Cells(seurat_subset)))
    
    atac <- LinkPeaks(
        object = seurat_subset,
        peak.assay = "ATAC",
        expression.assay = "RNA",
        genes.use = keep_genes,
        method = p_met,
        distance = as.numeric(w_size)             # Only consider peaks within x kb of gene TSS
    )
    
    message("Local link-peaks correlations completed!")
    ## inspect data
    head(Links(atac), n=3)

    ## prepare data to save cvs
    link_df <- as.data.frame(Links(atac))
    print(summary(link_df$score))
    
    f_name <- paste0(resolution_level, "_", seurat_cluster, "_local_link_peak_genes", f_sufix, ".csv")
    write.csv(
        link_df,
        file = here(cvsDir, f_name),
        row.names = FALSE
    )
    message("LinkPeaks saved: ", f_name)
    
}


message("Local link-peaks correlations completed!")


message("All done!!!")


# library("slurmjobs")
# job_single(
#   "00_link_peaks",
#   create_shell = TRUE,
#   partition = "katun",
#   memory = "30G",
#   cores = 2,
#   logdir = "logs",
#   command = "Rscript -e \"options(width = 120); sessioninfo::session_info()\"",
#   create_logdir = TRUE
# )

## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
