########################################################################
## EDA: Compute LinkPeaks() and filter High-Confident Peaks "PSEUDOBULK VERSION"
## INPUT: Peaks generated with CallPeaks() - MACS2
##
## CVS tables with links peaks "global" and "local" with
## - Spearman at 5e5 open-windows sized (check below details) 
##
## Authors. CSC
## Date. August 21, 2025
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
## Note. Seurat objects were created with module load conda_R/4.3.x
########################################################################

library("Seurat")
library("Signac")
library("BSgenome.Hsapiens.UCSC.hg38")  # full reference genome sequence / actual DNA bases (A/T/C/G) for each chromosome
library("GenomicRanges") 
library("tidyverse")
library("tidyr")
library("stringr")
library("here")

#===============================================================================
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
        "Processing job for resolution_level: ", resolution_level,
        "\nMethod: ", p_met,
        "\nWindow-size: ", w_size
    )
    f_sufix <- paste0(".", p_met, ".", w_size, ".cells_filtered_2perc")
} else {
    message("Input argument missed")
    stop()
}
f_sufix

# Check/create directories

## clusters renamed for Spatial-Registration on Visium project
inputRDS_Dir <- here(
  "processed-data",
  "05_Clustering_ARCr",
  "22_add_mid_level_clustering" # Seurat multiome final version with annotations
)
input_macs_file <- here(
    "processed-data",
    "06_peak_calling",
    "01_call_peaks_MACS2"
)

cvsDir <- here(
    "processed-data",
    "06_peak_calling",
    "02_link_peaks_pseudobulk_MACS2" # local peaks redo with Signac::CallPeaks()
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

# load link peak-gene csv

gene_peaks_csv <- here(input_macs_file, paste0("macs_peaks_Mid_resolution.csv"))

if (file.exists(gene_peaks_csv)) {
    peaks_df <- read.csv(gene_peaks_csv)
    message("File loaded!")
} else {
    stop(paste("File not found:", gene_peaks_csv))
}

colnames(peaks_df)
# [1] "seqnames"       "start"          "end"            "width"         
# [5] "strand"         "peak_called_in"

# CSV has columns like: seqnames, start, end, etc
peaks_gr <- makeGRangesFromDataFrame(
    peaks_df,
    keep.extra.columns = TRUE,       # keep additional metadata columns
    seqnames.field    = "seqnames",  # adjust if column name differs
    start.field       = "start",
    end.field         = "end"
)
# inspect
class(peaks_gr) # [1] "GenomicRanges"
head(peaks_gr)
# GRanges object with 6 ranges and 1 metadata column:
#     seqnames        ranges strand |         peak_called_in
# <Rle>     <IRanges>  <Rle> |            <character>
# [1]     chr1 181329-181534      * |             Inhib.Thal
# [2]     chr1 191217-191619      * | OPC,Oligo,Inhib.Thal..
# [3]     chr1 629146-629354      * |              Astrocyte
# [4]     chr1 629811-630032      * | LHb.7,Astrocyte,Olig..
# [5]     chr1 630189-630389      * |             Oligo,Endowhy a
# [6]     chr1 632189-632410      * |              Astrocyte
# -------
#     seqinfo: 34 sequences from an unspecified genome; no seqlengths
length(peaks_gr) # 355127
peaks_gr[1]
str(peaks_gr)

# After run CallPeaks() per cluster need to unify/re-quantify, so each cluster will have same peak set; And AggregateExpression() across the original assay will reflect those new peaks.
# reduce method will align the ranges and merge overlapping ranges to produce a simplified set.

# I have one GRanges that contains peaks from all clusters
union_peaks <- GenomicRanges::reduce(peaks_gr) 
length(union_peaks) # [1] 355127
head(union_peaks)

length(peaks_gr)            # 355127 - original count
length(union_peaks)         # 355127 - unified count (should be <= original)
any(width(union_peaks) <= 0)  # should be FALSE

# Note. No overlapping genomic ranges within your peaks_gr object that could be merged or "reduced" 

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

message("Pre-processing ready ...")


##==============================================================================

message("make a list to store all subsets: cluster level")

clusters <- levels(SeuratOBJ)

seurat_subsets <- list()

for (clust in clusters) {
    message("Subsetting cluster: ", clust)
    
    seurat_subsets[[clust]] <- subset(
        SeuratOBJ,
        idents = clust
    )
}

message("Seurat subsets by cell-type arrenged: ", length(seurat_subsets))


##==============================================================================
## Aggregate all cells in a cluster (or cell type) to get one accessibility profile per cluster

# Note. CallPeaks() run per cluster returns different GRanges per group.
# Aggregation across clusters needs a common peak matrix (same rows), otherwise per‑cluster sums aren’t comparable. 
# So, AggregateExpression() cannot be used directatly in the peaks dataset generated by CallPeaks()

# pb <- AggregateExpression(
#     SeuratOBJ,
#     assays =  c("RNA", "ATAC"),
#     group.by = "cluster_ann",
#     # layer = c("counts", "counts"), 
#     # return.seurat = TRUE
#     verbose = TRUE
# )   # matrix: peaks x clusters
# 
# #  Aggregated values are placed in the 'counts' layer of the returned object. 
# #  the data is then normalized by running NormalizeData on the aggregated counts. ScaleData is then run on the default assay before returning the object.
# Seurat_psedo_atac
# 
# pb_rna_counts  <- pb$RNA   # genes x clusters
# pb_atac_counts <- pb$ATAC  # peaks x clusters

##==============================================================================


message("Computing local link-peaks correlations ...")

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
    
    f_name <- paste0(resolution_level, "_pseudo_", seurat_cluster, "_local_link_peak_genes", f_sufix, ".csv")
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
