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
library("purrr")
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
input_macs_file <- here(
    "processed-data",
    "06_peak_calling",
    "01_call_peaks_MACS2"
)
cvsDir <- here(
    "processed-data",
    "06_peak_calling",
    "02_link_peaks_MACS2" # local peaks redo with Signac::CallPeaks()
)
output_RDS <- here(
    "processed-data",
    "06_peak_calling",
    "01_call_peaks_MACS2"
)

gene_peaks_csv <- here(input_macs_file, paste0("macs_peaks_Mid_resolution.csv"))

if (!dir.exists(cvsDir)) {
    dir.create(cvsDir)
}
if (!dir.exists(output_RDS)) {
    dir.create(output_RDS)
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

clusters <- levels(SeuratOBJ)

message("Seurat loaded and ready!")

clusters


##==============================================================================
## filter genes to those expressed in 2% of cells

rna_counts <- GetAssayData(SeuratOBJ, assay="RNA", layer="data")
length(rownames(rna_counts)) # [1] 36601
#length(rownames(rna_counts)[Matrix::rowSums(rna_counts > 0)]) # 34738
#length(rownames(rna_counts)[Matrix::rowSums(rna_counts > 0) > 0.02]) # 34738

# filter rna count expressed in at least 2% of the cells
keep_genes <- rownames(rna_counts)[Matrix::rowSums(rna_counts > 0) > 0.02 * ncol(rna_counts)]
length(keep_genes) # in count: [1] 15896
# Note. By default behavior on Seurat v5, after remove cells, ATAC assay was removed, which will be attached later
#SeuratOBJ_subset <- SeuratOBJ
#SeuratOBJ_subset <- SeuratOBJ_subset[keep_genes, ]  
# Convert GRanges to data frame and save for further analysis
head(keep_genes)

## save genes to avoid redundant step on other scripts
f_name <- paste0("rna_filtered_genes_2perc_cells.csv")
write.csv(
    keep_genes, 
    here(cvsDir, f_name),
    row.names = FALSE
)

message("Filtered gene names saved ...")

## make a list to store all Seurat subsets: cluster level

seurat_subsets <- list()

for (clust in clusters) {
    message("Subsetting cluster: ", clust)
    
    seurat_subsets[[clust]] <- subset(
        SeuratOBJ,
        idents = clust
    )
}

message("Seurat subsets by cell-type arranged: ", length(seurat_subsets))

##==============================================================================

# load link peak-gene csv

message("Loading ", basename(gene_peaks_csv))

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
class(peaks_gr)  # [1] "GenomicRanges"
head(peaks_gr)
length(peaks_gr) # 355127 -> mid resolution
head(peaks_gr)
# GRanges object with 6 ranges and 1 metadata column:
# seqnames        ranges strand |         peak_called_in
# <Rle>     <IRanges>  <Rle> |            <character>
# [1]     chr1 181329-181534      * |             Inhib.Thal
# [2]     chr1 191217-191619      * | OPC,Oligo,Inhib.Thal..
# [3]     chr1 629146-629354      * |              Astrocyte
# [4]     chr1 629811-630032      * | LHb.7,Astrocyte,Olig..
# [5]     chr1 630189-630389      * |             Oligo,Endo
# [6]     chr1 632189-632410      * |              Astrocyte


##==============================================================================

## Pre-processing: create new ATAC object from CallPeaks output fragments

## Set ATAC assay 
Assays(SeuratOBJ)
DefaultAssay(SeuratOBJ) <- "ATAC"

# Use *all* fragment files: meaning from all the samples 
## pick-up smallest cell-type size to run test:
# tmp <- purrr::map(seurat_subsets, ~ length(Cells(.x)))
# tmp
# # get the name of the object with minimum number of cells
# names(tmp)[which.min(tmp)]
# small_idx_to_test <- which.min(tmp)
# as.integer(small_idx_to_test)

frags_list <- Fragments(SeuratOBJ)
length(frags_list) # 10

## verification
stopifnot(length(frags_list) >= 1)
all_cells_frag <- unique(do.call(c, lapply(frags_list, Cells)))
length(all_cells_frag) # 55516
barcodes_in_frags <- unique(unlist(all_cells_frag, Cells(SeuratOBJ)))
mean(colnames(SeuratOBJ) %in% barcodes_in_frags)  # should be ~1.0

## Quantify peaks detected by Signac::CallPeaks() - not unified peaks
mat_macs2_peaks <- FeatureMatrix(
    fragments = frags_list,              # list of Fragment objects
    features  = peaks_gr,                # GRanges from CallPeaks
    cells     = colnames(SeuratOBJ),     # all barcodes across samples
    verbose   = TRUE
)

saveRDS(mat_macs2_peaks, file = "mat_macs2_peaks.rds")
# To load it back:
# loaded_features <- readRDS("mat_macs2_peaks.rds")

message("Chromatin counts saved ...")

# identical(colnames(mat_macs2_peaks), colnames(SeuratOBJ))
# head(mat_macs2_peaks)
# # > head(mat_macs2_peaks) # ge. Mid_level: 18 x 55,516 cells
# # 6 x 55516 sparse Matrix of class "dgCMatrix"
# # [[ suppressing 34 column names ‘S04_AAACAGCCAGAATGAC-1’, ‘S04_AAACAGCCAGCAAGGC-1’, ‘S04_AAACATGCACCTGGTG-1’ ... ]]

# create new ATAC assay with macs2 Granges
atac_macs2 <- CreateChromatinAssay(
    counts     = mat_macs2_peaks,
    ranges     = peaks_gr,
    annotation = tryCatch(Annotation(SeuratOBJ), error = function(e) NULL)
)
atac_macs2

## add the new chromatin object to existing Seurat 
SeuratOBJ[["ATAC_macs2"]] <- atac_macs2
SeuratOBJ
DefaultAssay(SeuratOBJ) <- "ATAC_macs2"

## milestone
f_name <- paste0("seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_", sufix_name, ".rds")
rds_name <- here(output_RDS, f_name)
saveRDS(SeuratOBJ, file = rds_name)

message("Seurat with new chromatin saved succesfully!")


# ## Calculate GC content for each peak and add it to the feature metadata
# SeuratOBJ <- RegionStats(
#     object = SeuratOBJ,
#     genome = genome,
#     assay = "ATAC"
# )
# colnames(SeuratOBJ)
# 
# message("GC content correction done!")
# 
# ## set a subset of genes to test
# #hb_cannonical_genes <- c("GPR151",  "POU4F1", "TAC3")
# 
# 
# ##==============================================================================
# ## Process "global" Link peak-genes
# 
# message("Computing global link-peaks correlations ...")
# 
# 
# ##==============================================================================
# ## Process "local" Link peak-genes (by cluster)
# 
# ## find peaks by cluster correlated with the expression of nearby genes 
# # access each subset by name
# 
# set.seed(22082025)
# 
# 
# for (seurat_cluster in names(seurat_subsets)) {
#     # seurat_cluster = "Endo"
#     
#     print(paste("Processing seurat cluster: ", seurat_cluster))
#     
#     seurat_subset <- seurat_subsets[[seurat_cluster]]
#     seurat_subset
#     
#     print(paste("Processing Peaks for ", seurat_cluster, "\n", 
#                 "Total cells found:", length(Cells(seurat_subset))))
#     
#     atac <- LinkPeaks(
#         object = seurat_subset,
#         peak.assay = "ATAC",
#         expression.assay = "RNA",
#         genes.use = keep_genes,
#         method = p_met,
#         distance = as.numeric(w_size)             # Only consider peaks within x kb of gene TSS
#     )
#     
#     print("Local link-peaks correlations completed!")
#     
#     ## inspect data
#     print(head(Links(atac), n=3))
# 
#     ## prepare data to save cvs
#     link_df <- as.data.frame(Links(atac))
#     print(summary(link_df$score))
#     
#     f_name <- paste0(resolution_level, "_", seurat_cluster, "_local_link_peak_genes", f_sufix, ".csv")
#     write.csv(
#         link_df,
#         file = here(cvsDir, f_name),
#         row.names = FALSE
#     )
#     print(paste("LinkPeaks saved: ", f_name))
#     
# }
# 
# 
# message("Local link-peaks correlations completed!")


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
