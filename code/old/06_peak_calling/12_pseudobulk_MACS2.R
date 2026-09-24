########################################################################
## Pseudobulk multiome rna+atac on both regular peaks and merged peaks
## - (1) Chromatin assay with none merged peaks
## - (2) Chromatin assay with merged peaks with  GenomicRanges::reduce()
## 
## Authors. CSC
## Date. Sep 04, 2025
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
## Note. Seurat objects were created with module load conda_R/4.3.x
########################################################################

library("Seurat")
library("Signac")
## I use BSgenome.Hsapiens.UCSC.hg38 for extracting DNA motifs, k-mers, sequence-based features and compute Tn5 bias correction
library("BSgenome.Hsapiens.UCSC.hg38")  # full reference genome sequence
library("GenomicRanges") 
library("purrr")
library("tidyverse")
library("tidyr")
library("stringr")
library("here")


## ATAC function's helper used globally
source(here("code", "06_peak_calling", "multiome_custom_functions", "multiome_idents_normalization_helper.R"))
# ls()

#===============================================================================
# resolution_level = "Mid"      # 18 cell-types
#===============================================================================

## read input arguments
args = commandArgs(trailingOnly = TRUE)
peaks_ds <- args[2]
# peak_ds=(regular merged) 
if (is.na(peaks_ds) || !nzchar(peaks_ds)) stop("Missing peak_dataset argument")

p_met = "spearman"
w_size_num <- 5e5 
w_size = "5e5" # for filenames only 
resolution_level = "Mid" 
# peaks_ds = "merged"

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

## Check/create directories

# inputRDS_Dir <- here(
#     "processed-data",
#     "06_peak_calling",
#     "02_link_peaks_MACS2",
#     "Seurat_subsets_links_rds"
# )
## setup the correct seurat data
if (peaks_ds=="merged") {
    inputRDS_Dir <- here(
        "processed-data",
        "06_peak_calling",
        "11_peaks_merge_MACS2"
    )
    Seurat_base_name <- "Seurat_peaks_merged_cell_level_Mid_resolution.rds"
    ATAC_assay_name = "ATAC_macs2_merged"
} else {
    inputRDS_Dir <- here(
        "processed-data",
        "06_peak_calling",
        "12_pseudobulk_MACS2"
    )
    Seurat_base_name <- "Seurat_peaks_macs2_cell_level_Mid_resolution.rds"
    ATAC_assay_name = "ATAC_macs2"
}

input_genes_filtered_file <- here(
    "processed-data",
    "06_peak_calling",
    "02_link_peaks_MACS2",
    "rna_filtered_genes_2perc_cells.csv"
)
output_Dir <- here(
    "processed-data",
    "06_peak_calling",
    "12_pseudobulk_MACS2"
)

if (!dir.exists(output_Dir)) {
    dir.create(output_Dir)
}

##==============================================================================
## Load Seurat / macs peaks / filtered genes. And make verification

seurat_name <- here(inputRDS_Dir, Seurat_base_name)
SeuratOBJ <- readRDS(here(inputRDS_Dir, Seurat_base_name))

DefaultAssay(SeuratOBJ) <- "RNA"
levels(SeuratOBJ)
message("Existing Assays:")
Assays(SeuratOBJ)
# [1] "RNA"        "ATAC"       "ATAC_macs2"

clusters <- levels(SeuratOBJ)
message("Current clustets at ", resolution_level, " level")
clusters



##==============================================================================
## load filtered genes expressed at least in 2% of cells & peaks from CallPeaks()

if (file.exists(input_genes_filtered_file)) {
    keep_genes <- read.csv(input_genes_filtered_file)$x
    head(keep_genes)
    message("Loaded ", length(keep_genes), " pre-filtered genes ...")
    # Loaded 15896 pre-filtered genes ...
} else {
    stop(paste("File not found:", input_genes_filtered_file))
}


##==============================================================================
# AggregateExpression() by cluster on the non-unified peaks dataset 

## inspect data
colnames(SeuratOBJ@meta.data)
df <- unique(SeuratOBJ[["mid_cluster"]])
rownames(df) <- NULL
df
# mid_cluster
# 1   Excit.Thal
# 2        LHb.4
# 3   Inhib.Thal
# ...

df2 <- unique(SeuratOBJ[["orig.ident"]])
rownames(df2) <- NULL
df2
# orig.ident
# 1    S04_Hb_r
# 2    S05_Hb_r
# 3    S06_Hb_r
# ...

## extract peak annotation from original ds
peaks_gr <- granges(SeuratOBJ[[ATAC_assay_name]]) 
head(peaks_gr)
length(peaks_gr)
# [1] 355127 / [1] 351037

## categories to pseudobulk data
grp_by_variables <- c("mid_cluster", "orig.ident")

# Subset peaks from the ATAC assay — get peak names
keep_peaks <- rownames(SeuratOBJ[[ATAC_assay_name]])

# Build the features list dynamically
feature_list <- list(
    RNA = keep_genes
)
feature_list[[ATAC_assay_name]] <- keep_peaks

Seurat_pb <- AggregateExpression(
    SeuratOBJ,
    features = feature_list,
    assays   = c("RNA", ATAC_assay_name),
    group.by = grp_by_variables,
    return.seurat = FALSE,
    verbose = TRUE
)   # matrix: peaks x clusters
# Names of identity class contain underscores ('_'), replacing with dashes ('-')

## check assays and features 
head(Seurat_pb)
# $ATAC_macs2
# 355127 x 169 sparse Matrix of class "dgCMatrix"
# [[ suppressing 31 column names ‘Astrocyte_S03-Hb-r’, ‘Astrocyte_S04-Hb-r’, ‘Astrocyte_S05-Hb-r’ ... ]]
# [[ suppressing 31 column names ‘Astrocyte_S03-Hb-r’, ‘Astrocyte_S04-Hb-r’, ‘Astrocyte_S05-Hb-r’ ... ]]
# 
# chr1-181329-181534   .  1    1   .   3    2   .    1    .   1  .  1   .   .   .   .   1   .  .
# chr1-191217-191619  13  .    1   2  13    8   7   13   13   6  1  .   1   1   .   2   .   .  .
# chr1-629146-629354   .  . 1129   1   1    2   .    .    1   2  .  . 136   .   .   .   .   .  .
# chr1-629811-630032  75 15  193  14  96  315  31  261  387 281  9  9  27  28  17  29  22  23  7


## Get the matrices directly, then add the ATAC matrix manually to the pseudobulk object
pb_rna_counts  <- Seurat_pb$RNA   # genes x clusters
pb_atac_counts <- Seurat_pb[[ATAC_assay_name]]  # peaks x clusters
all(colnames(pb_rna_counts) == colnames(pb_atac_counts))   # clusters align

## verify
pb_rna_counts_df  <- as.data.frame(pb_rna_counts)   # genes x clusters
dim(pb_rna_counts_df) # [1] 15896   169
#head(pb_rna_counts_df)

pb_atac_counts_df <- as.data.frame(pb_atac_counts)  # peaks x clusters
dim(pb_atac_counts_df) # [1] 355127    169
#head(pb_atac_counts_df)

## build pseudobulk RNA Seurat and compute normalization on pseudobulk
pb_obj <- CreateSeuratObject(counts = pb_rna_counts, assay = "RNA")
pb_obj <- global_rebuild_rna_normalization(pb_obj, "RNA")
pb_obj
# An object of class Seurat 
# 15896 features across 169 samples within 1 assay 
# Active assay: RNA (15896 features, 2000 variable features)
# 3 layers present: counts, data, scale.data

## Create ChromatinAssay with correctly ordered peaks

## pull original GRanges from ATAC_macs assay
gr_peaks <- granges(SeuratOBJ[[ATAC_assay_name]])
peak_ids <- Signac::GRangesToString(gr_peaks)

## ensure peaks order & ranges match. Avoids errors if any peaks were filtered out
missing_pk <- setdiff(peak_ids, rownames(pb_atac_counts))
# character(0)
if (length(missing_pk) > 0) {
    message("Warning: ", length(missing_pk), " peaks in gr_peaks not in pb_atac_counts; dropping those peaks.")
    keep_ids <- intersect(peak_ids, rownames(pb_atac_counts))
    gr_peaks  <- gr_peaks[match(keep_ids, peak_ids)]
    peak_ids  <- keep_ids
}
pb_atac_counts <- pb_atac_counts[peak_ids, , drop = FALSE]  # reorder to ranges

## Get the original gene annotation to reattach it later
orig_annot <- tryCatch(Annotation(SeuratOBJ[[ATAC_assay_name]]), error = function(e) NULL)
if (is.null(orig_annot)) {
    message("Original annotation not found. Building hg38 gene annotations via EnsDb...")
    suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))  # hg38-compatible EnsDb
    genes_ens <- genes(EnsDb.Hsapiens.v86)
    GenomeInfoDb::seqlevelsStyle(genes_ens) <- "UCSC"
    GenomeInfoDb::genome(genes_ens) <- "hg38"
    orig_annot <- genes_ens
}

## build chromatin and run normalization -- I will later redo by subset =======
## - This is not strictly necessary, but it keeps consistency in saved objects

PSEUDO_ASSAY <- paste0(ATAC_assay_name, "_pseudo")

pb_atac <- CreateChromatinAssay(
    counts     = pb_atac_counts,
    ranges     = gr_peaks,
    annotation = orig_annot 
)
pb_atac
# ChromatinAssay data with 355127 features for 169 cells
# Variable features: 0 
# Genome: 
#     Annotation present: TRUE 
# Motifs present: FALSE 
# Fragment files: 0 

## assign new chromatin to seurat
pb_obj[[PSEUDO_ASSAY]] <- pb_atac         
DefaultAssay(pb_obj) <- PSEUDO_ASSAY  

# Set genome to compute GC correction 
genome <- BSgenome.Hsapiens.UCSC.hg38

# Compute GC content for each peak
pb_obj <- RegionStats(
    object = pb_obj,
    assay = PSEUDO_ASSAY,
    genome = genome
)

message("GC content correction and normalization done!")

# ATAC normalization on pseudobulk: does TF-IDF/top-features/SVD
pb_obj <- global_rebuild_atac_normalization(pb_obj, PSEUDO_ASSAY)

## =======/

# verification
head(pb_obj@meta.data)
# orig.ident nCount_RNA nFeature_RNA nCount_ATAC_macs2_pseudo
# Astrocyte_S03-Hb-r  Astrocyte    2347221        15577                  2256439
# Astrocyte_S04-Hb-r  Astrocyte     316714        14726                   420938
# Astrocyte_S05-Hb-r  Astrocyte    1395057        15813                  1302169
# Astrocyte_S06-Hb-r  Astrocyte     177741        13875                   508489
# Astrocyte_S07-Hb-r  Astrocyte    1997327        15859                  2126439
# Astrocyte_S08-Hb-r  Astrocyte    2330802        15848                  1835445
# nFeature_ATAC_macs2_pseudo
# Astrocyte_S03-Hb-r                     334869
# Astrocyte_S04-Hb-r                     179509
# Astrocyte_S05-Hb-r                     319322
# Astrocyte_S06-Hb-r                     194464
# Astrocyte_S07-Hb-r                     327433
# Astrocyte_S08-Hb-r                     324172

## ============================================================================/

## Milestone: save Seurat pseudobulk
if (peaks_ds=="merged") {
    f_name <- paste0(resolution_level, "_pseudobulk.", p_met, ".", w_size, "_merged_peaks.rds")
} else {
    f_name <- paste0(resolution_level, "_pseudobulk.", p_met, ".", w_size, ".rds")
}
rds_name <- here(output_Dir, f_name)
saveRDS(pb_obj, file = rds_name)

message("Pseudobulk done!!!")


## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

