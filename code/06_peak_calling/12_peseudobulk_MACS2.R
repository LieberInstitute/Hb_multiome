########################################################################
## Workflow for pseudobulk peaks and rna-count
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
source(here("code", "06_peak_calling", "atac_custom_functions", "atac_normalization_helpers.R"))
# ls()

#===============================================================================
# resolution_level = "Mid"      # 18 cell-types
#===============================================================================

p_met = "spearman"
w_size_num <- 5e5 
w_size = "5e5" # for filenames only 
resolution_level = "Mid" 

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

inputRDS_Dir <- here(
    "processed-data",
    "06_peak_calling",
    "02_link_peaks_MACS2",
    "Seurat_subsets_links_rds"
)
input_genes_filtered_file <- here(
    "processed-data",
    "06_peak_calling",
    "02_link_peaks_MACS2",
    "rna_filtered_genes_2perc_cells.csv"
)
output_Dir <- here(
    "processed-data",
    "06_peak_calling",
    "12_peseudobulk_MACS2"
)

if (!dir.exists(output_Dir)) {
    dir.create(output_Dir)
}

##==============================================================================
## Load Seurat / macs peaks / filtered genes. And make verification
 
Seurat_base_name <- "Seurat_peaks_macs2_cell_level_Mid_resolution.rds"
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
# 4    Astrocyte
# 5      MHb.1.2
# 6        LHb.1
# 7          OPC
# 8        Oligo
# 9    Microglia
# 10     LHb.2.7
# 11        Endo
# 12   LHb.1.3.4
# 13       MHb.1
# 14       MHb.2
# 15       MHb.3
# 16        Thal
# 17     LHb.1.3
# 18       LHb.7

df2 <- unique(SeuratOBJ[["orig.ident"]])
rownames(df2) <- NULL
df2
# orig.ident
# 1    S04_Hb_r
# 2    S05_Hb_r
# 3    S06_Hb_r
# 4    S10_Hb_r
# 5    S11_Hb_r
# 6    S12_Hb_r
# 7    S03_Hb_r
# 8    S07_Hb_r
# 9    S08_Hb_r
# 10   S09_Hb_r


## categories to pseudobulk data
grp_by_variables <- c("mid_cluster", "orig.ident")

# Subset peaks from the ATAC assay — get peak names
keep_peaks <- rownames(SeuratOBJ[["ATAC_macs2"]])

Seurat_pb <- AggregateExpression(
    SeuratOBJ,
    features = list(RNA = keep_genes, ATAC_macs2 = keep_peaks),
    assays =  c("RNA", "ATAC_macs2"),
    group.by = grp_by_variables,
    return.seurat = FALSE,
    verbose = TRUE
)   # matrix: peaks x clusters

#  Aggregated values are placed in the 'counts' layer of the returned object.
#  the data is then normalized by running NormalizeData on the aggregated counts. ScaleData is then run on the default assay before returning the object.

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
pb_atac_counts <- Seurat_pb$ATAC_macs2  # peaks x clusters
all(colnames(pb_rna_counts) == colnames(pb_atac_counts))   # clusters align

## verify
pb_rna_counts_df  <- as.data.frame(pb_rna_counts)   # genes x clusters
dim(pb_rna_counts_df) # [1] 15896   169
head(pb_rna_counts_df)

pb_atac_counts_df <- as.data.frame(pb_atac_counts)  # peaks x clusters
dim(pb_atac_counts_df) # [1] 355127    169
head(pb_atac_counts_df)

## build pseudobulk RNA Seurat and compute normalization on pseudobulk
pb_obj <- CreateSeuratObject(counts = pb_rna_counts, assay = "RNA")
pb_obj <- global_rebuild_rna_normalization(pb_obj, "RNA")
pb_obj
# An object of class Seurat 
# 15896 features across 169 samples within 1 assay 
# Active assay: RNA (15896 features, 2000 variable features)
# 3 layers present: counts, data, scale.data

## Create ChromatinAssay with correctly ordered peaks
# ensure peaks order & ranges match
gr_peaks <- granges(SeuratOBJ[["ATAC_macs2"]])
peak_ids <- Signac::GRangesToString(gr_peaks)
pb_atac_counts <- pb_atac_counts[peak_ids, , drop = FALSE]  # reorder to ranges

# Get the gene annotation from the original object's assay to reattach it later to the subset 
# orig_annot <- Annotation(SeuratOBJ[["ATAC_macs2"]])
orig_annot <- tryCatch(Annotation(SeuratOBJ[["ATAC_macs2"]]), error = function(e) NULL)
if (is.null(orig_annot)) {
    message("Original annotation not found. Building hg38 gene annotations via EnsDb...")
    suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))  # hg38-compatible EnsDb
    genes_ens <- genes(EnsDb.Hsapiens.v86)
    GenomeInfoDb::seqlevelsStyle(genes_ens) <- "UCSC"
    GenomeInfoDb::genome(genes_ens) <- "hg38"
    orig_annot <- genes_ens
}

pb_atac <- CreateChromatinAssay(
    counts     = pb_atac_counts,
    ranges     = gr_peaks,
    annotation = orig_annot 
)
pb_atac
# ChromatinAssay data with 355127 features for 169 cells
# Variable features: 0 
# Genome: 
#     Annotation present: FALSE 
# Motifs present: FALSE 
# Fragment files: 0 

## assign new chromatin to seurat
pb_obj[["ATAC_macs2_pseudo"]] <- pb_atac
Assays(pb_obj)
# [1] "RNA"               "ATAC_macs2_pseudo"
DefaultAssay(pb_obj) <- "ATAC_macs2_pseudo"

# Set genome
genome <- BSgenome.Hsapiens.UCSC.hg38

# Compute GC content for each peak
pb_obj <- RegionStats(
    object = pb_obj,
    assay = "ATAC_macs2_pseudo",
    genome = genome
)

message("GC content correction and normalization done!")

# ATAC normalization on pseudobulk: does TF-IDF/top-features/SVD
pb_obj <- global_rebuild_atac_normalization(pb_obj, "ATAC_macs2_pseudo")

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

# ## Milestone: save Seurat pseudobulk
# f_name <- paste0("Seurat_pseudobulk_macs2_peaks_no_merged_", resolution_level, "_resolution.rds")
# rds_name <- here(output_Dir, f_name)
# saveRDS(pb_obj, file = rds_name)

## ============================================================================/

message("Computing peak-gene correlations on pseudobulk ... ")

DefaultAssay(pb_obj) <- "RNA"  # not strictly required but conventional

# Default min.cells = 10 works fine for large clusters (>2,000 cells), but it’s too strict for tiny clusters (<200 cells)
# I scaled min.cells with cluster size (n_cells); require that at least 5% of cells in that cluster support the peak
# and never drop below 3 cells minimum, so the calculation always has some robustness

n_cells <- ncol(pb_obj)
# [1] 169
message("Processing ", n_cells, " pseudobulk groups")

min_cells_lp <- max(3, floor(0.05 * n_cells))  # 5% or at least 3

keep_genes <- intersect(keep_genes, rownames(pb_obj[["RNA"]]))

pb_obj <- LinkPeaks(
    object = pb_obj,
    peak.assay = "ATAC_macs2_pseudo",
    expression.assay = "RNA",
    genes.use = keep_genes,
    distance = w_size_num,
    min.cells = min_cells_lp,
    method = p_met 
)

message("Links found:")
nrow(Links(pb_obj[["ATAC_macs2_pseudo"]])) > 0

message("Pseudobulk link-peak-genes for  ", resolution_level, " resolution level done!")

## inspect data
head(Links(atac), n=3)

## Extract and save links
links_df <- as.data.frame(Links(pb_obj[["ATAC_macs2_pseudo"]]))
summary(link_df$score)

f_name <- paste0(resolution_level, "_peseudobulk_link_peak_genes.csv")
write.csv(
    link_df,
    file = here(output_cvsDir, f_name),
    row.names = FALSE
)


message("LinkPeaks saved: ", f_name)

## Milestone: save Seurat pseudobulk
f_name <- paste0("Seurat_pseudobulk_macs2_peaks_no_merged_", resolution_level, "_resolution.rds")
rds_name <- here(output_Dir, f_name)
saveRDS(pb_obj, file = rds_name)


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

