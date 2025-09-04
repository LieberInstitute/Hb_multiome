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

# p_met = "spearman"
# w_size = "5e5"
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


# Check/create directories

## clusters renamed for Spatial-Registration on Visium project
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


Seurat_pb <- AggregateExpression(
    SeuratOBJ,
    features = keep_genes,
    assays =  c("RNA", "ATAC_macs2"),
    group.by = grp_by_variables,
    return.seurat = FALSE,
    verbose = TRUE
)   # matrix: peaks x clusters

#  Aggregated values are placed in the 'counts' layer of the returned object.
#  the data is then normalized by running NormalizeData on the aggregated counts. ScaleData is then run on the default assay before returning the object.

## check assays and features 
Seurat_pb
Assays(Seurat_pb)

## Get the matrices directly, then add the ATAC matrix manually to the pseudobulk object
pb_rna_counts  <- Seurat_pb$RNA   # genes x clusters
pb_atac_counts <- Seurat_pb$ATAC_macs2  # peaks x clusters
all(colnames(pb_rna_counts) == colnames(pb_atac_counts))   # clusters align


pb_rna_counts_df  <- as.data.frame(pb_rna_counts)   # genes x clusters
pb_atac_counts_df <- as.data.frame(pb_atac_counts)  # peaks x clusters
dim(pb_rna_counts_df) # [1] 36601    18
head(pb_rna_counts_df)

head(pb_atac_counts_df)

# build a cluster-level object
pb_obj <- CreateSeuratObject(counts = pb_rna_counts, assay = "RNA")
pb_obj <- NormalizeData(pb_obj) # RNA log-normalize per cluster
pb_obj <- FindVariableFeatures(pb_obj, selection.method = "vst") 
pb_obj <- ScaleData(pb_obj, features = rownames(pb_obj))

# ensure peaks order & ranges match
gr_peaks <- granges(SeuratOBJ[["ATAC_macs2"]])
peak_ids <- Signac::GRangesToString(gr_peaks)
pb_atac_counts <- pb_atac_counts[peak_ids, , drop = FALSE]  # reorder to ranges
# Error in methods::slot(object = object, name = layer) : 
# no slot of name "chr1-181329-181534" for this object of class "Assay"

pb_atac <- CreateChromatinAssay(
    counts     = pb_atac_counts,
    ranges     = gr_peaks,
    annotation = tryCatch(Annotation(SeuratOBJ), error = function(e) NULL)
)

SeuratOBJ[["ATAC_macs2_pseudo"]] <- pb_atac

DefaultAssay(Seurat_pb) <- "ATAC_macs2_pseudo"

# ATAC normalization & bias covariates on pseudobulk
Seurat_pb <- global_rebuild_atac_normalization(SeuratOBJ, "ATAC_macs2_pseudo")

# verification
head(Seurat_pbj@meta.data)
# orig.ident nCount_RNA nFeature_RNA nCount_ATAC_unified
# Astrocyte  SeuratProject   16100320        29950            14066550
# Endo       SeuratProject    1435841        22262             1340147
# Excit.Thal SeuratProject  165408249        32995           121029625
# Inhib.Thal SeuratProject   55039199        30951            49965552
# LHb.1      SeuratProject   19120834        28831            13670681
# LHb.1.3    SeuratProject    2960861        22701             2011885
# nFeature_ATAC_unified
# Astrocyte                 351031
# Endo                      307849
# Excit.Thal                351037
# Inhib.Thal                351037
# LHb.1                     351015
# LHb.1.3                   330675

# f_name <- paste0(resolution_level, "_pb_rna_counts_df.csv")
# write.csv(
#     pb_rna_counts_df,
#     file = here(cvs_ouputDir, f_name),
#     row.names = FALSE
# )
# f_name <- paste0(resolution_level, "_pb_atac_counts.csv")
# write.csv(
#     pb_atac_counts_df,
#     file = here(cvs_ouputDir, f_name),
#     row.names = FALSE
# )

message("Pseudobulk Done!")


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

