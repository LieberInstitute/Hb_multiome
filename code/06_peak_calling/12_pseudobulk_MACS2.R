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
source(here("code", "06_peak_calling", "multiome_custom_functions", "multiome_idents_normalization_helper.R"))
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
    "12_pseudobulk_MACS2"
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
peaks_gr <- granges(SeuratOBJ[["ATAC_macs2"]]) 
head(peaks_gr)
length(peaks_gr)
# [1] 355127

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
pb_atac_counts <- Seurat_pb$ATAC_macs2  # peaks x clusters
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
gr_peaks <- granges(SeuratOBJ[["ATAC_macs2"]])
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
orig_annot <- tryCatch(Annotation(SeuratOBJ[["ATAC_macs2"]]), error = function(e) NULL)
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
pb_obj[["ATAC_macs2_pseudo"]] <- pb_atac
Assays(pb_obj)
# [1] "RNA"               "ATAC_macs2_pseudo"
DefaultAssay(pb_obj) <- "ATAC_macs2_pseudo"

# Set genome to compute GC correction 
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
f_name <- paste0(resolution_level, "_pseudobulk.", p_met, ".", w_size, ".rds")
rds_name <- here(output_Dir, f_name)
saveRDS(pb_obj, file = rds_name)

message("Pseudobulk done!!!")


## ============================================================================/
## Next chunk with peak-gene correlations on pseudobulk works
## but as it requires a lot of mem resources I split the process by ct on the script:
## - 12_pseudobulk_LinkPeaks_MACS2_split_ct.R
## ============================================================================/

# message("Computing peak-gene correlations on pseudobulk ... ")
# 
# DefaultAssay(pb_obj) <- "RNA"  # not strictly required but conventional
# 
# # Default min.cells = 10 works fine for large clusters (>2,000 cells), but it’s too strict for tiny clusters (<200 cells)
# # I scaled min.cells with cluster size (n_samples); require that at least 5% of cells in that cluster support the peak
# # and never drop below 3 cells minimum, so the calculation always has some robustness
# 
# n_samples <- ncol(pb_obj)
# message("Processing ", n_samples, " pseudobulk groups")
# # Processing 169 pseudobulk groups
# 
# ## ensure genes.use exist in RNA assay and in annotation gene_name
# genes_in_rna <- rownames(pb_obj[["RNA"]])
# genes_in_annot <- if ("gene_name" %in% colnames(mcols(orig_annot))) unique(orig_annot$gene_name) else genes_in_rna
# genes_for_lp <- intersect(keep_genes, intersect(genes_in_rna, genes_in_annot))
# stopifnot(length(genes_for_lp) > 0)
# 
# 
# ##==============================================================================
# ## Compute cell-type specific (local) Link peak-genes
# 
# message("Making list of Seurat subsets: ", resolution_level)
# 
# # Parse cluster label from pb; like "Astrocyte_S03-Hb-r"
# if (!"cluster_pb" %in% colnames(pb_obj@meta.data)) {
#     # keep everything before the last "-" or "_"
#     pb_obj$cluster_pb <- sub("[-_].*$", "", colnames(pb_obj))
#     # "LHb.1.3.4_S07-Hb-r" - "LHb.1.3.4"
# }
# pb_obj$cluster_pb <- trimws(pb_obj$cluster_pb)
# stopifnot(!any(is.na(pb_obj$cluster_pb)), length(pb_obj$cluster_pb) == ncol(pb_obj))
# 
# # set idents
# Idents(pb_obj) <- pb_obj$cluster_pb
# clusters_pb <- levels(Idents(pb_obj))
# message("Clusters in pseudobulk: ", paste(clusters_pb, collapse = ", "))
# 
# # build per-cluster subsets
# seurat_subsets <- setNames(vector("list", length(clusters_pb)), clusters_pb)
# for (clust in clusters_pb) {
#     seurat_subsets[[clust]] <- subset(
#         pb_obj,
#         idents = clust
#     )
# }
# print(table(Idents(pb_obj)))
# 
# # Tag peaks with their calling clusters from the original per-cell GRanges
# peaks_gr$peak_id <- Signac::GRangesToString(peaks_gr)   # "chr-start-end"
# 
# peaks_by_cluster <- function(cluster) {
#     #finds the indices of rows where the cluster name is present / returns a vector of peak IDs 
#     hits <- grepl(paste0("(^|,)", cluster, "(,|$)"), peaks_gr$peak_called_in)
#     peaks_gr$peak_id[hits]
# }
# 
# message("Seurat subsets by cell-type: ", length(seurat_subsets))
# message("Starting LinkPeaks by cluster ... ")
# 
# set.seed(09062025)
# 
# links_list <- vector("list", length(seurat_subsets))
# names(links_list) <- names(seurat_subsets)
# 
# for (seurat_cluster in names(seurat_subsets)) {
#     # seurat_cluster = names(seurat_subsets)[2]
#     
#     seurat_subset <- seurat_subsets[[seurat_cluster]]
#     DefaultAssay(seurat_subset) <- "ATAC_macs2_pseudo"
#     
#     n_samples <- ncol(seurat_subset)
#     message("Processing ", unique(Idents(seurat_subset)), " (", n_samples, " samples)")
# 
#     if (n_samples < 3) {
#         message("Skipping ", seurat_cluster, " (only ", n_samples, " pseudobulk samples).")
#         next
#     }
#     
#     # Peaks called in this cluster AND present in the pseudobulk ATAC assay
#     cluster_peaks <- intersect(
#         peaks_by_cluster(seurat_cluster),
#         rownames(seurat_subset[["ATAC_macs2_pseudo"]])
#     )
#     if (length(cluster_peaks) == 0) {
#         message("No MACS2 peaks found in pseudobulk assay for ", seurat_cluster, "; skipping.")
#         next
#     }
#     
#     # Subset to peaks with enough support within this cluster’s samples
#     counts_data <- GetAssayData(
#         object = seurat_subset, 
#         assay = "ATAC_macs2_pseudo", 
#         slot = "counts"
#     )[cluster_peaks, , drop = FALSE]
#     
#     min_cells_sub <- max(3, floor(0.05 * n_samples))
#     
#     keep_peaks_sub <- rownames(counts_data)[Matrix::rowSums(counts_data > 0) >= min_cells_sub]
#     if (length(keep_peaks_sub) == 0) { 
#         message("No peaks pass support filter in ", seurat_cluster, "; skipping.")
#         next 
#     }    
# 
#     # Rebuild ChromatinAssay for the kept peaks with Ann (safer way)
#     peak_ranges <- Signac::StringToGRanges(keep_peaks_sub)
#     counts_keep <- counts_data[keep_peaks_sub, , drop = FALSE]
#     
#     # ensure rownames <-> ranges order match
#     stopifnot(identical(rownames(seurat_subset[["ATAC_macs2_pseudo"]]),
#                         Signac::GRangesToString(granges(seurat_subset[["ATAC_macs2_pseudo"]]))))
#     # ensure annotation exists
#     stopifnot(!is.null(Annotation(seurat_subset[["ATAC_macs2_pseudo"]])))
#     
#     new_assay <- CreateChromatinAssay(
#         counts = counts_keep, 
#         ranges = peak_ranges,
#         annotation = orig_annot 
#     )
#     
#     # Add the new assay back to the Seurat object
#     seurat_subset[["ATAC_macs2_pseudo"]] <- new_assay
# 
#     # Set the default assay again since this was replaced it
#     DefaultAssay(seurat_subset) <- "ATAC_macs2_pseudo"
#     
#     # Compute GC content for each peak
#     seurat_subset <- RegionStats(
#         object = seurat_subset,
#         assay = "ATAC_macs2_pseudo",
#         genome = genome
#     )
#     
#     message("GC content correction and normalization done!")
#     
#     seurat_subset <- global_rebuild_atac_normalization(seurat_subset, "ATAC_macs2_pseudo")
# 
#     # Genes present in this subset (use per-subset genes, not global)
#     genes_sub <- intersect(genes_for_lp, rownames(seurat_subset[["RNA"]]))
#     if (length(genes_sub) == 0) {
#         message("No genes left after filtering for ", seurat_cluster, "; skipping.")
#         next
#     }
#     
#     message(" Cell-type=", seurat_cluster,
#             "\n Peaks=", length(keep_peaks_sub),
#             "\n Genes=", length(genes_sub),
#             "\n min.cells=", min_cells_sub,
#             "\n distance=", w_size_num)
#     
#     #=====/
#     
#     atac <- Signac::LinkPeaks(
#         object = seurat_subset,
#         peak.assay = "ATAC_macs2_pseudo",
#         expression.assay = "RNA",
#         genes.use = genes_sub,
#         distance = w_size_num,
#         min.cells = min_cells_sub,
#         method = p_met 
#     )
#  
#     lk <- Links(seurat_subset[["ATAC_macs2_pseudo"]])
#     
#     if (length(lk) > 0) {
#         # Save per-cluster links
#         mcols(lk)$cluster <- seurat_cluster
#         links_list[[seurat_cluster]] <- lk
#         f_name <- here(
#             output_Dir, 
#             paste0(resolution_level, "_", seurat_cluster, "_pseudobulk_link_peak_genes.", p_met, ".", w_size, ".csv")
#             )
#         write.csv(as.data.frame(lk), f_name, row.names = FALSE)
#         message("links: ", length(lk), " (saved: ", basename(f_name), ")")
#     } else {
#         message("no links for ", seurat_cluster)
#     }
#     
#     # Save the subset if you like (kept from your script)
#     f_name <- paste0(resolution_level, "_", seurat_cluster, "_seurat_subset.rds")
#     saveRDS(seurat_subset, file = here(output_Dir, f_name))
#     
#     message("Seurat subset and LinkPeaks saved!")
# 
# }
# 
# ##=============================================================================/
# 
# message("Pseudobulk LinkPeaks for  ", resolution_level, " resolution level done!")
# 
# ## Extract and save ALL links
# links_all <- do.call(c, links_list[ lengths(links_list) > 0 ])
# 
# if (length(links_all) > 0) {
#     
#     Links(pb_obj[["ATAC_macs2_pseudo"]]) <- unique(links_all)
#     f_name <- here(output_Dir, paste0(resolution_level, "_pseudobulk_LinkPeaks_perCluster.ALL.", p_met, ".", w_size, ".csv"))
#     write.csv(as.data.frame(Links(pb_obj[["ATAC_macs2_pseudo"]])), f_name, row.names = FALSE)
#     message("Combined per-cluster links: ", length(Links(pb_obj[["ATAC_macs2_pseudo"]])))
#     message("ALL LinkPeaks saved: ", f_name)
#     
# }
# 
#
# message("All done!!!")


## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

