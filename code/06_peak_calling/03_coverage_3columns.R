########################################################################
## Plot Coverage Plots on WNN clusters
##
## Authors. CSC
## Date. March 24, 2025
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
## Note. Seurat objects were created with module load conda_R/4.3.x
########################################################################

library("Seurat")
library("Signac")
library("ggplot2")
library("tidyverse")
library("stringr")
library("here")

## input directories

here()

# Check/create directories

## clusters renamed for Spatial-Registration on Visium project
inputRDS_Dir <- here(
    "processed-data",
    "05_Clustering_ARCr",
    "08_wnn_gene_expression_plts_renamed_idents"
)
inputCVS_Dir <- here(
    "processed-data",
    "05_Clustering_ARCr",
    "02_Hb_celltypes_from_seurat_reanalyze_v3",
    "cvs_files_markers"
)
plotDir <- here(
    "plots",
    "06_peak_calling"
)

## Check directories
if (!dir.exists(plotDir)) {
    dir.create(plotDir)
}

## Load Seurat
# Use Seurat with clusters renamed for Spatial-Registration on Visium project
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium.rds"
seurat_name <- here(inputRDS_Dir, Seurat_base_name)
SeuratOBJ <- readRDS(here(inputRDS_Dir, Seurat_base_name))
levels(SeuratOBJ)

DefaultAssay(SeuratOBJ) <- "ATAC"
class(SeuratOBJ[["ATAC"]])
Seurat_base_name <- str_extract(seurat_name, regex("C\\.\\w+"))
# C.leiden_lsi_r2_renamed_visium


## Read DEG to prepare coverage plots of the top 5 genes highly expressed

DEG_file_name <- "WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_cellTypes_integrated_top50.csv"
DEG_file_name <- here(inputCVS_Dir, DEG_file_name)
dge_cluster_names <- read.csv(DEG_file_name)
dge_cluster_names <- dge_cluster_names |> drop_na(cell_type)
head(dge_cluster_names)
# p_val avg_log2FC pct.1 pct.2 p_val_adj cluster     gene            cell_type
# 1     0   4.336764 0.938 0.100         0       1 OTX2-AS1        DD_Inhib.Thal
# 2     0   3.849508 0.900 0.083         0       1      KIT        DD_Inhib.Thal
# 3     0   3.682782 0.927 0.122         0       1    MEIS2      LB_Thalamus/MDm
unique(dge_cluster_names$cluster)
# [1]  1  2  3  4  5  6  7  8  9 10 11 12 14 15 16 17 18 19 20 21 22 23 24 25 26
# [26] 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41

## Identified and subset clusters annotated for `habenula`.

message("Cluster-IDs from `WNN`")

SeuOBJ_clusters <- Idents(SeuratOBJ)
all_clusters <- unlist(levels(SeuOBJ_clusters))
unique(SeuratOBJ$seurat_clusters)

# Create 3 list with clusters to plot: LHb, MHb, None

lhb_clusters <- hb_clusters[grepl("LHb", all_clusters)]
lhb_cluster_ids <- str_extract(lhb_clusters, "(?<=C\\.)\\d+")

mhb_clusters <- hb_clusters[grepl("MHb", all_clusters)]
mhb_cluster_ids <- str_extract(mhb_clusters, "(?<=C\\.)\\d+")

noneHb_clusters <- hb_clusters[!grepl("MHb|LHb", all_clusters)]
noneHb_clusters_ids <- str_extract(noneHb_clusters, "(?<=C\\.)\\d+")


# Create 3 subsets of Seurat objects to plot: LHb, MHb, None

levels(Idents(SeuratOBJ))
# [1] "C.05 DD_LHb" "C.07 DD_MHb" "C.10 DD_MHb" "C.11 DD_MHb" "C.14 DD_MHb"
# [6] "C.16 DD_MHb" "C.18 DD_LHb" "C.23 DD_LHb" "C.24 DD_LHb" "C.30 DD_LHb"
# [11] "C.33 DD_LHb" "C.36 DD_MHb" "C.40 DD_LHb" "C.01"        "C.02"       
# [16] "C.03"        "C.04"        "C.06"        "C.08"        "C.09"       
# [21] "C.12"        "C.13"        "C.15"        "C.17"        "C.19"       
# [26] "C.20"        "C.21"        "C.22"        "C.25"        "C.26"       
# [31] "C.27"        "C.28"        "C.29"        "C.31"        "C.32"       
# [36] "C.34"        "C.35"        "C.37"        "C.38"        "C.39"       
# [41] "C.41"        "C.42"   

SeuratOBJ_LHb <- subset(SeuratOBJ, idents = lhb_clusters)
SeuratOBJ_MHb <- subset(SeuratOBJ, idents = mhb_clusters)
SeuratOBJ_None <- subset(SeuratOBJ, idents = noneHb_clusters)


## canonical genes to plot

features <- c("POU4F1")
# features <- c("POU4F1", "GPR151", "TAC3")

##### Build plt for column 1

plt1 <- CoveragePlot(
    object = SeuratOBJ_MHb,
    region = features,
    features = features,
    extend.upstream = 500,
    extend.downstream = 500,
    peaks = TRUE,
    links = TRUE
)
plt1 <- plt1 +
    labs(title = paste0("Clusters from WNN: ", seurat_name)) +
    theme(
        text = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(hjust = 0.5)
    )
print(plt1)



##### Build plt for column 1

## make first coverage plt for LHb clusters

## filter the top 5
# unique(dge_cluster_names$cluster)
# top5 <- dge_cluster_names |>
#     filter(cluster %in% hb_clusters_idx) |>
#     group_by(cluster) |>
#     top_n(n = 5, wt = avg_log2FC)
# head(top5)
# #     p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene    cell_type
# # <dbl>      <dbl> <dbl> <dbl>     <dbl>   <int> <chr>   <chr>
# # 1     0       3.04 0.768 0.145         0       5 RFTN1   DD_LHb
# # 2     0       3.03 0.798 0.194         0       5 CBLN2   DD_LHb
# # 3     0       3.37 0.73  0.129         0       5 GALR1   DD_LHb
# # 4     0       3.15 0.669 0.124         0       5 HTR4    DD_LHb
# # 5     0       3.38 0.956 0.457         0       5 COL25A1 DD_LHb


# top5_by_clust <- top5 |>
#     filter(cluster == 5) |>
#     group_by(cluster) |>
#     top_n(n = 5, wt = avg_log2FC)


# ## Prepare and save coverage plot
# 
# 
# for (clus in unique(top5$cluster)) {
#     # testing: clus = 5
#     tmp_name <- paste0(
#         Seurat_base_name,
#         "_PEAKS_hb-cluster-",
#         clus,
#         ".pdf"
#     )
#     
#     message("Processing habenula cluster: ", clus, "; Saved as: ", tmp_name)
#     
#     top5_cluster <- top5 |>
#         filter(cluster == clus)
#     
# #    pdf(file = here(plotDir, tmp_name))
#     
#     walk(
#         seq_along(top5_cluster$gene),
#         ~ {
#             tryCatch(
#                 {
#                     message(paste0("Processing gene ", top5_cluster$gene[.x]))
#                     
#                     features <- top5_cluster$gene[.x]
#                     plt1 <- CoveragePlot(
#                         object = SeuratOBJ,
#                         region = features,
#                         features = features,
#                         extend.upstream = 500,
#                         extend.downstream = 500,
#                         peaks = TRUE,
#                         links = TRUE
#                     )
#                     plt1 <- plt1 +
#                         labs(title = paste0("Clusters from WNN: ", seurat_name)) +
#                         theme(
#                             text = element_text(size = 8),
#                             axis.text.x = element_text(size = 7),
#                             axis.text.y = element_text(size = 7),
#                             plot.title = element_text(hjust = 0.5)
#                         )
#                     print(plt1)
#                 },
#                 error = function(e) {
#                     message(paste0(
#                         "Error occurred while processing gene ",
#                         top5_cluster$gene[.x],
#                         ": ",
#                         e$message
#                     ))
#                 }
#             )
#         }
#     )
#     
#     dev.off()
# }


message("Coverage plots completed")

# library("slurmjobs")
# job_single(
#   "01_coverage_basic",
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
