########################################################################
## Search Peaks on ATAC modality
##
## Authors. CSC
## Date. June 13, 2025
## Recommended resources on interactive mode: srun --pty --mem=80GB --x11 bash
## Note. Seurat objects were created with module load conda_R/4.3.x
########################################################################

library("Seurat")
library("Signac")
library("ggplot2")
#library("tidyverse")
library("here")

# Check/create directories
## clusters renamed for Spatial-Registration on Visium project

inputRDS_Dir <- here(
    "processed-data",
    "05_Clustering_ARCr",
    #"08_wnn_gene_expression_plts_renamed_idents"
    "05_rename_idents"
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

#Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium.rds"
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2.rds"
seurat_name <- here(inputRDS_Dir, Seurat_base_name)
SeuratOBJ <- readRDS(here(inputRDS_Dir, Seurat_base_name))


Seurat_subset <- subset(SeuratOBJ, idents = grep("MHb|LHb", Idents(SeuratOBJ), value = TRUE))


DefaultAssay(SeuratOBJ) <- "ATAC"
class(SeuratOBJ[["ATAC"]])
SeuratOBJ[["ATAC"]]
# ChromatinAssay data with 262951 features for 55702 cells
# Variable features: 249866 
# Genome: 
#     Annotation present: TRUE 
# Motifs present: FALSE 
# Fragment files: 10 

## extract base_name to save plots 
Seurat_base_name <- str_extract(seurat_name, regex("C\\.\\w+")) 
Seurat_base_name
# C.leiden_lsi_r2_renamed_visium


##==============================================================================
message("Searching peaks in ATAC modality ...")
message("Searching on Habenula ATAC clusters ...")

# Store peak sets per cluster
cluster_peaks <- list()

for (clus in clusters) {
    message("Processing cluster: ", clus)
    
    # Subset cells from one cluster
    cells_in_cluster <- WhichCells(SeuratOBJ, idents = clus)
    
    # Subset ATAC counts matrix
    atac_counts <- GetAssayData(SeuratOBJ, assay = "ATAC", slot = "counts")[, cells_in_cluster]
    
    # Find peaks with any signal in this cluster
    peaks_present <- rownames(atac_counts)[Matrix::rowSums(atac_counts) > 0]
    
    # Store
    cluster_peaks[[clus]] <- peaks_present
}


##==============================================================================
## Preparing to compute pair peaks comparison in the Habenula clusters

all_cluster_IDs <- levels(SeuratOBJ)
## extract Hb clusters
Hb_cluster_IDs <- cluster_IDs[grepl("MHb|LHb", all_cluster_IDs)]
Hb_cluster_IDs
# [1] "C.05.DD_LHb" "C.07.DD_MHb" "C.10.DD_MHb" "C.11.DD_MHb" "C.14.DD_MHb"
# [6] "C.16.DD_MHb" "C.18.DD_LHb" "C.23.DD_LHb" "C.24.DD_LHb" "C.30.DD_LHb"
# [11] "C.33.DD_LHb" "C.36.DD_MHb" "C.40.DD_LHb"
pairwise_combinations <- combn(Hb_cluster_IDs, 2, simplify = FALSE)
pairwise_combinations
pairwise_df <- do.call(rbind, pairwise_combinations)
colnames(pairwise_df) <- c("cluster_1", "cluster_2")
pairwise_df <- as.data.frame(pairwise_df)
pairwise_df


for (clust_p in pairwise_df)
    da_peaks <- FindMarkers(
        object = seurat_atac,
        ident.1 = clust_p["cluster1"],
        ident.2 = clust_p["cluster2"],
        test.use = 'LR',
        min.pct = 0.05
    )
}


