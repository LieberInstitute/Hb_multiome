########################################################################
## Merge WNN Clusters
## INPUT:
##      (1) Seurat object with WNN renamed clusters
## OUPUT:
##      (1) Seurat with clusters annotated and MERGED meta-data added
##      (2) Visualization: DimPlot
##      (3) Visualization: Heatmap
##          Ref.https://github.com/LieberInstitute/LFF_spatial_ERC/blob/7622879b82b6ee5cd8f21a2c5ac18006af7a7aec/code/04_snRNA-seq/32_sn_subcluster_hierarchical_cluster.R#L22-L27)
## Authors. CSC 
## Date. Jun, 2025
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
########################################################################

library("Seurat")
library("purrr")
library("tibble")
library("ggplot2")
library("dplyr")
library("stringr")
library("here")

## input directories

if (packageVersion("Seurat") != "5.0.1") stop("Incompatible Seurat version")

inputRDS_Dir <- here("processed-data", "05_Clustering_ARCr", "01_clustering_std_method") 
outputRDS_Dir <- here("processed-data", "05_Clustering_ARCr", "14_merge_renamed_clusters")
plotDir <- here("plots", "05_Clustering_ARCr", "14_merge_renamed_clusters")

## Check directories
if (!dir.exists(plotDir)) {dir.create(plotDir)}
if (!dir.exists(outputRDS_Dir)) {dir.create(outputRDS_Dir)}

Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium.rds"

message("Processing ", Seurat_base_name)


##### (1) Load Seurat with WNN idents given by default 

seurat_RDSname <- here(inputRDS_Dir, paste0(Seurat_base_name, ".rds"))

SeuratOBJ <- readRDS(seurat_RDSname)

total_cells <- length(Cells(x = SeuratOBJ))

message("Renaming ", nrow(unique(SeuratOBJ[["seurat_clusters"]])), " clusters for ", Seurat_base_name)

message("Clusters before merge")
levels(SeuratOBJ)
# [1] "C.04.LHb.4"         "C.05.LHb.2.7"       "C.07.MHb.2"        
# [4] "C.08.LHb.4"         "C.09.LHb.4"         "C.10.MHb.1"        
# [7] "C.11.MHb.1.2"       "C.13.LHb.4"         "C.14.MHb.1"        
# [10] "C.16.LHb.6"         "C.18.LHb.1.3.4"     "C.23.LHb.1"        
# [13] "C.24.DD_LHb"        "C.30.DD_LHb"        "C.33.LHb.1.3"      
# [16] "C.36.MHb.3"         "C.40.LHb.4"         "C.01.Inhib.Thal"   
# [19] "C.02.Oligo"         "C.03.Excit.Thal"    "C.06.DD_Excit.Thal"
# [22] "C.12.Excit.Thal"    "C.15.Excit.Thal"    "C.17.Excit.Thal"   
# [25] "C.19.Inhib.Thal"    "C.20.Astrocyte"     "C.21.Astrocyte"    
# [28] "C.22.Oligo"         "C.25.Excit.Thal"    "C.26.OPC"          
# [31] "C.27.Microglia"     "C.28.Inhib.Thal"    "C.29.Endo"         
# [34] "C.31.DD_Excit.Thal" "C.32.Excit.Thal"    "C.34.DD_Oligo"     
# [37] "C.35.Excit.Thal"    "C.37.Thal"          "C.38.Inhib.Thal"   
# [40] "C.39.Inhib.Thal"    "C.41.Microglia"

