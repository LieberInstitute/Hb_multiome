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

message("Processing coverage plots for `POU4F1` and `GPR151` genes")
features <- c("POU4F1", "GPR151")

## Coverage plot with canonical Habenula genes 

f_name <- paste0(Seurat_base_name,"_peaks_GPR151.png")
features <- "GPR151"
plt1 <- CoveragePlot(
  object = SeuratOBJ,
  region = features,
  features = features,
  extend.upstream = 500,
  extend.downstream = 500,
  peaks = TRUE,
  links = TRUE
)  
plt1 <- plt1 +
  labs(title = paste0("Clusters from WNN: ", seurat_name)) +
  theme(text = element_text(size = 8),
        axis.text.x= element_text(size = 7), axis.text.y= element_text(size = 7),
        plot.title=element_text(hjust=0.5))
ggsave(plt1, filename = here(plotDir, f_name), height = 12, width = 6)


f_name <- paste0(Seurat_base_name,"_peaks_POU4F1.png")
features <- "POU4F1"
plt1 <- CoveragePlot(
  object = SeuratOBJ,
  region = features,
  features = features,
  extend.upstream = 500,
  extend.downstream = 500,
  peaks = TRUE,
  links = TRUE
)  
plt1 <- plt1 +
  labs(title = paste0("Clusters from WNN: ", seurat_name)) +
  theme(text = element_text(size = 8),
        axis.text.x= element_text(size = 7), axis.text.y= element_text(size = 7),
        plot.title=element_text(hjust=0.5))
ggsave(plt1, filename = here(plotDir, f_name), height = 12, width = 6)


f_name <- paste0(Seurat_base_name,"_peaks_TAC3.png")
features <- "TAC3"
plt1 <- CoveragePlot(
  object = SeuratOBJ,
  region = features,
  features = features,
  extend.upstream = 500,
  extend.downstream = 500,
  peaks = TRUE,
  links = TRUE
)  
plt1 <- plt1 +
  labs(title = paste0("Clusters from WNN: ", seurat_name)) +
  theme(text = element_text(size = 8),
        axis.text.x= element_text(size = 7), axis.text.y= element_text(size = 7),
        plot.title=element_text(hjust=0.5))
ggsave(plt1, filename = here(plotDir, f_name), height = 12, width = 6)



## ========================================================================== ##

## Read DEG to prepare coverage plots of the top 5 genes highly expressed

# All DEG

DEG_file_name <- "WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_cellTypes_integrated_top50.csv"
DEG_file_name <- here(inputCVS_Dir, DEG_file_name)
df_cluster_names <- read.csv(DEG_file_name)
df_cluster_names <- df_cluster_names |> drop_na(cell_type)
head(df_cluster_names)
# p_val avg_log2FC pct.1 pct.2 p_val_adj cluster     gene            cell_type
# 1     0   4.336764 0.938 0.100         0       1 OTX2-AS1        DD_Inhib.Thal
# 2     0   3.849508 0.900 0.083         0       1      KIT        DD_Inhib.Thal
# 3     0   3.682782 0.927 0.122         0       1    MEIS2      LB_Thalamus/MDm

## Identified and subset clusters annotated for `habenula`.

message("Cluster-IDs from `WNN`")

SeuOBJ_clusters <- Idents(SeuratOBJ)
hb_clusters <- unlist(levels(SeuOBJ_clusters))
## Get top 5. Filter habenula clusters only
no_hb_clust = list()
for (idx in seq_along(hb_clusters)) {
  if (nchar(hb_clusters[idx]) <= 4) {
    no_hb_clust <- append(no_hb_clust, hb_clusters[idx])
  }
}
no_hb_clust <- c(unlist(no_hb_clust))
hb_clusters <- hb_clusters[!hb_clusters %in% c(no_hb_clust)]
as.vector(hb_clusters)
# [1] "C.05 DD_LHb" "C.07 DD_MHb" "C.10 DD_MHb" "C.11 DD_MHb" "C.14 DD_MHb"
# [6] "C.16 DD_MHb" "C.18 DD_LHb" "C.23 DD_LHb" "C.24 DD_LHb" "C.30 DD_LHb"
# [11] "C.33 DD_LHb" "C.36 DD_MHb" "C.40 DD_LHb"
# length(hb_clusters)

##  Use length of clustersto extract clusters IDs
hb_clusters <- as.integer(substr(hb_clusters, 3, 4))
# [1]  5  7 10 11 14 16 18 23 24 30 33 36 40

## filter the top 5
unique(df_cluster_names$cluster)
top5 <- df_cluster_names |>
  filter(cluster %in% hb_clusters) |>
  group_by(cluster) |>
  top_n(n = 5, wt = avg_log2FC)
head(top5)
#     p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene    cell_type
# <dbl>      <dbl> <dbl> <dbl>     <dbl>   <int> <chr>   <chr>    
# 1     0       3.04 0.768 0.145         0       5 RFTN1   DD_LHb   
# 2     0       3.03 0.798 0.194         0       5 CBLN2   DD_LHb   
# 3     0       3.37 0.73  0.129         0       5 GALR1   DD_LHb   
# 4     0       3.15 0.669 0.124         0       5 HTR4    DD_LHb   
# 5     0       3.38 0.956 0.457         0       5 COL25A1 DD_LHb  


## Prepare and save coverage plot

for (clus in unique(top5$cluster)) {
  # testing: clus = 5
  tmp_name <- paste0(
    Seurat_base_name,
    "_PEAKS_hb-cluster-",
    clus,
    ".pdf"
  )
  
  message("Processing habenula cluster: ", clus, "; Saved as: ", tmp_name)
  
  top5_cluster <- top5 |>
    filter(cluster == clus)
  
  pdf(file = here(plotDir, tmp_name))
  
  walk(seq_along(top5_cluster$gene), ~ {

    tryCatch({
      
      message(paste0("Processing gene ", top5_cluster$gene[.x]))
      
      features <- top5_cluster$gene[.x]
      plt1 <- CoveragePlot(
        object = SeuratOBJ,
        region = features,
        features = features,
        extend.upstream = 500,
        extend.downstream = 500,
        peaks = TRUE,
        links = TRUE
      )
      plt1 <- plt1 +
        labs(title = paste0("Clusters from WNN: ", seurat_name)) +
        theme(text = element_text(size = 8),
              axis.text.x= element_text(size = 7), axis.text.y= element_text(size = 7),
              plot.title=element_text(hjust=0.5))
      print(plt1)
      
    }, error = function(e) {
      
      message(paste0("Error occurred while processing gene ", top5_cluster$gene[.x], ": ", e$message))
      
    })
    
  })
  
  dev.off()
  
}
# Error occurred while processing gene AC109466.1: Gene not found
# Error occurred while processing gene LINC02143: Gene not found

message("Coverage plots completed")



## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


