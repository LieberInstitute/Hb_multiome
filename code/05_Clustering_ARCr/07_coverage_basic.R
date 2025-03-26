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
# library("patchwork")
# library("purrr")
library("tidyverse")
# library("stringr")
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
  "05_Clustering_ARCr",
  "07_coverage_basic"
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

## Coverage plot with cannonical habenula genes 

# f_name <- paste0(
#   Seurat_base_name,
#   "_peaks_POU4F1_GPR151.png"
# )
# plt1 <- CoveragePlot(
#   object = SeuratOBJ,
#   region = features,
#   extend.upstream = 500,
#   extend.downstream = 00,
#   peaks = TRUE,
#   links = TRUE
# )  

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






## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


