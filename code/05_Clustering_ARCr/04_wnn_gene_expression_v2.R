########################################################################
## Produce Violin Plots (s) for GEX on selected WNN clustering results
## 
## Input:Seurat object with WNN Leiden res=1 and knn=30. Both assay Harmonized
## Output: Violoin plots with cannonical Hb gene-markers
##
## Authors. CSC 
## Date. Jan 24, 2024
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
########################################################################

library("Seurat")
library("Signac")
library("ggplot2")
library("patchwork")
library("purrr")
library("stringr")
library("here")

## input directories

# Check/create directories
inputCVS_Dir <- here("code", "05_Clustering_ARCr")
inputRDS_Dir <- here("processed-data", "05_Clustering_ARCr", "01_clustering_std_method")
plotDir <- here("plots", "05_Clustering_ARCr", "04_wnn_gene_expression_v2")

## Check directories
if (!dir.exists(plotDir)) {dir.create(plotDir)}
if (!dir.exists(outputCVS_Dir)) {dir.create(outputCVS_Dir)}

## Read all the clustering results (RDS wnn name) selected to build the Violin Plots 
# tmp_dir <- here(inputCVS_Dir, "input_wnn_rds_names.txt")
# wnn_file_names = readLines(tmp_dir)[1:4] # only resolution r1, excluded r2
# wnn_file_names_lst <- here(inputRDS_Dir, paste0(wnn_file_names, ".rds"))
wnn_file_names_lst <- here(inputRDS_Dir,"seurat.norm_counts_CRr_WNN_rnaHarm_atacLSI_k30_C.leiden_lsi_r1.5.rds")
wnn_file_names_lst <- c(wnn_file_names_lst, here(inputRDS_Dir,"seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r1.5.rds"))

## function to build Vplots for the features selected

f_plt_violin <- function(seurat_name){
  
  SeuratOBJ <- readRDS(seurat_name)
  DefaultAssay(SeuratOBJ) <- "RNA"
  Seurat_base_name <- str_extract(seurat_name, regex("C\\.\\w+")) 
  # features <- c("POU4F1", "GPR151", "TAC3")
  features <- c("POU4F1", "GPR151")
  
  plt1 <- VlnPlot(object = SeuratOBJ, layer = "data",
                  features = features,
                  pt.size = 0) +
    labs(x = paste0("**Clusters from WNN: ", Seurat_base_name)) &
    theme(text = element_text(size = 8), 
          axis.text.x= element_text(size = 7), axis.text.y= element_text(size = 7),
          plot.title=element_text(hjust=0.5)) 
  
  suffix <- str_extract(basename(seurat_name), "WNN+.+k30")
  tmp_name <- paste0(Seurat_base_name, "_POU4F1_GPR151_VlnPlot_", suffix ,".pdf")
  ggsave(plt1, filename = here(plotDir, tmp_name), height = 4, width = 17)
  
  return(plt1)
  
}


## Load the Seurat with WNN clusters

message("Reading WNN to evalute gene expression of `POU4F1` and `GPR151`")

# test: Vplot_lst <- f_plt_violin(wnn_file_names_lst[1])

# we have 4 clustering results selected 
#tmp_file_name <- str_extract(basename(wnn_file_names_lst[[1]]), "WNN+.+k30")
Vplot_lst <- map(wnn_file_names_lst, ~ f_plt_violin(.x))
length(Vplot_lst)
# Vplot_lst[[1]]

tmp_name <- paste0(Seurat_base_name, "_ALL_POU4F1_GPR151_VlnPlot_HarmonyRNA-ATAC.pdf")
pdf(file = here(plotDir, tmp_name))

plt_Vplots_cols <- Reduce("/", Vplot_lst)
print(plt_Vplots_cols)

dev.off()


message('\nPlots saved `', plotDir, '`')


# library("slurmjobs")
# slurmjobs::job_single(
#   name = "04_wnn_gene_expression_v2", memory = "60G", cores = 2, create_shell = TRUE,
#   task_num = 8
# )


## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


