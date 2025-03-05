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
library("tidyr")
library("here")

## input directories

# Check/create directories
# inputCVS_Dir <- here("code", "05_Clustering_ARCr")
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


## Read DEG to plot the top 5 genes

top_deg_file <- here(cvsDirOUT, "cvs_files_markers", suffix_clust_names)
basename(top_deg_file) # WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r1_cellTypes_integrated_top50.csv

df_cluster_names <- read.csv(top_deg_file)
df_cluster_names <- df_cluster_names |> drop_na(cell_type)
head(df_cluster_names)
# p_val avg_log2FC pct.1 pct.2 p_val_adj cluster       gene     cell_type
# 1     0       1.66 0.626  0.28         0       1     GABRG3 DD_Excit.Thal
# 2     0      -3.92 0.034  0.34         0       1    RASGRP1        DD_MHb
# 3     0       0.61 0.935  0.63         0       1     RBFOX1 DD_Excit.Thal
# 4     0       2.37 0.673  0.24         0       2 AC007614.1        DD_LHb
# 5     0       2.18 0.612  0.23         0       2       CDH4        DD_LHb
# 6     0       2.89 0.540  0.19         0       3       GNG8        DD_MHb

# lst_clusters <- unique(df_cluster_names$cluster)
# lst_clust <- list()
# 
# for (clust in lst_clusters) {
#   # testing:  clust=1
#   message("Cluster: ", clust)
#   lst_top5 <- df_cluster_names$gene[df_cluster_names$cluster==clust][1:5]
#   lst_top5 <- lst_top5[!is.na(lst_top5)]
#   print(lst_top5)
#   
#   ## Build the list of lists with top 5 genes by cluster
#   lst_name <- paste0("C_",as.character(clust))
#   if (!is_null(lst_top5)) { lst_clust <- append(lst_clust, list(names(lst_name)[1] <- lst_top5)) }
# }
# 
# lst_clust[[1]]
# # [1] "GABRG3"  "RASGRP1" "RBFOX1" 

## f_plt_violin_top5 <- function(seurat_name, l_clus){
  l_clus <- 1
  seurat_name <- wnn_file_names_lst[1]
  
  SeuratOBJ <- readRDS(seurat_name)
  DefaultAssay(SeuratOBJ) <- "RNA"
  Seurat_base_name <- str_extract(seurat_name, regex("C\\.\\w+")) 
  # Assign the top5 features
  features <- l_clus[[1]]
  
  plt1 <- VlnPlot(object = SeuratOBJ, layer = "data",
                  features = features,
                  pt.size = 0) +
    labs(x = paste0("**Clusters from WNN: ", Seurat_base_name)) &
    theme(text = element_text(size = 8), 
          axis.text.x= element_text(size = 7), axis.text.y= element_text(size = 7),
          plot.title=element_text(hjust=0.5)) 
  
  suffix <- str_extract(basename(seurat_name), "WNN+.+k30")
  tmp_name <- paste0(Seurat_base_name, "_C", l_clus, "_VlnPlot_", suffix ,".pdf")
  ggsave(plt1, filename = here(plotDir, tmp_name), height = 4, width = 17)
  
  return(plt1)
  
}

seurat_name <- basename(wnn_file_names_lst[[2]])
Seurat_base_name <- str_extract(seurat_name, regex("C\\.\\w+")) 
suffix <- str_extract(basename(seurat_name), "WNN+.+k30")
tmp_name <- paste0(Seurat_base_name, "_Top5markers_", suffix ,".pdf")

## Get top 5 
head(df_cluster_names)
hb_clust <- c(2,4,7,9,15,25)
top5 <- df_cluster_names |>
  filter(cluster %in% hb_clust) |>
  group_by(cluster) |>
  top_n(n = 5, wt = avg_log2FC)
head(top5)
unique(top5$cluster)

message(" Ploting gene-expression of the top 5 DEG: ", length(top5))

pdf(file=here(plotDir, tmp_name))

walk(seq_along(top5$gene), ~ {
  # message(paste0("Processing gene ",.x, " ", top5$gene[.x]))
  tryCatch({
    
    # message(paste0("Generating Vplot for gene ", top5$gene[.x]))
    
    plt1 <- VlnPlot(object = SeuratOBJ, layer = "data",
                    #features = "COL1A2",
                    features = top5$gene[.x],
                    pt.size = 0) +
      labs(x = paste0("**Clusters from WNN: ", Seurat_base_name)) &
      theme(text = element_text(size = 8), 
            axis.text.x= element_text(size = 7), axis.text.y= element_text(size = 7),
            plot.title=element_text(hjust=0.5))     
    print(plt1)
    message(paste0("Gene ", .x , " ", top5$gene[.x], " completed successfully"))
    
  }, error = function(e) {
    
    message(paste0("Error occurred while processing gene ", top5$gene[.x], ": ", e$message))
    
  })
  
})

dev.off()

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


