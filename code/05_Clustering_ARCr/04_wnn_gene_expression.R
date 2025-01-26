########################################################################
## Plot GEX on selected WNN clustering results
##
## Authors. CSC 
## Date. Jan 24, 2024
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
########################################################################

library("Seurat")
library("ggplot2")
library("patchwork")
library("purrr")
library("stringr")
library("here")

## input directories

# Check/create directories
inputCVS_Dir <- here("code", "05_Clustering_ARCr")
inputRDS_Dir <- here("processed-data", "05_Clustering_ARCr", "01_clustering_std_method")
plotDir <- here("plots", "05_Clustering_ARCr", "04_wnn_gene_expression.R")

## Check directories
if (!dir.exists(plotDir)) {dir.create(plotDir)}
if (!dir.exists(outputCVS_Dir)) {dir.create(outputCVS_Dir)}


## Load input with RDS wnn to compare

## read input arguments ( name of RDS Seurat file with wnn clustering to parse )
# Seurat_base_name <- commandArgs(trailingOnly = TRUE)
# ## Some WNN clustering results of interest
# ## For testing:
Seurat_base_name <- "seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.louvain_lsi_r1"

tmp_dir <- here(inputCVS_Dir, "input_wnn_rds_names.txt")
wnn_file_names = readLines(tmp_dir)[1:4] # only resolution r1, excluded r2

# any("r1" %in% wnn_file_names)

wnn_file_names_lst <- here(inputRDS_Dir, paste0(wnn_file_names, ".rds"))

f_plt_violin <- function(seurat_name){
  
  SeuratOBJ <- readRDS(seurat_name)
  DefaultAssay(SeuratOBJ) <- "RNA"
  Seurat_base_name <- str_extract(seurat_name, regex("C\\.\\w+")) 
  features <- c("POU4F1", "GPR151")
  
  plt1 <- VlnPlot(object = SeuratOBJ, layer = "data",
                  features = features,
                  pt.size = 0) +
    labs(x = paste0("WNN: ", Seurat_base_name)) &
    theme(text = element_text(size = 8), 
          axis.text.x= element_text(size = 7), axis.text.y= element_text(size = 7),
          plot.title=element_text(hjust=0.5)) 
    
  #plt1 <- plt1 + plot_annotation(paste0("WNN: ", Seurat_base_name), 
  #                               caption = 'Cell Ranger ARC reanalize',theme=theme(plot.title=element_text(hjust=0.5)))
  
  tmp_name <- paste0("Vplot_POU4F1_GPR151_", Seurat_base_name,".pdf")
  ggsave(plt1, filename = here(plotDir, tmp_name), height = 4, width = 16)
  
  return(plt1)
  
}


## Load the Seurat with WNN clusters

message("Reading WNN to evalute gene expression of `POU4F1` and `GPR151` on: ", str_extract(Seurat_base_name, regex("C\\.\\w+")))

# Vplot_lst <- f_plt_violin(wnn_file_names_lst[1])

Vplot_lst <- map(as.vector(wnn_file_names_lst), ~ f_plt_violin(.x))
length(Vplot_lst)

pdf_name <-  here(plotDir, "Vplot_WNN.POU4F1_GPR151.pdf")
pdf(file = pdf_name)
par(mfrow=c(2,2))
Vplot_lst[[1]] / Vplot_lst[[2]] / Vplot_lst[[3]] / Vplot_lst[[4]]
dev.off()

# walk(seq_along(Vplot_lst), ~ {
# 
#   tryCatch({
#     #print(.x)
#     Vplot_lst[[.x]]
#   }, error = function(e) {
#     message("Error occurred while plotting ", .x)
#   })
# 
# })
# dev.off()


    
    
message('\nPlots saved `', plotDir, '`')

# pdf(file=file_name_pdf)
# par(mfrow = c(2,1)) # 4 rows x 1 column 
# 
# for (i in seq_len(4)) {
#   print(Vplot_lst[i])
# }
# 
# dev.off()


## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


