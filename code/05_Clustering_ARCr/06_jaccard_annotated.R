########################################################################
## Compares GEX clustering results against pre-selected WNN clustering results
## Based on: https://bioconductor.org/books/3.14/OSCA.advanced/clustering-redux.html
## INPUT:
##      (1) First Seurat with WNN to compare
##      (2) Second Seurat with WNN to compare
## OUPUT:
##      1) Jaccard index heat-map  and approximate silhouette (plot and cvs files) 
##
## Authors. CSC 
## Date. Jan 29, 2025
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
## Note: module conda_R/4.3.x is required to load Seurat Objects processed with previous R version
########################################################################

library("Seurat")
library("Signac")
library("bluster")
library("pheatmap")
library("dplyr")
library("ggplot2")
library("stringr")
library("grid")
library("here")
library("viridisLite")

## input directories

inputRDS_Dir <- here("processed-data", "05_Clustering_ARCr", "05_rename_idents")
outputCVS_Dir <- here("processed-data", "05_Clustering_ARCr", "06_jaccard_annotated")
plotDir <- here("plots", "05_Clustering_ARCr", "06_jaccard_annotated")

## Check directories

if (!dir.exists(plotDir)) {dir.create(plotDir)}
if (!dir.exists(outputCVS_Dir)) {dir.create(outputCVS_Dir)}


## Load input with RDS wnn to compare

## read input arguments ( name of RDS Seurat file with wnn clustering to parse )

Seurat_base_name <- commandArgs(trailingOnly = TRUE)
## Some WNN clustering results of interest
## For testing:
# Old names: 
# Seurat_base_name <- "seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.louvain_lsi_r1,seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.louvainM_lsi_r1"
# New names:
# Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k40_C.leiden_lsi_r2,seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2"



message("Processing: \n", str_split(Seurat_base_name, ",")[[1]][1], "\n", str_split(Seurat_base_name, ",")[[1]][2])

## Prepare RDS seurat names and file names 
Seurat_base_name_1 <- trimws(strsplit(Seurat_base_name, ",")[[1]][1])
Seurat_base_name_2 <- trimws(strsplit(Seurat_base_name, ",")[[1]][2])
seurat_RDSname_1 <- here(inputRDS_Dir, paste0(Seurat_base_name_1, ".rds"))
seurat_RDSname_2 <- here(inputRDS_Dir, paste0(Seurat_base_name_2, ".rds"))
str_extract(Seurat_base_name_1, regex("C\\.\\w+"))
# Names used when we were only comparing methods, now we have to extend the variable name to add the knn and the resolution
# Seurat_base_name_1 <- str_extract(Seurat_base_name_1, regex("C\\.\\w+")) #C.louvain_lsi_r1
# Seurat_base_name_2 <- str_extract(Seurat_base_name_2, regex("C\\.\\w+")) #C.louvainM_lsi_r1
Seurat_base_name_1 <- str_extract(Seurat_base_name_1, regex("k[3|4]0\\_C\\.\\w+")) #k30_C.leiden_lsi_r1
Seurat_base_name_2 <- str_extract(Seurat_base_name_2, regex("k[3|4]0\\_C\\.\\w+")) #k40_C.leiden_lsi_r1
file_name_all  <- paste0(Seurat_base_name_1, "-", Seurat_base_name_2)

## Load the Seurats with WNN clusters

message("Reading Seurat(s) to evalute WNN-Clusters and compute Jaccard Index:\nWNN.1: ",
        Seurat_base_name_1, "\nWNN.2: ", Seurat_base_name_2)


## Prepare data based on Seurat(s) with WNN to compare clustering 

f_prepare_data_to_plot <- function(seurat_name){
  
  ## Load Seurat 
  SeuratOBJ <- readRDS(seurat_name)
  # Reductions(SeuratOBJ) # "integrated.harmony" "umap.lovain" "wnn.umap"
  message("WNN clustering loaded!\nCells: ", length(Cells(x = SeuratOBJ)))
  message("\nWNN `", seurat_name, "` containing ", nrow(unique(SeuratOBJ[["seurat_clusters"]])), " clusters")
  # table(SeuratOBJ[["seurat_clusters"]])
  # 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
  # 5516 5240 3041 3023 2378 2330 2295 2261 2123 2077 2041 1971 1809 1758 1743 1700 ...

  ## extract cells and clusters
  message("Cells-IDs from `integrated.harmony reduction`")
  SeuOBJ_cellEmbeddings <- Embeddings(SeuratOBJ, reduction = "integrated.harmony")
  message("Cluster-IDs from `WNN`")
  SeuOBJ_clusters <- Idents(SeuratOBJ)
  
  return(list(cellsEmb=SeuOBJ_cellEmbeddings, clust=SeuOBJ_clusters))
  
}

## Prepare function to:
## (1) Plot approximate silhouette for evaluating cluster separation
## (2) Identified and save closest neighboring cluster for each cell in each cluster 

f_plot_approxSilhouette <- function(cellsID, clustID, fn){

  sil.approx <- approxSilhouette(cellsID, clustID)
  #sil.approx
  # DataFrame with 84177 rows and 3 columns
  # cluster    other      width
  # <factor> <factor>  <numeric>
  # 10C_AAACAGCCAATCATGT-1       5        13  0.2307765
  # 10C_AAACAGCCACTTCACT-1       16       2   0.4228183
  # 10C_AAACAGCCAGGACCTT-1       0        6   0.1671839
  
  sil.data <- as.data.frame(sil.approx)

  sil.data$closest <- factor(ifelse(sil.data$width > 0, clustID, sil.data$other))
  
  sil.data$cluster <- clustID

    ## identified the closest neighboring cluster for each cell in each cluster
  tbl_aprox_sil <- table(Cluster = clustID, sil.data$closest)
  cvs_file <- paste0(fn, "_Silhouette_WNN.cvs")
  cvs_file <- here(outputCVS_Dir, cvs_file)
  write.csv(tbl_aprox_sil, cvs_file)
  message("approximate silhouette cvs saved!")
  
  plt_title <- paste0("Silhouette.Approx for WNN ", fn)
  plt1 <- ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) +
    ggbeeswarm::geom_quasirandom(method="smiley", alpha=.4) + labs(title = plt_title) + labs(x='') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(legend.position = 'none')
  ggsave(plt1, filename = here(plotDir, paste0(fn,"_Silhouette_WNN.png")), height = 6, width = 10)
  
  message("Approximate-silhouette cvs files and plot saved!")
  
  return(plt1)
  
}


## extract cells and clusters from WNN-1 and WNN-2

if ( !length(list.files(inputRDS_Dir, basename(seurat_RDSname_1))==1) ) { message("Seurat 1 object missed!");  stop() }

# testing
# seurat_RDSname_1 = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/05_Clustering_ARCr/05_rename_idents/seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.louvain_lsi_r1.rds"
get_cell_info <- f_prepare_data_to_plot(seurat_RDSname_1) # first argument
names(get_cell_info)

message("cell info from ", Seurat_base_name_1, " processed!")

if ( !length(list.files(inputRDS_Dir, basename(seurat_RDSname_2))==1) ) { message("Seurat 1 object missed!");  stop() }

get_cell_info2 <- f_prepare_data_to_plot(seurat_RDSname_2) # second argument
names(get_cell_info2)

message("cell info from ", Seurat_base_name_2, " processed!")


## Process Approximate-silhouette, save csv files and plots 

message("Starting approximate-silhouette for evaluating cluster separation ...")

f_plot_approxSilhouette(get_cell_info$cellsEmb,
                        trimws(gsub("\\(\\d*\\.\\d*\\%\\)", "", get_cell_info$clust)), # remove % symbol and values to clean labels on plot
                        Seurat_base_name_1)
f_plot_approxSilhouette(get_cell_info2$cellsEmb,
                        trimws(gsub("\\(\\d*\\.\\d*\\%\\)", "", get_cell_info2$clust)),
                        Seurat_base_name_2)


## Comparing the clustering data sets 

message("Starting Jaccard Index Processing ...")

clust.wnn1 <- get_cell_info$clust # 1: louvain
levels(clust.wnn1)
clust.wnn2 <- get_cell_info2$clust # 2: louvainM
levels(clust.wnn2)


## Running Jaccard

# compute Jaccard
jacc.mat <- linkClustersMatrix(clust.wnn1, clust.wnn2)

# rename clusters
# gsub("C\\.|\\_lsi_r1","", Seurat_base_name_1)
# gsub("C\\.|\\_lsi_r1","", Seurat_base_name_2)
rownames(jacc.mat) <- trimws(gsub("\\(\\d*\\.\\d*\\%\\)", "", rownames(jacc.mat))) #rownames(jacc.mat) # remove % values to clean plot
colnames(jacc.mat) <- trimws(gsub("\\(\\d*\\.\\d*\\%\\)", "", colnames(jacc.mat))) #colnames(jacc.mat)
substr(rownames(jacc.mat), 1, 12)


## Save Jacquard plot

tmp_png <- here(plotDir, paste0(file_name_all, "_Jaccard_WNN.png"))
tmp_title <- unlist(strsplit(file_name_all,"-",fixed=T))
plt1<-pheatmap(jacc.mat, color=viridis::viridis(100), 
                 cluster_cols=FALSE, 
                 cluster_rows=FALSE, # show hierarchical clust
                 angle_col = 90,
                 main = paste0("Jaccard WNN\n (x-axis)", tmp_title[2], " (y-axis) ", tmp_title[1]),
                 fontsize = 10,
                 legend = TRUE)
                 #display_numbers=T,
                 #number_format="%.1f",
                 #fontsize_number=7) 

png(tmp_png,width=12,height=10,units="in",res=1200)
plt1
dev.off()

message("Process completed!")


## Reproducibility information

library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


