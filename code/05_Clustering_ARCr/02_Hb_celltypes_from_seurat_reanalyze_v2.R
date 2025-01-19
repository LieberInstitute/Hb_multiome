########################################################################
## Read Seurat clusters to search cell-types based on a custom marker gene list
## INPUT:
##      A Seurat Harmony corrected dataset
##      A csv file with Gene-marker list. INTEGRATE ALL GENE MARKERS IN ONE LIST
##      A csv file with DGE from a Seurat Harmony (not pseudo-bulked)
## 
## OUPUT:
##      1) A csv files with clusters cell-type identification
##
## NOTE. Identify cell-types on WNN clusters from CellRangerARC-reanalyze filtered datasets
##
## Authors. CSC 
## Date. Dec 09th, 2024
########################################################################

## load libraries
library("Seurat")
library("Signac")
library("tidyverse")
library("dplyr")
library("stringr")
library("data.table")
library("magrittr")
library("here")

here::here()

## read input arguments ( name of RDS Seurat file with wnn clustering to parse )
Seurat_base_name <- commandArgs(trailingOnly = TRUE)

# testing:
# Seurat_base_name <- "seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.leiden_lsi_r1"


message(" Reading: ", Seurat_base_name)

## input directories

# Check/create directories
inputRDS_Dir <- here("processed-data", "05_Clustering_ARCr", "01_clustering_std_method")
inputDir_cvs <- here(inputRDS_Dir, "cvs_files")
processedDir <- here("processed-data", "05_Clustering_ARCr", "02_Hb_celltypes_from_seurat_reanalyze_v2")
cvsDir <- here(processedDir, "cvs_files_markers")

## Check directories
if (!dir.exists(processedDir)) {dir.create(processedDir)}
if (!dir.exists(cvsDir)) {dir.create(cvsDir)}

## Contains marker lists 
source(here("code", "04_DiffExpr_Clustering_seurat", "remote_DGE_marker_gene_lists.R"))       # Call functions to read paths


#############################           Initials        ################################

message("Reading files to annotate cell-types in WNN clusters")

## Some WNN clustering results of interest. Testing:
# Seurat_base_name <- "seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.leiden_lsi_r1"
# Seurat_base_name <- "seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.louvain_lsi_r1"
# Seurat_base_name <- "seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.louvainM_lsi_r1"
# Seurat_base_name <- "seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.SLM_lsi_r1"

seurat_RDSname = paste0(Seurat_base_name, ".rds")

# Validate seurat exists
if (length(list.files(inputRDS_Dir, pattern = seurat_RDSname)==1)) {
  message("Processing: ", Seurat_base_name)
} else {
  message("Input seurat object missed!")
  stop()
}

seurat_RDSname <- here(inputRDS_Dir, seurat_RDSname)

SeuratOBJ <- readRDS(seurat_RDSname)
## verification
length(Cells(x = SeuratOBJ))
nrow(unique(SeuratOBJ[["seurat_clusters"]]))

message("Seurat loaded! \nStarting cell-type identification ...")

## Gene markers lists. New function to join LB and DD gene markers lists

markers.custom <- get_multiple_markers_genes_lst()
tmp <- names(markers.custom)
tmp <- paste(tmp, collapse=', ')
message("Processing ", length(markers.custom), " categories of gene-markers list \n *****(", tmp, ")*****")

## sub gene markers list
# print(names(markers.custom))

## Check unique and number of duplicate marker genes
# x <- markers.custom
# unlist(x)
# table(unname(unlist(x)))
# x1 <- duplicated(unname(unlist(x)))
# sum(x1, na.rm=TRUE)
# summary(table(unname(unlist(x))))

prefix_name <- 'all_gm'                                    # prefix to save matched markers found in the clusters

## set the number of top DGE genes to pick up

n_slice <- 50  


#############################  Set the DGE list to parse  ################################

## Extract cluster data
message("Samples to process: ", paste(unique(SeuratOBJ@meta.data$orig.ident), collapse = ", "))

md <- SeuratOBJ@meta.data %>% as.data.table

## Apply vertical format to unique cluster with number of UMIs, arranged by sample and cluster number
mdT <- md[, .N, by = c("orig.ident", "seurat_clusters")] %>%
  arrange(., orig.ident, seurat_clusters, .by_group = FALSE)
df_mdT <- as.data.frame(mdT)
head(df_mdT)
sum(df_mdT$N)

## Save cluster information

# Extract base name to easily identify files  

tmp <- str_extract(Seurat_base_name, "ARC+.+")
# Ex. ARCr_QCed_WNN_k30_C.leiden_lsi_r1

cvs_name <- here(processedDir, paste0(tmp, '_cluster_info.csv'))
write.csv(df_mdT, cvs_name)

## extract unique clusters in ascending order
clusters <- unique(df_mdT$seurat_clusters)
clusters <- as.integer(levels(clusters)[as.integer(clusters)])

## Read All markers CVS file for all clusters
DGE_cvs_name <- here(inputDir_cvs, paste0(Seurat_base_name, "_markers.csv"))

message('/nSaved cluster info. \nIdentifying cell types for ', length(clusters),' clusters using `', basename(DGE_cvs_name), "`")
# Identifying cell types for 33 clusters using `seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.leiden_lsi_r1_markers.csv`

seurat_clust_DEG <- read.csv(DGE_cvs_name, header = TRUE) #seurat_clust_DEG <- as.data.frame(read.csv(DGE_cvs_name, header = TRUE))
seurat_clust_DEG <- seurat_clust_DEG[,-1]
head(seurat_clust_DEG, n=3)
# p_val avg_log2FC pct.1 pct.2 p_val_adj cluster       gene
# 1     0  -3.229510 0.222 0.693         0       1      CPNE4
# 2     0  -3.692433 0.184 0.640         0       1       VAV3
# 3     0  -4.193924 0.029 0.460         0       1 AC119673.2

## Save Top 50 DEG by sample before cell-type annotation for reference

top_DGE_clust <- seurat_clust_DEG |> 
  group_by(cluster) |>
  filter(p_val_adj < 0.05) |> 
  slice_head(n = n_slice) |> 
  arrange(cluster, p_val_adj)
tail(top_DGE_clust)

all_markers_cvs_name <- here(processedDir, paste0(tmp, "_all_DEG_Top50_WNN.csv"))
write.csv(top_DGE_clust, all_markers_cvs_name)

message("\nSaved CVS file with Top50 DEG from: ", tmp)


####### Parse the 10/20 DGE genes from GEX cluster against the marker genes list provided ####### 

## Build df to save cell-types that match with the gene-marker-list
# names(markers.custom) #[1] "literature_base" and "data_driven" in the same list
# [1] "DD_Astrocyte"                 "DD_Endo"                     
# [3] "DD_Excit.Thal"                "DD_Inhib.Thal"               
# [5] "DD_LHb"                       "DD_MHb"                      
# [7] "DD_Microglia"                 "DD_Oligo"                    
# [9] "DD_OPC"                       "LB_neuron"                   
# [11] "LB_excitatory_neuron"         "LB_inhibitory_neuron"        
# [13] "LB_Hb neuron specific"        "LB_MHB neuron specific"      
# [15] "LB_LHB neuron specific"       "LB_oligodendrocyte"          
# [17] "LB_oligodendrocyte_precursor" "LB_microglia"                
# [19] "LB_astrocyte"                 "LB_Endo/CP"                  
# [21] "LB_Thalamus broad"            "LB_Thalamus/MDm/Endo"        
# [23] "LB_Thalamus/MDm/Endo-"        "LB_Thalamus/MDm"  

markers.lst <- markers.custom
# get a vector with all cell-types with their marker genes
gm_lst <- as.vector(as.list(markers.lst))

## tbl to save top 50 genes by cluster
all_gene_match <- setNames(data.frame(matrix(ncol = 5, nrow = 0)),
                           c("Feature.ID", "Feature.Name", "p_val_adj", "cell-type", "cluster")) #Cluster.Adjusted.p.value

message("\nSearching cell-types for all gene markers lists")

## Compose file name with cell-types identified, for every group of clusters defined above, for every marker-list reference provided
prefix_name <- paste0('DD_LB_top', n_slice)

## Annotate cell types based on the reference of gene markers DD+LB

for (clust in clusters) {
  # # Testing:
  # # clust <- 0
  # message("Parsing cluster ", as.character(clust))
  # top_DGE_clust <- seurat_clust_DEG |>
  #   dplyr::filter(cluster == clust, p_val_adj < 0.05) |> slice_head(n = n_slice)
  # top_DGE_clust
  # 
  # i_pos <- 0      # reset gene-marker list position
  # 
  # for ( gm in gm_lst ) {
  #   #for testing: gm <- gm_lst[[3]]
  #   i_pos <- i_pos+1                        # control cell type position
  #   # Match top N genes with the marker genes for the cell-type x
  #   gene_match <- top_DGE_clust |> filter_all(any_vars(. %in% gm))
  #   # add matched genes to a dataframe
  #   if ( nrow(gene_match) > 0 ) {
  #     names(gene_match)[names(gene_match) == clust ] <- "Cluster.Adjusted.p.value" # rename cols to rbind
  #     gene_match['cell-type']  <- names(gm_lst[i_pos])
  #     gene_match['cluster']  <- clust
  #     all_gene_match <- rbind(all_gene_match, gene_match)
  #   }
  # }
}

# habenula_markers_cvs_name <- here(cvsDir, 
#                                   paste0(Seurat_base_name, '_cellTypes_', prefix_name, ".csv"))
# print(paste("Printing results in ", habenula_markers_cvs_name))
# write.csv(all_gene_match, habenula_markers_cvs_name, row.names=FALSE)

message(' Cell type identification done!')


library("sessioninfo")
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

