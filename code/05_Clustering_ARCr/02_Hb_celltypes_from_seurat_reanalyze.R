########################################################################
## Read Seurat clusters to search cell-types based on a custom marker gene list
## INPUT:
##      A Seurat Harmony corrected dataset
##      A csv file with Gene-marker list
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
library("data.table")
library("magrittr")
library("here")

here::here()

## read input arguments ( name of RDS Seurat file with wnn clustering to parse )
args = commandArgs(trailingOnly=TRUE)
Seurat_base_name <- args[2]

## input directories

# Check/create directories
inputRDS_Dir <- here("processed-data", "05_Clustering_ARCr", "01_clustering_std_method")
inputDir_cvs <- here("processed-data", "05_Clustering_ARCr", "01_clustering_std_method", "cvs_files")
processedDir <- here("processed-data", "05_Clustering_ARCr", "02_Hb_celltypes_from_seurat_reanalyze")
cvsDir <- here("processed-data", "05_Clustering_ARCr", "02_Hb_celltypes_from_seurat_reanalyze", "cvs_files_markers")

## Check directories
if (!dir.exists(processedDir)) {dir.create(processedDir)}
if (!dir.exists(cvsDir)) {dir.create(cvsDir)}

## Contains marker lists 
source(here("code", "04_DiffExpr_Clustering_seurat", "remote_DGE_marker_gene_lists.R"))       # Call functions to read paths


#############################           Initials        ################################

message("Reading files to annotate cell-types in WNN clusters")

## Some WNN clustering results of interest
# seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.leiden_lsi_r1.rds
# seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.leiden_lsi_r2.rds
# seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.louvain_lsi_r1.rds
# seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.louvain_lsi_r2.rds
# seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.louvainM_lsi_r1.rds
# seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.louvainM_lsi_r2.rds
# seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.SLM_lsi_r1.rds
# seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.SLM_lsi_r2.rds

# testing
# Seurat_base_name <- "seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.louvain_lsi_r1"

seurat_RDSname = paste0(Seurat_base_name, ".rds")

## Validate seurat exists
if (length(list.files(inputRDS_Dir, pattern = seurat_RDSname)==1)) {
  message("Processing ", Seurat_base_name)
} else {
  message("Input seurat object missed!")
  stop()
}
seurat_RDSname <- here(inputRDS_Dir, seurat_RDSname)

SeuratOBJ <- readRDS(seurat_RDSname)
## verification
length(Cells(x = SeuratOBJ))
nrow(unique(SeuratOBJ[["seurat_clusters"]]))

message("Seurat loaded! Starting cell-type identification")

## Select gene markers lists. We have 3.
markers.custom = list()

## Gene markers lists
markers.custom[["data_driven"]] <- get_Top50r_markers_genes_Hb()
markers.custom[["literature_base"]] <- get_erik_and_Hb_markers_genes()  

## sub-population list
#names(markers.custom$literature_base)
#names(markers.custom)

## Check unique marker genes
# x <- markers.custom[["literature_base"]]
# unlist(x)
# table(unname(unlist(x)))
# duplicated(unname(unlist(x)))
# summary(table(unname(unlist(x))))

prefix_name <- 'all_gm'                                    # prefix to save matched markers found in the clusters

## set the number of top DGE genes to pick up
n_slice <- 20  


#############################  Set the DGE list to parse  ################################

## Extract cluster data
unique(SeuratOBJ@meta.data$orig.ident)

md <- SeuratOBJ@meta.data %>% as.data.table

## Apply vertical format to unique cluster with number of UMIs, arranged by sample and cluster number
mdT <- md[, .N, by = c("orig.ident", "seurat_clusters")] %>%
    arrange(., orig.ident, seurat_clusters, .by_group = FALSE)
df_mdT <- as.data.frame(mdT)
#head(df_mdT)
#sum(df_mdT$N)

## Save cluster information
cvs_name <- here(cvsDir, paste0(Seurat_base_name, '_cluster_info.csv'))
## e.g: seurat.data_counts_Harmony_cluster_info.csv
write.csv(df_mdT, cvs_name)

## extract unique clusters in ascending order
clusters <- unique(df_mdT$seurat_clusters)
clusters <- as.integer(levels(clusters)[as.integer(clusters)])

message('Identifing cell types for ', length(clusters),' clusters from CR_reanalyze QCed dataset')

## Read All markers CVS file for all clusters
DGE_cvs_name <- here(inputDir_cvs, paste0(Seurat_base_name, "_markers.csv"))

seurat_clust <- as.data.frame(read.csv(DGE_cvs_name, header = TRUE))
head(seurat_clust, n=3)

message('Parsing ', length(markers.custom), ' gene-markers lists on ', length(clusters) ,' clusters in CR_reanalyze QCed dataset ')


####### Parse the 10/20 DGE genes from GEX cluster against the marker genes list provided ####### 

## Build df to save cell-types that match with the gene-marker-list
names(markers.custom) #[1] "literature_base" "data_driven" 
idx_lst <- 0 

for (markers.lst in markers.custom) {
  # for testing: markers.lst <- markers.custom$literature_base
  # for testing: markers.lst <- markers.custom$data_driven
  
  all_gene_match <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("Feature.ID", "Feature.Name", "Cluster.Adjusted.p.value", "cell-type", "cluster"))
  
  idx_lst <- idx_lst + 1
  
  message("\nSearching cell-types for ")
  names(markers.lst)
  
  ## Compose file name with cell-types identified, for every group of clusters defined above, for every marker-list reference provided 
  prefix_name <- paste0(names(markers.custom[idx_lst]), '_top', n_slice)
  
  # Parse every cluster and extract the top <n_slice> genes
  for (clust in clusters) {
    # Testing: clust <- 0
    message("Parsing cluster ", as.character(clust))
    top_DGE_clust <- seurat_clust |> 
      dplyr::filter(cluster == clust, p_val_adj < 0.05) |> slice_head(n = n_slice)
    
    # get a vector with all cell-types with their marker genes
    gm_lst <- as.vector(as.list(markers.lst))
    i_pos <- 0      # reset gene-marker list position
    
    for ( gm in gm_lst ) {
      #for testing: gm <- gm_lst[[2]]
      i_pos <- i_pos+1                        # control cell type position
      # Match top10genes with the marker genes for the cell-type x 
      gene_match <- top_DGE_clust |> filter_all(any_vars(. %in% gm))
      # add matched genes to a dataframe
      if ( nrow(gene_match) > 0 ) {
        names(gene_match)[names(gene_match) == clust ] <- "Cluster.Adjusted.p.value" # rename cols to rbind
        gene_match['cell-type']  <- names(gm_lst[i_pos])
        gene_match['cluster']  <- clust
        all_gene_match <- rbind(all_gene_match, gene_match)
      }
      
    }
  }
  
  habenula_markers_cvs_name <- here(cvsDir, 
                                    paste0(Seurat_base_name, '_cellTypes_', prefix_name, ".csv"))
  print(paste("Printing results in ", habenula_markers_cvs_name))
  write.csv(all_gene_match, habenula_markers_cvs_name, row.names=FALSE)
  rm("gene_match", "all_gene_match")
}

message(' Cell type identification done!')


# slurmjobs::job_loop(
#   loops = list(cellranger_pipe = c("CR_arc_reanalyze", "CR_arc_reanalyze_outliers", "CR_arc_reanalyze_outliers_ATAC")),
#   name = "01_Hb_celltypes_from_seurat_reanalyze_v4",
#   cores = 2,
#   create_shell = TRUE
# )

library("sessioninfo")
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

