########################################################################
## Read Seurat clusters to search cell-types based on a custom marker gene list
## INPUT:
##      A Seurat batch corrected data with `min_cells` by cluster filtered
##      A csv file with Gene-marker list.
##      A csv file with DGE from a Seurat Harmony/CCA data (not pseudo-bulked)
## 
## OUPUT:
##      1) A csv files with clusters cell-type identification
##
## NOTE. This pipeline only identify cell-types for cellRangerARC-reanalyze datasets
##
## Authors. CSC 
## Date. Feb 27th, 2024
########################################################################

## load libraries
library("Seurat")
library("tidyverse")
library("dplyr")
library("data.table")
library("magrittr")
library("here")

here::here()

## read input arguments
args = commandArgs(trailingOnly=TRUE)
cellranger_pipe <- args[2]
## For testing:
# cellranger_pipe <- "CR_arc_reanalyze" # This is multiome not QCed
# cellranger_pipe <- "CR_arc_reanalyze_outliers" # This is for only GEX
# cellranger_pipe <- "CR_arc_reanalyze_outliers_ATAC" # This is for only ATAC

## input directories
if (length(cellranger_pipe)) {
  
  ## Avoid to re-run data processed before
  if (cellranger_pipe=="CR_arc_reanalyze_outliers" || cellranger_pipe=="CR_arc_reanalyze" ) { stop() }
  
  message("CellRanger ARC input: ", cellranger_pipe)
  # Check/create directories
  
  if (cellranger_pipe=="CR_arc_reanalyze") {
    # stop()
    inputDir <- here("processed-data", "03_pseudobulking", "cellrangerARC_reanalyze")
    inputDir_cvs <- here("processed-data", "03_pseudobulking", "cellrangerARC_reanalyze", "cvs_files_markers")
    processedDir <- here("processed-data", "04_DiffExpr_Clustering_seurat", "cellrangerARC_reanalyze")
    cvsDir <- here("processed-data", "04_DiffExpr_Clustering_seurat", "cellrangerARC_reanalyze", "cvs_files_markers")
    
  } else { # CR_arc_reanalyze_outliers (GEX) or CR_arc_reanalyze_outliers_ATAC 
    
    inputDir <- here("processed-data", "03_pseudobulking", "cellrangerARC_reanalyze_outliers") 
    inputDir_cvs <- here("processed-data", "03_pseudobulking", "cellrangerARC_reanalyze_outliers", "cvs_files_markers")
    processedDir <- here("processed-data", "04_DiffExpr_Clustering_seurat", "cellrangerARC_reanalyze_outliers")
    cvsDir <- here("processed-data", "04_DiffExpr_Clustering_seurat", "cellrangerARC_reanalyze_outliers", "cvs_files_markers")
  }  
  
} else {
  
  message("Input argument missed")
  message("CellRanger input: ", cellranger_pipe)
  stop()
  
}
  
## Check directories
if (!dir.exists(processedDir)) {dir.create(processedDir)}
if (!dir.exists(cvsDir)) {dir.create(cvsDir)}

## Contains marker lists 
source(here("code", "04_DiffExpr_Clustering_seurat", "remote_DGE_marker_gene_lists.R"))       # Call functions to read paths


#############################           Initials        ################################

## Set count-mtx type and integration model (CCA or Harmony)

count_mtx_type <- 'norm_counts' 
Seurat_reduction <- 'Harmony' 
minCells <- 1
# Seurat_reduction <- 'CCA'
# count_mtx_type <- 'data_counts'
## Minimum cells by cluster

if (count_mtx_type=='data_counts') { Seurat_base_name <- 'seurat.data_counts' } else { Seurat_base_name <- 'seurat.norm_counts' }

## Build Seurat object name. `subset` suffix means clusters with fewer cells than `minCells` had been filtered. 
ifelse (Seurat_reduction=='CCA', Seurat_base_name <- paste0(Seurat_base_name, '_CCA_All'), Seurat_base_name <- paste0(Seurat_base_name, '_Harmony_All'))
                                                             
## Validate seurat exists
if (length(list.files(inputDir, pattern = Seurat_base_name)==1)) {
  message("Processing ", cellranger_pipe)
} else {
  message("Input seurat object missed!")
  stop()
}

message("Starting cell-type identification for `", cellranger_pipe, "`")

## Load Seurat Integrated with cluster information
seurat_RDSname <- case_when(
  cellranger_pipe == "CR_arc_reanalyze_outliers" ~ here(inputDir, paste0(Seurat_base_name, "_GEX_subset_Outliers.rds")),
  cellranger_pipe == "CR_arc_reanalyze_outliers_ATAC" ~ here(inputDir, paste0(Seurat_base_name, "_ATAC_subset_Outliers.rds")),
  .default = as.character(here(inputDir, paste0(Seurat_base_name, ".rds")))
)
# if (cellranger_pipe == "CR_arc_reanalyze_outliers" || cellranger_pipe == "CR_arc_reanalyze_outliers") {
#   seurat_RDSname <- here(inputDir, paste0(Seurat_base_name, "_GEX_subset_Outliers.rds"))
# } else {
#   seurat_RDSname <-here(inputDir, paste0(Seurat_base_name, ".rds"))
# }
basename(seurat_RDSname)
SeuratOBJ <- readRDS(seurat_RDSname)
## verification
length(Cells(x = SeuratOBJ))

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
cvs_name <- case_when(
  cellranger_pipe == "CR_arc_reanalyze_outliers" ~ here(cvsDir, paste0(Seurat_base_name, '_cluster_info_GEX.csv')),
  cellranger_pipe == "CR_arc_reanalyze_outliers_ATAC" ~ here(cvsDir, paste0(Seurat_base_name, '_cluster_info_ATAC.csv')),
  .default = as.character( here(cvsDir, paste0(Seurat_base_name, '_cluster_info.csv')))
)
## e.g: seurat.data_counts_Harmony_cluster_info.csv
write.csv(df_mdT, cvs_name)

## extract unique clusters in ascending order
clusters <- unique(df_mdT$seurat_clusters)
clusters <- as.integer(levels(clusters)[as.integer(clusters)])

message('Identifing cell types for ', length(clusters),' clusters from `', cellranger_pipe, '` dataset')

## Read All markers CVS file for all clusters
DGE_cvs_name <- case_when(
  cellranger_pipe == "CR_arc_reanalyze_outliers" ~ here(inputDir_cvs, paste0(Seurat_base_name, "markers_GEX.csv")) ,
  cellranger_pipe == "CR_arc_reanalyze_outliers_ATAC" ~ here(inputDir_cvs, paste0(Seurat_base_name, "markers_ATAC.csv")) ,
  .default = as.character(here(inputDir_cvs, paste0(Seurat_base_name, "markers.csv")) )
)
#DGE_cvs_name <- here(inputDir_cvs, paste0(Seurat_base_name, "markers_GEX.csv")) 
#seurat.norm_counts_Harmony_Allmarkers_GEX.csv

seurat_clust <- as.data.frame(read.csv(DGE_cvs_name, header = TRUE))
head(seurat_clust, n=3)

message('Parsing ', length(markers.custom), ' gene-markers lists on ', length(clusters) ,' clusters in `', cellranger_pipe, '` dataset ', 
        Seurat_reduction, ' reduction')


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
  prefix_name <- case_when(
    cellranger_pipe == "CR_arc_reanalyze_outliers" ~ paste0(names(markers.custom[idx_lst]), '_top', n_slice, "_GEX"),
    cellranger_pipe == "CR_arc_reanalyze_outliers_ATAC" ~ paste0(names(markers.custom[idx_lst]), '_top', n_slice, "_ATAC"),
    .default = as.character(paste0(names(markers.custom[idx_lst]), '_top', n_slice))
  )
  #prefix_name <- paste0(names(markers.custom[idx_lst]), '_top', n_slice, "_GEX")
  
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

