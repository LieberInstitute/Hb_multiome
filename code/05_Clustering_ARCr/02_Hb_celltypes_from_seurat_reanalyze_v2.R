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

tmp_file_name <- str_extract(Seurat_base_name, "ARC+.+")
# Ex. ARCr_QCed_WNN_k30_C.leiden_lsi_r1

cvs_name <- here(processedDir, paste0(tmp_file_name, '_cluster_info.csv'))
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

all_markers_cvs_name <- here(cvsDir, paste0(tmp_file_name, "_DEG_Top50_WNN.csv"))
write.csv(top_DGE_clust, all_markers_cvs_name)

message("\nSaved CVS file with Top50 DEG from: ", basename(all_markers_cvs_name))


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

## tbl to save top 50 genes by cluster
all_gene_match <- setNames(data.frame(matrix(ncol = 5, nrow = 0)),
                           c("Feature.ID", "Feature.Name", "p_val_adj", "cell_type", "cluster")) #Cluster.Adjusted.p.value

message("\nSearching cell-types for all gene markers lists")


## Annotate cell types based on the reference of gene markers DD+LB

for (clust in clusters) {
  # Test: clust <- 3
  
  message("Parsing cluster ", clust)

  DGE_by_clust <- top_DGE_clust |> filter(cluster == clust)

  if (!nrow(DGE_by_clust) > 0) {
    
    message("No information available for the cluster: ", clust)
    
  } else {
    
    for (gm_idx in seq_along(markers.lst) ) {
      # Test: gm_idx = 3
      gm_lst = markers.lst[[gm_idx]]
      gm_cell_type = names(markers.lst[gm_idx])
      # Match top N genes with the marker genes for the cell-type x
      gene_match <- DGE_by_clust |> filter_all(any_vars(. %in% gm_lst))
      gene_match
      # add matched genes to a dataframe
      if (nrow(gene_match) > 0 ) {
        names(gene_match)[names(gene_match) == clust ] <- "Cluster.Adjusted.p.value" # rename cols to rbind
        gene_match['cell_type']  <- gm_cell_type # md.csc 'cell_type' by 'cell-type'
        gene_match['cluster']  <- clust
        all_gene_match <- rbind(all_gene_match, gene_match)
      }
      
    }
    
  }
  
}

message(nrow(all_gene_match), " total matches.")

all_markers_cvs_name <- here(cvsDir, paste0(tmp_file_name, "_DEG_Top50_WNN_DD_LB_matching_markers.csv"))

write.csv(all_gene_match, all_markers_cvs_name, row.names=FALSE)

message("\nSaved CVS file with matching genes on: ", basename(all_markers_cvs_name))


## Joint Top50 DEG and add matching genes - Annotate cell types based on the reference of gene markers DD+LB

nrow(top_DGE_clust)
# [1] 1650
nrow(all_gene_match)
# [1] 375
# delete columns with redundant data
gene_match_subset <- all_gene_match |> select(gene, cell_type)

# keeps all observations in the top50 DEG and add `cell_type` column of matching genes
integrate_tbl <- left_join(top_DGE_clust, gene_match_subset, by = c("cluster", "gene"))
nrow(integrate_tbl)
print(integrate_tbl, n=50)
#    p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene       cell_type
#     <dbl>      <dbl> <dbl> <dbl>     <dbl>   <int> <chr>      <chr>    
# 1     0      -3.23 0.222 0.693         0       1 CPNE4      NA       
# 2     0      -3.69 0.184 0.64          0       1 VAV3       NA  
# 19     0     -2.07  0.399 0.752         0       1 FAT3       NA           
# 20     0     -2.48  0.09  0.441         0       1 PDE3A      NA           
# 21     0     -2.08  0.268 0.618         0       1 SNCA       NA           
# 22     0     -1.87  0.293 0.639         0       1 ST6GALNAC3 NA           
# 23     0      1.66  0.628 0.285         0       1 GABRG3     DD_Excit.Thal

message(nrow(integrate_tbl), " total matches.")

all_markers_cvs_name <- here(processedDir, paste0(tmp_file_name, "_DEG_Top50_WNN_DD_LB_matching_markers_integrated.csv"))

write.csv(integrate_tbl, all_markers_cvs_name, row.names=FALSE)

message("\nSaved CVS file with Top50 matching genes integrated on: ", basename(all_markers_cvs_name))
message(' Cell type identification completed!')



library("sessioninfo")
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

