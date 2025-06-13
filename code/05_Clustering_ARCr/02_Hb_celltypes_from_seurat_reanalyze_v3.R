########################################################################
## Read Seurat clusters to search cell-types based on a custom marker gene list. 
##      Remove redundant (duplicated) marker genes from DD gene markers reference
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
library("purrr")
library("data.table")
library("magrittr")
library("here")

here::here()

## read input arguments ( name of RDS Seurat file with wnn clustering to parse )
Seurat_base_name <- commandArgs(trailingOnly = TRUE)
# Debug with our best cluster representation is: 
# Seurat_base_name = "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2"

## Some WNN clustering results of interest. Testing:
# Seurat_base_name <- "seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.leiden_lsi_r1"
# Seurat_base_name <- "seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.louvain_lsi_r1"
# Seurat_base_name <- "seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.louvainM_lsi_r1"
# Seurat_base_name <- "seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.SLM_lsi_r1"

## new option with normalized data in both rna and atac
# Seurat_base_name <-  "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.louvainM_lsi_r1"
# Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r1"
# Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacLSI_k30_C.leiden_lsi_r1"
# Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2"


message(" Reading: ", Seurat_base_name)

## input directories

# Check/create directories
inputRDS_Dir <- here("processed-data", "05_Clustering_ARCr", "01_clustering_std_method")
inputDir_cvs <- here(inputRDS_Dir, "cvs_files")
processedDir <- here("processed-data", "05_Clustering_ARCr", "02_Hb_celltypes_from_seurat_reanalyze_v3")
cvsDir <- here(processedDir, "cvs_files_markers")

## Check directories
if (!dir.exists(processedDir)) {dir.create(processedDir)}
if (!dir.exists(cvsDir)) {dir.create(cvsDir)}

## Contains marker lists 
source(here("code", "04_DiffExpr_Clustering_seurat", "remote_DGE_marker_gene_lists.R"))   


#############################           Initials        ################################

message("Reading files to annotate cell-types in WNN clusters")
# seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r1
# seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2
# seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k40_C.leiden_lsi_r1
# seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k40_C.leiden_lsi_r2

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
## make slim Seurat 
SeuratOBJ <- DietSeurat(SeuratOBJ, 
                        assays = "RNA")

## verification
SeuratOBJ
length(Cells(x = SeuratOBJ))
message(nrow(unique(SeuratOBJ[["seurat_clusters"]])), " clusters found")

message("Seurat loaded!")
message("Starting cell-type identification ...")

## Gene markers lists. New function to join LB and DD gene markers lists

message("Load mixed gene reference: Literature-Based + Data-Driven")

markers.custom <- get_multiple_markers_genes_lst()
tmp <- names(markers.custom)
tmp <- paste(tmp, collapse=', ')
message("Processing ", length(markers.custom), " categories of gene-markers list \n *****(", tmp, ")*****")

## Check and remove duplicated genes markers

names(markers.custom)
x <- markers.custom
length(unlist(x)) # 481
# table(unname(unlist(x)))
v_dup <- duplicated(unname(unlist(x)))
dup_genes <- unname(unlist(x))[v_dup]

## delete duplicated genes before annotate

if (length(dup_genes>0)) {
  
  message("Eliminating ", length(dup_genes), " duplicated marker genes from DD gene markers reference")
  remove_lst <- dup_genes
  
  ## get DD genes and index list 
  DD <- c("DD_Astrocyte","DD_Endo","DD_Excit.Thal", "DD_Inhib.Thal", "DD_LHb","DD_MHb","DD_Microglia","DD_Oligo","DD_OPC")
  DD_idx <- which(names(markers.custom) %in% DD)
  
  ## remove duplicates from general DD+LB gene markers list
  f_remove_duplicates <- function(tmp, DD_idx, dup_genes){
    for (gen in dup_genes) { 
      # dd_idex = 1
      tmp <- map(DD_idx, ~ tmp[[.x]][tmp[[.x]] != gen]) 
    }
    return(tmp)
  }
  
  tmp <- f_remove_duplicates(markers.custom, DD_idx, dup_genes)
  for (i in DD_idx) { markers.custom[[i]] <- tmp[[i]] }
  
  ## verify
  x <- markers.custom
  v_dup <- duplicated(unname(unlist(x)))
  dup_genes <- unname(unlist(x))[v_dup]
  message("Removed from DD gene-markers list ", length(remove_lst), " duplicated genes!\nTotal kept it ", length(unlist(x)))
  
}


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
# total_sample <- df_mdT |>
#   group_by(orig.ident) |>
#   summarise(mean = mean(N))
# total_sample

## Save cluster information

# Extract base name to easily identify files  

tmp_file_name <- str_extract(Seurat_base_name, "WNN+.+")
# Ex. ARCr_QCed_WNN_k30_C.leiden_lsi_r1

cvs_name <- here(cvsDir, paste0(tmp_file_name, '_cluster_info.csv'))
write.csv(df_mdT, cvs_name)

## extract unique clusters in ascending order
clusters <- unique(df_mdT$seurat_clusters)
clusters <- as.integer(levels(clusters)[as.integer(clusters)])

## Read All markers CVS file for all clusters
DGE_cvs_name <- here(inputDir_cvs, paste0(Seurat_base_name, "_markers.csv"))

message('Saved cluster info. \nAnnotating cell types for ', length(clusters),' clusters using `', basename(DGE_cvs_name), "`")
# Annotating cell types for 42 clusters using `seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_markers.csv`

seurat_clust_DEG <- read.csv(DGE_cvs_name, header = TRUE)
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

message("Saved CVS file with Top50 DEG from: ", basename(all_markers_cvs_name))


####### Parse the 10/20 DGE genes from GEX cluster against the gene markers list provided ####### 

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

message("Searching cell-types for all gene markers in **literature_base** and **data_driven**")


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

# all_markers_cvs_name <- here(cvsDir, paste0(tmp_file_name, "_DEG_Top50_WNN_DD_LB_matching_markers.csv"))

# write.csv(all_gene_match, all_markers_cvs_name, row.names=FALSE)

# message("\nSaved CVS file with matching genes on: ", basename(all_markers_cvs_name))


## Joint Top50 DEG and add matching genes - Annotate cell types based on the reference of gene markers DD+LB

nrow(top_DGE_clust)
nrow(all_gene_match)
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
# 2050 total matches.

all_markers_cvs_name <- here(cvsDir, paste0(tmp_file_name, "_cellTypes_integrated_top50.csv"))
                                                                  
write.csv(integrate_tbl, all_markers_cvs_name, row.names=FALSE)

message("Saved CVS file with Top50 matching genes integrated on: ", basename(all_markers_cvs_name))
message('Cell type identification completed!')



library("sessioninfo")
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# [1] "Reproducibility information:"
# > Sys.time()
# [1] "2025-06-12 00:26:36 EDT"
# > proc.time()
# user   system  elapsed 
# 186.020   36.698 5671.495 
# > options(width = 120)
# > session_info()
# ─ Session info ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.2 Patched (2024-02-08 r85876)
# os       Rocky Linux 9.4 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2025-06-12
# pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3.x/bin/pandoc
# quarto   NA
# 
# ─ Packages ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package          * version    date (UTC) lib source
# abind              1.4-8      2024-09-12 [1] CRAN (R 4.3.2)
# BiocGenerics       0.48.1     2023-11-01 [2] Bioconductor
# BiocParallel       1.36.0     2023-10-24 [2] Bioconductor
# Biostrings         2.70.3     2024-03-13 [1] Bioconductor 3.18 (R 4.3.2)
# bitops             1.0-9      2024-10-03 [1] CRAN (R 4.3.2)
# cellranger         1.1.0      2016-07-27 [2] CRAN (R 4.3.2)
# cli                3.6.5      2025-04-23 [1] CRAN (R 4.3.2)
# cluster            2.1.6      2023-12-01 [3] CRAN (R 4.3.2)
# codetools          0.2-19     2023-02-01 [3] CRAN (R 4.3.2)
# colorspace         2.1-1      2024-07-26 [1] CRAN (R 4.3.2)
# cowplot            1.1.3      2024-01-22 [2] CRAN (R 4.3.2)
# crayon             1.5.3      2024-06-20 [1] CRAN (R 4.3.2)
# data.table       * 1.17.2     2025-05-12 [1] CRAN (R 4.3.2)
# deldir             2.0-2      2023-11-23 [2] CRAN (R 4.3.2)
# dichromat          2.0-0.1    2022-05-02 [2] CRAN (R 4.3.2)
# digest             0.6.37     2024-08-19 [1] CRAN (R 4.3.2)
# dotCall64          1.1-1      2023-11-28 [2] CRAN (R 4.3.2)
# dplyr            * 1.1.4      2023-11-17 [2] CRAN (R 4.3.2)
# farver             2.1.2      2024-05-13 [1] CRAN (R 4.3.2)
# fastDummies        1.7.3      2023-07-06 [2] CRAN (R 4.3.2)
# fastmap            1.2.0      2024-05-15 [1] CRAN (R 4.3.2)
# fastmatch          1.1-4      2023-08-18 [2] CRAN (R 4.3.2)
# fitdistrplus       1.1-11     2023-04-25 [2] CRAN (R 4.3.2)
# forcats          * 1.0.0      2023-01-29 [2] CRAN (R 4.3.2)
# future             1.33.1     2023-12-22 [2] CRAN (R 4.3.2)
# future.apply       1.11.1     2023-12-21 [2] CRAN (R 4.3.2)
# generics           0.1.4      2025-05-09 [1] CRAN (R 4.3.2)
# GenomeInfoDb       1.38.8     2024-03-15 [1] Bioconductor 3.18 (R 4.3.2)
# GenomeInfoDbData   1.2.11     2024-02-09 [2] Bioconductor
# GenomicRanges      1.54.1     2023-10-29 [2] Bioconductor
# ggplot2          * 3.5.2      2025-04-09 [1] CRAN (R 4.3.2)
# ggrepel            0.9.6      2024-09-07 [1] CRAN (R 4.3.2)
# ggridges           0.5.6      2024-01-23 [2] CRAN (R 4.3.2)
# globals            0.16.2     2022-11-21 [2] CRAN (R 4.3.2)
# glue               1.8.0      2024-09-30 [1] CRAN (R 4.3.2)
# goftest            1.2-3      2021-10-07 [2] CRAN (R 4.3.2)
# gridExtra          2.3        2017-09-09 [2] CRAN (R 4.3.2)
# gtable             0.3.6      2024-10-25 [1] CRAN (R 4.3.2)
# here             * 1.0.1      2020-12-13 [2] CRAN (R 4.3.2)
# hms                1.1.3      2023-03-21 [2] CRAN (R 4.3.2)
# htmltools          0.5.8.1    2024-04-04 [1] CRAN (R 4.3.2)
# htmlwidgets        1.6.4      2023-12-06 [2] CRAN (R 4.3.2)
# httpuv             1.6.16     2025-04-16 [1] CRAN (R 4.3.2)
# httr               1.4.7      2023-08-15 [2] CRAN (R 4.3.2)
# ica                1.0-3      2022-07-08 [2] CRAN (R 4.3.2)
# igraph             2.0.1.9008 2024-02-09 [2] Github (igraph/rigraph@39158c6)
# IRanges            2.36.0     2023-10-24 [2] Bioconductor
# irlba              2.3.5.1    2022-10-03 [2] CRAN (R 4.3.2)
# jsonlite           2.0.0      2025-03-27 [1] CRAN (R 4.3.2)
# KernSmooth         2.23-22    2023-07-10 [3] CRAN (R 4.3.2)
# later              1.4.2      2025-04-08 [1] CRAN (R 4.3.2)
# lattice            0.22-5     2023-10-24 [3] CRAN (R 4.3.2)
# lazyeval           0.2.2      2019-03-15 [2] CRAN (R 4.3.2)
# leiden             0.4.3.1    2023-11-17 [1] CRAN (R 4.3.2)
# lifecycle          1.0.4      2023-11-07 [2] CRAN (R 4.3.2)
# listenv            0.9.1      2024-01-29 [2] CRAN (R 4.3.2)
# lmtest             0.9-40     2022-03-21 [2] CRAN (R 4.3.2)
# lubridate        * 1.9.3      2023-09-27 [2] CRAN (R 4.3.2)
# magrittr         * 2.0.3      2022-03-30 [2] CRAN (R 4.3.2)
# MASS               7.3-60.0.1 2024-01-13 [3] CRAN (R 4.3.2)
# Matrix             1.6-5      2024-01-11 [3] CRAN (R 4.3.2)
# matrixStats        1.5.0      2025-01-07 [1] CRAN (R 4.3.2)
# mime               0.13       2025-03-17 [1] CRAN (R 4.3.2)
# miniUI             0.1.1.1    2018-05-18 [2] CRAN (R 4.3.2)
# nlme               3.1-164    2023-11-27 [3] CRAN (R 4.3.2)
# parallelly         1.36.0     2023-05-26 [2] CRAN (R 4.3.2)
# patchwork          1.2.0      2024-01-08 [2] CRAN (R 4.3.2)
# pbapply            1.7-2      2023-06-27 [2] CRAN (R 4.3.2)
# pillar             1.10.2     2025-04-05 [1] CRAN (R 4.3.2)
# pkgconfig          2.0.3      2019-09-22 [2] CRAN (R 4.3.2)
# plotly             4.10.4     2024-01-13 [2] CRAN (R 4.3.2)
# plyr               1.8.9      2023-10-02 [2] CRAN (R 4.3.2)
# png                0.1-8      2022-11-29 [2] CRAN (R 4.3.2)
# polyclip           1.10-6     2023-09-27 [2] CRAN (R 4.3.2)
# progressr          0.14.0     2023-08-10 [2] CRAN (R 4.3.2)
# promises           1.3.2      2024-11-28 [1] CRAN (R 4.3.2)
# purrr            * 1.0.4      2025-02-05 [1] CRAN (R 4.3.2)
# R6                 2.6.1      2025-02-15 [1] CRAN (R 4.3.2)
# RANN               2.6.1      2019-01-08 [2] CRAN (R 4.3.2)
# RColorBrewer       1.1-3      2022-04-03 [2] CRAN (R 4.3.2)
# Rcpp               1.0.14     2025-01-12 [1] CRAN (R 4.3.2)
# RcppAnnoy          0.0.22     2024-01-23 [2] CRAN (R 4.3.2)
# RcppHNSW           0.6.0      2024-02-04 [2] CRAN (R 4.3.2)
# RcppRoll           0.3.0      2018-06-05 [1] CRAN (R 4.3.2)
# RCurl              1.98-1.17  2025-03-22 [1] CRAN (R 4.3.2)
# readr            * 2.1.5      2024-01-10 [2] CRAN (R 4.3.2)
# readxl           * 1.4.3      2023-07-06 [2] CRAN (R 4.3.2)
# reshape2           1.4.4      2020-04-09 [2] CRAN (R 4.3.2)
# reticulate         1.35.0     2024-01-31 [2] CRAN (R 4.3.2)
# rlang              1.1.6      2025-04-11 [1] CRAN (R 4.3.2)
# ROCR               1.0-11     2020-05-02 [2] CRAN (R 4.3.2)
# rprojroot          2.0.4      2023-11-05 [2] CRAN (R 4.3.2)
# Rsamtools          2.18.0     2023-10-24 [2] Bioconductor
# RSpectra           0.16-2     2024-07-18 [1] CRAN (R 4.3.2)
# Rtsne              0.17       2023-12-07 [2] CRAN (R 4.3.2)
# S4Vectors          0.40.2     2023-11-23 [2] Bioconductor 3.18 (R 4.3.2)
# scales             1.4.0      2025-04-24 [1] CRAN (R 4.3.2)
# scattermore        1.2        2023-06-12 [2] CRAN (R 4.3.2)
# sctransform        0.4.1      2023-10-19 [2] CRAN (R 4.3.2)
# sessioninfo      * 1.2.3      2025-02-05 [1] CRAN (R 4.3.2)
# Seurat           * 5.0.1      2023-11-17 [2] CRAN (R 4.3.2)
# SeuratObject     * 5.0.1      2023-11-17 [2] CRAN (R 4.3.2)
# shiny              1.10.0     2024-12-14 [1] CRAN (R 4.3.2)
# Signac           * 1.13.0     2024-04-04 [1] CRAN (R 4.3.2)
# sp               * 2.1-3      2024-01-30 [2] CRAN (R 4.3.2)
# spam               2.10-0     2023-10-23 [2] CRAN (R 4.3.2)
# spatstat.data      3.0-4      2024-01-15 [2] CRAN (R 4.3.2)
# spatstat.explore   3.2-6      2024-02-01 [2] CRAN (R 4.3.2)
# spatstat.geom      3.2-8      2024-01-26 [2] CRAN (R 4.3.2)
# spatstat.random    3.2-2      2023-11-29 [2] CRAN (R 4.3.2)
# spatstat.sparse    3.0-3      2023-10-24 [2] CRAN (R 4.3.2)
# spatstat.utils     3.0-5      2024-06-17 [1] CRAN (R 4.3.2)
# stringi            1.8.7      2025-03-27 [1] CRAN (R 4.3.2)
# stringr          * 1.5.1      2023-11-14 [2] CRAN (R 4.3.2)
# survival           3.5-7      2023-08-14 [3] CRAN (R 4.3.2)
# tensor             1.5        2012-05-05 [2] CRAN (R 4.3.2)
# tibble           * 3.2.1      2023-03-20 [2] CRAN (R 4.3.2)
# tidyr            * 1.3.1      2024-01-24 [2] CRAN (R 4.3.2)
# tidyselect         1.2.1      2024-03-11 [1] CRAN (R 4.3.2)
# tidyverse        * 2.0.0      2023-02-22 [2] CRAN (R 4.3.2)
# timechange         0.3.0      2024-01-18 [2] CRAN (R 4.3.2)
# tzdb               0.4.0      2023-05-12 [2] CRAN (R 4.3.2)
# utf8               1.2.5      2025-05-01 [1] CRAN (R 4.3.2)
# uwot               0.2.3      2025-02-24 [1] CRAN (R 4.3.2)
# vctrs              0.6.5      2023-12-01 [2] CRAN (R 4.3.2)
# viridisLite        0.4.2      2023-05-02 [2] CRAN (R 4.3.2)
# withr              3.0.2      2024-10-28 [1] CRAN (R 4.3.2)
# xtable             1.8-4      2019-04-21 [2] CRAN (R 4.3.2)
# XVector            0.42.0     2023-10-24 [2] Bioconductor
# zlibbioc           1.48.2     2024-03-13 [1] Bioconductor 3.18 (R 4.3.2)
# zoo                1.8-12     2023-04-13 [2] CRAN (R 4.3.2)
# 
# [1] /users/csoto/R/4.3.x
# [2] /jhpce/shared/community/core/conda_R/4.3.x/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.3.x/R/lib64/R/library
# * ── Packages attached to the search path.
# 
