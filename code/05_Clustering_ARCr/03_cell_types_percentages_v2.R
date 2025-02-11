###############################################################################
##
##  Calculate and summarize percentage of pre-selected cell types in Leiden r=1 knn=30
##
###############################################################################

library("stringr")
library("tidyverse")
library("janitor")
library("scales")
library("data.table")
library("here")

here::here()
options(digits=2)

## read input arguments
# args = commandArgs(trailingOnly=TRUE)
# cellranger_pipe <- args[2]
# marker_lst <- args[4]
cellranger_pipe = "CR_arc_reanalyze" # processed and QCed multiome data

# ## We are only using the norm count with Harmony
# count_mtx_type <- 'norm_counts'
# Seurat_reduction <- 'Harmony'

## Selected manually the clusters based on the cell-type identification gene-marker lists
marker_lst <- "integrated"
# marker_lst="literature_base"
# marker_lst="data_driven"

## input validations
if (length(cellranger_pipe) == 0 || length(marker_lst) == 0) {
  message("Cellranger pipeline or input marker list missed or not valid!")
  stop()
} else {
  message("CellRanger dataset input: ", cellranger_pipe)
  message("Marker list input: ", (marker_lst))
}


message("Processing dataset ", cellranger_pipe)


## path to input Directory and suffix of cluster data

cvsDirIN <- here("processed-data", "05_Clustering_ARCr", "02_Hb_celltypes_from_seurat_reanalyze_v2", "cvs_files_markers")
cvsDirOUT <- here("processed-data", "05_Clustering_ARCr", "02_Hb_celltypes_from_seurat_reanalyze_v2")
message("Processing cell types generated with: `", basename(cvsDirIN), "` script")

suffix <- '_Harmony_All_cluster'


##########These are CellRangerARC or CellRangerARC-reanalyze datasets

## LEIDEN resolution=r1

if (cellranger_pipe=="CR_arc_reanalyze") {
  hb <- c(2,3,5,7,9,10,12,17)
  thal <- c(1,8,11,19,21,22,25,26,29,30,31)
  endo <- c(13,32) 
  glia <- c(4,15,23,24,27,28)
  undeterminated <- c(6,14,16,18,20,33)
} 

# ## LEIDEN r2
# if (cellranger_pipe=="CR_arc_reanalyze") {
#   hb <- c(1,3,6,8,9,11,13,15,21,25,41)
#   thal <- c(7,16,19,20,22,23,27,29,34,35,37,39,42)
#   endo <- c(18,40,43) 
#   glia <- c(4,24,26,28,31,32,33,36)
#   undeterminated <- c(2,5,10,12,14,17,20,38)
# } 


message(cellranger_pipe, " groups of clusters assigned!")


## Read cluster info. Set count-mtx type and integration model (CCA or Harmony)

count_mtx_type <- "norm_counts"
Seurat_base_name <- "ARCr_QCed_WNN_k30_C.leiden_lsi_r1"

## Prepare file name for retrieve cluster info from corresponding pipeline

input_csv <- here(cvsDirIN, paste0(Seurat_base_name, "_cluster_info.csv"))
basename(input_csv)

## Prepare file name for retrieve cell-types from corresponding pipeline

suffix_clust_names <- paste0(Seurat_base_name, "_cellTypes_integrated_top50.csv")
basename(suffix_clust_names)
# [1] "ARCr_QCed_WNN_k30_C.leiden_lsi_r1_cellTypes_integrated_top50.csv"

message("Calculating percentage of neurons, habenula, thalamus and glia cells")
message("Active gene marker: ", marker_lst)

###### Load clusters by sample from processed_data directory ######

## Format sample ID names and reorder
df_mdT <- read.csv(input_csv, row.names = 1)
df_mdT$orig.ident <- sprintf("S%02d_Hb_r", parse_number(df_mdT$orig.ident))
ref_sort <- sort(unique(df_mdT$orig.ident))
df_mdT <- df_mdT[order(sapply(df_mdT$orig.ident, function(x) which(x == ref_sort))), ]

unique(df_mdT$orig.ident)

### Preparing list of clusters to summarize
allT <- c(sort(unique(df_mdT[["seurat_clusters"]])))

## Build the list of lists with grouo of cluster to quantify (%)
lst_clust <- list(allTypes = allT)
# [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
# [26] 26 27 28 29 30 31 32 33
if (!is_null(hb)) { lst_clust <- append(lst_clust, list(hb = hb)) }
if (!is_null(thal)) { lst_clust <- append(lst_clust, list(thal= thal)) }
if (!is_null(neu)) { lst_clust <- append(lst_clust, list(neu = neu)) }
if (!is_null(glia)) { lst_clust <- append(lst_clust, list(glia = glia)) }
if (!is_null(undeterminated)) { lst_clust <- append(lst_clust, list(undeterminated = undeterminated)) }
if (!is_null(endo)) { lst_clust <- append(lst_clust, list(endo = endo)) }

#lst_clust <- list(hb = hb, neu = neu, thal= thal, allTypes = allT)
message("Group of clusters ready!")
names(lst_clust)
lst_clust$allTypes

## Parse total cells by cluster and calculate percentages
## for testing: i <- 1
for (i in seq_along(lst_clust)) {

  clust <- lst_clust[[i]]
  clust_types <- names(lst_clust)[[i]]
  message(length(clust), " clusters to compute in ", clust_types, " group")
  ref_sort <- sort(unique(df_mdT$orig.ident))
  
  ## Filtering clusters of interest in the given samples
  grp_clusters <- df_mdT |>
    dplyr::filter(seurat_clusters %in% clust) |>
    pivot_wider(names_from = orig.ident, values_from = N) |>  # make wider format the table
    left_join(df_mdT |>
                group_by(seurat_clusters) |>
                summarise(total_clust = sum(N))) |> # sum total cells by cluster
    mutate(seurat_clusters = as.character(seurat_clusters))
  ## replace NAs
  grp_clusters[is.na(grp_clusters)] <- 0
  colnames(grp_clusters)

  ########  Calculate perceptual values by CLUSTER ########

  totals_grp_clusters <- grp_clusters
  
  ref_sort <- intersect(colnames(grp_clusters), ref_sort)
  
  ## Add total cells by cluster in the group
  totals_grp_clusters <- totals_grp_clusters |> select(all_of(c(ref_sort,"seurat_clusters"))) |> pivot_longer(!seurat_clusters)
  totals_grp_clusters <- totals_grp_clusters |> group_by(seurat_clusters) |> summarise(total_clust = sum(value))

  ## Add percent cells by cluster in the group
  total_cells <- sum(df_mdT$N)
  totals_grp_clusters <- totals_grp_clusters |> mutate(Perc.Cluster = totals_grp_clusters$total_clust * 100 / total_cells)

  ## Load cluster info and cell-types to collapse names in `cell.type` column (description)

  top_deg_file <- here(cvsDirOUT, "cvs_files_markers", suffix_clust_names)
  #top_deg_file <- here(cvsDirOUT, "cvs_files_markers", paste0(Seurat_base_name, suffix_clust_names))
  df_cluster_names <- read.csv(top_deg_file)
  df_cluster_names <- df_cluster_names[c("cluster", "cell_type")] |>
    group_by(cluster) |> summarise(cell_types = paste(cell_type, collapse=",")) |>
    rename(seurat_clusters = cluster) |> dplyr::filter(seurat_clusters %in% clust)

  ## Update in nice-readable format the `cell.types` column with `IDs+counts by marker` by cluster
  list_ct <- df_cluster_names[["cell_types"]]
  nc <- list()
  
  ## add number of markers that match each cell-type. Ex. "inhibitory_neuron (2) Thalamus/MDm (1)" 
  for (x in seq_along(list_ct)) {
    ct_cluster <- list_ct[x]
    x1 <- sapply(ct_cluster, function(x) strsplit(x, ","))
    for (ct in x1) {
      ids <- unique(ct)
      num_rep <- table(ct)
      col_new <- noquote(c(rbind(ids, paste0("(", num_rep, ")"))))
    }
    if (length(col_new)==2) {col_new <- paste0("***", col_new[1], col_new[2])}
    nc <- append(nc, paste0(col_new, collapse = " "))
  }
  
  df_cluster_names$cell_types <- noquote(unlist(nc))

  ## Save detail counts and percentages by grp of clusters
  totals_grp_clusters <- merge(grp_clusters, totals_grp_clusters, by = "seurat_clusters", all.x = TRUE, sort = FALSE)
  totals_grp_clusters <- merge(totals_grp_clusters, df_cluster_names, by = "seurat_clusters", all.x = TRUE, sort = FALSE)
  totals_grp_clusters["total_clust.y"] <- list(NULL) ## Delete column

  # cvs_name <- paste0("DETAIL_", marker_lst, "_" ,count_mtx_type, suffix, "_", clust_types ,"_v3.csv")    
  # write.csv(totals_grp_clusters, here(cvsDirOUT, "cvs_files_markers", cvs_name))

  ########  Calculate perceptual values by SAMPLE ########

  ## Add total cells by sample in the group
  df_summary <- totals_grp_clusters |> select(all_of(ref_sort)) |> pivot_longer(cols = everything(), names_to = "sampleID") |>  
    group_by(sampleID) |> reframe(total_sample = sum(value))

  ## Calculate `total cells` by sample
  df_summary <- full_join(df_summary, (df_mdT |> 
    group_by(orig.ident) |> select(orig.ident, N) |>
      summarise(total_sample = sum(N))), join_by(sampleID == orig.ident)) |>
    mutate(Percentage.Sample = total_sample.x * 100 / total_sample.y)
  df_summary <- subset(df_summary, sampleID %in% c(ref_sort))
  
  df_summary <- rename(df_summary, all_of(c(Total.Cells.Sample = "total_sample.x", Total.Cells.ALL = "total_sample.y") ))
  df_summary <- as.data.frame(t(df_summary)) |> row_to_names(1)
  ## replace NAs
  df_summary[is.na(df_summary)] <- 0
  df_summary

  # cvs_name <- paste0('SUMMARY_', marker_lst, "_" ,count_mtx_type, suffix, "_", clust_types,"_v3.csv")
  # write.csv(df_summary, here(cvsDirOUT, "cvs_files_markers", cvs_name), row.names=TRUE)

  ## bind detail with totals into the main table
  df_summary <- data.frame(seurat_clusters = c("Total.Cells.Sample", "Total.Cells.ALL", "Percentage.Sample"), df_summary)
  new_columns <- list(total_clust.x = NA, Perc.Cluster = NA, cell_types = NA)
  df_summary <- df_summary |> mutate(!!!new_columns)
  colnames(totals_grp_clusters)
  colnames(df_summary)
  df_summary <- rbind(totals_grp_clusters, df_summary)

  ## add suffix to file name to identify the corresponding report
  cvs_name <- paste0('FULL_SUMMARY_LEIDENr1_knn30_', marker_lst, "_" , count_mtx_type, suffix, "_", clust_types,"_v1.csv")
  
  write.csv(df_summary, here(cvsDirOUT, cvs_name), row.names=FALSE)

}

message("Done!")


# # slurm script reproducibility
# library("slurmjobs")
# slurmjobs::job_loop(
#   loops = list(cellranger_pipe = c("CR_crossBarcodes", "CR_complementBarcodes", "CR_arc_reanalyze", "CR_arc_reanalyze_outliers", "CR_arc_reanalyze_outliers_ATAC"),
#                marker_lst = c("literature_base", "data_drive")),
#   name = "02_cell_types_percentages_v4",
#   cores = 2,
#   create_shell = TRUE
# )


