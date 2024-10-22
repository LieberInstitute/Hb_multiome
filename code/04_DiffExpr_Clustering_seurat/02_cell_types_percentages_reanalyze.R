###############################################################################
##
##  Calculate and summarize percentage statistics of pre-selected cell types in different categories
## 
###############################################################################

library("stringr")
library("tidyverse")
library("scales")
library("here")
library("ggplot2")

options(digits=2)

## Data from cellrangerARC-reanalyze
cvsDir <- here("processed-data", "04_DiffExpr_Clustering_seurat", "cellrangerARC_reanalyze")
plotDir <- here("plots", "04_DiffExpr_Clustering_seurat", "cellrangerARC_reanalyze")

## Data from cellranger-count
# cvsDir <- here("processed-data", "04_DiffExpr_Clustering_seurat", "cellranger_count")
# plotDir <- here("plots", "04_DiffExpr_Clustering_seurat", "cellranger_count")

if (!dir.exists(cvsDir)) { dir.create(cvsDir) }
if (!dir.exists(plotDir)) { dir.create(plotDir) }

## read input arguments
args = commandArgs(trailingOnly=TRUE)
marker_lst <- args[2]

## Selected manually the clusters based on the cell-type identification gene-marker lists
## args opt for testing: marker_lst="literature_base" 
##                       marker_lst="data_driven" 
if (marker_lst=="literature_base") {
  neu <- c(0,3,4,5,8,9,17,22,24) 
  hb <- c(3,4,5,9) 
  thal <- c(0,8,24)
} else if (marker_lst=="data_driven") {
  neu <- c(0,1,2,3,4,5,8,9,10,11,14,16,21,24,26) 
  hb <- c(2,4,5,9,10,14,16)
  thal <- c(0,1,3,8,11,21,24)
} else {
  message("Input marker list not valid!")
  stop()
}



# By now only available for norm count with Harmomy
#count_mtx_type <- 'data_counts'
count_mtx_type <- 'norm_counts' 
Seurat_reduction <- 'Harmony'

## Set count-mtx type and integration model (CCA or Harmony)
if (count_mtx_type=='data_counts') { Seurat_base_name <- 'seurat.data_counts' } else { Seurat_base_name <- 'seurat.norm_counts' }
## Compose Seurat object name processed before

if (Seurat_reduction=='CCA') {

    message("Calculating percentage of neurons, habenula and thalamus cell types for ", Seurat_reduction)
  # suffix <- '_CCA_cluster'
  # suffix_clust_names <- '_CCA_subset_cell_types_all_gm20.csv'
  # input_csv <- here(cvsDir, "cvs_files_markers", paste0(Seurat_base_name, suffix))
  
} else if (Seurat_reduction=='Harmony') {
  message("Calculating percentage of neurons, habenula and thalamus cell types for ", Seurat_reduction)
  suffix <- '_Harmony_All_cluster'
  input_csv <- here(cvsDir, "cvs_files_markers", paste0(Seurat_base_name, suffix))
  if (marker_lst=="literature_base") {
    suffix_clust_names <- '_Harmony_All_cellTypes_literature_base_top20.csv'  
    suffix_clust_plt <- '_Harmony_All_cellTypes_literature_base_top20'  
  } else if (marker_lst=="data_driven") {
    suffix_clust_names <- '_Harmony_All_cellTypes_data_driven_top20.csv'
    suffix_clust_plt <- '_Harmony_All_cellTypes_data_driven_top20'
  }
} else {
  message("Not file available!")
  stop()
}
  
input_csv <- paste0(input_csv, "_info.csv")

message("Loading information from ", basename(input_csv), ". Calculating percentages in progress ... ")

###### Load clusters by sample from processed_data directory ######

## Format sample ID names and reorder
df_mdT <- read.csv(input_csv, row.names = 1)
df_mdT$orig.ident <- sprintf("S%02d_Hb_r", extract_numeric(df_mdT$orig.ident))
ref_sort <- sort(unique(df_mdT$orig.ident))
df_mdT <- df_mdT[order(sapply(df_mdT$orig.ident, function(x) which(x == ref_sort))), ]
#print(unique(df_mdT$orig.ident))

### Preparing list of clusters to summarize 
## Add all cell-types to the end of the list of lists
allT <- c(sort(unique(df_mdT[["seurat_clusters"]])))
lst_clust <- list(hb = hb, neu = neu, thal= thal, allTypes = allT) 
names(lst_clust)


## Parse total cells by cluster and calculate percentages
## for testing: i <- 1
for (i in seq_along(lst_clust)) {  
  
  clust <- lst_clust[[i]]
  clust_types <- names(lst_clust)[[i]]
  message(length(clust), " clusters to summarize in ", clust_types)
  
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

  ########  Calculate perceptual values by CLUSTER ######## 
  
  totals_grp_clusters <- grp_clusters

  ## Add total cells by cluster in the group
  totals_grp_clusters <- totals_grp_clusters |> select(all_of(c(ref_sort,"seurat_clusters"))) |> pivot_longer(!seurat_clusters)
  totals_grp_clusters <- totals_grp_clusters |> group_by(seurat_clusters) |> summarise(total_clust = sum(value))
  
  ## Add percent cells by cluster in the group
  total_cells <- sum(df_mdT$N)
  totals_grp_clusters <- totals_grp_clusters |> mutate(Perc.Cluster = totals_grp_clusters$total_clust * 100 / total_cells)
  
  ## Load cluster info and cell-types to collapse names in `cell.type` column (description) 
  df_cluster_names <- read.csv(here(cvsDir, "cvs_files_markers", paste0(Seurat_base_name, suffix_clust_names)))
  df_cluster_names <- df_cluster_names[c("cluster", "cell.type")] |> 
    group_by(cluster) |> summarise(cell.types = paste(cell.type, collapse=",")) |>
    rename(seurat_clusters = cluster) |> dplyr::filter(seurat_clusters %in% clust)
  
  ## Update in nice-readable format the `cell.types` column with `IDs+counts by marker` by cluster
  list_ct <- df_cluster_names[["cell.types"]]
  nc <- list()
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
  df_cluster_names$cell.types <- noquote(unlist(nc))
  
  ## Save detail counts and percentages by grp of clusters
  totals_grp_clusters <- merge(grp_clusters, totals_grp_clusters, by = "seurat_clusters", all.x = TRUE, sort = FALSE)
  totals_grp_clusters <- merge(totals_grp_clusters, df_cluster_names, by = "seurat_clusters", all.x = TRUE, sort = FALSE)
  
  cvs_name <- paste0("DETAIL_", marker_lst, "_" ,count_mtx_type, suffix, "_", clust_types ,"_v2.csv")
  #write.csv(grp_clusters, here(cvsDir, cvs_name)) 
  write.csv(totals_grp_clusters, here(cvsDir, cvs_name))
  
  
  ########  Calculate perceptual values by SAMPLE ######## 
  
  ## Add total cells by sample in the group
  df_summary <- totals_grp_clusters |> select(all_of(ref_sort)) |> pivot_longer(cols = everything(), names_to = "sampleID") |> 
    group_by(sampleID) |> reframe(total_sample = sum(value))

  ## Calculate `total cells` by sample
  df_summary <- full_join(df_summary, (df_mdT |>
    group_by(orig.ident) |> select(orig.ident, N) |>
      summarise(total_sample = sum(N))), join_by(sampleID == orig.ident)) |>
    mutate(Percentage.Sample = total_sample.x * 100 / total_sample.y)
  
  df_summary <- rename(df_summary, all_of(c(Total.Cells.Sample = "total_sample.x", Total.Cells.ALL = "total_sample.y") ))
  df_summary <- as.data.frame(t(df_summary)) |> row_to_names(1)
  
  cvs_name <- paste0('SUMMARY_', marker_lst, "_" ,count_mtx_type, suffix, "_", clust_types,"_v2.csv")
  write.csv(df_summary, here(cvsDir, cvs_name), row.names=FALSE)
  
  ## bind rows to the main table
  df_summary <- data.frame(seurat_clusters = c("Total.Cells.Sample", "Total.Cells.ALL", "Percentage.Sample"), df_summary)
  new_columns <- list(total_clust = NA, Perc.Cluster = NA, cell.types = NA)
  df_summary <- df_summary |> mutate(!!!new_columns)
  df_summary <- rbind(totals_grp_clusters, df_summary)
  
  cvs_name <- paste0('SUMMARY_', marker_lst, "_" ,count_mtx_type, suffix, "_", clust_types,"_v2.csv")
  write.csv(df_summary, here(cvsDir, cvs_name), row.names=FALSE)

}

message("Done!")


## slurm script reproducibility

# slurmjobs::job_loop(
#   loops = list(marker_lst = c("literature_base", "data_drive")),
#   name = "02_cell_types_percentages_reanalyze",
#   cores = 2,
#   create_shell = TRUE
# )


