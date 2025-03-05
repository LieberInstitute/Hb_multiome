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

## read input arguments from input_wnn_rds_names_v2 
Seurat_base_name <- commandArgs(trailingOnly = TRUE)
# Seurat_base_name = "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2" (selected)
# Seurat_base_name = "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k40_C.leiden_lsi_r2"

## Selected manually the clusters based on the cell-type identification gene-marker lists
marker_lst <- "integrated"
# marker_lst="literature_base"
# marker_lst="data_driven"

## input validations
if (length(Seurat_base_name) == 0 || length(Seurat_base_name) == 0) {
  message("Cellranger pipeline or input marker list missed or not valid!")
  stop()
} else {
  message("CellRanger dataset input: ", Seurat_base_name)
  message("Marker list input: ", (marker_lst))
}


message("Processing dataset ", Seurat_base_name)


## path to input Directory and suffix of cluster data

cvsDirIN <- here("processed-data", "05_Clustering_ARCr", "02_Hb_celltypes_from_seurat_reanalyze_v3", "cvs_files_markers")
cvsDirOUT <- here("processed-data", "05_Clustering_ARCr", "02_Hb_celltypes_from_seurat_reanalyze_v3")
message("Processing cell types generated with: `", basename(cvsDirIN), "` script")


##########These are CellRangerARC or CellRangerARC-reanalyze datasets

hb <- c()
thal <- c()
neu <- c()
glia <- c()
undeterminated <- c()
endo <- c() 

## LEIDEN resolution=r1

if (Seurat_base_name=="WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r1") {
  ## rna harmonized and atac harmonized
  hb <- c(2,4,7,9,15,25)
  thal <- c(1,5,8,16,21,23,26)
  endo <- c(22) 
  glia <- c(3,13,17,18,19,20,24)
  undeterminated <- c(6,10,11,14)
} else { (Seurat_base_name=="WNN_rnaHarm_atacLSI_k30_C.leiden_lsi_r1") 
  ## rna harmonized and lsi
  hb <- c(2,3,4,7,9,10,11,14)
  thal <- c(1,6,12,19,20,21,23,25,28,29,30)
  neu <- (13)
  endo <- c(16,31) 
  glia <- c(5,22,24,26,27)
  undeterminated <- c(8,15,17,18)
}

message("Groups of clusters assigned!")

## Prepare file name for retrieve cluster info from corresponding pipeline

input_csv <- here(cvsDirIN, paste0(Seurat_base_name, "_cluster_info.csv"))
basename(input_csv)

## Prepare file name for retrieve cell-types from corresponding pipeline

suffix_clust_names <- paste0(Seurat_base_name, "_cellTypes_integrated_top50.csv")
# [1] "*_cellTypes_integrated_top50.csv"

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
# [26] 26 27 28 29 30 31
if (!is_null(hb)) { lst_clust <- append(lst_clust, list(hb = hb)) }
if (!is_null(thal)) { lst_clust <- append(lst_clust, list(thal= thal)) }
if (!is_null(neu)) { lst_clust <- append(lst_clust, list(neu = neu)) }
if (!is_null(glia)) { lst_clust <- append(lst_clust, list(glia = glia)) }
if (!is_null(undeterminated)) { lst_clust <- append(lst_clust, list(undeterminated = undeterminated)) }
if (!is_null(endo)) { lst_clust <- append(lst_clust, list(endo = endo)) }

#lst_clust <- list(hb = hb, neu = neu, thal= thal, allTypes = allT)
message("Processing ", length(lst_clust$allTypes) , " clusters")
names(lst_clust)
lst_clust$allTypes

## Parse total cells by cluster and calculate percentages
## for testing: i <- 7 
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
                summarise(total_clust = round(sum(N), digits = 0))) |> # sum total cells by cluster
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
  totals_grp_clusters <- totals_grp_clusters |> mutate(Perc.Cluster = round((totals_grp_clusters$total_clust * 100 / total_cells), digits = 3))

  ## Load cluster info and cell-types to collapse names in `cell.type` column (description)

  top_deg_file <- here(cvsDirOUT, "cvs_files_markers", suffix_clust_names)
  basename(top_deg_file) # WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r1_cellTypes_integrated_top50.csv
  
  df_cluster_names <- read.csv(top_deg_file)
  df_cluster_names <- df_cluster_names |> drop_na(cell_type)
  df_cluster_names <- df_cluster_names[c("cluster", "cell_type")] |>
    group_by(cluster) |> summarise(cell_types = paste(cell_type, collapse=",")) |>
    rename(seurat_clusters = cluster) |> dplyr::filter(seurat_clusters %in% clust)

  ## Update in nice-readable format the `cell.types` column with `IDs+counts by marker` by cluster
  list_ct <- df_cluster_names[["cell_types"]]
  nc <- list()
  
  ## parse clusters to sum number of markers that match each cell-type. Ex. "inhibitory_neuron (2) Thalamus/MDm (1)" 
  for (x in seq_along(list_ct)) {
    ## test: x=2
    ct_cluster <- list_ct[x]
    x1 <- sapply(ct_cluster, function(x) strsplit(x, ","))
    for (ct in x1) {
      ids <- unique(ct)
      num_rep <- table(ct)
      col_new <- noquote(c(rbind(ids, paste0("(", num_rep, ") "))))
    }
    if (length(col_new)==2) {col_new <- paste0("*** ", col_new[1], col_new[2])}
    nc <- append(nc, paste0(col_new, collapse = " "))
  }
  
  df_cluster_names$cell_types <- noquote(unlist(nc))

  ## Save detail counts and percentages by grp of clusters
  totals_grp_clusters <- merge(grp_clusters, totals_grp_clusters, by = "seurat_clusters", all.x = TRUE, sort = FALSE)
  totals_grp_clusters <- merge(totals_grp_clusters, df_cluster_names, by = "seurat_clusters", all.x = TRUE, sort = FALSE)
  totals_grp_clusters["total_clust.y"] <- list(NULL) ## Delete column

  #cvs_name <- paste0("DETAIL_", marker_lst, "_" ,count_mtx_type, suffix, "_", clust_types ,"_v3.csv")    
  cvs_name <- paste0("DETAIL_", Seurat_base_name, "_" , clust_types ,"_v3.csv")    
  write.csv(totals_grp_clusters, here(cvsDirOUT, "cvs_files_markers", cvs_name))

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
  # colnames(totals_grp_clusters)
  # colnames(df_summary)
  df_summary <- rbind(totals_grp_clusters, df_summary)

  ## add suffix to file name to identify the corresponding report
  # cvs_name <- paste0('FULL_SUMMARY_LEIDENr1_knn30_WNN_rnaHarm_atacHarm_', marker_lst, "_", clust_types,"_v3.csv")
  cvs_name <- paste0('FULL_SUMMARY_',  Seurat_base_name, "_" , clust_types ,"_v3.csv") 
  write.csv(df_summary, here(cvsDirOUT, cvs_name), row.names=FALSE)

}

message("Done!")


# ##  slurm script reproducibility
# library("slurmjobs")
# job_single(
#   name = "03_cell_types_perc_WNN_leidenR1_knn30_v2", memory = "60G", cores = 2, create_shell = TRUE
# )

