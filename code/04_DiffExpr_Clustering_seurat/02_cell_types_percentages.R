###############################################################################
##
##  Calculate and summarize percentage statistics of pre-selected cell types in different categories
##
###############################################################################

library("stringr")
library("tidyverse")
library("janitor")
library("scales")
library("here")

here::here()
options(digits=2)

## read input arguments
args = commandArgs(trailingOnly=TRUE)
cellranger_pipe <- args[2]
marker_lst <- args[4]
## We are only using the norm count with Harmony
count_mtx_type <- 'norm_counts'
Seurat_reduction <- 'Harmony'

## Selected manually the clusters based on the cell-type identification gene-marker lists
## args opt for testing:
# marker_lst="literature_base"
# marker_lst="data_driven"
# cellranger_pipe = "CR_crossBarcodes"
# cellranger_pipe = "CR_complementBarcodes"
# cellranger_pipe = "CR_arc_reanalyze"
# cellranger_pipe = "CR_arc_reanalyze_outliers" 
# cellranger_pipe = "CR_arc_reanalyze_outliers" # new argument for process ATAC outliers

## input validations
if (length(cellranger_pipe)== 0 || length(marker_lst)== 0) {
  message("Cellranger pipeline or input marker list missed or not valid!")
  stop()
} else {
  message("CellRanger input: ", cellranger_pipe)
  message("Marker list input: ", (marker_lst))
}

## Avoid to re-run data processed before
if (cellranger_pipe=="CR_crossBarcodes" || cellranger_pipe=="CR_complementBarcodes" || cellranger_pipe=="CR_arc_reanalyze" ) { stop() }

## path to input Directory and suffix of cluster data
cvsDir <- case_when(
  cellranger_pipe == "CR_crossBarcodes" ~ here("processed-data", "04_DiffExpr_Clustering_seurat", "cellranger_count"),
  cellranger_pipe == "CR_complementBarcodes" ~ here("processed-data", "04_DiffExpr_Clustering_seurat", "cellranger_count_complement"), 
  cellranger_pipe == "CR_arc_reanalyze" ~ here("processed-data", "04_DiffExpr_Clustering_seurat", "cellrangerARC_reanalyze"),
  cellranger_pipe == "CR_arc_reanalyze_outliers" ~ here("processed-data", "04_DiffExpr_Clustering_seurat", "cellrangerARC_reanalyze_outliers")
)
suffix <- case_when(
  cellranger_pipe == "CR_crossBarcodes" ~ '_Harmony_All_cluster',
  cellranger_pipe == "CR_complementBarcodes" ~ '_Harmony_All_subset_cluster',
  cellranger_pipe == "CR_arc_reanalyze" ~ '_Harmony_All_cluster',
  cellranger_pipe == "CR_arc_reanalyze_outliers" ~ '_Harmony_All_cluster'
)

## CellRanger-count pipelines
## Manually pre-selected clusters for all the available CellRanger datasets

if (cellranger_pipe=="CR_crossBarcodes" || cellranger_pipe=="CR_complementBarcodes") {
  
  if (marker_lst=="literature_base") {
    neu <- c(0,1,2,5,7,8,10,13,14,17,25,26,29)
    hb <- c(0,2,5,7,10,13,14)
    thal <- c()
    glia_other <- c()
  } else if (marker_lst=="data_driven") {
    neu <- c(0,1,2,5,7,8,9,10,13,14,15,16,26,27,28,29)
    hb <- c(2,5,7,10,13,14,16,27)
    thal <- c(1,8,15,26,28,29)
    glia_other <- c()
  }

} else if (cellranger_pipe=="CR_arc_reanalyze" || cellranger_pipe=="CR_arc_reanalyze_outliers") {
  
  if (cellranger_pipe=="CR_arc_reanalyze") {
    
      if (marker_lst=="literature_base") {
        neu <- c(0,3,4,5,8,9,17,22,24)
        hb <- c(3,4,5,9)
        thal <- c(8)
        glia_other <- c()
        
      } else if (marker_lst=="data_driven") {
        neu <- c(0,1,2,3,4,5,8,9,10,11,14,16,21,24,26)
        hb <- c(2,4,5,9,14,16)
        thal <- c(0,1,3,8,21,24)
        glia_other <- c()
      }
    
  } else { #CR_arc_reanalyze_outliers
    
    if (marker_lst=="literature_base") {
      hb <- c(2,3,5,9)
      thal <- c(8)
      neu <- c(0,4,15, hb, thal)
      glia <- c(7,12,17,19,20,22)
      undeterminated <- c()
      endo <- c(23, 27) #endo
      
    } else if (marker_lst=="data_driven") {
      hb <- c(2,4,5,9,14,16)
      thal <- c(0,1,8,24)
      neu <- c(10, hb, thal) 
      glia <- c(3,7,12,17,19,20,22,28)
      undeterminated <- c(6,21) # glia ^ neuron)
      endo <- c(23) # endo 
    } 
  }
}

# message("Removed previous summary reports from `", cellranger_pipe ,"` dataset!")

## Set count-mtx type and integration model (CCA or Harmony)
if (count_mtx_type=='data_counts') { Seurat_base_name <- 'seurat.data_counts' } else { Seurat_base_name <- 'seurat.norm_counts' }

## Build CSV file with DEG from corresponding Seurat processed before
input_csv <- here(cvsDir, "cvs_files_markers", paste0(Seurat_base_name, suffix, "_info"))
## add suffix to read the corresponding file
input_csv <- case_when(
  cellranger_pipe == "CR_crossBarcodes" || cellranger_pipe == "CR_complementBarcodes" || cellranger_pipe == "CR_arc_reanalyze" ~ input_csv,
  cellranger_pipe == "CR_arc_reanalyze_outliers" ~ paste0(input_csv, "_GEX"),
  cellranger_pipe == "CR_arc_reanalyze_outliers_ATAC" ~ paste0(input_csv, "_ATAC")
)
input_csv <- paste0(input_csv, ".csv")
basename(input_csv)
# seurat.norm_counts_Harmony_All_cluster_info_GEX.csv

## build file name for cell-types annotations
if (marker_lst=="literature_base") {
  suffix_clust_names <- '_Harmony_All_cellTypes_literature_base_top20'
} else if (marker_lst=="data_driven") {
  suffix_clust_names <- '_Harmony_All_cellTypes_data_driven_top20'
}
## add suffix to read the corresponding file
suffix_clust_names <- case_when(
  cellranger_pipe == "CR_crossBarcodes" || cellranger_pipe == "CR_complementBarcodes" || cellranger_pipe == "CR_arc_reanalyze" ~ suffix_clust_names,
  cellranger_pipe == "CR_arc_reanalyze_outliers" ~ paste0(suffix_clust_names, "_GEX"),
  cellranger_pipe == "CR_arc_reanalyze_outliers_ATAC" ~ paste0(suffix_clust_names, "_ATAC")
)
suffix_clust_names <- paste0(suffix_clust_names, ".csv")
basename(suffix_clust_names)
# "_Harmony_All_cellTypes_literature_base_top20_GEX.csv"

message("Calculating percentage of neurons, habenula and thalamus cell types")
message("Active gene marker list ref: ", marker_lst)

###### Load clusters by sample from processed_data directory ######

## Format sample ID names and reorder
df_mdT <- read.csv(input_csv, row.names = 1)
df_mdT$orig.ident <- sprintf("S%02d_Hb_r", parse_number(df_mdT$orig.ident))
ref_sort <- sort(unique(df_mdT$orig.ident))
df_mdT <- df_mdT[order(sapply(df_mdT$orig.ident, function(x) which(x == ref_sort))), ]
#print(unique(df_mdT$orig.ident))

### Preparing list of clusters to summarize
allT <- c(sort(unique(df_mdT[["seurat_clusters"]])))

## Build the list with cluster to parse
lst_clust <- list(allTypes = allT)
if (!is_null(hb)) { lst_clust <- append(lst_clust, list(hb = hb)) }
if (!is_null(thal)) { lst_clust <- append(lst_clust, list(thal= thal)) }
if (!is_null(neu)) { lst_clust <- append(lst_clust, list(neu = neu)) }
if (!is_null(glia)) { lst_clust <- append(lst_clust, list(glia = glia)) }
if (!is_null(undeterminated)) { lst_clust <- append(lst_clust, list(undeterminated = undeterminated)) }
if (!is_null(endo)) { lst_clust <- append(lst_clust, list(endo = endo)) }

#lst_clust <- list(hb = hb, neu = neu, thal= thal, allTypes = allT)
message("Group of clusters prepared: ")
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
  ## Note. Alternative can be used the re-run of DEG calculated in the complement Seurat from CR-count after removed the cross-barcodes
  ##                  <<Directory: top20_subset>>
  ##      Here I am using the DEG originally calculated in the CR-count harmnony data
  ##      df_cluster_names <- read.csv(here(cvsDir, "cvs_files_markers", "top20_subset", paste0(Seurat_base_name, suffix_clust_names)))
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
  totals_grp_clusters["total_clust.y"] <- list(NULL) ## Delete column

  # cvs_name <- paste0("DETAIL_", marker_lst, "_" ,count_mtx_type, suffix, "_", clust_types ,"_v3.csv")    
  # write.csv(totals_grp_clusters, here(cvsDir, "cvs_files_markers", cvs_name))

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

  # cvs_name <- paste0('SUMMARY_', marker_lst, "_" ,count_mtx_type, suffix, "_", clust_types,"_v3.csv")
  # write.csv(df_summary, here(cvsDir, "cvs_files_markers", cvs_name), row.names=TRUE)

  ## bind rows to the main table
  df_summary <- data.frame(seurat_clusters = c("Total.Cells.Sample", "Total.Cells.ALL", "Percentage.Sample"), df_summary)
  new_columns <- list(total_clust.x = NA, Perc.Cluster = NA, cell.types = NA)
  df_summary <- df_summary |> mutate(!!!new_columns)
  df_summary <- rbind(totals_grp_clusters, df_summary)

  ## add suffix to read the corresponding file
  if (cellranger_pipe == "CR_arc_reanalyze_outliers_ATAC") {
    cvs_name <- paste0('FULL_SUMMARY_', marker_lst, "_" ,count_mtx_type, suffix, "_", clust_types,"_ATAC_v3.csv")
  } else {
    cvs_name <- paste0('FULL_SUMMARY_', marker_lst, "_" ,count_mtx_type, suffix, "_", clust_types,"_v3.csv")
  }
  # cvs_name <- paste0('FULL_SUMMARY_', marker_lst, "_" ,count_mtx_type, suffix, "_", clust_types,"_v3.csv")
  write.csv(df_summary, here(cvsDir, cvs_name), row.names=FALSE)

}

message("Done!")


# # slurm script reproducibility
# library("slurmjobs")
# slurmjobs::job_loop(
#   loops = list(cellranger_pipe = c("CR_crossBarcodes", "CR_complementBarcodes", "CR_arc_reanalyze", "CR_arc_reanalyze_outliers"), 
#                marker_lst = c("literature_base", "data_drive")),
#   name = "02_cell_types_percentages_v3",
#   cores = 2,
#   create_shell = TRUE
# )


