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

cvsDir <- here("processed-data", "04_DiffExpr_Clustering_seurat", "cellranger_count")
plotDir <- here("plots", "04_DiffExpr_Clustering_seurat", "cellranger_count")

if (!dir.exists(cvsDir)) { dir.create(cvsDir) }
if (!dir.exists(plotDir)) { dir.create(plotDir) }

## read input arguments
args = commandArgs(trailingOnly=TRUE)
marker_lst <- args[2]

## Selected manually the clusters based on the cell-type identification gene-marker lists
## args opt for testing: marker_lst="literature_base" 
##                       marker_lst="data_driven" 
if (marker_lst=="literature_base") {
  neu <- c(0,1,2,5,7,8,10,13,14,17,25,26,29) 
  hb <- c(0,2,5,7,10,13,14) 
  thal <- c(1,8,26,29) # all clusters are mixed with other type of neurons
} else if (marker_lst=="data_driven") {
  neu <- c(0,1,2,5,7,8,9,10,13,14,15,16,26,27,28,29)
  hb <- c(2,5,7,9,10,13,14,16,27)
  thal <- c(0,1,8,15,26,28,29)
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
basename(input_csv)

message("Calculating percentages in progress .../n ")

###### Load clusters by sample from processed_data directory ######

## Format sample ID names and reorder
df_mdT <- read.csv(input_csv, row.names = 1)
df_mdT$orig.ident <- sprintf("S%02d_Hb_gex", extract_numeric(df_mdT$orig.ident))
ref_sort <- sort(unique(df_mdT$orig.ident))
df_mdT <- df_mdT[order(sapply(df_mdT$orig.ident, function(x) which(x == ref_sort))), ]
print(unique(df_mdT$orig.ident))

# df_mdT$orig.ident <- gsub("^([0-9]+)C", "S\\1", df_mdT$orig.ident)
# df_mdT$orig.ident <- gsub("^(S[0-9]+_Hb)", "\\1_gex", df_mdT$orig.ident)
# df_mdT$orig.ident <- str_extract(df_mdT$orig.ident, "^(S[0-9]+_Hb_gex)")      
# ref_sort <- c("S3_Hb_gex", "S4_Hb_gex", "S5_Hb_gex", "S6_Hb_gex", "S7_Hb_gex", "S8_Hb_gex", "S9_Hb_gex", "S10_Hb_gex", "S11_Hb_gex", "S12_Hb_gex")
# df_mdT[order(sapply(df_mdT$orig.ident, function(x) which(x == ref_sort))), ]
#unique(df_mdT$orig.ident)

### Preparing list of clusters to summarize 
## Add all cell-types to the end of the list of lists
allT <- c(unique(df_mdT[["seurat_clusters"]]))
lst_clust <- list(hb = hb, neu = neu, thal= thal, allTypes = allT) 
names(lst_clust)

i<-1
## Parse total cells and percentages from 3 cell types
for (clust in lst_clust) {
  ## for testing literature base: clust <- lst_clust[1]

  clust_types <- names(lst_clust[i])
  message(length(clust), " Clusters to summarize")
  i<-i+1
  
  ## Filtering clusters of interest in the given samples
  grp_clusters <- df_mdT |>
    dplyr::filter(seurat_clusters %in% clust) |>
    pivot_wider(names_from = orig.ident, values_from = N) |>  # make wider format the table
    left_join(df_mdT |>
                group_by(seurat_clusters) |> 
                summarise(total_clust = sum(N))) |> # sum total cells by cluster
    mutate(seurat_clusters = as.character(seurat_clusters)) |>
    bind_rows(df_mdT |>   # bind row with total cells by sample
                dplyr::filter(seurat_clusters %in% clust) |>
                group_by(orig.ident) |> 
                summarise(total_sample = sum(N)) |>
                mutate(seurat_clusters = "Total.Cells.Sample") |>
                pivot_wider(id_cols = seurat_clusters, names_from = orig.ident, values_from = total_sample) |>
                mutate(total_clust = rowSums(across(where(is.numeric))))) # sum total cells by sample
  #table(is.na.data.frame(grp_clusters))
  grp_clusters[is.na(grp_clusters)] <- 0
  
  ## save cluster information. g.e: SUMMARY_literature_base_norm_counts_Harmony_All_cluster_hb_v2.csv
  cvs_name <- paste0('SUMMARY_', marker_lst, "_" ,count_mtx_type, suffix, "_", clust_types ,"_v2.csv")
  write.csv(grp_clusters, here(cvsDir, cvs_name))  

  #### Calculate perceptual values by cluster and sample 
  df_summary_Hb <- read.csv(here(cvsDir, cvs_name), row.names = 1, check.names=FALSE)
  ## Extract total cells by specific cell-types by sample
  total_cells_by_cell_type_grp <- tail(df_summary_Hb, 1) |>
    pivot_longer(!seurat_clusters, names_to = "sample", values_to = "counts_grp") |>
    select(c(sample, counts_grp))
  
  ## calculate `total cells by sample by specific group of clusters`
  n <- dim(total_cells_by_cell_type_grp)[1]-1
  tbl_total_cells_grp <- total_cells_by_cell_type_grp[1:n,]
  colnames(tbl_total_cells_grp)[2] <- clust_types

  ## Calculate `total cells` by sample
  total_cellsS <- df_mdT |>   
    group_by(orig.ident) |> 
    summarise(total_sample = sum(N)) |>
    mutate(seurat_clusters = "Total.Cells.ALL") |> #"Total_sample"
    pivot_wider(id_cols = seurat_clusters, names_from = orig.ident, values_from = total_sample) |>
    mutate(total_clust = rowSums(across(where(is.numeric)))) 
  df_summary_Hb <- bind_rows(df_summary_Hb, total_cellsS)
  
  ## Calculate percentages by clusters
  total_cells <- sum(df_mdT$N)
  #sample_ids <- noquote(paste(ref_sort, sep="", collapse=","))
  df_summary_Hb <- df_summary_Hb |> rowwise() |> 
    #mutate(Perc.Cluster = sum(c(ref_sort) * 100 / total_cells))
    mutate(Perc.Cluster = sum(c(S03_Hb_gex,S04_Hb_gex,S05_Hb_gex,S06_Hb_gex,S07_Hb_gex,S08_Hb_gex,S09_Hb_gex,S10_Hb_gex,S11_Hb_gex,S12_Hb_gex) * 100 / total_cells))
  df_summary_Hb$Perc.Cluster <- round(df_summary_Hb$Perc.Cluster, digits = 2)
  ## ## Calculate percentages by sample
  total_cellsS <- total_cellsS |>   
    pivot_longer(!seurat_clusters, names_to = "sample", values_to = "countsT") |>
    select(c(countsT))
  
  vect1 <- df_summary_Hb[df_summary_Hb$seurat_clusters=="Total.Cells.Sample", ][0:length(ref_sort)+2] 
  vect2 <- df_summary_Hb[df_summary_Hb$seurat_clusters=="Total.Cells.ALL", ][0:length(ref_sort)+2] 
  percent.sample <- round((vect1*100 / vect2), digits = 2)
  percent.sample <- append(list(seurat_clusters = "Percentage.Sample"), percent.sample, 1)
  df_summary_Hb <- bind_rows(df_summary_Hb, percent.sample)

  ## Load cluster info and cell-types to collapse names in `cell.type` column (description) 
  df_cluster_names <- read.csv(here(cvsDir, "cvs_files_markers", paste0(Seurat_base_name, suffix_clust_names)))
  df_cluster_namesv2 <- df_cluster_names[c("cluster", "cell.type")] |> # df_cluster_names_long
    group_by(cluster) |> summarise(cell.types = paste(cell.type, collapse=",")) |>
    rename(seurat_clusters = cluster) |> dplyr::filter(seurat_clusters %in% clust)
  
  ## Update `cell.types` column with `IDs+counts by marker` by cluster
  list_ct <- df_cluster_namesv2[["cell.types"]]
  nc <- list()
  idx<-1
  for (ct_cluster in list_ct) {
    x1 <- sapply(ct_cluster, function(x) strsplit(x, ","))
    for (ct in x1) {
      ids <- unique(ct)
      num_rep <- table(ct)
      col_new <- noquote(c(rbind(ids, paste0("(", num_rep, ")"))))
    }
    if (length(col_new)==2) {col_new <- paste0("***", col_new[1], col_new[2])}
    nc <- append(nc, paste0(col_new, collapse = " "))
    idx <- idx+1
  }
  df_cluster_namesv2$cell.types <- noquote(unlist(nc))
  
  ## v2 has gene information
  df_summary_Hb <- merge(df_summary_Hb, df_cluster_namesv2, by = "seurat_clusters", all.x = TRUE, sort = FALSE)
  cvs_name <- paste0('SUMMARY_', marker_lst, "_" ,count_mtx_type, suffix, "_", clust_types,"_v2.csv")
  write.csv(df_summary_Hb, here(cvsDir, cvs_name))
  
}

message("Done!")


## slurm script reproducibility

# slurmjobs::job_loop(
#   loops = list(marker_lst = c("literature_base", "data_drive")),
#   name = "02_cell_types_percentages_CRcount",
#   cores = 2,
#   create_shell = TRUE
# )


