###############################################################################
##
##  Calculate and summarize percentage statistics of pre-selected cell types in different categories
## 
###############################################################################

library(tidyverse)
library(scales)
library(here)

options(digits=2)

## Set count-mtx type and integration model (CCA or Harmony)

cvsDir <- here("processed-data", "04_DiffExpr_Clustering_seurat")

count_mtx_type <- 'data_counts'
#count_mtx_type <- 'norm_counts' 
#Seurat_reduction <- 'CCA'
Seurat_reduction <- 'Harmony'

if (count_mtx_type=='data_counts') { Seurat_base_name <- 'seurat.combined.data_counts' } else { Seurat_base_name <- 'seurat.combined.norm_counts' }
## Compose Seurat object name processed before
if (Seurat_reduction=='CCA') {
  suffix <- '_CCA_subset_cluster'
  suffix_clust_names <- '_CCA_subset_cell_types_all_gm20.csv'
  input_csv <- here(cvsDir, "cvs_files_markers", "pilot_results", paste0(Seurat_base_name, suffix))
} else {
  suffix <- '_Harmony_subset_cluster'
  input_csv <- here(cvsDir, "cvs_files_markers", "pilot_results", paste0(Seurat_base_name, suffix))
  suffix_clust_names <- '_Harmony_subset_cell_types_all_gm20.csv'
}
input_csv <- paste0(input_csv, "_info.csv")
basename(input_csv)


## Load clusters by sample from processed_data directory
df_mdT <- read.csv(input_csv, row.names = 1)
head(df_mdT)

### Preparing list of clusters to summarize 
allT <- c(unique(df_mdT[["seurat_clusters"]]))
hb <- c(6, 12, 14, 15) ## Habenula clusters
neu <- c(6, 9, 10, 12, 14, 15) ## Neuron clusters
lst_clust <- list(hb = hb, neu = neu, allTypes = allT)
names(lst_clust)

i<-1

for (x in lst_clust) {
  ## for testing: x <-  c(0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15)
  clust <- x
  clust_types <- names(lst_clust[i])
  message(" Summarizing ", clust_types," clusters. ")
  message(length(clust), " Clusters to summarize")
  print(clust)
  i<-i+1
  
  ## Filtering clusters of interest
  hb_clusters <- df_mdT |>
    dplyr::arrange(seurat_clusters) |> 
    dplyr::filter(seurat_clusters %in% clust) |>
    pivot_wider(names_from = orig.ident, values_from = N) |>  # make wider format the table
    left_join(df_mdT |>
                group_by(seurat_clusters) |> 
                summarise(total_clust = sum(N))) |> # sum total cells by cluster
    mutate(seurat_clusters = as.character(seurat_clusters)) |>
    bind_rows(df_mdT |>   # bid row with total cells by sample
                dplyr::filter(seurat_clusters %in% clust) |>
                group_by(orig.ident) |> 
                summarise(total_sample = sum(N)) |>
                #mutate(seurat_clusters = "Hb_total_sample") |>
                mutate(seurat_clusters = "Total.Cells.Sample") |>
                pivot_wider(id_cols = seurat_clusters, names_from = orig.ident, values_from = total_sample) |>
                mutate(total_clust = rowSums(across(where(is.numeric))))) # sum total cells by sample
  
  #hb_clusters <- hb_clusters |> rename(S4_Hb_KDM = `4S_Hb_KDM`, S5_Hb_KDM = `5S_Hb_KDM`, S6_Hb_KDM = `6S_Hb_KDM`)
  #print(hb_clusters)

  cvs_name <- paste0('SUMMARY_', count_mtx_type, suffix, "_", clust_types ,"_pilot.csv")
  cvs_file <- here(cvsDir, "pilot_results", cvs_name)
  write.csv(hb_clusters, cvs_file)  

  
  ######## Load summary csv rpt in wider format
  df_summary_Hb <- read.csv(cvs_file, row.names = 1, check.names=FALSE)

  #### Calculate perceptual values by cluster and sample 
  
  ## Extract total cells by specific cell-types by sample
  total_cells_by_Hb_sample <- tail(df_summary_Hb, 1) |>
    pivot_longer(!seurat_clusters, names_to = "sample", values_to = "countsHb") |>
    select(c(sample, countsHb))
  
  ## Calculate Total cells by sample
  total_cellsS <- df_mdT |>   
    group_by(orig.ident) |> 
    summarise(total_sample = sum(N)) |>
    mutate(seurat_clusters = "Total.Cells.ALL") |> #"Total_sample"
    pivot_wider(id_cols = seurat_clusters, names_from = orig.ident, values_from = total_sample) |>
    mutate(total_clust = rowSums(across(where(is.numeric)))) #|> 
    #rename(S4_Hb_KDM = `4S_Hb_KDM`, S5_Hb_KDM = `5S_Hb_KDM`, S6_Hb_KDM = `6S_Hb_KDM`)
  df_summary_Hb <- bind_rows(df_summary_Hb, total_cellsS)
  
  ## Calculate percent cells by sample
  total_cells <- sum(df_mdT$N)
  df_summary_Hb <- df_summary_Hb |> rowwise() |> 
    mutate(Perc.Cluster = sum(c(S1_Hb_KDM, S2_Hb_KDM) * 100 / total_cells))

  ## Extract total cells by by sample
  total_cells_by_sample <- total_cellsS |> 
    pivot_longer(!seurat_clusters, names_to = "sample", values_to = "countsT") |>
    select(c(countsT))
  
  df_tmp <- bind_cols(total_cells_by_Hb_sample, total_cells_by_sample) 
  df_tmp <- df_tmp |> rowwise() |> mutate(Percentage_sample = round((countsHb*100 / countsT), digits = 2))
  df_tmp <- df_tmp[c("sample", "Percentage_sample")] |> pivot_wider(names_from = sample, values_from = Percentage_sample)
  df_tmp <- cbind(seurat_clusters = "Perc.Sample", df_tmp, "Perc.Cluster" = 0) 
  df_summary_Hb <- bind_rows(df_summary_Hb, df_tmp)
  
  ## Load markers used to label clusters and collapse names in cluster cell.type description
  df_cluster_names <- read.csv(here(cvsDir, "cvs_files_markers", "pilot_results", paste0(Seurat_base_name, suffix_clust_names)))
  df_cluster_names <- df_cluster_names[c("cluster", "cell.type")] |>
    group_by(cluster) |> summarise(cell.types = paste(cell.type, collapse=",")) |>
    rename(seurat_clusters = cluster) |> dplyr::filter(seurat_clusters %in% clust)
  df_summary_Hb <- merge(df_summary_Hb, df_cluster_names, by = "seurat_clusters", all.x = TRUE, sort = FALSE)
  print(df_summary_Hb)
  
  write.csv(df_summary_Hb, cvs_file)

}

