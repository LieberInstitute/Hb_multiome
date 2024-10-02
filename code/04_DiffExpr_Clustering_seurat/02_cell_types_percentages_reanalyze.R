###############################################################################
##
##  Calculate and summarize percentage statistics of pre-selected cell types in different categories
## 
###############################################################################

library("tidyverse")
library("scales")
library("here")
library("ggplot2")

options(digits=2)

cvsDir <- here("processed-data", "04_DiffExpr_Clustering_seurat", "cellrangerARC_reanalyze")
if (!dir.exists(cvsDir)) { dir.create(cvsDir) }

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
  neu <- c(0,1,3,4,5,8,9,10,11,14,16,21,24,2,26) 
  hb <- c(4,5,9,10,14,16,2)
  thal <- c(0,1,3,11,21,24,26)
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

## Load clusters by sample from processed_data directory
df_mdT <- read.csv(input_csv, row.names = 1)
df_mdT$orig.ident <- gsub("^([0-9])S", "S\\1", df_mdT$orig.ident)
#text <- gsub("^([S0-9*]_Hb)", "\\1_r", df_mdT$orig.ident)
table(df_mdT)
#str_remove(unlist(plt_grp$sample), "_KDM_reanalysis")

## rename sample IDs starting with numeric character 
ref_sort <- c("S3_Hb_KDM_reanalysis", "S4_Hb_KDM_reanalysis", "S5_Hb_KDM_reanalysis",
  "S6_Hb_KDM_reanalysis", "S7_Hb_KDM_reanalysis", "S8_Hb_KDM_reanalysis",
  "S9_Hb_KDM_reanalysis", "S10_Hb_KDM_reanalysis", "S11_Hb_KDM_reanalysis", "S12_Hb_KDM_reanalysis")
df_mdT[order(sapply(df_mdT$orig.ident, function(x) which(x == ref_sort))), ]
#unique(df_mdT$orig.ident)

### Preparing list of clusters to summarize 
## Add all cell-types to the end of the list of lists
allT <- c(unique(df_mdT[["seurat_clusters"]]))
lst_clust <- list(hb = hb, neu = neu, thal= thal, allTypes = allT) 
names(lst_clust)

i<-1

## Parse total cells and percentages from 3 cell types
for (x in lst_clust) {
  ## for testing literature base: x <- c(0,3,4,5,9,17,22,24)
  ## for testing data-driven base: x <- c(0,1,3,4,5,8,9,10,11,14,16,21,24,2,26) 
  clust <- x
  clust_types <- names(lst_clust[i])
  message(length(clust), " Clusters to summarize")
  i<-i+1
  
  ## Filtering clusters of interest in the given samples
  hb_clusters <- df_mdT |>
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
                mutate(seurat_clusters = "Total.Cells.Sample") |>
                pivot_wider(id_cols = seurat_clusters, names_from = orig.ident, values_from = total_sample) |>
                mutate(total_clust = rowSums(across(where(is.numeric))))) # sum total cells by sample
  #table(is.na.data.frame(hb_clusters))
  hb_clusters[is.na(hb_clusters)] <- 0
  
  ## save cluster information. g.e: SUMMARY_literature_base_norm_counts_Harmony_All_cluster_hb.csv
  cvs_name <- paste0('SUMMARY_', marker_lst, "_" ,count_mtx_type, suffix, "_", clust_types ,".csv")
  # Reorder column by sample name
  hb_clusters <- hb_clusters[,  c("seurat_clusters", ref_sort, "total_clust")]
  write.csv(hb_clusters, here(cvsDir, cvs_name))  

  #### Calculate perceptual values by cluster and sample 
  df_summary_Hb <- read.csv(here(cvsDir, cvs_name), row.names = 1, check.names=FALSE)
  ## Extract total cells by specific cell-types by sample
  total_cells_by_cell_type_grp <- tail(df_summary_Hb, 1) |>
    pivot_longer(!seurat_clusters, names_to = "sample", values_to = "counts_grp") |>
    select(c(sample, counts_grp))
  
  ## plot total cells for the specific group of clusters
  n <- dim(total_cells_by_cell_type_grp)[1]-1
  plt_grp <- total_cells_by_cell_type_grp[1:n,]
  plt_grp$sample <- str_remove(unlist(plt_grp$sample), "_KDM_reanalysis")
  plt1 <- ggplot(plt_grp, aes(x=sample, y=counts_grp, group=1)) +
    geom_line() + geom_point() + ggtitle(paste("Cell-type for" , clust_types , "group")) + 
    theme(axis.text.x = element_text())  # csc. it reorder the samples.to fix (angle = 45, hjust=0.3)
  ggsave(here(cvsDir, paste0("plt", suffix_clust_plt,"_", clust_types, ".png")))
  
  ## Calculate Total cells by sample
  total_cellsS <- df_mdT |>   
    group_by(orig.ident) |> 
    summarise(total_sample = sum(N)) |>
    mutate(seurat_clusters = "Total.Cells.ALL") |> #"Total_sample"
    pivot_wider(id_cols = seurat_clusters, names_from = orig.ident, values_from = total_sample) |>
    mutate(total_clust = rowSums(across(where(is.numeric)))) 
  df_summary_Hb <- bind_rows(df_summary_Hb, total_cellsS)
  
  ## Calculate percent cells by sample
  total_cells <- sum(df_mdT$N)
  #sum_col <- noquote(paste(ref_sort, sep="", collapse=","))
  df_summary_Hb <- df_summary_Hb |> rowwise() |> 
    mutate(Perc.Cluster = sum(c(S3_Hb_KDM_reanalysis,S4_Hb_KDM_reanalysis,S5_Hb_KDM_reanalysis,S6_Hb_KDM_reanalysis,S7_Hb_KDM_reanalysis,S8_Hb_KDM_reanalysis,S9_Hb_KDM_reanalysis,S10_Hb_KDM_reanalysis,S11_Hb_KDM_reanalysis,S12_Hb_KDM_reanalysis) * 100 / total_cells))
    #mutate(Perc.Cluster = sum(c(S4_Hb_KDM_reanalysis, S5_Hb_KDM_reanalysis, S6_Hb_KDM_reanalysis) * 100 / total_cells))

  ## Extract total cells by sample
  #is.na.data.frame(total_cellsS)
  total_cells_by_sample <- total_cellsS |> 
    pivot_longer(!seurat_clusters, names_to = "sample", values_to = "countsT") |>
    select(c(countsT))
  
  df_tmp <- bind_cols(total_cells_by_cell_type_grp, total_cells_by_sample) 
  df_tmp <- df_tmp |> rowwise() |> mutate(Percentage_sample = round((counts_grp*100 / countsT), digits = 2))
  df_tmp <- df_tmp[c("sample", "Percentage_sample")] |> pivot_wider(names_from = sample, values_from = Percentage_sample)
  df_tmp <- cbind(seurat_clusters = "Perc.Sample", df_tmp, "Perc.Cluster" = 0) 
  df_summary_Hb <- bind_rows(df_summary_Hb, df_tmp)
  
  ## Load markers used to label clusters and collapse names in cluster cell.type description
  # seurat.combined.data_counts_Harmony_cell_types_all_gm20.csv
  df_cluster_names <- read.csv(here(cvsDir, "cvs_files_markers", paste0(Seurat_base_name, suffix_clust_names)))
  df_cluster_names <- df_cluster_names[c("cluster", "cell.type")] |>
    group_by(cluster) |> summarise(cell.types = paste(cell.type, collapse=",")) |>
    rename(seurat_clusters = cluster) |> dplyr::filter(seurat_clusters %in% clust)
  df_summary_Hb <- merge(df_summary_Hb, df_cluster_names, by = "seurat_clusters", all.x = TRUE, sort = FALSE)
  print(df_summary_Hb)
  
  write.csv(df_summary_Hb, here(cvsDir, cvs_name))

}


message("Done!")


## slurm script reproducibility

# slurmjobs::job_loop(
#   loops = list(marker_lst = c("literature_base", "data_drive")),
#   name = "02_cell_types_percentages_reanalyze",
#   cores = 2,
#   create_shell = TRUE
# )


