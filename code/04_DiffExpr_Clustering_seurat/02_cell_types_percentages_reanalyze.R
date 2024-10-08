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

cvsDir <- here("processed-data", "04_DiffExpr_Clustering_seurat", "cellrangerARC_reanalyze")
plotDir <- here("plots", "04_DiffExpr_Clustering_seurat", "cellrangerARC_reanalyze")

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

###### Load clusters by sample from processed_data directory ######

## Format sample ID names and reorder
df_mdT <- read.csv(input_csv, row.names = 1)
df_mdT$orig.ident <- gsub("^([0-9])S", "S\\1", df_mdT$orig.ident)
df_mdT$orig.ident <- gsub("^(S[0-9]+_Hb)", "\\1_r", df_mdT$orig.ident)
df_mdT$orig.ident <- str_extract(df_mdT$orig.ident, "^(S[0-9]+_Hb_r)")      
#table(df_mdT)
ref_sort <- c("S3_Hb_r", "S4_Hb_r", "S5_Hb_r", "S6_Hb_r", "S7_Hb_r", "S8_Hb_r", "S9_Hb_r", "S10_Hb_r", "S11_Hb_r", "S12_Hb_r")
df_mdT[order(sapply(df_mdT$orig.ident, function(x) which(x == ref_sort))), ]
#unique(df_mdT$orig.ident)

### Preparing list of clusters to summarize 
## Add all cell-types to the end of the list of lists
allT <- c(unique(df_mdT[["seurat_clusters"]]))
lst_clust <- list(hb = hb, neu = neu, thal= thal, allTypes = allT) 
names(lst_clust)

i<-1
lst_plt <- list()

## Parse total cells and percentages from 3 cell types
for (x in lst_clust) {
  ## for testing literature base: x <- c(0,3,4,5,9,17,22,24)
  ## for testing data-driven base: x <- c(0,1,3,4,5,8,9,10,11,14,16,21,24,2,26) 
  clust <- x
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
    bind_rows(df_mdT |>   # bid row with total cells by sample
                dplyr::filter(seurat_clusters %in% clust) |>
                group_by(orig.ident) |> 
                summarise(total_sample = sum(N)) |>
                mutate(seurat_clusters = "Total.Cells.Sample") |>
                pivot_wider(id_cols = seurat_clusters, names_from = orig.ident, values_from = total_sample) |>
                mutate(total_clust = rowSums(across(where(is.numeric))))) # sum total cells by sample
  #table(is.na.data.frame(grp_clusters))
  grp_clusters[is.na(grp_clusters)] <- 0
  
  ## save cluster information. g.e: SUMMARY_literature_base_norm_counts_Harmony_All_cluster_hb.csv
  cvs_name <- paste0('SUMMARY_', marker_lst, "_" ,count_mtx_type, suffix, "_", clust_types ,".csv")
  # Reorder column by sample name
  grp_clusters <- grp_clusters[,  c("seurat_clusters", ref_sort, "total_clust")]
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
  message("Adding total cells for ", clust_types," cluster group.")
  lst_new <-list(tbl_total_cells_grp[2])
  names(lst_new[i]) <- names(lst_clust[i])
  lst_plt <- append(lst_plt, lst_new)
  
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
    mutate(Perc.Cluster = sum(c(S3_Hb_r,S4_Hb_r,S5_Hb_r,S6_Hb_r,S7_Hb_r,S8_Hb_r,S9_Hb_r,S10_Hb_r,S11_Hb_r,S12_Hb_r) * 100 / total_cells))
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

  ## Load markers used to label clusters and collapse names in `cell.type` description 
  # seurat.combined.data_counts_Harmony_cell_types_all_gm20.csv
  df_cluster_names <- read.csv(here(cvsDir, "cvs_files_markers", paste0(Seurat_base_name, suffix_clust_names)))
  df_cluster_names_long <- df_cluster_names[c("cluster", "cell.type")] |>
    group_by(cluster) |> summarise(cell.types = paste(cell.type, collapse=",")) |>
    rename(seurat_clusters = cluster) |> dplyr::filter(seurat_clusters %in% clust)
  df_summary_Hb_long <- merge(df_summary_Hb, df_cluster_names_long, by = "seurat_clusters", all.x = TRUE, sort = FALSE)
  cvs_name_long <- paste0('SUMMARY_', marker_lst, "_" ,count_mtx_type, suffix, "_", clust_types ,"_long_CT.csv")
  write.csv(df_summary_Hb_long, here(cvsDir, cvs_name_long))
  ## Second version drop duplicate cell.types
  df_cluster_names <- df_cluster_names[c("cluster", "cell.type")] |>
    group_by(cluster) |> distinct() |> summarise(cell.types = paste(cell.type, collapse=",")) |>
    rename(seurat_clusters = cluster) |> distinct() |> dplyr::filter(seurat_clusters %in% clust)
  df_summary_Hb <- merge(df_summary_Hb, df_cluster_names, by = "seurat_clusters", all.x = TRUE, sort = FALSE)
  write.csv(df_summary_Hb, here(cvsDir, cvs_name))
  #print(df_summary_Hb)  
}

# print(lst_plt)
# plot(unlist(lst_plt[10]))
# message(" List for plot")
# plt_total_cells <- ggplot(lst_plt, aes(x=sample, y=counts_grp, group=1)) +
#   geom_line() + geom_point() + ggtitle(paste("Cell-type for" , clust_types , "group")) +
#   theme(axis.text.x = element_text())  # csc. it reorder the samples.to fix (angle = 45, hjust=0.3)
# ggsave(here(cvsDir, paste0("plt", suffix_clust_plt,"_", clust_types, ".png")))


message("Done!")


## slurm script reproducibility

# slurmjobs::job_loop(
#   loops = list(marker_lst = c("literature_base", "data_drive")),
#   name = "02_cell_types_percentages_reanalyze",
#   cores = 2,
#   create_shell = TRUE
# )


