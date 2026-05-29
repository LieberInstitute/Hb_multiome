#Using deconvobuddies to get Mean Ratio and 1vs all cell-type markers


library("SingleCellExperiment")
library("dplyr")
library("ggplot2")
library("DeconvoBuddies")
library(qs2)
library(here)

here::here()


#Path to save any generated data
new_data_path = here('processed-data', '05_03_annotation_adjustments', '14_deconvoBuddies_markers')
#Path to plot directory
plot_path = here('plots','05_03_annotation_adjustments', '14_deconvoBuddies_markers')

if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)


multiome_path = here('processed-data', '05_03_annotation_adjustments','06_refined_annotations')
multiome_sce = qs_read(paste0(multiome_path, '/refined_annotation_multiomeHab_SCE.qs2'))

#Drop ATAC assay
altExps(multiome_sce) <- NULL
gc()


#MeanRatio markers
marker_stats_MeanRatio <- get_mean_ratio(
    sce = multiome_sce, # sce is the SingleCellExperiment with our data
    assay_name = "logcounts", ## assay to use, we recommend logcounts [default]
    cellType_col = "refined_mid_cluster", # column in colData with cell type info
)


#1vsAll markers
marker_stats_1vAll <- findMarkers_1vAll(
     sce = multiome_sce, # sce is the SingleCellExperiment with our data
     assay_name = "logcounts",
     cellType_col = "refined_mid_cluster", # column in colData with cell type info
     mod = "~orig.ident" # Control for donor stored in "orig.ident" with mod
 )



## join the two marker_stats tables
marker_stats <- marker_stats_MeanRatio |>
    left_join(marker_stats_1vAll, by = join_by(gene, cellType.target))


#Save results

saveRDS(marker_stats_MeanRatio, file = paste0(new_data_path, '/marker_stats_MeanRatio.rds'))
saveRDS(marker_stats_1vAll , file = paste0(new_data_path, '/marker_stats_1vAll.rds'))
saveRDS(marker_stats, file = paste0(new_data_path, '/marker_stats_combo.rds'))



## Check stats for our top genes
#marker_stats |>
#    filter(MeanRatio.rank == 1) |>
#    select(gene, cellType.target, MeanRatio, MeanRatio.rank, std.logFC, std.logFC.rank)

# create hockey stick plots to compare MeanRatio and standard logFC values.
hockey_plot = marker_stats |>
    ggplot(aes(MeanRatio, std.logFC)) +
    geom_point() +
    facet_wrap(~cellType.target, scales = 'free') + theme_bw()


ggsave(plot = hockey_plot, filename = 'hab_hockey_plots.pdf', path = plot_path, width = 12, height = 8, device = 'pdf')


#And plot top 10 markers
for(celltype in unique(multiome_sce$refined_mid_cluster)){
  p_markers = plot_marker_express(
      sce = multiome_sce,
      stats = marker_stats,
      cell_type = celltype,
      n_genes = 10,
      cellType_col = "refined_mid_cluster"
  )

  ggsave(plot = p_markers, filename = paste0('hab_top10_markers_', celltype, '.pdf'), path = plot_path, width = 12, height = 8, device = 'pdf')

}

