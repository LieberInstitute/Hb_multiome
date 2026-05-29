#Checking out the marker geners from MetaMarkers and the MeanRatio markers
#Get GO enrichment results for these marker sets too


library(ggplot2)
library(MetaMarkers)
library(dplyr)
library(here)
library(MetaNeighbor)
library(org.Hs.eg.db)
library(GO.db)
library(qs2)

here::here()


#Path to save any generated data
new_data_path = here('processed-data', '05_03_annotation_adjustments', '15_marker_GO_enrich')
#Path to plot directory
plot_path = here('plots','05_03_annotation_adjustments', '15_marker_GO_enrich')

if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)


deconvo_marker_path = here('processed-data','05_03_annotation_adjustments','14_deconvoBuddies_markers')
metamarker_path = here('processed-data', '05_03_annotation_adjustments', '08_metamarkers')
go_path = here('processed-data','05_03_annotation_adjustments', '10_celltype_func_annot')

marker_stats_MeanRatio = readRDS(file = paste0(deconvo_marker_path, '/marker_stats_MeanRatio.rds'))
marker_stats_1vAll = readRDS(file = paste0(deconvo_marker_path, '/marker_stats_1vAll.rds'))
marker_stats = readRDS(file = paste0(deconvo_marker_path, '/marker_stats_combo.rds'))


multiome_mid_metaMarkers = read_meta_markers(paste0(metamarker_path, '/multiome_refined_mid_meta_markers.csv.gz'))

multiome_path = here('processed-data', '05_03_annotation_adjustments','06_refined_annotations')
multiome_sce = qs_read(paste0(multiome_path, '/refined_annotation_multiomeHab_SCE.qs2'))

gene_dataset_universe = rownames(multiome_sce)


#Set up the GO terms. Filter for the genes present in the dataset and then filter gene sets down to max 100 genes
go_sets = readRDS(file.path(go_path, "go_human.rds"))

go_sets = lapply(go_sets, function(gene_set) {gene_set[gene_set %in% gene_dataset_universe]})
min_size = 10
max_size = 100
go_set_size = sapply(go_sets, length)
go_sets = go_sets[go_set_size >= min_size & go_set_size <= max_size]
length(go_sets)


top_50_metamarkers = multiome_mid_metaMarkers |> filter(rank <= 100) |> select(cell_type, gene) |>
  mutate(type = 'metamarker')

top_50_meanratio = marker_stats_MeanRatio |> filter(MeanRatio.rank <= 100) |> select(cellType.target, gene) |>
  mutate(type = 'meanRatio')
colnames(top_50_meanratio) = c('cell_type', 'gene', 'type')

top_50_1vall = marker_stats_1vAll |> filter(std.logFC.rank <= 100) |> select(cellType.target, gene) |>
  mutate(type = '1vall')
colnames(top_50_1vall) = c('cell_type', 'gene', 'type')


all_top_markers = rbind(top_50_metamarkers, top_50_meanratio, top_50_1vall)



overlap_pct <- all_top_markers |>
  group_by(cell_type) |>
  summarize(
    # Get unique genes per marker type
    metamarker_genes = list(unique(gene[type == "metamarker"])),
    meanratio_genes = list(unique(gene[type == "meanRatio"])),
    onevall_genes = list(unique(gene[type == "1vall"])),  
    
    # Calculate pairwise overlaps
    overlap_meta_mr = length(intersect(metamarker_genes[[1]], meanratio_genes[[1]])),
    overlap_meta_ova = length(intersect(metamarker_genes[[1]], onevall_genes[[1]])),
    overlap_mr_ova = length(intersect(meanratio_genes[[1]], onevall_genes[[1]])),
    
    # Convert to percentage of max possible
    n_meta = length(metamarker_genes[[1]]),
    n_mr = length(meanratio_genes[[1]]),
    n_ova = length(onevall_genes[[1]]),
    
    pct_overlap_meta_mr = (overlap_meta_mr / pmin(n_meta, n_mr)) * 100,
    pct_overlap_meta_ova = (overlap_meta_ova / pmin(n_meta, n_ova)) * 100,
    pct_overlap_mr_ova = (overlap_mr_ova / pmin(n_mr, n_ova)) * 100,
    
    .groups = "drop"
  ) |>
  select(cell_type, starts_with("pct_overlap"))

overlap_pct




test_1 = all_top_markers |> filter(cell_type == 'Astrocyte' & type == 'metamarker') |> pull(gene)
test_2 = all_top_markers |> filter(cell_type == 'Astrocyte' & type == 'meanRatio') |> pull(gene)
test_3 = all_top_markers |> filter(cell_type == 'Astrocyte' & type == '1vall') |> pull(gene)


length(intersect(test_1, test_2)) / length(test_1)



#Test out Go enrichment for the intersect of the genes in the top 100 for the MeanRatio and the 1vall

conf_markers = all_top_markers |> group_by(cell_type, gene) |> filter(type %in% c('meanRatio','1vall')) |>
  summarise(n = n()) |> filter(n == 2)


#conf_markers = all_top_markers |> filter(type == 'metamarker') 

conf_markers = all_top_markers |> filter(type == 'meanRatio') 

#conf_markers = all_top_markers |> filter(type == '1vall') 


enrichment_results = list()
cell_types = unique(conf_markers$cell_type)

for(i in 1:length(cell_types)){

  test_celltype_markers = conf_markers |> filter(cell_type == cell_types[i]) |> pull(gene)
  celltype_df <- data.frame(
    GO_term = names(go_sets),
    overlap = sapply(go_sets, \(x) length(intersect(test_celltype_markers, x))),
    go_length = sapply(go_sets, length),
    cell_type = cell_types[i]
  ) |>
    mutate(
      celltype_length = length(test_celltype_markers),
      universe_term = length(gene_dataset_universe) - go_length,
      p_value = phyper(overlap - 1, go_length, universe_term, celltype_length, lower.tail = FALSE)
    ) |>
    dplyr::select(GO_term, overlap, p_value, cell_type)

    enrichment_results[[i]] = celltype_df

}


enrichment_df = do.call('rbind', enrichment_results)
  
enrichment_df$fdr_adjusted_p = p.adjust(enrichment_df$p_value, method = 'BH')

enrichment_df |> filter(fdr_adjusted_p < 0.05) |> View() 



set_of_interest = "GO:0097061|dendritic spine organization|BP"
plotDotPlot(dat = multiome_sce,
experiment_labels = multiome_sce$orig.ident,
celltype_labels = multiome_sce$refined_mid_cluster,
gene_set = go_sets[[set_of_interest]]) + ggtitle(set_of_interest)
