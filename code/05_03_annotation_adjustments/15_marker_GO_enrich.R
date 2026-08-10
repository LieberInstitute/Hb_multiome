#Checking out the marker geners from MetaMarkers and the MeanRatio markers
#Get GO enrichment results for these marker sets too


library(ggplot2)
library(MetaMarkers)
library(SingleCellExperiment)
library(ComplexHeatmap)
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

#MetaNeighbor GO terms
go_aurocs = read.table(here(go_path, 'multiomeHab_functional_aurocs.txt'))
go_aurocs$average = rowMeans(go_aurocs)

marker_stats_MeanRatio = readRDS(file = paste0(deconvo_marker_path, '/marker_stats_MeanRatio.rds'))
marker_stats_1vAll = readRDS(file = paste0(deconvo_marker_path, '/marker_stats_1vAll.rds'))
marker_stats = readRDS(file = paste0(deconvo_marker_path, '/marker_stats_combo.rds'))


multiome_mid_metaMarkers = read_meta_markers(paste0(metamarker_path, '/multiome_refined_mid_meta_markers.csv.gz'))

multiome_path = here('processed-data', '05_03_annotation_adjustments','06_refined_annotations')
multiome_sce = qs_read(paste0(multiome_path, '/refined_annotation_multiomeHab_SCE.qs2'))

#Filter for just the habenula cell-types
hab_celltypes = c('GABA_LHb_C.1','GABA_LHb_C.2','LHb_C','LHb_A','LHb_B','MHb_A','MHb_B','MHb_C','MHb_D', 'Excit.Thal', 'Inhib.Thal')
multiome_nonHab_sce = multiome_sce[ ,!multiome_sce$refined_mid_cluster %in% hab_celltypes]
multiome_sce = multiome_sce[ ,multiome_sce$refined_mid_cluster %in% hab_celltypes]
gc()


gene_dataset_universe = rownames(multiome_sce)


#Set up the GO terms. Filter for the genes present in the dataset and then filter gene sets down to max 100 genes
go_sets = readRDS(file.path(go_path, "go_human.rds"))

go_sets = lapply(go_sets, function(gene_set) {gene_set[gene_set %in% gene_dataset_universe]})
min_size = 10
max_size = 100
go_set_size = sapply(go_sets, length)
go_sets = go_sets[go_set_size >= min_size & go_set_size <= max_size]
length(go_sets)


gs_size = sapply(go_sets, length)
go_aurocs$n_genes = gs_size[rownames(go_aurocs)]



top_50_metamarkers = multiome_mid_metaMarkers |> filter(rank <= 50) |> dplyr::select(cell_type, gene) |>
  mutate(type = 'metamarker')

top_50_meanratio = marker_stats_MeanRatio |> filter(MeanRatio.rank <= 50) |> dplyr::select(cellType.target, gene) |>
  mutate(type = 'meanRatio')
colnames(top_50_meanratio) = c('cell_type', 'gene', 'type')

top_50_1vall = marker_stats_1vAll |> filter(std.logFC.rank <= 50) |> dplyr::select(cellType.target, gene) |>
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
  dplyr::select(cell_type, starts_with("pct_overlap"))

overlap_pct




test_1 = all_top_markers |> filter(cell_type == 'Astrocyte' & type == 'metamarker') |> pull(gene)
test_2 = all_top_markers |> filter(cell_type == 'Astrocyte' & type == 'meanRatio') |> pull(gene)
test_3 = all_top_markers |> filter(cell_type == 'Astrocyte' & type == '1vall') |> pull(gene)


length(intersect(test_1, test_2)) / length(test_1)



#Test out Go enrichment for the intersect of the genes in the top 100 for the MeanRatio and the 1vall

#conf_markers = all_top_markers |> group_by(cell_type, gene) |> filter(type %in% c('meanRatio','1vall')) |>
#  summarise(n = n()) |> filter(n == 2)


#conf_markers = all_top_markers |> filter(type == 'metamarker') 

#conf_markers = all_top_markers |> filter(type == 'meanRatio') 

conf_markers = all_top_markers |> filter(type == '1vall') 


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

#Filter down the terms to those with also high MN GO term scores
index = match(enrichment_df$GO_term, rownames(go_aurocs))
enrichment_df$MN_go_score = go_aurocs$average[index]
enrichment_df$n_go_genes= go_aurocs$n_genes[index]

enrichment_df$alpha_label = ifelse(enrichment_df$n_go_genes <= 30 & enrichment_df$fdr_adjusted_p <= 0.05 & enrichment_df$MN_go_score >= 0.75, 'Gene sets of interest', 'All other gene sets')

gene_set_scatter_p = ggplot(filter(enrichment_df, cell_type %in% hab_celltypes), 
aes(x = MN_go_score, y = -log10(fdr_adjusted_p), alpha = alpha_label, color = alpha_label)) +
  geom_point() +
  scale_alpha_manual(values = c('Gene sets of interest' = 1, 'All other gene sets' = .25), name = 'Gene set size') +
  scale_color_manual(values = c('Gene sets of interest' = 'red', 'All other gene sets' = 'black'), name = 'Gene set size') +
  geom_hline(yintercept = -log10(.05), color = 'red') +
  geom_vline(xintercept = .75, color = 'red') +
  theme_bw() + ggtitle('1vall markers: GO term enrichment and MN AUROC scores') +
  ylab('GO term enrichment for top markers (-log10(FDR-adjusted p-value))') + xlab('Cross-cell type predictability (MetaNeighbor AUROC)')

gene_set_scatter_p 

sig_go_terms = enrichment_df |> filter(fdr_adjusted_p <= .05 & cell_type %in% hab_celltypes & MN_go_score >= .75 & n_go_genes <= 30  ) |> pull(GO_term) |> unique()
length(sig_go_terms)

for(i in 1:length(sig_go_terms)){

  set_of_interest = sig_go_terms[i]
  p = plotDotPlot(dat = multiome_sce,
  experiment_labels = multiome_sce$orig.ident,
  celltype_labels = multiome_sce$refined_mid_cluster,
  gene_set = go_sets[[set_of_interest]]) + ggtitle(set_of_interest)

  print(p)

}

#Inhibitory terms
#"GO:0007214|gamma-aminobutyric acid signaling pathway|BP" 
#"GO:1904862|inhibitory synapse assembly|BP"              
#"GO:0016917|GABA receptor activity|MF"                    
#"GO:0060077|inhibitory synapse|CC" 

unique(unlist(go_sets[sig_go_terms[grepl('GABA', sig_go_terms) | grepl('inhibitory', sig_go_terms) | grepl('gamma',sig_go_terms)]]))
#[1] "ATF4"    "CACNB4"  "GABBR1"  "GABRA1"  "GABRA2"  "GABRA3"  "GABRA4"  "GABRA5"  "GABRA6"  "GABRB1"  "GABRB2"  "GABRB3" 
#[13] "GABRD"   "GABRE"   "GABRG1"  "GABRG2"  "GABRG3"  "GABRR1"  "GABRR2"  "GNAI2"   "HTR1A"   "PLCL1"   "SLC12A2" "GABBR2" 
#[25] "PLCL2"   "PHF24"   "GPR156"  "PIANP"   "GABRR3"  "SHISA7"  "CBLN1"   "FGF13"   "PLXNB1"  "WNT5A"   "CLSTN3"  "SEMA4D" 
#[37] "SRGAP2"  "LGI2"    "NLGN2"   "CLSTN2"  "SEMA4A"  "CBLN4"   "MDGA1"   "NPAS4"   "LHFPL4"  "HAPLN4"  "SRGAP2C" "GABRP"  
#[49] "GRID1"   "GABRQ"   "DTNB"    "GAD1"    "GAD2"    "GLRA1"   "NLGN4Y"  "SYT11"   "NLGN3"   "NLGN4X"  "IGSF9"   "IGSF21" 
#[61] "SLC32A1" "CEP112"  "IQSEC3"


#Excitatory terms
#"GO:0008066|glutamate receptor activity|MF"       
#"GO:0004970|glutamate-gated receptor activity|MF"
unique(unlist(go_sets[sig_go_terms[grepl('glutamate', sig_go_terms)]]))
#[1] "GRIA1"  "GRIA2"  "GRIA3"  "GRIA4"  "GRID1"  "GRID2"  "GRIK1"  "GRIK2"  "GRIK3"  "GRIK4"  "GRIK5"  "GRIN1"  "GRIN2A"
#[14] "GRIN2B" "GRIN2C" "GRIN2D" "GRM1"   "GRM2"   "GRM3"   "GRM4"   "GRM5"   "GRM6"   "GRM7"   "GRM8"   "GRIN3A" "GRIN3B"

#Save scatter plot
ggsave(plot = gene_set_scatter_p, filename = here(plot_path, '1vall_markers_GO_enrichment_MN_AUROC_scatter.pdf'), 
width = 8, height = 3, device = 'pdf')





#And get heatmaps of expression for the top 10 1vall markers for all cell-types

top_10_1vall_neurons = marker_stats_1vAll |> filter(std.logFC.rank <= 10 & cellType.target %in% hab_celltypes) |> 
  dplyr::select(cellType.target, gene) |> 
  mutate(type = '1vall')
colnames(top_10_1vall_neurons) = c('cell_type', 'gene', 'type')

dim(top_10_1vall_neurons)

top_10_1vall_neurons = top_10_1vall_neurons[!duplicated(top_10_1vall_neurons$gene), ]

gene_index = sapply(strsplit(split = '.', rownames(multiome_sce), fixed = TRUE), '[[', 1)  %in% top_10_1vall_neurons$gene
marker_exp = assay(multiome_sce, 'logcounts')[gene_index, ]


cell_annot_matrix <- Matrix::sparse.model.matrix(~ 0 + refined_mid_cluster, data = colData(multiome_sce))
colnames(cell_annot_matrix) <- gsub("refined_mid_cluster", "", colnames(cell_annot_matrix))

cell_counts = colSums(cell_annot_matrix)
pseudobulk_marker_exp = marker_exp %*% cell_annot_matrix

table(colnames(pseudobulk_marker_exp) == names(cell_counts))

avg_marker_exp <- sweep(pseudobulk_marker_exp, 2, cell_counts, "/")

scaled_avg_marker_exp = t(scale(t(avg_marker_exp)))

celltype_order = c('Inhib.Thal','Excit.Thal','GABA_LHb_C.1','GABA_LHb_C.2','LHb_C','LHb_B','LHb_A','MHb_A','MHb_C','MHb_B','MHb_D')
col_order = match(celltype_order, colnames(scaled_avg_marker_exp))

top_10_1vall_neurons$cell_type = factor(top_10_1vall_neurons$cell_type, levels = celltype_order)
gene_order = c(top_10_1vall_neurons[order(top_10_1vall_neurons$cell_type), 'gene' ])
row_order = match(gene_order$gene, sapply(strsplit(split = '.',rownames(scaled_avg_marker_exp), fixed = TRUE), '[[', 1))


scaled_avg_marker_exp = scaled_avg_marker_exp[row_order, col_order]

scaled_avg_marker_exp = t(scaled_avg_marker_exp) 
colnames(scaled_avg_marker_exp) = sapply(strsplit(split = '.', colnames(scaled_avg_marker_exp), fixed = TRUE), '[[', 1)



index = match(colnames(scaled_avg_marker_exp), top_10_1vall_neurons$gene)
column_groups = top_10_1vall_neurons$cell_type[index]

col_func = auroc_cols <- rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
marker_heatmap = Heatmap(scaled_avg_marker_exp, name = 'Scaled average expression', show_row_names = TRUE, show_column_names = TRUE,
        row_title = 'Cell types', column_title = 'Top 10 1vall markers',
        cluster_rows = FALSE, cluster_columns = FALSE, col = col_func,
      height = unit(30, "mm"),
      row_names_gp = gpar(fontsize = 5),      # Row text size
      column_names_gp = gpar(fontsize = 5),
    column_split = column_groups )

marker_heatmap = draw(marker_heatmap)


pdf(here(plot_path, 'top_10_1vall_markers_heatmap.pdf'), width = 10, height = 4)
draw(marker_heatmap)
dev.off()

