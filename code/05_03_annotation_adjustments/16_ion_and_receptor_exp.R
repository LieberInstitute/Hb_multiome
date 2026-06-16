#Looking at the expression of ion channels and various receptors across the habenula cell-types
#These may not necessarily be the most differentially expressed genes, but they are of interest to us


library(SingleCellExperiment)
library(MetaMarkers)
library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(qs2)
library(here)
library(sessioninfo)

plot_path = here('plots','05_03_annotation_adjustments', '16_ion_and_receptor_exp')
if (!dir.exists(plot_path)) dir.create(plot_path)

#colors
source(here('code','05_03_annotation_adjustments','celltype_colors.R'))


message('Loading data...')

multiome_path = here('processed-data', '05_03_annotation_adjustments','06_refined_annotations')
multiome_sce = qs_read(paste0(multiome_path, '/refined_annotation_multiomeHab_SCE.qs2'))

assay(multiome_sce, 'cpm') = MetaMarkers::convert_to_cpm(assay(multiome_sce, 'counts'))

#Drop ATAC assay
altExps(multiome_sce) <- NULL
gc()

#Filter for just the habenula cell-types
hab_celltypes = c('Inhib_LHb_4.1','Inhib_LHb_4.2','LHb.4','LHb.2.7','LHb.1.3.4','MHb.1','MHb.1.2','MHb.2','MHb.3')
multiome_nonHab_sce = multiome_sce[ ,!multiome_sce$refined_mid_cluster %in% hab_celltypes]
multiome_sce = multiome_sce[ ,multiome_sce$refined_mid_cluster %in% hab_celltypes]


go_path = here('processed-data','05_03_annotation_adjustments', '10_celltype_func_annot')
#GO terms
go_sets = readRDS(file.path(go_path, "go_human.rds"))

#Filter for genes present in the dataset
known_genes = rownames(multiome_sce)
go_sets = lapply(go_sets, function(gene_set) {gene_set[gene_set %in% known_genes]})


#Identify the ion channel go terms
go_sets[1:10]


table(grepl('histamine',names(go_sets)))
go_sets[grepl('histamine',names(go_sets))]

#ion channel genes will largely start with these
#KCN HCN SCN CACNA
#HTR genes for serotonin
#CHRN amd CHRM for acetylcholine


potassium_genes = rownames(multiome_sce)[grepl('KCN', rownames(multiome_sce))]

hcn_genes = rownames(multiome_sce)[grepl('HCN', rownames(multiome_sce))]

sodium_genes = rownames(multiome_sce)[grepl('SCN', rownames(multiome_sce)) & !grepl('FSCN', rownames(multiome_sce)) & !grepl('OBSCN', rownames(multiome_sce))]

calcium_genes = rownames(multiome_sce)[grepl('CACNA', rownames(multiome_sce))]

dopamine_genes = c("DRD1", "DRD2", "DRD3", "DRD4", "DRD5", 'TH')

serotonin_genes = rownames(multiome_sce)[grepl('HTR', rownames(multiome_sce))]

acetylcholine_genes = rownames(multiome_sce)[grepl('CHRN|CHRM', rownames(multiome_sce))]

endocannabinoid_genes = unique(unlist(unname(go_sets[grepl('endocannabinoid',names(go_sets))])))

opioid_genes = unique(unlist(unname(go_sets[grepl('opioid',names(go_sets))])))

glp_genes = c('GLP1R', 'GLP2R','ZGLP1')

table(multiome_sce$merged_cluster, multiome_sce$refined_cluster_ann)

get_gene_exp_heatmap = function(sce, assay_type, gene_list, gene_list_name, exp_filt){

  #Grab the expression of the genes of interest
  gene_index = rownames(sce) %in% gene_list
  marker_exp = assay(sce, assay_type)[gene_index, ]
  
  #Get the average expression of each gene by cell-type
  #This sums counts across cells per gene
  cell_annot_matrix <- Matrix::sparse.model.matrix(~ 0 + refined_mid_cluster, data = colData(sce))
  colnames(cell_annot_matrix) <- gsub("refined_mid_cluster", "", colnames(cell_annot_matrix))

  cell_counts = colSums(cell_annot_matrix)
  pseudobulk_marker_exp = marker_exp %*% cell_annot_matrix

  #Compute average and zscored expression
  avg_marker_exp <- sweep(pseudobulk_marker_exp, 2, cell_counts, "/")

  # Filter out genes with zero variance (they cause NaN in scaling)
  gene_vars = apply(avg_marker_exp, 1, var)
  avg_marker_exp = avg_marker_exp[gene_vars > 0, ]

  #Filter genes on expression level
  if (!missing(exp_filt)) {

    num_passing = sum(rowMeans(avg_marker_exp) >= exp_filt)
    avg_marker_exp = avg_marker_exp[rowMeans(avg_marker_exp) >= exp_filt, ]
  
    if (num_passing == 0) {
      stop("No genes passed the expression filter.")
    }
  }

  #Get scaled expression
  scaled_avg_marker_exp = t(scale(t(avg_marker_exp)))

  # Create row annotation with average expression across cell types as barplot
  gene_avg_exp = rowMeans(avg_marker_exp)
  
  row_anno = HeatmapAnnotation(
    avg_exp = anno_barplot(gene_avg_exp, baseline = 0, 
                           gp = gpar(col = NA, fill = "#A50026"),
                           height = unit(3, "cm")),
    annotation_name_side = "top",
    annotation_name_gp = gpar(fontsize = 8),
    which = "row",
    show_annotation_name = TRUE,
    show_legend = FALSE
  )

  # Create column annotation based on merged_cluster
  # Get the merged_cluster for each cell type (refined_mid_cluster)
  celltype_to_merged <- tapply(sce$merged_cluster, sce$refined_mid_cluster, function(x) x[1])
  col_anno_data <- data.frame(merged_cluster = celltype_to_merged[colnames(avg_marker_exp)])
  rownames(col_anno_data) <- colnames(avg_marker_exp)
  
  # Define colors for merged_cluster values
  merged_cluster_colors <- c("LHb" = my_colors_class[['LHb']], "MHb" = my_colors_class[['MHb']])
  
  col_anno = HeatmapAnnotation(
    df = col_anno_data,
    col = list(merged_cluster = merged_cluster_colors),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 8),
    which = "column",
    show_annotation_name = TRUE,
    show_legend = TRUE,
    simple_anno_size = unit(4, "mm")
  )

  col_func = rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
  zscore_heatmap = Heatmap(scaled_avg_marker_exp, name = 'Scaled average expression', 
  show_row_names = TRUE, show_column_names = TRUE,
  row_title = 'Cell types', column_title = gene_list_name,
  cluster_rows = TRUE, cluster_columns = TRUE, 
  clustering_method_rows = 'ward.D2', clustering_method_columns = 'ward.D2',
  col = col_func,
  #height = unit(30, "mm"),
  row_names_gp = gpar(fontsize = 6),      # Row text size
  column_names_gp = gpar(fontsize = 8),
  top_annotation = col_anno,
  right_annotation = row_anno)

  #zscore_heatmap = draw(zscore_heatmap)

  #And now at the finer resolution
  cell_annot_matrix <- Matrix::sparse.model.matrix(~ 0 + refined_cluster_ann, data = colData(sce))
  colnames(cell_annot_matrix) <- gsub("refined_cluster_ann", "", colnames(cell_annot_matrix))

  cell_counts = colSums(cell_annot_matrix)
  pseudobulk_marker_exp = marker_exp %*% cell_annot_matrix

  #Compute average and zscored expression
  avg_marker_exp <- sweep(pseudobulk_marker_exp, 2, cell_counts, "/")

  # Filter out genes with zero variance (they cause NaN in scaling)
  gene_vars = apply(avg_marker_exp, 1, var)
  avg_marker_exp = avg_marker_exp[gene_vars > 0, ]
  
  #Filter genes on expression level
  if (!missing(exp_filt)) {

    num_passing = sum(rowMeans(avg_marker_exp) >= exp_filt)
    avg_marker_exp = avg_marker_exp[rowMeans(avg_marker_exp) >= exp_filt, ]
  
    if (num_passing == 0) {
      return(zscore_heatmap)  # Return the previous heatmap if no genes pass the filter
    }
  }

  scaled_avg_marker_exp = t(scale(t(avg_marker_exp)))

  # Create row annotation with average expression across cell types as barplot
  gene_avg_exp = rowMeans(avg_marker_exp)
  
  row_anno = HeatmapAnnotation(
    avg_exp = anno_barplot(gene_avg_exp, baseline = 0, 
                           gp = gpar(col = NA, fill = "#A50026"),
                           height = unit(3, "cm")),
    annotation_name_side = "top",
    annotation_name_gp = gpar(fontsize = 8),
    which = "row",
    show_annotation_name = TRUE,
    show_legend = FALSE
  )

  #Plot with complex heatmap, cluster rows and columns, with avewrage expression and z-scored expression

  # Create column annotation based on merged_cluster
  # Get the merged_cluster for each cell type (refined_cluster_ann)
  celltype_to_merged_fine <- tapply(sce$merged_cluster, sce$refined_cluster_ann, function(x) x[1])
  col_anno_data_fine <- data.frame(merged_cluster = celltype_to_merged_fine[colnames(avg_marker_exp)])
  rownames(col_anno_data_fine) <- colnames(avg_marker_exp)
  
  col_anno_fine = HeatmapAnnotation(
    df = col_anno_data_fine,
    col = list(merged_cluster = merged_cluster_colors),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 8),
    which = "column",
    show_annotation_name = TRUE,
    show_legend = TRUE,
    simple_anno_size = unit(4, "mm")
  )

  col_func = rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
  fine_zscore_heatmap = Heatmap(scaled_avg_marker_exp, name = 'Scaled average expression', 
  show_row_names = TRUE, show_column_names = TRUE,
  row_title = 'Cell types', column_title = gene_list_name,
  cluster_rows = TRUE, cluster_columns = TRUE, 
  clustering_method_rows = 'ward.D2', clustering_method_columns = 'ward.D2',
  col = col_func,
  #height = unit(30, "mm"),
  row_names_gp = gpar(fontsize = 6),      # Row text size
  column_names_gp = gpar(fontsize = 8),
  top_annotation = col_anno_fine,
  right_annotation = row_anno)


  return(list(zscore_heatmap = zscore_heatmap, fine_zscore_heatmap = fine_zscore_heatmap))

}

plot_assay = 'logcounts'
expression_filter = log10(1.5)
#hcn_heatmaps = get_gene_exp_heatmap(multiome_sce, plot_assay, hcn_genes, 'HCN genes', exp_filt = log10(2))
#draw(hcn_heatmaps[[1]])
#draw(hcn_heatmaps[[2]])

calcium_heatmaps = get_gene_exp_heatmap(multiome_sce, plot_assay, calcium_genes, 'Calcium genes', exp_filt = expression_filter)
draw(calcium_heatmaps[[1]])
draw(calcium_heatmaps[[2]])

sodium_heatmaps = get_gene_exp_heatmap(multiome_sce, plot_assay, sodium_genes, 'Sodium genes', exp_filt = expression_filter)
draw(sodium_heatmaps[[1]])
draw(sodium_heatmaps[[2]])

potassium_heatmaps = get_gene_exp_heatmap(multiome_sce, plot_assay, potassium_genes, 'Potassium genes', exp_filt = expression_filter)
draw(potassium_heatmaps[[1]])
draw(potassium_heatmaps[[2]])

#dopamine_heatmaps = get_gene_exp_heatmap(multiome_sce, plot_assay, dopamine_genes, 'Dopamine genes', exp_filt = expression_filter)
#draw(dopamine_heatmaps[[1]])
#draw(dopamine_heatmaps[[2]])

serotonin_heatmaps = get_gene_exp_heatmap(multiome_sce, plot_assay, serotonin_genes, 'Serotonin genes', exp_filt = expression_filter)
draw(serotonin_heatmaps[[1]])
draw(serotonin_heatmaps[[2]])

acetylcholine_heatmaps = get_gene_exp_heatmap(multiome_sce, plot_assay, acetylcholine_genes, 'Acetylcholine genes', exp_filt = expression_filter)
draw(acetylcholine_heatmaps[[1]])
draw(acetylcholine_heatmaps[[2]])

endocannabinoid_heatmaps = get_gene_exp_heatmap(multiome_sce, plot_assay, endocannabinoid_genes, 'Endocannabinoid genes', exp_filt = expression_filter)
draw(endocannabinoid_heatmaps[[1]])
draw(endocannabinoid_heatmaps[[2]])

opioid_heatmaps = get_gene_exp_heatmap(multiome_sce, plot_assay, opioid_genes, 'Opioid genes', exp_filt = expression_filter)
draw(opioid_heatmaps[[1]])
draw(opioid_heatmaps[[2]])

#glp_heatmaps = get_gene_exp_heatmap(multiome_sce, plot_assay, glp_genes, 'GLP genes', exp_filt = expression_filter)
#draw(glp_heatmaps[[1]])
#draw(glp_heatmaps[[2]])



#Save plots

pdf(file.path(plot_path, 'calcium_genes_heatmap.pdf'), width = 8, height = 4)
draw(calcium_heatmaps[[1]])
dev.off()

pdf(file.path(plot_path, 'sodium_genes_heatmap.pdf'), width = 8, height = 4)
draw(sodium_heatmaps[[1]])
dev.off()

pdf(file.path(plot_path, 'potassium_genes_heatmap.pdf'), width = 8, height = 4)
draw(potassium_heatmaps[[1]])
dev.off()

pdf(file.path(plot_path, 'serotonin_genes_heatmap.pdf'), width = 8, height = 4)
draw(serotonin_heatmaps[[1]])
dev.off()

pdf(file.path(plot_path, 'acetylcholine_genes_heatmap.pdf'), width = 8, height = 4)
draw(acetylcholine_heatmaps[[1]])
dev.off()

pdf(file.path(plot_path, 'endocannabinoid_genes_heatmap.pdf'), width = 8, height = 4)
draw(endocannabinoid_heatmaps[[1]])
dev.off()

pdf(file.path(plot_path, 'opioid_genes_heatmap.pdf'), width = 8, height = 4)
draw(opioid_heatmaps[[1]])
dev.off()




