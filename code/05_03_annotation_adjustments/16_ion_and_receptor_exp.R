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
new_data_path = here('processed-data','05_03_annotation_adjustments', '16_ion_and_receptor_exp')
if (!dir.exists(plot_path)) dir.create(plot_path)
if (!dir.exists(new_data_path)) dir.create(new_data_path)

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
hab_celltypes = c('GABA_LHb_C.1','GABA_LHb_C.2','LHb_C','LHb_A','LHb_B','MHb_A','MHb_C','MHb_B','MHb_D')
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

#These are particularly interesting because they determine a tonic firing pattern at rest, these are cation channels that are always open. 
hcn_genes = rownames(multiome_sce)[grepl('HCN', rownames(multiome_sce))]

sodium_genes = rownames(multiome_sce)[grepl('SCN', rownames(multiome_sce)) & !grepl('FSCN', rownames(multiome_sce)) & !grepl('OBSCN', rownames(multiome_sce))]

calcium_genes = rownames(multiome_sce)[grepl('CACNA', rownames(multiome_sce))]

#dopamine_genes = c("DRD1", "DRD2", "DRD3", "DRD4", "DRD5", 'TH')

serotonin_genes = rownames(multiome_sce)[grepl('HTR', rownames(multiome_sce))]
serotonin_genes = c(serotonin_genes, 'SLC6A4') #Including the serotonin transporter protein, primary target of SSRIs

acetylcholine_genes = rownames(multiome_sce)[grepl('CHRN|CHRM', rownames(multiome_sce))]

endocannabinoid_genes = unique(unlist(unname(go_sets[grepl('endocannabinoid',names(go_sets))])))

opioid_genes = c(unique(unlist(unname(go_sets[grepl('opioid',names(go_sets))]))))

#glp_genes = c('GLP1R', 'GLP2R','ZGLP1')

#And the GPR orphan receptors
gpr_genes = rownames(multiome_sce)[grepl('GPR', rownames(multiome_sce)) & !grepl('GPRASP', rownames(multiome_sce)) & !grepl('GPRIN', rownames(multiome_sce))]

#And the associated channels with the non-canonical inhibition of habenula neurons
#Supposedly, high levels of SLC12A1 and SLC12A2 (Na-K-Cl cotransporters) and a lack of expression of SLC12A5 (KCC2, K and Cl pump out of cells)
# Leads to a high level of intracellular Cl- and thus GABA-A receptor activation leads to depolarization instead of hyperpolarization
#In combo, calcium activated chloride channels are the main driver of inactivation of habenula neurons
#In mice, the primary CaCC is TMEM16A in the habenula, also ANO1, ANO2 is TMEM16B
non_canon_inhib_genes = c('SLC12A1', 'SLC12A2', 'SLC12A5', 'ANO1', 'ANO2')

#And add in the GABA and Glut receptors from the gene set enrichment analysis
gaba_genes = c("ATF4", "CACNB4", "GABBR1", "GABRA1", "GABRA2", "GABRA3", "GABRA4", "GABRA5", "GABRA6", "GABRB1", "GABRB2", "GABRB3",
               "GABRD", "GABRE", "GABRG1", "GABRG2", "GABRG3", "GABRR1", "GABRR2", "GNAI2", "HTR1A", "PLCL1", "SLC12A2", "GABBR2",
               "PLCL2", "PHF24", "GPR156", "PIANP", "GABRR3", "SHISA7", "CBLN1", "FGF13", "PLXNB1", "WNT5A", "CLSTN3", "SEMA4D",
               "SRGAP2", "LGI2", "NLGN2", "CLSTN2", "SEMA4A", "CBLN4", "MDGA1", "NPAS4", "LHFPL4", "HAPLN4", "SRGAP2C", "GABRP",
               "GRID1", "GABRQ", "DTNB", "GAD1", "GAD2", "GLRA1", "NLGN4Y", "SYT11", "NLGN3", "NLGN4X", "IGSF9", "IGSF21",
               "SLC32A1", "CEP112", "IQSEC3")


glut_genes = c("GRIA1", "GRIA2", "GRIA3", "GRIA4", "GRID1", "GRID2", "GRIK1", "GRIK2", "GRIK3", "GRIK4", "GRIK5", "GRIN1", "GRIN2A",
               "GRIN2B", "GRIN2C", "GRIN2D", "GRM1", "GRM2", "GRM3", "GRM4", "GRM5", "GRM6", "GRM7", "GRM8", "GRIN3A", "GRIN3B")


#Save a table with the genes of interest
genes_of_interest_df = data.frame(gene = c(hcn_genes, potassium_genes, sodium_genes, calcium_genes, serotonin_genes, acetylcholine_genes, endocannabinoid_genes, opioid_genes, gpr_genes, non_canon_inhib_genes, gaba_genes, glut_genes), 
category = c(rep('HCN', length(hcn_genes)), rep('Potassium', length(potassium_genes)), rep('Sodium', length(sodium_genes)), rep('Calcium', length(calcium_genes)), rep('Serotonin', length(serotonin_genes)), rep('Acetylcholine', length(acetylcholine_genes)), 
rep('Endocannabinoid', length(endocannabinoid_genes)), rep('Opioid', length(opioid_genes)), rep('GPR', length(gpr_genes)), rep('Non-canonical inhibitory', length(non_canon_inhib_genes)), rep('GABA', length(gaba_genes)), rep('Glutamate', length(glut_genes))))

genes_of_interest_df
write.csv(genes_of_interest_df, file.path(new_data_path, 'genes_of_interest.csv'), row.names = FALSE)




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

    num_passing = sum(matrixStats::rowMaxs(as.matrix(avg_marker_exp)) >= exp_filt)
    avg_marker_exp = avg_marker_exp[matrixStats::rowMaxs(as.matrix(avg_marker_exp)) >= exp_filt, ]
  
    if (num_passing == 0) {
      stop("No genes passed the expression filter.")
    }
  }

  #Get scaled expression
  scaled_avg_marker_exp = t(scale(t(avg_marker_exp)))

  #Get max-normalized expression
  avg_marker_exp <- as.matrix(avg_marker_exp)
  max_normalized_exp = sweep(avg_marker_exp, 1, matrixStats::rowMaxs(avg_marker_exp), "/")

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
  merged_cluster_colors <- c("LHb" = my_colors_class[['LHb']], "MHb" = my_colors_class[['MHb']], 
  'Excit_Thal' = my_colors_class[['Thalamus']], 'Inhib_Thal' = my_colors_class[['Thalamus']])
  
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

  col_func = grDevices::colorRampPalette(c("white", "#FFFFBF", "#FEE08B", "#FDAE61", "#F46D43", "#D73027", "#A50026"))(100)
  zscore_heatmap = Heatmap(max_normalized_exp , name = 'Max norm. avg. exp.', 
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

    num_passing = sum(matrixStats::rowMaxs(as.matrix(avg_marker_exp)) >= exp_filt)
    avg_marker_exp = avg_marker_exp[matrixStats::rowMaxs(as.matrix(avg_marker_exp)) >= exp_filt, ]
  
    if (num_passing == 0) {
      return(zscore_heatmap)  # Return the previous heatmap if no genes pass the filter
    }
  }

  scaled_avg_marker_exp = t(scale(t(avg_marker_exp)))

  #Get max-normalized expression
  avg_marker_exp <- as.matrix(avg_marker_exp)
  max_normalized_exp = sweep(avg_marker_exp, 1, matrixStats::rowMaxs(avg_marker_exp), "/")

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

  col_func = grDevices::colorRampPalette(c("white", "#FFFFBF", "#FEE08B", "#FDAE61", "#F46D43", "#D73027", "#A50026"))(100)
  fine_zscore_heatmap = Heatmap(max_normalized_exp , name = 'Max norm. avg. exp.', 
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
hcn_heatmaps = get_gene_exp_heatmap(multiome_sce, plot_assay, hcn_genes, 'HCN genes', exp_filt = expression_filter)
draw(hcn_heatmaps[[1]])
draw(hcn_heatmaps[[2]])

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

gpr_heatmaps = get_gene_exp_heatmap(multiome_sce, plot_assay, gpr_genes, 'GPR genes', exp_filt = expression_filter)
draw(gpr_heatmaps[[1]])
draw(gpr_heatmaps[[2]])

#glp_heatmaps = get_gene_exp_heatmap(multiome_sce, plot_assay, glp_genes, 'GLP genes', exp_filt = expression_filter)
#draw(glp_heatmaps[[1]])
#draw(glp_heatmaps[[2]])

glut_heatmaps = get_gene_exp_heatmap(multiome_sce, plot_assay, glut_genes, 'Glutamate genes', exp_filt = expression_filter)
draw(glut_heatmaps[[1]])
draw(glut_heatmaps[[2]])

gaba_heatmaps = get_gene_exp_heatmap(multiome_sce, plot_assay, gaba_genes, 'GABA genes', exp_filt = expression_filter)
draw(gaba_heatmaps[[1]])
draw(gaba_heatmaps[[2]])

#Adjust the minimum expression filter for this one
non_canon_inhib_heatmaps = get_gene_exp_heatmap(multiome_sce, plot_assay, non_canon_inhib_genes, 
  'Non-canonical inhibitory genes', 
exp_filt = 0)
draw(non_canon_inhib_heatmaps[[1]])
draw(non_canon_inhib_heatmaps[[2]])


#Save plots

pdf(file.path(plot_path, 'HCN_genes_heatmap.pdf'), width = 8, height = 4)
draw(hcn_heatmaps[[1]])
dev.off()

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

pdf(file.path(plot_path, 'GPR_genes_heatmap.pdf'), width = 8, height = 4)
draw(gpr_heatmaps[[1]])
dev.off()

pdf(file.path(plot_path, 'glutamate_genes_heatmap.pdf'), width = 8, height = 4)
draw(glut_heatmaps[[1]])
dev.off()


pdf(file.path(plot_path, 'GABA_genes_heatmap.pdf'), width = 8, height = 4)
draw(gaba_heatmaps[[1]])
dev.off()

pdf(file.path(plot_path, 'non_canonical_inhibitory_genes_heatmap.pdf'), width = 8, height = 4)
draw(non_canon_inhib_heatmaps[[1]])
dev.off()


#And check out expression in the UMAP too

#Reload the full dataset
multiome_sce = qs_read(paste0(multiome_path, '/refined_annotation_multiomeHab_SCE.qs2'))
assay(multiome_sce, 'cpm') = MetaMarkers::convert_to_cpm(assay(multiome_sce, 'counts'))

plt2 <- scater::plotReducedDim(multiome_sce, 
                dimred = "umap.integrated",
                colour_by = "refined_mid_cluster", 
              point_size = .5) +
  scale_color_manual(values = my_colors_mid) +
  guides(colour = guide_legend(override.aes = list(size = 3), reverse = TRUE)) +
  labs(
    title = "RNA UMAP",
    x = "UMAP 1",
    y = "UMAP 2",
    colour = "Cell type"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )

plt2



plot_gene_on_umap <- function(sce, gene_name, reduction = "umap.integrated", assay_type = "logcounts") {
  
  # Check that gene exists
  if (!gene_name %in% rownames(sce)) {
    stop(paste0("Gene '", gene_name, "' not found in object"))
  }
  
  # Extract UMAP coordinates
  umap_coords <- reducedDim(sce, reduction)
  
  # Extract gene expression
  expr <- assay(sce, assay_type)[gene_name, ]
  
  # Create dataframe
  plot_df <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    expression = expr
  )
  
  # Create plot
  p <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = expression)) +
    geom_point(size = 1) +
    scale_color_gradient(low = "white", high = "red", name = paste0(gene_name, "\nExpression")) +
    labs(title = paste0(gene_name, " expression on UMAP (integrated)"),
         x = "UMAP 1", y = "UMAP 2") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      aspect.ratio = 1
    )
  
  return(p)
}

# Test with Oprm1
plot_gene_on_umap(multiome_sce, "OPRM1")

plot_gene_on_umap(multiome_sce, "SLC6A4")
plot_gene_on_umap(multiome_sce, "MME")
plot_gene_on_umap(multiome_sce, "GPR149")
plot_gene_on_umap(multiome_sce, "GSDME")
plot_gene_on_umap(multiome_sce, "MIR548XHG")

plot_gene_on_umap(multiome_sce, "CHRNA6")
plot_gene_on_umap(multiome_sce, "CHRNB3")
plot_gene_on_umap(multiome_sce, "CHRNA3")
plot_gene_on_umap(multiome_sce, "CHRNA5")
plot_gene_on_umap(multiome_sce, "CHRNB4")

plot_gene_on_umap(multiome_sce, "CHRM2")
plot_gene_on_umap(multiome_sce, "CHRM3")


plot_gene_on_umap(multiome_sce, "MGLL")
plot_gene_on_umap(multiome_sce, "FABP5")
plot_gene_on_umap(multiome_sce, "CNR1")
plot_gene_on_umap(multiome_sce, "CNR2")

plot_gene_on_umap(multiome_sce, "GRM5")
plot_gene_on_umap(multiome_sce, "PLCB1")


plot_gene_on_umap(multiome_sce, "CACNA1C")
plot_gene_on_umap(multiome_sce, "CACNA2D2")
plot_gene_on_umap(multiome_sce, "CACNA1E")

plot_gene_on_umap(multiome_sce, "CACNA1D")
plot_gene_on_umap(multiome_sce, "CACNA2D3")
plot_gene_on_umap(multiome_sce, "CACNA1B")
plot_gene_on_umap(multiome_sce, "CACNA1A")
plot_gene_on_umap(multiome_sce, "CACNA1A")

#Non-canonical inhibitory associated genes
plot_gene_on_umap(multiome_sce, "SLC12A1")
plot_gene_on_umap(multiome_sce, "SLC12A2")
plot_gene_on_umap(multiome_sce, "SLC12A5")
plot_gene_on_umap(multiome_sce, "ANO1")
plot_gene_on_umap(multiome_sce, "ANO2")

#Some mouse evidence that PVALB+ GABA+ neurons are the locally projecting inhibitory neurons in the lateral habenula
plot_gene_on_umap(multiome_sce, "PVALB")
plot_gene_on_umap(multiome_sce, "SLC32A1")
plot_gene_on_umap(multiome_sce, "SST")

plot_gene_on_umap(multiome_sce, "HCN1")
plot_gene_on_umap(multiome_sce, "HCN2")
plot_gene_on_umap(multiome_sce, "HCN3")

plot_gene_on_umap(multiome_sce, "CACNA1G")
plot_gene_on_umap(multiome_sce, "CACNA1I")

plot_gene_on_umap(multiome_sce, "NALCN")

plot_gene_on_umap(multiome_sce, "HTR2A")
plot_gene_on_umap(multiome_sce, "HTR2C")
plot_gene_on_umap(multiome_sce, "HTR4")
plot_gene_on_umap(multiome_sce, "HTR7")

plot_gene_on_umap(multiome_sce, "GNAO1")
plot_gene_on_umap(multiome_sce, "ADCY8")
plot_gene_on_umap(multiome_sce, "OGFRL1")
plot_gene_on_umap(multiome_sce, "OPRM1")
plot_gene_on_umap(multiome_sce, "GNAS")
plot_gene_on_umap(multiome_sce, "SYP")
plot_gene_on_umap(multiome_sce, "OGFR")
plot_gene_on_umap(multiome_sce, "SIGMAR1")
plot_gene_on_umap(multiome_sce, "PENK")
plot_gene_on_umap(multiome_sce, "PNOC")

plot_gene_on_umap(multiome_sce, "PVALB")

source(here('code','05_03_annotation_adjustments','celltype_colors.R'))




plot_violin_maxnorm <- function(sce, gene, assay_name = "logcounts", celltype_col, study_label) {
  expr <- assay(sce, assay_name)[gene, ]
  celltype <- colData(sce)[[celltype_col]]
  
  df <- data.frame(
    expression = as.numeric(expr),
    celltype = celltype
  )
  
  # Max-normalize
  max_exp <- max(df$expression)
  
  if (max_exp == 0) {
    warning("Gene '", gene, "' has zero expression across all cell types.")
    df$norm_expr <- 0
  } else {
    df$norm_expr <- df$expression / max_exp
  }
  
  # Compute medians per celltype for the dot overlay
  medians <- df |>
    summarise(median_expr = median(expression), .by = celltype)
  
  ggplot(df, aes(x = celltype, y = expression, fill = celltype)) +
    geom_violin(
      scale = "width",
      width = 0.9,
      color = NA,
      trim = TRUE
    ) +
    scale_fill_manual(values = my_colors_mid, guide = 'none') +
    geom_point(
      data = medians,
      aes(x = celltype, y = median_expr),
      size = 1.5,
      color = "black"
    ) +
    scale_y_continuous(limits = c(0, 6), breaks = c(0, 2, 4, 6)) +
    labs(x = NULL, y = "Logcounts", title = sprintf('%s - %s', study_label, gene)) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
      axis.ticks.x = element_blank(),
      plot.title = element_text(face = "italic", size = 11),
      panel.grid = element_blank(),
      # Condensed vertical height — control via coord_cartesian or plot sizing
      aspect.ratio = 0.15
    )
}


#Load up the Yalcinbas data too to look at KCC2 expression
yalcinbas_path = '/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/sce_objects'
load(paste0(yalcinbas_path, '/official_final_sce.RDATA'))
yalcinbas_sce = sce
table(yalcinbas_sce$final_Annotations) 

yalcinbas_sce$consensus_annot = yalcinbas_sce$final_Annotations
yalcinbas_sce$consensus_annot[yalcinbas_sce$final_Annotations == 'MHb.1'] = 'MHb_A'
yalcinbas_sce$consensus_annot[yalcinbas_sce$final_Annotations == 'MHb.2'] = 'MHb_B'
yalcinbas_sce$consensus_annot[yalcinbas_sce$final_Annotations == 'MHb.3'] = 'MHb_D'
yalcinbas_sce$consensus_annot[yalcinbas_sce$final_Annotations %in% c('LHb.2', 'LHb.7')] = 'LHb_A'
yalcinbas_sce$consensus_annot[yalcinbas_sce$final_Annotations %in% c('LHb.1', 'LHb.5')] = 'LHb_B'
yalcinbas_sce$consensus_annot[yalcinbas_sce$final_Annotations %in% c('LHb.3', 'LHb.4')] = 'LHb_C'

table(yalcinbas_sce$consensus_annot, yalcinbas_sce$final_Annotations)

yalcinbas_neuron = yalcinbas_sce[, yalcinbas_sce$consensus_annot %in% c('MHb_A', 'MHb_B', 'LHb_A', 'LHb_B', 'LHb_C')]
yalcinbas_neuron$consensus_annot = factor(yalcinbas_neuron$consensus_annot, levels = c('MHb_B', 'MHb_A', 'LHb_A', 'LHb_B', 'LHb_C'))

yalcinbas_kcc2_p = plot_violin_maxnorm(yalcinbas_neuron, "SLC12A5", assay_name = "logcounts", 
celltype_col = "consensus_annot", study_label = "Yalcinbas")

yalcinbas_slc6a4_p = plot_violin_maxnorm(yalcinbas_neuron, "SLC6A4", assay_name = "logcounts", 
celltype_col = "consensus_annot", study_label = "Yalcinbas")


multiome_sce$consensus_annot = multiome_sce$refined_mid_cluster
multiome_sce$consensus_annot[multiome_sce$refined_mid_cluster %in% c('GABA_LHb_C.1','GABA_LHb_C.2')] = 'LHb_C'

multiome_neuron = multiome_sce[, multiome_sce$consensus_annot %in% c('MHb_B', 'MHb_A', 'LHb_A', 'LHb_B', 'LHb_C')]
multiome_neuron$consensus_annot = factor(multiome_neuron$consensus_annot, levels = c('MHb_B', 'MHb_A', 'LHb_A', 'LHb_B', 'LHb_C'))

multiome_kcc2_p = plot_violin_maxnorm(multiome_neuron, "SLC12A5", assay_name = "logcounts", 
celltype_col = "consensus_annot", study_label = "Multiome")

multiome_slc6a4_p = plot_violin_maxnorm(multiome_neuron, "SLC6A4", assay_name = "logcounts", 
celltype_col = "consensus_annot", study_label = "Multiome")

yalcinbas_kcc2_p
multiome_kcc2_p
ggsave(plot = yalcinbas_kcc2_p, filename = file.path(plot_path, 'Yalcinbas_SLC12A5_violin.pdf'), width = 6, height = 2)
ggsave(plot = multiome_kcc2_p, filename = file.path(plot_path, 'Multiome_SLC12A5_violin.pdf'), width = 6, height = 2)


#Adding regression lines to the violin plots
# Extract SLC12A5 expression and cell type for both datasets
# Yalcinbas
yalcinbas_df <- data.frame(
  expression = as.numeric(logcounts(yalcinbas_neuron["SLC12A5", ])),
  celltype = yalcinbas_neuron$consensus_annot
)
yalcinbas_df$celltype_numeric <- as.numeric(yalcinbas_df$celltype)  # 1=MHb_B, 2=MHb_A, ..., 5=LHb_C

# Multiome
multiome_df <- data.frame(
  expression = as.numeric(logcounts(multiome_neuron["SLC12A5", ])),
  celltype = multiome_neuron$consensus_annot
)
multiome_df$celltype_numeric <- as.numeric(multiome_df$celltype)

# Fit linear trend models
yalcinbas_lm <- lm(expression ~ celltype_numeric, data = yalcinbas_df)
multiome_lm <- lm(expression ~ celltype_numeric, data = multiome_df)

summary(yalcinbas_lm)
summary(multiome_lm)$coefficients[2, 4]


# Add regression line and slope annotation to the violin plots
add_trend_line <- function(p, lm_fit, df) {
  slope <- coef(lm_fit)[2]
  pval <- summary(lm_fit)$coefficients[2, 4]
  r2 <- summary(lm_fit)$r.squared
  
  if(pval < 2.2e-16) {
    pval_text <- "p < 2.2e-16"
  } else {
    pval_text <- sprintf("p = %.3f", pval)
  }
  label <- sprintf("slope = %.3f, R² = %.3f, %s", slope, r2, pval_text)
  
  # Get predicted values at each level for the line
  pred_df <- data.frame(
    celltype_numeric = 1:5,
    celltype = levels(df$celltype)
  )
  pred_df$predicted <- predict(lm_fit, newdata = pred_df)
  
  p + 
    geom_line(data = pred_df, aes(x = celltype, y = predicted, group = 1), 
              color = "red", linewidth = 0.8) +
    geom_point(data = pred_df, aes(x = celltype, y = predicted), 
               color = "red", size = 2) +
    annotate("text", x = 3, y = max(df$expression) * 0.9, 
             label = label, size = 3, color = "red")
}

yalcinbas_kcc2_trend <- add_trend_line(yalcinbas_kcc2_p, yalcinbas_lm, yalcinbas_df)
yalcinbas_kcc2_trend
multiome_kcc2_trend <- add_trend_line(multiome_kcc2_p, multiome_lm, multiome_df)
multiome_kcc2_trend

ggsave(plot = yalcinbas_kcc2_trend, filename = file.path(plot_path, 'Yalcinbas_SLC12A5_violin_trend.pdf'), width = 6, height = 2)
ggsave(plot = multiome_kcc2_trend, filename = file.path(plot_path, 'Multiome_SLC12A5_violin_trend.pdf'), width = 6, height = 2)


