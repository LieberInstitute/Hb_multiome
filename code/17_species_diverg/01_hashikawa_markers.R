#Generating the cell-type marker stats of the consensus cell-types in the hashikawa dataset

library(here)
library("DeconvoBuddies")
library(sessioninfo)
library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(tidyverse)
library(Matrix)
library(edgeR)
library(MetaMarkers)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(ggplot2)



hashikawa_path = '/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/09_cross_species_analysis/Hashikawa_data'
hashikawa_data_path = here('processed-data','05_02_external_Hb_comparisons','04_wallace_hashikawa_mouse')

wallace_path = here('processed-data', '05_02_external_Hb_comparisons', '02_qc_and_clust_wallace_2019')
wallace_data_path = here('processed-data', '05_02_external_Hb_comparisons', '03_metaMarkers_wallace_2019')


new_data_path = here('processed-data','17_species_diverg','01_hashikawa_markers')
plot_path = here('plots','17_species_diverg','01_hashikawa_markers')
dir.create(new_data_path, showWarnings = FALSE)
dir.create(plot_path, showWarnings = FALSE)
#colors
source(here('code','05_03_annotation_adjustments','celltype_colors.R'))


#Path to orthologs
path_to_orthologs = here('processed-data', '05_02_external_Hb_comparisons', 'human_mouse_zebrafish_orthologs.txt.gz')


#Wallace mouse data
wallace_sce = readRDS(paste0(wallace_path, '/all_donor_sce_with_denovo_clusters.rds'))

# Add the mouse metadata
current_mouse_metadata = readRDS(paste0(wallace_data_path, '/wallace_mouse_metaclust_celltype_annot_metadata.rds'))
colData(wallace_sce) = S4Vectors::DataFrame(current_mouse_metadata)




#Loads the sce_mouse_sub list, has two sce objects. We have the consensus annotations for just the neuronal subset object
# Add the consensus annotations, and then grab the non-neurons from the full object to piece back together

load(paste0(hashikawa_path, '/sce_mouse_habenula.Rdata'))
hashikawa_sce_all = sce_mouse_sub$all #Using just the neuron subset
hashikawa_sce_neuron = sce_mouse_sub$neuron

assay(hashikawa_sce_all, 'cpm') = MetaMarkers::convert_to_cpm(assay(hashikawa_sce_all, 'counts'))
assay(hashikawa_sce_neuron, 'cpm') = MetaMarkers::convert_to_cpm(assay(hashikawa_sce_neuron, 'counts'))

# Add in the updated mouse annotations we generated for this subset
neuron_hashikawa_metadata = readRDS(paste0(hashikawa_data_path, '/hashikawa_mouse_neuron_metaclust_celltype_annot_metadata.rds'))
#Double checking this aligns with the neuronal subset cells
table(rownames(neuron_hashikawa_metadata) == colnames(hashikawa_sce_neuron))
#Pass in the updated cell-type annotations
colData(hashikawa_sce_neuron) = neuron_hashikawa_metadata 

#Grab the non-neurons from the full dataset.
nonN_index = !grepl('Neuron',hashikawa_sce_all$celltype )
temp_nonNs_sce = hashikawa_sce_all[ , nonN_index]

#Fix up the metadata for the two objects, just swap the celltype column of the neurons for the meta_clust_celltyp_annot
hashikawa_sce_neuron$celltype = hashikawa_sce_neuron$meta_clust_celltype_annot

col_index = colnames(colData(hashikawa_sce_neuron)) != 'meta_clust_celltype_annot'
colData(hashikawa_sce_neuron) = colData(hashikawa_sce_neuron)[ , col_index]

#Fix up misaligned genes too, the neuronal subset has fewer genes
shared_genes = intersect(rownames(temp_nonNs_sce), rownames(hashikawa_sce_neuron))
temp_nonNs_sce = temp_nonNs_sce[rownames(temp_nonNs_sce) %in% shared_genes, ]
hashikawa_sce_neuron = hashikawa_sce_neuron[rownames(hashikawa_sce_neuron) %in% shared_genes, ]

#And check the order
table(rownames(temp_nonNs_sce) == rownames(hashikawa_sce_neuron))
#And combine
hashikawa_final_sce = cbind(temp_nonNs_sce, hashikawa_sce_neuron)
table(hashikawa_final_sce$celltype)


#This is the consensus cross-species mapping. Notably, we're leaving the mouse MHb_subP_cholinergic cell-type as is, as did not find a strong match in humans
#Also group the non-neuronal subtypes
hashikawa_final_sce$consensus_annot = hashikawa_final_sce$celltype
hashikawa_final_sce$consensus_annot[hashikawa_final_sce$celltype %in% c('MHb_subP')] = 'MHb_A'
hashikawa_final_sce$consensus_annot[hashikawa_final_sce$celltype %in% c(paste0('MHb_cholinergic_', 1:4))] = 'MHb_B'
hashikawa_final_sce$consensus_annot[hashikawa_final_sce$celltype %in% c('LHb_2', 'LHb_1_5')] = 'LHb_A'
hashikawa_final_sce$consensus_annot[hashikawa_final_sce$celltype %in% c('LHb_1_2', 'LHb_1_4')] = 'LHb_B'
hashikawa_final_sce$consensus_annot[hashikawa_final_sce$celltype %in% c('LHb_1_1', 'LHb_1_3')] = 'LHb_C'
hashikawa_final_sce$consensus_annot[grepl('Astrocyte', hashikawa_final_sce$celltype) ] = 'Astrocyte'
hashikawa_final_sce$consensus_annot[grepl('Oligo', hashikawa_final_sce$celltype) ] = 'Oligo'
hashikawa_final_sce$consensus_annot[grepl('OPC', hashikawa_final_sce$celltype) ] = 'OPC'

table(hashikawa_final_sce$consensus_annot)

#wallace annotations
wallace_sce$consensus_annot = wallace_sce$meta_clust_celltype_annot
wallace_sce$consensus_annot[wallace_sce$meta_clust_celltype_annot %in% c('MHb_subP')] = 'MHb_A'
wallace_sce$consensus_annot[wallace_sce$meta_clust_celltype_annot %in% c('MHb_cholinergic')] = 'MHb_B'
wallace_sce$consensus_annot[wallace_sce$meta_clust_celltype_annot %in% c('LHb_2')] = 'LHb_A'
wallace_sce$consensus_annot[wallace_sce$meta_clust_celltype_annot %in% c('LHb_1')] = 'LHb_C'
wallace_sce$consensus_annot[wallace_sce$meta_clust_celltype_annot %in% c('Astrocytes')] = 'Astrocyte'
wallace_sce$consensus_annot[wallace_sce$meta_clust_celltype_annot %in% c('Oligodendrocytes')] = 'Oligo'
wallace_sce$consensus_annot[wallace_sce$meta_clust_celltype_annot %in% c('Differentiating Oligodendrocytes')] = 'OPC'

table(wallace_sce$consensus_annot)


#Gene symbols as the rownames
rownames(hashikawa_final_sce) = rowData(hashikawa_final_sce)$Symbol
#rownames(wallace_sce) = rowData(wallace_sce)$Symbol

#Maybe clear up some memory, dont need the reduced dims for the t-stats
reducedDims(hashikawa_final_sce) <- list()
rm(temp_nonNs_sce,hashikawa_sce_neuron,hashikawa_sce_all , sce_mouse_sub  )
gc()


#And swap in the human gene name 
hu_mu_zf_ortholog_df = data.table::fread(path_to_orthologs)
table(hu_mu_zf_ortholog_df$`Mouse homology type`)

#Get the 1to1 orthologs for the mouse
hu_mu_zf_ortholog_df <- hu_mu_zf_ortholog_df %>%
  filter(
    `Mouse homology type` == 'ortholog_one2one'
  ) %>%
  filter(!duplicated(`Gene name`))
dim(hu_mu_zf_ortholog_df)

index = match(rownames(wallace_sce), hu_mu_zf_ortholog_df$`Mouse gene name`)
rownames(wallace_sce) = hu_mu_zf_ortholog_df$`Gene name`[index]

index = match(rownames(hashikawa_final_sce), hu_mu_zf_ortholog_df$`Mouse gene name`)
rownames(hashikawa_final_sce) = hu_mu_zf_ortholog_df$`Gene name`[index]

#And filter all the NAs
na_gene_index = !is.na(rownames(hashikawa_final_sce))
hashikawa_final_sce = hashikawa_final_sce[na_gene_index, ]
table(is.na(rownames(hashikawa_final_sce)))

na_gene_index = !is.na(rownames(wallace_sce))
wallace_sce = wallace_sce[na_gene_index, ]
table(is.na(rownames(wallace_sce)))


#logcounts normalization
hashikawa_final_sce = logNormCounts(hashikawa_final_sce)
wallace_sce = logNormCounts(wallace_sce)

#Now get the marker stats
#MeanRatio markers
marker_stats_MeanRatio <- get_mean_ratio(
    sce = hashikawa_final_sce, # sce is the SingleCellExperiment with our data
    assay_name = "logcounts", ## assay to use, we recommend logcounts [default]
    cellType_col = "consensus_annot", # column in colData with cell type info
)


#1vsAll markers
marker_stats_1vAll <- findMarkers_1vAll(
     sce = hashikawa_final_sce, # sce is the SingleCellExperiment with our data
     assay_name = "logcounts",
     cellType_col = "consensus_annot" # column in colData with cell type info
 )



## join the two marker_stats tables
marker_stats <- marker_stats_MeanRatio |>
    left_join(marker_stats_1vAll, by = join_by(gene, cellType.target))


#Save results

saveRDS(marker_stats_MeanRatio, file = paste0(new_data_path, '/marker_stats_MeanRatio.rds'))
saveRDS(marker_stats_1vAll , file = paste0(new_data_path, '/marker_stats_1vAll.rds'))
saveRDS(marker_stats, file = paste0(new_data_path, '/marker_stats_combo.rds'))


#And now generate metaMarkers too, based on the stim vs cntl labels from hashikawa, the only sample info we have
table(hashikawa_final_sce$stim)
hashikawa_cntl = hashikawa_final_sce[ , hashikawa_final_sce$stim == 'cntl' ]
hashikawa_stim = hashikawa_final_sce[ , hashikawa_final_sce$stim == 'stim' ]


#Also add in the Wallace data here too
colnames(colData(wallace_sce))
table(wallace_sce$putative_donor, wallace_sce$consensus_annot)
wallace_d1 = wallace_sce[ , wallace_sce$putative_donor == '160822' ]
wallace_d2 = wallace_sce[ , wallace_sce$putative_donor == '161102' ]
wallace_d3 = wallace_sce[ , wallace_sce$putative_donor == '161103' ]
wallace_d4 = wallace_sce[ , wallace_sce$putative_donor == '161105' ]

mouse_marker_list = list(
  cntl = compute_markers(assay(hashikawa_cntl, 'cpm'), hashikawa_cntl$consensus_annot),
  stim = compute_markers(assay(hashikawa_stim, 'cpm'), hashikawa_stim$consensus_annot),
  wallace_d1 = compute_markers(assay(wallace_d1, 'cpm'), wallace_d1$consensus_annot),
  wallace_d2 = compute_markers(assay(wallace_d2, 'cpm'), wallace_d2$consensus_annot),
  wallace_d3 = compute_markers(assay(wallace_d3, 'cpm'), wallace_d3$consensus_annot),
  wallace_d4 = compute_markers(assay(wallace_d4, 'cpm'), wallace_d4$consensus_annot)
  )


mouse_metaMarkers = make_meta_markers(mouse_marker_list , detailed_stats = TRUE)

#Save the metamarkers
export_meta_markers(mouse_metaMarkers, 
  paste0(new_data_path, '/mouse_meta_markers.csv'), 
  names(mouse_metaMarkers))




# Function to create condensed violin plots similar to Tasic et al. 2018 Fig 4C
# Y-axis is max-normalized expression across cell types
plot_violin_maxnorm <- function(sce, gene, assay_name = "logcounts", celltype_col, study_label, y_max = 6) {
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
    scale_y_continuous(limits = c(0, y_max), breaks = seq(0, y_max, by = ifelse(y_max <= 1, 1, 2))) +
    labs(x = NULL, y = "Max-normalized expression", title = sprintf('%s - %s', study_label, gene)) +
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


hashikawa_neuron = hashikawa_final_sce[ , hashikawa_final_sce$consensus_annot %in% c('MHb_A', 'MHb_B', 'LHb_A', 'LHb_B', 'LHb_C')]
wallace_neuron = wallace_sce[ , wallace_sce$consensus_annot %in% c('MHb_A', 'MHb_B', 'LHb_A', 'LHb_B', 'LHb_C')]

hashikawa_neuron$consensus_annot = factor(hashikawa_neuron$consensus_annot, levels = c('MHb_B', 'MHb_A', 'LHb_A', 'LHb_B', 'LHb_C'))
wallace_neuron$consensus_annot = factor(wallace_neuron$consensus_annot, levels = c('MHb_B', 'MHb_A', 'LHb_A', 'LHb_C'))

p1 = plot_violin_maxnorm(hashikawa_neuron, gene = "SLC12A5", 
assay_name = "logcounts", celltype_col = "consensus_annot", study_label = "Hashikawa")

p2 = plot_violin_maxnorm(wallace_neuron, gene = "SLC12A5", 
assay_name = "logcounts", celltype_col = "consensus_annot", study_label = "Wallace")

p1
p2

ggsave(plot = p1, filename = paste0(plot_path, '/hashikawa_SLC12A5_violin.pdf'), 
width = 6, height = 2, device = 'pdf')
ggsave(plot = p2, filename = paste0(plot_path, '/wallace_SLC12A5_violin.pdf'), 
width = 6, height = 2, device = 'pdf')

# Add the regression lines

#Adding regression lines to the violin plots
# Extract SLC12A5 expression and cell type for both datasets
# Yalcinbas
wallace_df <- data.frame(
  expression = as.numeric(logcounts(wallace_neuron["SLC12A5", ])),
  celltype = wallace_neuron$consensus_annot
)
wallace_df$celltype_numeric <- as.numeric(wallace_df$celltype)  # 1=MHb_B, 2=MHb_A, ..., 5=LHb_C

# Multiome
hashikawa_df <- data.frame(
  expression = as.numeric(logcounts(hashikawa_neuron["SLC12A5", ])),
  celltype = hashikawa_neuron$consensus_annot
)
hashikawa_df$celltype_numeric <- as.numeric(hashikawa_df$celltype)

# Fit linear trend models
wallace_lm <- lm(expression ~ celltype_numeric, data = wallace_df)
hashikawa_lm <- lm(expression ~ celltype_numeric, data = hashikawa_df)

summary(wallace_lm)
summary(hashikawa_lm)


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
  num_celltypes = max(df$celltype_numeric)
  pred_df <- data.frame(
    celltype_numeric = 1:num_celltypes,
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

wallace_trend <- add_trend_line(p2, wallace_lm, wallace_df)
wallace_trend
hashikawa_trend <- add_trend_line(p1, hashikawa_lm, hashikawa_df)
hashikawa_trend

ggsave(plot = hashikawa_trend, filename = paste0(plot_path, '/hashikawa_SLC12A5_violin_trend.pdf'), 
width = 6, height = 2, device = 'pdf')
ggsave(plot = wallace_trend, filename = paste0(plot_path, '/wallace_SLC12A5_violin_trend.pdf'), 
width = 6, height = 2, device = 'pdf')


#And other genes of interest
hashikawa_cnr1_p = plot_violin_maxnorm(hashikawa_neuron, "CNR1", assay_name = "logcounts",
celltype_col = "consensus_annot", study_label = "Hashikawa", y_max = 2)

hashikawa_oprm1_p = plot_violin_maxnorm(hashikawa_neuron, "OPRM1", assay_name = "logcounts", 
celltype_col = "consensus_annot", study_label = "Hashikawa", y_max = 1)

hashikawa_htr2c_p = plot_violin_maxnorm(hashikawa_neuron, "HTR2C", assay_name = "logcounts",
celltype_col = "consensus_annot", study_label = "Hashikawa", y_max = 4)

hashikawa_htr4_p = plot_violin_maxnorm(hashikawa_neuron, "HTR4", assay_name = "logcounts", 
celltype_col = "consensus_annot", study_label = "Hashikawa", y_max = 4)

hashikawa_htr7_p = plot_violin_maxnorm(hashikawa_neuron, "HTR7", assay_name = "logcounts", 
celltype_col = "consensus_annot", study_label = "Hashikawa", y_max = 4)

hashikawa_htr2a_p = plot_violin_maxnorm(hashikawa_neuron, "HTR2A", assay_name = "logcounts", 
celltype_col = "consensus_annot", study_label = "Hashikawa", y_max = 4)

#Gene not present
#hashikawa_htr1e_p = plot_violin_maxnorm(hashikawa_neuron, "HTR1E", assay_name = "logcounts", 
#celltype_col = "consensus_annot", study_label = "Hashikawa")

hashikawa_htra1_p = plot_violin_maxnorm(hashikawa_neuron, "HTRA1", assay_name = "logcounts", 
celltype_col = "consensus_annot", study_label = "Hashikawa", y_max = 4)

hashikawa_chrna3_p = plot_violin_maxnorm(hashikawa_neuron, "CHRNA3", assay_name = "logcounts", 
celltype_col = "consensus_annot", study_label = "Hashikawa", y_max = 4)

hashikawa_chrnb4_p = plot_violin_maxnorm(hashikawa_neuron, "CHRNB4", assay_name = "logcounts", 
celltype_col = "consensus_annot", study_label = "Hashikawa", y_max = 4)

hashikawa_chrna5_p = plot_violin_maxnorm(hashikawa_neuron, "CHRNA5", assay_name = "logcounts",
celltype_col = "consensus_annot", study_label = "Hashikawa")

#Wallace

wallace_cnr1_p = plot_violin_maxnorm(wallace_neuron, "CNR1", assay_name = "logcounts",
celltype_col = "consensus_annot", study_label = "Wallace", y_max = 2)

wallace_oprm1_p = plot_violin_maxnorm(wallace_neuron, "OPRM1", assay_name = "logcounts", 
celltype_col = "consensus_annot", study_label = "Wallace", y_max = 1)

wallace_htr2c_p = plot_violin_maxnorm(wallace_neuron, "HTR2C", assay_name = "logcounts",
celltype_col = "consensus_annot", study_label = "Wallace", y_max = 4)

wallace_htr4_p = plot_violin_maxnorm(wallace_neuron, "HTR4", assay_name = "logcounts", 
celltype_col = "consensus_annot", study_label = "Wallace", y_max = 4)

wallace_htr7_p = plot_violin_maxnorm(wallace_neuron, "HTR7", assay_name = "logcounts", 
celltype_col = "consensus_annot", study_label = "Wallace", y_max = 4)

wallace_htr2a_p = plot_violin_maxnorm(wallace_neuron, "HTR2A", assay_name = "logcounts", 
celltype_col = "consensus_annot", study_label = "Wallace", y_max = 4)

#Gene not present
#wallace_htr1e_p = plot_violin_maxnorm(wallace_neuron, "HTR1E", assay_name = "logcounts", 
#celltype_col = "consensus_annot", study_label = "Wallace")

wallace_htra1_p = plot_violin_maxnorm(wallace_neuron, "HTRA1", assay_name = "logcounts", 
celltype_col = "consensus_annot", study_label = "Wallace", y_max = 4)

wallace_chrna3_p = plot_violin_maxnorm(wallace_neuron, "CHRNA3", assay_name = "logcounts", 
celltype_col = "consensus_annot", study_label = "Wallace", y_max = 4)

wallace_chrnb4_p = plot_violin_maxnorm(wallace_neuron, "CHRNB4", assay_name = "logcounts", 
celltype_col = "consensus_annot", study_label = "Wallace", y_max = 4)

wallace_chrna5_p = plot_violin_maxnorm(wallace_neuron, "CHRNA5", assay_name = "logcounts",
celltype_col = "consensus_annot", study_label = "Wallace", y_max = 4)


ggsave(plot = wallace_cnr1_p, filename = file.path(plot_path, 'Wallace_CNR1_violin.pdf'), width = 6, height = 2)
ggsave(plot = hashikawa_cnr1_p, filename = file.path(plot_path, 'Hashikawa_CNR1_violin.pdf'), width = 6, height = 2)

ggsave(plot = wallace_oprm1_p, filename = file.path(plot_path, 'Wallace_OPRM1_violin.pdf'), width = 6, height = 2)
ggsave(plot = hashikawa_oprm1_p, filename = file.path(plot_path, 'Hashikawa_OPRM1_violin.pdf'), width = 6, height = 2)

ggsave(plot = wallace_htr2c_p, filename = file.path(plot_path, 'Wallace_HTR2C_violin.pdf'), width = 6, height = 2)
ggsave(plot = hashikawa_htr2c_p, filename = file.path(plot_path, 'Hashikawa_HTR2C_violin.pdf'), width = 6, height = 2)

ggsave(plot = wallace_htr4_p, filename = file.path(plot_path, 'Wallace_HTR4_violin.pdf'), width = 6, height = 2)
ggsave(plot = hashikawa_htr4_p, filename = file.path(plot_path, 'Hashikawa_HTR4_violin.pdf'), width = 6, height = 2)

ggsave(plot = wallace_htr7_p, filename = file.path(plot_path, 'Wallace_HTR7_violin.pdf'), width = 6, height = 2)
ggsave(plot = hashikawa_htr7_p, filename = file.path(plot_path, 'Hashikawa_HTR7_violin.pdf'), width = 6, height = 2)

ggsave(plot = wallace_htr2a_p, filename = file.path(plot_path, 'Wallace_HTR2A_violin.pdf'), width = 6, height = 2)
ggsave(plot = hashikawa_htr2a_p, filename = file.path(plot_path, 'Hashikawa_HTR2A_violin.pdf'), width = 6, height = 2)

ggsave(plot = wallace_htra1_p, filename = file.path(plot_path, 'Wallace_HTRA1_violin.pdf'), width = 6, height = 2)
ggsave(plot = hashikawa_htra1_p, filename = file.path(plot_path, 'Hashikawa_HTRA1_violin.pdf'), width = 6, height = 2)

ggsave(plot = wallace_chrna3_p, filename = file.path(plot_path, 'Wallace_CHRNA3_violin.pdf'), width = 6, height = 2)
ggsave(plot = hashikawa_chrna3_p, filename = file.path(plot_path, 'Hashikawa_CHRNA3_violin.pdf'), width = 6, height = 2)

ggsave(plot = wallace_chrnb4_p, filename = file.path(plot_path, 'Wallace_CHRNB4_violin.pdf'), width = 6, height = 2)
ggsave(plot = hashikawa_chrnb4_p, filename = file.path(plot_path, 'Hashikawa_CHRNB4_violin.pdf'), width = 6, height = 2)

ggsave(plot = wallace_chrna5_p, filename = file.path(plot_path, 'Wallace_CHRNA5_violin.pdf'), width = 6, height = 2)
ggsave(plot = hashikawa_chrna5_p, filename = file.path(plot_path, 'Hashikawa_CHRNA5_violin.pdf'), width = 6, height = 2)


session_info()
