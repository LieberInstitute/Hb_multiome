#exploring and plotting the gene set annotation results

library(MetaNeighbor)
library(MetaMarkers)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(qs2)
library(here)
library(sessioninfo)

message('Loading data...')

multiome_path = here('processed-data', '05_03_annotation_adjustments','06_refined_annotations')
multiome_sce = qs_read(paste0(multiome_path, '/refined_annotation_multiomeHab_SCE.qs2'))

assay(multiome_sce, 'cpm') = MetaMarkers::convert_to_cpm(assay(multiome_sce, 'counts'))

#Filter for just the habenula cell-types
hab_celltypes = c('Inhib_LHb_4.1','Inhib_LHb_4.2','LHb.4','LHb.2.7','LHb.1.3.4','MHb.1','MHb.1.2','MHb.2','MHb.3', 'Excit.Thal', 'Inhib.Thal')
multiome_nonHab_sce = multiome_sce[ ,!multiome_sce$refined_mid_cluster %in% hab_celltypes]
multiome_sce = multiome_sce[ ,multiome_sce$refined_mid_cluster %in% hab_celltypes]


#Drop ATAC assay
altExps(multiome_sce) <- NULL
gc()

kegg_path = here('processed-data','05_03_annotation_adjustments', '11_kegg_term_annot')
go_path = here('processed-data','05_03_annotation_adjustments', '10_celltype_func_annot')

#GO and KEGG results
go_aurocs = read.table(here(go_path, 'multiomeHab_functional_aurocs.txt'))
kegg_aurocs = read.table(here(kegg_path, 'multiomeHab_KEGG_functional_aurocs.txt'))

kegg_nonHab_aurocs = read.table(here(kegg_path, 'multiomeNonHab_KEGG_functional_aurocs.txt'))
go_nonHab_aurocs = read.table(here(go_path, 'multiomeNonHab_GO_functional_aurocs.txt'))


#GO and KEGG terms
kegg_sets = readRDS(file.path(kegg_path, "kegg_human.rds"))
go_sets = readRDS(file.path(go_path, "go_human.rds"))

#Filter for genes present in the dataset
known_genes = rownames(multiome_sce)
kegg_sets = lapply(kegg_sets, function(gene_set) {gene_set[gene_set %in% known_genes]})
go_sets = lapply(go_sets, function(gene_set) {gene_set[gene_set %in% known_genes]})

plot_path = here('plots','05_03_annotation_adjustments', '13_viz_celltype_func_annots')
if (!dir.exists(plot_path)) dir.create(plot_path)


#Set up kegg set stats
gs_size = sapply(kegg_sets, length)
kegg_aurocs_df = data.frame(kegg_term = rownames(kegg_aurocs), kegg_aurocs)
kegg_aurocs_df$average = rowMeans(kegg_aurocs)
kegg_aurocs_df$n_genes = gs_size[rownames(kegg_aurocs)]
kegg_aurocs_df$ontology = 'KEGG'

kegg_aurocs_df

#Top KEGG gene sets 
head(kegg_aurocs_df[order(kegg_aurocs_df$average, decreasing = TRUE),], 10)
#Worst KEGG gene sets 
head(kegg_aurocs_df[order(kegg_aurocs_df$average, decreasing = FALSE),], 10)

small_kegg_sets = kegg_aurocs_df[kegg_aurocs_df$n_genes <= 50,]
head(small_kegg_sets[order(small_kegg_sets$average, decreasing = TRUE),], 20) 

large_kegg_sets = kegg_aurocs_df[kegg_aurocs_df$n_genes >= 50 & kegg_aurocs_df$n_genes <= 100,]
head(large_kegg_sets[order(large_kegg_sets$average, decreasing = TRUE),], 20) 

set_of_interest = 'Nicotine addiction'
plotDotPlot(dat = multiome_sce,
experiment_labels = multiome_sce$orig.ident,
celltype_labels = multiome_sce$refined_mid_cluster,
gene_set = kegg_sets[[set_of_interest]]) + ggtitle(set_of_interest)


#Highlight the high scores of the neuropsychiatric terms
terms_to_label = c('Amphetamine addiction','Long-term depression','Nicotine addiction','Morphine addiction','Cocaine addiction')
kegg_aurocs_df$label_terms = ifelse(kegg_aurocs_df$kegg_term %in% terms_to_label, kegg_aurocs_df$kegg_term, NA)
kegg_aurocs_df$point_color = ifelse(kegg_aurocs_df$kegg_term %in% terms_to_label, 'red', 'black')
pos <- position_jitter(width = 0.1, seed = 2)

ggplot(kegg_aurocs_df, aes(x = ontology, y = average, label = label_terms)) + 
  geom_boxplot(fill = 'grey') + 
  ggrepel::geom_label_repel(min.segment.length = 0, max.overlaps = Inf,
    box.padding = 1, position = pos,
    xlim  = c(NA ,1) ) +
  geom_point(position = pos, aes(color = point_color), show.legend = FALSE) + 
  scale_color_manual(values = c('red' = 'red', 'black' = 'black')) +
  theme_bw() + ylab('Predictive strength (avg. AUROC)') +
  ggtitle('KEGG gene sets predicting Habenula cell-types')



#Compare performance to average expression of the gene set and gene set size
ggplot(kegg_aurocs_df, aes(x = n_genes, y = average)) +
  geom_point() + theme_bw() + xlab('Number of genes in gene set') + ylab('Predictive strength (avg. AUROC)') +
  ggtitle('KEGG gene set size vs. performance')

#And now gene set average expression
kegg_sets
all_genes = rownames(multiome_sce)

onehot_mat <- do.call(cbind, lapply(kegg_sets, \(gene_set) {
  as.integer(all_genes %in% gene_set)
}))

rownames(onehot_mat) <- all_genes
colnames(onehot_mat) <- names(kegg_sets)

#Gives you the sum CPMs of each gene set per cell
sum_gene_cpms = t(assay(multiome_sce, 'cpm')) %*% onehot_mat
#Gives the average CPM per gene set across all cells
gene_set_avg = colSums(sum_gene_cpms) / dim(sum_gene_cpms)[1]
names(gene_set_avg) = colnames(sum_gene_cpms)


index = match(kegg_aurocs_df$kegg_term, names(gene_set_avg))
kegg_aurocs_df$avg_expression = gene_set_avg[index]

ggplot(kegg_aurocs_df, aes(x = log10(avg_expression), y = average)) +
  geom_point() + theme_bw() + xlab('log10 Avg. CPM of gene set') + ylab('Predictive strength (avg. AUROC)') +
  ggtitle('KEGG gene set expression vs. performance')

#What about average expression of gene and then average expression of gene set?
gene_avg = rowMeans(assay(multiome_sce, 'cpm'))

geneset_avg_sum = t(onehot_mat) %*% gene_avg
gene_avg_geneset_avg = geneset_avg_sum / colSums(onehot_mat)

index = match(kegg_aurocs_df$kegg_term, rownames(gene_avg_geneset_avg))
kegg_aurocs_df$gene_avg_geneset_avg = gene_avg_geneset_avg[index]

ggplot(kegg_aurocs_df, aes(x = log10(gene_avg_geneset_avg), y = average)) +
  geom_point() + theme_bw() + xlab('log10 Avg. Gene -> Avg. Gene set') + ylab('Predictive strength (avg. AUROC)') +
  ggtitle('KEGG gene set expression vs. performance')


#####################
# Kegg results in the non-neurons
#####################

#Set up kegg set stats
gs_size = sapply(kegg_sets, length)
kegg_nonHab_aurocs_df = data.frame(kegg_term = rownames(kegg_nonHab_aurocs ), kegg_nonHab_aurocs )
kegg_nonHab_aurocs_df$average = rowMeans(kegg_nonHab_aurocs )
kegg_nonHab_aurocs_df$n_genes = gs_size[rownames(kegg_nonHab_aurocs )]
kegg_nonHab_aurocs_df$ontology = 'KEGG'

#Top KEGG gene sets 
head(kegg_nonHab_aurocs_df[order(kegg_nonHab_aurocs_df$average, decreasing = TRUE),], 10)
#Worst KEGG gene sets 
head(kegg_nonHab_aurocs_df[order(kegg_nonHab_aurocs_df$average, decreasing = FALSE),], 10)

small_kegg_sets = kegg_nonHab_aurocs_df[kegg_nonHab_aurocs_df$n_genes <= 50,]
head(small_kegg_sets[order(small_kegg_sets$average, decreasing = TRUE),], 20) 

large_kegg_sets = kegg_nonHab_aurocs_df[kegg_nonHab_aurocs_df$n_genes >= 50 & kegg_nonHab_aurocs_df$n_genes <= 100,]
head(large_kegg_sets[order(large_kegg_sets$average, decreasing = TRUE),], 20) 

set_of_interest = 'Nicotine addiction'
plotDotPlot(dat = multiome_nonHab_sce,
experiment_labels = multiome_nonHab_sce$orig.ident,
celltype_labels = multiome_nonHab_sce$refined_mid_cluster,
gene_set = kegg_sets[[set_of_interest]]) + ggtitle(set_of_interest)


#Highlight the high scores of the neuropsychiatric terms
terms_to_label = c('Amphetamine addiction','Long-term depression','Nicotine addiction','Morphine addiction','Cocaine addiction')
kegg_nonHab_aurocs_df$label_terms = ifelse(kegg_nonHab_aurocs_df$kegg_term %in% terms_to_label, kegg_nonHab_aurocs_df$kegg_term, NA)
kegg_nonHab_aurocs_df$point_color = ifelse(kegg_nonHab_aurocs_df$kegg_term %in% terms_to_label, 'red', 'black')
pos <- position_jitter(width = 0.1, seed = 2)

ggplot(kegg_nonHab_aurocs_df, aes(x = ontology, y = average, label = label_terms)) + 
  geom_boxplot(fill = 'grey') + 
  ggrepel::geom_label_repel(min.segment.length = 0, max.overlaps = Inf,
    box.padding = 1, position = pos,
    xlim  = c(NA ,1) ) +
  geom_point(position = pos, aes(color = point_color), show.legend = FALSE) + 
  scale_color_manual(values = c('red' = 'red', 'black' = 'black')) +
  theme_bw() + ylab('Predictive strength (avg. AUROC)') +
  ggtitle('KEGG gene sets predicting non-Habenula cell-types')


#Compare between the Hab and non-Hab neurons
table(rownames(kegg_nonHab_aurocs_df) == rownames(kegg_aurocs_df))

kegg_aurocs_df$nonHab_average = kegg_nonHab_aurocs_df$average

ggplot(kegg_aurocs_df, aes(x = average, y = nonHab_average, label = label_terms)) +
  geom_point() + theme_bw() + xlab('Habenula AUROC') + ylab('non-Habenula AUROC') + 
  geom_abline( intercept = 0, slope = 1, color = 'red') +
  coord_cartesian(clip = "off") +
  ggrepel::geom_label_repel(min.segment.length = 0, max.overlaps = Inf,
    box.padding = 1,
    fill = "white", xlim = c(-Inf, Inf), ylim = c(-Inf, Inf)) +
  ggtitle('KEGG gene set performance in Hab vs. non-Hab neurons')

##########################
#Go set results
##########################



#Set up go set stats
gs_size = sapply(go_sets, length)
go_aurocs_df = data.frame(go_term = rownames(go_aurocs), go_aurocs)
go_aurocs_df$average = rowMeans(go_aurocs)
go_aurocs_df$n_genes = gs_size[rownames(go_aurocs)]
go_aurocs_df$ontology = 'GO'


head(go_aurocs_df[order(go_aurocs_df$average, decreasing = TRUE),], 20) %>% View()
#Poor performing gene sets
head(go_aurocs_df[order(go_aurocs_df$average, decreasing = FALSE),], 20) %>% View()

small_go_sets = go_aurocs_df[go_aurocs_df$n_genes <= 30,]
head(small_go_sets[order(small_go_sets$average, decreasing = TRUE),], 20)

set_of_interest = "GO:0016917|GABA receptor activity|MF"
plotDotPlot(dat = multiome_sce,
experiment_labels = multiome_sce$orig.ident,
celltype_labels = multiome_sce$refined_mid_cluster,
gene_set = go_sets[[set_of_interest]]) + ggtitle(set_of_interest)



#Highlight the high scores of the top terms
terms_to_label = go_aurocs_df %>% filter(n_genes <= 30) %>% slice_max(order_by = average, n = 10) %>% pull(go_term)
go_aurocs_df$label_terms = ifelse(go_aurocs_df$go_term %in% terms_to_label, go_aurocs_df$go_term, NA)
go_aurocs_df$point_color = ifelse(go_aurocs_df$go_term %in% terms_to_label, 'red', 'black')
pos <- position_jitter(width = 0.1, seed = 2)

ggplot(go_aurocs_df, aes(x = ontology, y = average, label = label_terms)) + 
  geom_boxplot(fill = 'grey') + 
  ggrepel::geom_label_repel(min.segment.length = 0, max.overlaps = Inf,
    box.padding = 1, position = pos,
    xlim  = c(NA ,1) ) +
  geom_point(position = pos, aes(color = point_color), show.legend = FALSE) + 
  scale_color_manual(values = c('red' = 'red', 'black' = 'black')) +
  theme_bw() + ylab('Predictive strength (avg. AUROC)') +
  ggtitle('go gene sets predicting Habenula cell-types')




#Compare performance to average expression of the gene set and gene set size
ggplot(go_aurocs_df, aes(x = n_genes, y = average)) +
  geom_point() + theme_bw() + xlab('Number of genes in gene set') + ylab('Predictive strength (avg. AUROC)') +
  ggtitle('GO gene set size vs. performance')

#And now gene set average expression
all_genes = rownames(multiome_sce)

onehot_mat <- do.call(cbind, lapply(go_sets, \(gene_set) {
  as.integer(all_genes %in% gene_set)
}))

rownames(onehot_mat) <- all_genes
colnames(onehot_mat) <- names(go_sets)

#Gives you the sum CPMs of each gene set per cell
#sum_gene_cpms = t(assay(multiome_sce, 'cpm')) %*% onehot_mat
#Gives the average CPM per gene set across all cells
#gene_set_avg = colSums(sum_gene_cpms) / dim(sum_gene_cpms)[1]
#names(gene_set_avg) = colnames(sum_gene_cpms)


#index = match(go_aurocs_df$go_term, names(gene_set_avg))
#go_aurocs_df$avg_expression = gene_set_avg[index]

#ggplot(go_aurocs_df, aes(x = log10(avg_expression), y = average)) +
#  geom_point() + theme_bw() + xlab('log10 Avg. CPM of gene set') + ylab('Predictive strength (avg. AUROC)') +
#  ggtitle('GO gene set expression vs. performance')

#What about average expression of gene and then average expression of gene set?
gene_avg = rowMeans(assay(multiome_sce, 'cpm'))

geneset_avg_sum = t(onehot_mat) %*% gene_avg
gene_avg_geneset_avg = geneset_avg_sum / colSums(onehot_mat)

index = match(go_aurocs_df$go_term, rownames(gene_avg_geneset_avg))
go_aurocs_df$gene_avg_geneset_avg = gene_avg_geneset_avg[index]

ggplot(go_aurocs_df, aes(x = log10(gene_avg_geneset_avg), y = average)) +
  geom_point() + theme_bw() + xlab('log10 Avg. Gene -> Avg. Gene set') + ylab('Predictive strength (avg. AUROC)') +
  ggtitle('GO gene set expression vs. performance')



#####################
# go results in the non-neurons
#####################

#Set up go set stats
gs_size = sapply(go_sets, length)
go_nonHab_aurocs_df = data.frame(go_term = rownames(go_nonHab_aurocs ), go_nonHab_aurocs )
go_nonHab_aurocs_df$average = rowMeans(go_nonHab_aurocs )
go_nonHab_aurocs_df$n_genes = gs_size[rownames(go_nonHab_aurocs )]
go_nonHab_aurocs_df$ontology = 'go'

#Top go gene sets 
head(go_nonHab_aurocs_df[order(go_nonHab_aurocs_df$average, decreasing = TRUE),], 10)
#Worst go gene sets 
head(go_nonHab_aurocs_df[order(go_nonHab_aurocs_df$average, decreasing = FALSE),], 10)

small_go_sets = go_nonHab_aurocs_df[go_nonHab_aurocs_df$n_genes <= 30,]
head(small_go_sets[order(small_go_sets$average, decreasing = TRUE),], 20) 


set_of_interest = 'GO:1904862|inhibitory synapse assembly|BP'
plotDotPlot(dat = multiome_nonHab_sce,
experiment_labels = multiome_nonHab_sce$orig.ident,
celltype_labels = multiome_nonHab_sce$refined_mid_cluster,
gene_set = go_sets[[set_of_interest]]) + ggtitle(set_of_interest)


#Compare between the Hab and non-Hab neurons
table(rownames(go_nonHab_aurocs_df) == rownames(go_aurocs_df))

go_aurocs_df$nonHab_average = go_nonHab_aurocs_df$average

ggplot(go_aurocs_df, aes(x = average, y = nonHab_average, color = point_color)) +
  geom_point() + theme_bw() + xlab('Habenula AUROC') + ylab('non-Habenula AUROC') + 
  geom_abline( intercept = 0, slope = 1, color = 'red') +
  ggtitle('GO gene set performance in Hab vs. non-Hab neurons')

