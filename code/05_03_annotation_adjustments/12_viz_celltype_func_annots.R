#exploring and plotting the gene set annotation results

library(MetaNeighbor)
library(MetaMarkers)
library(SingleCellExperiment)
library(qs2)
library(here)
library(sessioninfo)

message('Loading data...')

multiome_path = here('processed-data', '05_03_annotation_adjustments','06_refined_annotations')
multiome_sce = qs_read(paste0(multiome_path, '/refined_annotation_multiomeHab_SCE.qs2'))

assay(multiome_sce, 'cpm') = MetaMarkers::convert_to_cpm(assay(multiome_sce, 'counts'))

#Drop ATAC assay
altExps(multiome_sce) <- NULL
gc()

kegg_path = here('processed-data','05_03_annotation_adjustments', '11_kegg_term_annot')
go_path = here('processed-data','05_03_annotation_adjustments', '10_celltype_func_annot')

#GO and KEGG results
go_aurocs = read.table(here(go_path, 'multiomeHab_functional_aurocs.txt'))
kegg_aurocs = read.table(here(kegg_path, 'multiomeHab_KEGG_functional_aurocs.txt'))

#GO and KEGG terms
kegg_sets = readRDS(file.path(kegg_path, "kegg_human.rds"))
go_sets = readRDS(file.path(go_path, "go_human.rds"))

plot_path = here('plots','05_03_annotation_adjustments', '12_viz_celltype_func_annots')
if (!dir.exists(plot_path)) dir.create(plot_path)


#Set up kegg set stats
gs_size = sapply(kegg_sets, length)
kegg_aurocs_df = data.frame(kegg_term = rownames(kegg_aurocs), kegg_aurocs)
kegg_aurocs_df$average = rowMeans(kegg_aurocs)
kegg_aurocs_df$n_genes = gs_size[rownames(kegg_aurocs)]

kegg_aurocs_df


head(kegg_aurocs_df[order(kegg_aurocs_df$average, decreasing = TRUE),], 20)


plotDotPlot(dat = multiome_sce,
experiment_labels = multiome_sce$orig.ident,
celltype_labels = multiome_sce$refined_mid_cluster,
gene_set = kegg_sets[['Morphine addiction']])



#Set up go set stats
gs_size = sapply(go_sets, length)
go_aurocs_df = data.frame(go_term = rownames(go_aurocs), go_aurocs)
go_aurocs_df$average = rowMeans(go_aurocs)
go_aurocs_df$n_genes = gs_size[rownames(go_aurocs)]

head(go_aurocs_df[order(go_aurocs_df$average, decreasing = TRUE),], 20)

small_go_sets = go_aurocs_df[go_aurocs_df$n_genes < 20,]
head(small_go_sets[order(small_go_sets$average, decreasing = TRUE),], 20)


plotDotPlot(dat = multiome_sce,
experiment_labels = multiome_sce$orig.ident,
celltype_labels = multiome_sce$refined_mid_cluster,
gene_set = go_sets[['GO:0005001|transmembrane receptor protein tyrosine phosphatase activity|MF']])





