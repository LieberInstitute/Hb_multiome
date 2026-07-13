#Testing out some lateral habenual specific dimension reduction and or clustering to see what
#variability might exist related to the described tonic vs burst firing patterns observed in the lateral habenula

library(SingleCellExperiment)
library(harmony)
library(scater)
library(scran)
library(here)
library(ggplot2)
library(dplyr)
library(qs2)


#Path to save any generated data
#new_data_path = here('processed-data', '05_03_annotation_adjustments', '18_LHb_ephys_subclust')
#Path to plot directory
plot_path = here('plots','05_03_annotation_adjustments', '18_LHb_ephys_subclust')

#if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)


multiome_path = here('processed-data', '05_03_annotation_adjustments','06_refined_annotations')
multiome_sce = qs_read(paste0(multiome_path, '/refined_annotation_multiomeHab_SCE.qs2'))

assay(multiome_sce, 'cpm') = MetaMarkers::convert_to_cpm(assay(multiome_sce, 'counts'))
colnames(colData(multiome_sce))

SingleCellExperiment::altExp(multiome_sce, "ATAC") <- NULL

#colors
source(here('code','05_03_annotation_adjustments','celltype_colors.R'))

go_path = here('processed-data','05_03_annotation_adjustments', '10_celltype_func_annot')
#GO terms
go_sets = readRDS(file.path(go_path, "go_human.rds"))

#Filter for genes present in the dataset
known_genes = rownames(multiome_sce)
go_sets = lapply(go_sets, function(gene_set) {gene_set[gene_set %in% known_genes]})

s1 = names(go_sets)[grepl('potassium',names(go_sets)) & grepl('channel',names(go_sets))]
s2 = names(go_sets)[grepl('sodium',names(go_sets)) & grepl('channel',names(go_sets))]
s3 = names(go_sets)[grepl('calcium',names(go_sets)) & grepl('channel',names(go_sets))]
s4 = names(go_sets)[grepl('chloride',names(go_sets)) & grepl('channel',names(go_sets))]
s5 = names(go_sets)[grepl('GABA',names(go_sets))]
s6 = names(go_sets)[grepl('gamma-aminobutyric acid',names(go_sets))]
s7 = names(go_sets)[grepl('glutamate',names(go_sets)) & grepl('receptor',names(go_sets))]

ion_go_names = c(s1, s2, s3, s4, s5, s6, s7)
ion_go_terms = go_sets[ion_go_names]
ion_go_terms = unique(unlist(unname(ion_go_terms)))

#Just dim reduce the lateral habenula neurons
lateral_sce_sub = multiome_sce[, multiome_sce$refined_mid_cluster %in% c('LHb_A', 'LHb_B', 'LHb_C', 'GABA_LHb_C.1', 'GABA_LHb_C.2')]
lateral_sce_sub = logNormCounts(lateral_sce_sub,name = "logcounts")

#And now subset for the ion channel and receptor genes
gene_index = rownames(lateral_sce_sub) %in% ion_go_terms
lateral_sce_sub = lateral_sce_sub[gene_index,]

#dec.lateral_sce_sub <- modelGeneVar(lateral_sce_sub, assay.type = "logcounts")
#hvgs <- getTopHVGs(dec.lateral_sce_sub, n = 2000)

#Use all genes in PCA
hvgs = rownames(lateral_sce_sub)

lateral_sce_sub <- runPCA(lateral_sce_sub,
    subset_row = hvgs,
    ncomponents = 30,
    name = "PCA"
)


#Run Harmony on the donors

## Run Harmony
message("Running Harmony - ", Sys.time())
lateral_sce_sub <- RunHarmony(lateral_sce_sub, group.by.vars = 'orig.ident', verbose = TRUE)



set.seed(123)
message("Running TSNE - ", Sys.time())
lateral_sce_sub <- runTSNE(lateral_sce_sub, dimred = "HARMONY", name = "TSNE.HARMONY")
colnames(reducedDim(lateral_sce_sub, "TSNE.HARMONY")) <- c("TSNE1", "TSNE2")

set.seed(123)
message("Running UMAP - ", Sys.time())
lateral_sce_sub <- runUMAP(lateral_sce_sub, dimred = "HARMONY", name = "UMAP.HARMONY")
colnames(reducedDim(lateral_sce_sub, "UMAP.HARMONY")) <- c("UMAP1", "UMAP2")


plotReducedDim(
  lateral_sce_sub,
  dimred = "UMAP.HARMONY", colour_by = "orig.ident"
)

plotReducedDim(
  lateral_sce_sub,
  dimred = "UMAP.HARMONY", colour_by = "refined_mid_cluster"
)

plotReducedDim(
  lateral_sce_sub,
  dimred = "TSNE.HARMONY", colour_by = "refined_mid_cluster"
)

plotReducedDim(
  lateral_sce_sub,
  dimred = "PCA", colour_by = "refined_mid_cluster"
)


plot_gene_on_umap <- function(sce, gene_name, reduction = "UMAP.HARMONY", assay_type = "logcounts") {
  
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

plot_gene_on_umap(lateral_sce_sub, "HCN1")
plot_gene_on_umap(lateral_sce_sub, "KCNJ10")

plot_gene_on_umap(lateral_sce_sub, "CACNA1G")
plot_gene_on_umap(lateral_sce_sub, "CACNA1I")

plot_gene_on_umap(lateral_sce_sub, "GRIN1")
plot_gene_on_umap(lateral_sce_sub, "GRIN2A")
plot_gene_on_umap(lateral_sce_sub, "GRIN2B")
plot_gene_on_umap(lateral_sce_sub, "GRIN2D")
plot_gene_on_umap(lateral_sce_sub, "GRIN3A")

plot_gene_on_umap(lateral_sce_sub, "SLC12A1")
plot_gene_on_umap(lateral_sce_sub, "SLC12A2")
plot_gene_on_umap(lateral_sce_sub, "SLC12A5")
plot_gene_on_umap(lateral_sce_sub, "ANO1")
plot_gene_on_umap(lateral_sce_sub, "ANO2")
