
#Computing doublet scores and checking the doublet proportions of each fine resolution cluster
#Drop any doublet enriched clusters
#Redo the dimensionality reduction
#Save for downstream analyses

library(SingleCellExperiment)
library(scDblFinder)
library(Seurat)
library(Signac)
library(qs2)
library(harmony)
library(dplyr)
library(ggplot2)
library(here)

here::here()

#Path for new data generated
new_data_path = here('processed-data', '05_5_drop_doublets','01_drop_doublets_and_reDimReduce')
#Path to plot directory
plot_path = here('plots', '05_5_drop_doublets', '01_drop_doublets_and_reDimReduce')

if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)


## Seurat object with the mid-level cluster annots and dim reductions saved
midSeurat_Dir <- here(
    "processed-data",
    "05_Clustering_ARCr",
    "22_add_mid_level_clustering"
)
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_HD.rds"

mid_file_name <- here(midSeurat_Dir, Seurat_base_name)

midSeurat = readRDS(mid_file_name)


#Run scDblFinder to compute doublet scores
sce <- scDblFinder(GetAssayData(midSeurat , slot="counts"), samples = midSeurat$orig.ident)
# port the resulting scores back to the Seurat object:
midSeurat$scDblFinder.score <- sce$scDblFinder.score
midSeurat$scDblFinder.class <- sce$scDblFinder.class

table(midSeurat$scDblFinder.class)
#singlet doublet 
#  50184    5332 

#Fine resolution clusters
cluster_order <- midSeurat@meta.data |>
  dplyr::count(cluster_ann, scDblFinder.class) |>
  dplyr::group_by(cluster_ann) |>
  dplyr::mutate(prop = n / sum(n)) |>
  dplyr::filter(scDblFinder.class == "doublet") |>
  dplyr::arrange(prop) |>
  dplyr::pull(cluster_ann)

midSeurat@meta.data |>
  dplyr::mutate(cluster_ann = factor(cluster_ann, levels = cluster_order)) |>
  ggplot(aes(x = cluster_ann, fill = scDblFinder.class)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Fine clusters", y = "Proportion", fill = "scDblFinder") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

midSeurat@meta.data |>
  dplyr::mutate(cluster_ann = factor(cluster_ann, levels = cluster_order)) |>
  ggplot(aes(x = cluster_ann, y = scDblFinder.score)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = "Fine clusters") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


#Check out current UMAP and the clusters we will drop
## =============================================================================
## Picked up Hex-color codes similar across cell-type

my_colors <- c(
    LHb = "#1f78b4",
    MHb = "#ad1d8c",
    Oligo = "#384a08",
    Astrocyte = "#532222", 
    OPC = "#829454",
    Microglia = "#141b02",
    Endo = "#d95f02",
    Inhib_Thal = "#9a9fe7",
    Excit_Thal = "#42467b",
    Thal = "#4d55b7"
)

## assign color gradients to mid resolution clusters based on Broad cell-types

# extract LHb and MHb clusters
cluster_levels <- levels(midSeurat)
LHb_clusters <- grep("LHb", cluster_levels, value = TRUE)
MHb_clusters <- grep("MHb", cluster_levels, value = TRUE)

# Create tonal gradients for LHb and MHb
LHb_colors <- colorspace::sequential_hcl(length(LHb_clusters), h = 210, c = 80, l = c(30, 80))
MHb_colors <- colorspace::sequential_hcl(length(MHb_clusters), h = 320, c = 80, l = c(30, 80))

# Build full cluster color map
my_colors_mid <- setNames(rep("#bdbdbd", length(cluster_levels)), cluster_levels)
my_colors_mid[LHb_clusters] <- LHb_colors
my_colors_mid[MHb_clusters] <- MHb_colors

# assign base color for other types from your existing palette
for (category in c("Oligo", "Astrocyte", "OPC", "Microglia", "Endo", "Inhib.Thal", "Excit.Thal", "Thal")) {
    matched <- grep(category, cluster_levels, value = TRUE)
    my_colors_mid[matched] <- my_colors[[gsub("\\.", "_", category)]]
}

## =============================================================================

plt1 <- DimPlot(midSeurat, 
                label = TRUE, 
                reduction = "wnn.umap",
                group.by = "mid_cluster",
                label.size = 3,
                cols = my_colors_mid) + 
    NoLegend() +
    labs(title = "WNN cell types (Mid-resolution)")

plt1

#Based on cluster proportions of doublets, these are the ones we would exclude
#c('C.37.Thal','C.39.Inhib.Thal','C.24.LHb.4','C.30.LHb.7','C.40.LHb.4')
midSeurat$doublet_exclude = midSeurat$mid_cluster
midSeurat$doublet_exclude[midSeurat$cluster_ann %in% c('C.37.Thal','C.39.Inhib.Thal','C.24.LHb.4','C.30.LHb.7','C.40.LHb.4')] = 'Exclude'

my_colors_mid["Exclude"] <- "#ff3b3bff"

#Verify original UMAP
plt2 <- DimPlot(midSeurat, 
                label = TRUE, 
                reduction = "wnn.umap",
                group.by = "doublet_exclude",
                label.size = 3,
                cols = my_colors_mid) + 
    NoLegend() +
    labs(title = "WNN cell types (Mid-resolution)")


plt2


ggsave(plt1, filename = 'original_mid_res_WNN_umap.pdf', path = plot_path, device = 'pdf',
width = 6, height = 5, useDingbats = FALSE)
ggsave(plt2, filename = 'doubletRed_mid_res_WNN_umap.pdf', path = plot_path, device = 'pdf',
width = 6, height = 5, useDingbats = FALSE)


#Drop the doublet enriched clusters
midSeurat = subset(midSeurat, doublet_exclude != 'Exclude')

table(midSeurat$cluster_ann)


#Redo essential steps, normalization, variable feature selection, scaling, integrating, and dimensionality reductions
#Just start from counts again to simplify, drop the existing dimensionality reductions.

drop_cols <- c("ATAC.weight", "RNA.weight", "doublet_exclude", "RNA_snn_res.1", 
"C.leiden", "C.leiden_atac", "C.leiden_wnn", "unintegrated_clusters" ,  "seurat_clusters" )

# Keep everything else
midSeurat@meta.data <- midSeurat@meta.data |>
  dplyr::select(-dplyr::any_of(drop_cols))


colnames(midSeurat[[]])

#Drop all the current reductions
# See current reductions
Reductions(midSeurat)

# Drop all dimensional reductions
midSeurat@reductions <- list()

# Confirm
Reductions(midSeurat)



###############################
#From here, going to be using the same code as in the original reductions. Some of these steps are spread out across scripts
#But mostly getting code from 03_pseudobulking/08_harmony_CR_ARCr.R
###############################

#Rescale the data and find variable RNA features
midSeurat <- FindVariableFeatures(midSeurat, selection.method = "vst") 

all.genes <- rownames(midSeurat)
midSeurat<- ScaleData(midSeurat, features = all.genes)

## Run PCA
midSeurat <- RunPCA(midSeurat)

## Re run UMAP on QCed dataset
midSeurat<- RunUMAP(midSeurat, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

p1 <- DimPlot(midSeurat, reduction = 'umap.unintegrated', group.by = "orig.ident") + ggtitle("UMAP on RNA")

p1.umap.rna.sample <- DimPlot(midSeurat, reduction = 'umap.unintegrated', group.by = "orig.ident",
                              label = T, label.size = 2.5, repel = T) + theme(legend.position="none") + ggtitle("UMAP on RNA QCed (sample)")

p1
p1.umap.rna.sample

#Re-prep the ATAC data

DefaultAssay(midSeurat) <- "ATAC"

# Run term frequency inverse document frequency (TF-IDF) normalization on a matrix.
midSeurat <- RunTFIDF(midSeurat,
                        method = 1,  # computes log(𝑇𝐹×𝐼𝐷𝐹).
                        scale.factor = 10000)
midSeurat <- FindTopFeatures(midSeurat,
                               min.cutoff = 'q5', # 95% most common features coverage as VariableFeatures
                               verbose = TRUE)
midSeurat <- RunSVD(midSeurat)

#Unintegrated UMAP on the ATAC data
midSeurat <- RunUMAP(midSeurat, dims = 2:20, reduction = "lsi", reduction.name = "umap.lsi.unintegrated")

p1.umap.atac.sample <- DimPlot(midSeurat, reduction = 'umap.lsi.unintegrated', group.by = "orig.ident",
                                label = T, label.size = 2.5, repel = T)  + theme(legend.position="none") + ggtitle("UMAP on ATAC QCed (sample)")

p1.umap.atac.sample



#Rerun harmony on the RNA
# dd/mm/yyyy. #Original seed used
set.seed(12022025)

# use k-means centroids initialization
midSeurat  <- midSeurat |>
  RunHarmony(group.by.vars = "orig.ident",
             reduction = "pca",
             assay.use = "RNA",
             reduction.save = "integrated.harmony",
             plot_convergence = TRUE,
             #nclust = 50,                     # Number of clusters in model. nclust=1 equivalent to simple linear regression
             #max.iter = 10,                   # One round of Harmony involves one clustering and one correction step
             #max.iter.cluster = 20,          # Maximum number of rounds to run clustering at each round of Harmony
             early_stop = T)
## rewrite harmony assay
Reductions(midSeurat )

p1.harm.rna.sample <- DimPlot(midSeurat , reduction = 'integrated.harmony', group.by = "orig.ident") + 
  ggtitle("Harmony on RNA QCed")
p1.harm.rna.sample

## Run UMAP on harmonized RNA data
midSeurat  <- RunUMAP(midSeurat , dims = 1:30, reduction = "integrated.harmony", reduction.name = "umap.integrated")

p1.umap.harm.rna.sample <- DimPlot(midSeurat , reduction = 'umap.integrated', group.by = "orig.ident") +
  ggtitle("UMAP on Harmony RNA QCed")
p1.umap.harm.rna.sample 

## Run TSNE
midSeurat  <- RunTSNE(midSeurat , dims = 1:30, reduction = "integrated.harmony", reduction.name = "tsne.integrated")

p1.tsne.harm.rna.sample <- DimPlot(midSeurat , reduction = 'tsne.integrated', group.by = "orig.ident") +
  ggtitle("TSNE on Harmony RNA QCed")
p1.tsne.harm.rna.sample 


#And now Harmony on the ATAC data
midSeurat <- midSeurat|> 
  RunHarmony(group.by.vars = "orig.ident", 
             reduction.save = "integrated.lsi.harmony",
             assay.use = "ATAC",
             reduction.use= 'lsi',
             plot_convergence = TRUE,
             #max.iter = 10,  # To avoid warning() message: Quick-TRANSfer stage steps exceeded maximum (= 2785100)  
             #                  left empty as it could need more or less lters to complete
             early_stop = T,
             project.dim = F) 

## Re-join layers after RNA integration
# Assays(SeuratOBJ.1)
# [1] "RNA"  "ATAC"
midSeurat[["RNA"]] <- JoinLayers(midSeurat[["RNA"]])

## more plots after correction for comparison purposes 
p1.harm.atac.sample <- DimPlot(midSeurat, reduction = 'integrated.lsi.harmony', group.by = "orig.ident") +
  ggtitle("Harmony on ATAC QCed")
p1.harm.atac.sample

## Run UMAP on harmonized RNA data
midSeurat <- RunUMAP(midSeurat, dims = 2:30, reduction = "integrated.lsi.harmony", reduction.name = "umap.lsi.integrated")

p1.umap.harm.atac.sample <- DimPlot(midSeurat, reduction = 'umap.lsi.integrated', group.by = "orig.ident") +
  ggtitle("UMAP on Harmony ATAC QCed")
p1.umap.harm.atac.sample 

## Run TSNE
midSeurat<- RunTSNE(midSeurat, dims = 2:30, reduction = "integrated.lsi.harmony", reduction.name = "tsne.lsi.integrated")

p1.tsne.harm.atac.sample <- DimPlot(midSeurat, reduction = 'tsne.lsi.integrated', group.by = "orig.ident") +
  ggtitle("TSNE on Harmony ATAC QCed")
p1.tsne.harm.atac.sample 



#Save plots and seurat object up till now
#unintegrated plots
p1.umap.rna.sample
p1.umap.atac.sample

ggsave(p1.umap.rna.sample, filename = 'unint_donor_RNA_umap.pdf', path = plot_path, device = 'pdf',
width = 6, height = 5, useDingbats = FALSE)
ggsave(p1.umap.atac.sample, filename = 'unint_donor_ATAC_umap.pdf', path = plot_path, device = 'pdf',
width = 6, height = 5, useDingbats = FALSE)

#Integrated RNA PCA, UMAP, TSNE
p1.harm.rna.sample
p1.umap.harm.rna.sample 
p1.tsne.harm.rna.sample 

ggsave(p1.harm.rna.sample, filename = 'RNA_Harmony_pca.pdf', path = plot_path, device = 'pdf',
width = 6, height = 5, useDingbats = FALSE)
ggsave(p1.umap.harm.rna.sample , filename = 'RNA_Harmony_umap.pdf', path = plot_path, device = 'pdf',
width = 6, height = 5, useDingbats = FALSE)
ggsave(p1.tsne.harm.rna.sample , filename = 'RNA_Harmony_tsne.pdf', path = plot_path, device = 'pdf',
width = 6, height = 5, useDingbats = FALSE)

#Integrated ATAC PCA, UMAP, TSNE
p1.harm.atac.sample
p1.umap.harm.atac.sample 
p1.tsne.harm.atac.sample 

ggsave(p1.harm.atac.sample, filename = 'ATAC_Harmony_pca.pdf', path = plot_path, device = 'pdf',
width = 6, height = 5, useDingbats = FALSE)
ggsave(p1.umap.harm.atac.sample , filename = 'ATAC_Harmony_umap.pdf', path = plot_path, device = 'pdf',
width = 6, height = 5, useDingbats = FALSE)
ggsave(p1.tsne.harm.atac.sample , filename = 'ATAC_Harmony_tsne.pdf', path = plot_path, device = 'pdf',
width = 6, height = 5, useDingbats = FALSE)

###########################################
#And now do the WNN approach using both assays
#
#This code comes from 05_Clustering_ARCr/01_clustering_std_method.R
###########################################

Reductions(midSeurat)

DefaultAssay(midSeurat) <- "RNA"

midSeurat <- FindMultiModalNeighbors(
  midSeurat,
  k.nn = 30,
  reduction.list = list("integrated.harmony", "integrated.lsi.harmony"),
  dims.list = list(1:30, 2:20)
)

midSeurat <- RunUMAP(
  midSeurat,
  n.neighbors = 30, 
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_"
)

plt3 <- DimPlot(
  midSeurat,
  reduction = "wnn.umap",
  group.by = 'mid_cluster', 
  label = TRUE,
  label.size = 2.5,
  cols = my_colors_mid
) +
  ggtitle("WNN cell types (Mid-resolution) after doublet removal")
plt3


plt4 <- DimPlot(
  midSeurat,
  reduction = "wnn.umap",
  group.by = 'orig.ident', 
  label = TRUE,
  label.size = 2.5
) +
  ggtitle("WNN cell types (Mid-resolution) after doublet removal")
plt4

ggsave(plt3, filename = 'doubletRemoved_mid_res_WNN_umap.pdf', path = plot_path, device = 'pdf',
width = 6, height = 5, useDingbats = FALSE)
ggsave(plt4, filename = 'doubletRemoved_donor_WNN_umap.pdf', path = plot_path, device = 'pdf',
width = 6, height = 5, useDingbats = FALSE)


#And now save the seurat object for downstream analyses
#Go with the qs2 version

#Save as seurat
qs_save(midSeurat, paste0(new_data_path, '/reprocessed_doubletRemoved_multiomeHab_seurat.qs2'))

#Save as SingleCellExperiment
# Keep both RNA + ATAC
sce <- as.SingleCellExperiment(
  midSeurat,
  assay = c("RNA", "ATAC")
)

# Re-add all Seurat reductions explicitly
for (red in Reductions(midSeurat)) {
  emb <- Embeddings(midSeurat, reduction = red)
  emb <- emb[colnames(sce), , drop = FALSE]
  reducedDim(sce, red) <- emb
}

# Check what was kept
assayNames(sce)
altExpNames(sce)
reducedDimNames(sce)

qs_save(sce, paste0(new_data_path, '/reprocessed_doubletRemoved_multiomeHab_SCE.qs2'))


