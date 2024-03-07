########################################################################
## Aggregateexpression for conserved gene markers
## Authors. CSC
## Date. March 5th, 2024
## Last.Adaptation: xxx
##
## Input: Seurat integrated object with samples S1 and S2 after CCA correction
## Output:  
##
## NOTES:
## For slurm env: runsrun --x11 --pty --partition=interactive bash
########################################################################

library('Seurat')                                 # 4.9.9.9045 2023-05-17 [1] Github (satijalab/seurat@7d1094c)
## Additional required packages for aggregation
library('multtest')
library('metap')
## Additional required packages for customize clusters
library('scCustomize')
library('magrittr')
library('tidyverse')
## Packages to plot
library('ggplot2')
library('patchwork')
library('cowplot')
theme_set(theme_cowplot())

#library(SeuratDisk)                             
library(here)

here::here()

if (!packageVersion("Seurat")=='4.9.9.9060') {
    stop
    message('This pipeline was implemented with Seurat v5 and Signac v1.11+ ')
    message('You need the laterst Seurat v5 (‘4.9.9.9060’)')
    message('Current available repository on: https://satijalab.org/seurat/articles/install.html  ') }

# Check if processed_data directory exists, if not create it
if (!dir.exists(here("processed-data/06_pseudobulking/"))) {
    dir.create(here("processed-data/06_pseudobulking/"))
}
# Check if plot directory exists, if not create it
if (!dir.exists(here("plots/06_pseudobulking/"))) {
    dir.create(here("plots/06_pseudobulking/"))
}

#source(here("code/functions_custom", "remote_plot_functions.R"))    # Call to plot GEX assay
#source(here("code/functions_custom", "remote_filtering_functions.R"))   # Call functions to subset the Seurat object

########################    Initials ########################  

## select the count-mtx to merge (raw or normalized data)
count_mtx_type <- 'data_counts'      
#count_mtx_type <- 'norm_counts' 

## load pre-existing Seurat
get_seurat <- function(name) {

        #Ex. "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/04_merge_seurats/seurat.combined.data_counts_PCA_CCA.rds"
    sobj <- readRDS(name)
    # verification of the integration
    print(table(sobj$orig.ident))
    # S1_Hb_KDM S2_Hb_KDM
    # 8178      9816
    print(head(sobj, n=2))
    return(sobj)
    
}

# get QC violin plots
plot_violinQC <- function(seuratOBJ, sfeature, stitle) {
    p1 <- VlnPlot(object = seuratOBJ, features = sfeature, 
                  group.by = 'orig.ident', pt.size = 0) & geom_boxplot() &
        theme(legend.position = 'none',
              axis.text.x = element_text(angle=0, hjust=1, size=8),  #10
              axis.text.y = element_text(size=8), 
              axis.title.x = element_blank(),
              axis.title.y = element_blank()) #&
    #labs(title = "", x = 'Samples', y ="")
    ggtitle(stitle)
    return(p1)
}



# load pre-existing seurat objects
if (count_mtx_type=='data_counts') { s_sample <- 'seurat.combined.data_counts_PCA_CCA' } else { s_sample <- 'seurat.combined.norm_counts_PCA_CCA' }
rds_name <- here('processed-data/04_merge_seurats', paste0(s_sample, '.rds'))
# ~/seurat.combined.data_counts_PCA_CCA.rds"
SeuratOBJ <- get_seurat(rds_name)
SeuratOBJ2 <- SeuratOBJ
# An object of class Seurat 
# 36601 features across 17994 samples within 1 assay 
# Active assay: RNA (36601 features, 2000 variable features)
# 3 layers present: data, counts, scale.data
# 4 dimensional reductions calculated: pca, umap.unintegrated, integrated.cca, umap

# verification of the integration
table(SeuratOBJ$orig.ident)


# S1_Hb_KDM S2_Hb_KDM 
# 8178      9816 
head(colnames(SeuratOBJ))
tail(colnames(SeuratOBJ))
colnames(SeuratOBJ@meta.data)
# [1] "orig.ident"                 "atac_peak_region_fragments"
# [3] "atac_fragments"             "nCount_RNA"                
# [5] "nFeature_RNA"               "log10GenesPerUMI"          
# [7] "percent.mt"                 "percent.ribo"              
# [9] "MTRatio"                    "unintegrated_clusters"     
# [11] "seurat_clusters"            "RNA_snn_res.1"   



######################  Customize clusters  ######################


oldIdent <- levels(Idents(SeuratOBJ))
# [1] "0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15" "16" "17" "18"

newIdent <- paste("C", 0:18, sep = "_")
# [1] "C_0"  "C_1"  "C_2"  "C_3"  "C_4"  "C_5"  "C_6"  "C_7"  "C_8"  "C_9"  "C_10" "C_11" "C_12" "C_13" "C_14" "C_15"
# [17] "C_16" "C_17" "C_18

# count cells by clusters
table(Idents(SeuratOBJ))
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18 
# 4122 1625 1615 1330 1311 1250 1100 1026  821  787  683  565  564  515  273  163  136   75   33 

# rename clusters to make them more readable
# require scCustomize/Wrapper funtion to rename clusters
SeuratOBJ <- Rename_Clusters(SeuratOBJ, new_idents = newIdent,
                         meta_col_name = "seurat_clusters.renamed")
table(Idents(SeuratOBJ))
# C_0  C_1  C_2  C_3  C_4  C_5  C_6  C_7  C_8  C_9 C_10 C_11 C_12 C_13 C_14 C_15 C_16 C_17 C_18 
# 4122 1625 1615 1330 1311 1250 1100 1026  821  787  683  565  564  515  273  163  136   75   33 
#View(table(Idents(SeuratOBJ)))

head(SeuratOBJ)


# ## Trying to change orig.ident meta.data
# library(stringr)
# 
# ## Confirm how many samples do we have
# # suffixes <- str_extract(string = colnames(SeuratOBJ), pattern = "[:digit:]$")
# # unique(suffixes)
# 
# unique(SeuratOBJ@meta.data$orig.ident)
# 
# # Create dataframe by sample that contains matching orig.ident code
# meta_by_sample <- tibble::tribble(
#   ~orig.ident,  ~sample_name,
#   1, "S1_Hb",
#   2, "S2_Hb" 
# )
# 
# # Change orig.ident column to factor so that it can be joined later
# meta_by_sample$orig.ident <- as.factor(meta_by_sample$orig.ident)
# 
# # Pull existing meta data where samples are specified by orig.ident and remove everything but orig.ident
# OBJ_meta <- SeuratOBJ@meta.data %>% 
#   select(orig.ident) %>% 
#   rownames_to_column("barcodes")
# 
# # Use full join with object meta data in x position so that by sample meta dataframe is propagated across the by cell meta dataframe from the object.  And then remove orig.ident because it's already present in object meta data.
# full_new_meta <- full_join(x = OBJ_meta, y = meta_by_sample) %>% 
#   column_to_rownames("barcodes") %>% 
#   select(-orig.ident)
# 
# # Use AddMetaData to add new meta data to object
# OBJ <- AddMetaData(object = OBJ, metadata = full_new_meta)



######################  Identify conserved cell type markers ######################
## Implementation from: https://satijalab.org/seurat/articles/integration_introduction.html#identify-conserved-cell-type-markers 

## Run in an integrated Seurat
SeuratOBJ[["RNA"]] <- JoinLayers(SeuratOBJ[["RNA"]])

## unique(Idents(SeuratOBJ))
## Hb.markers <- FindConservedMarkers(SeuratOBJ, ident.1 = "Clust_0", grouping.var = "orig.ident", verbose = FALSE)
## head(nk.markers)

## To avoid issue when having few cells need to adjust the minimum number of cells
## For example, if there is a cluster "15" that has 0 cells, the function will skip that cluster with a warning (that's perfect). Also if the number of cells is between min.cells.groups (default = 3) and 0, an error is thrown and it stops working. That is why I previously remove from the Seurat Object the cells of the clusters with 3 or less cells for each condition/sample. 
few_cells_samples <- unique(SeuratOBJ@meta.data$orig.ident)
few_cells <- vector()

for (i in 1:length(few_cells_samples)) {   # remove cellstype w/ less than 3 cells in each sample/condition
  few_cells_tmp <- table(SeuratOBJ@meta.data$seurat_clusters.renamed[SeuratOBJ@meta.data$orig.ident == few_cells_samples[i]]) <= 3
  few_cells_tmp <- names(few_cells_tmp)[few_cells_tmp == "TRUE"]
  few_cells <- c(few_cells,few_cells_tmp)
}

# > few_cells
# [1] "18"

clusters <- sort(unique(SeuratOBJ@meta.data$seurat_clusters.renamed))
clusters <- clusters[clusters %!in% few_cells]  # need to check CSC


# ## Determine the number of clusters
## https://github.com/satijalab/seurat/issues/6076

# num_clusters <- max(as.numeric(as.character(
#   SeuratOBJ@meta.data$seurat_clusters)))
# 
# ## Cycle through each cluster finding the conserved markers**
# ## Store each dataframe of markers in the misc slot**
# for (i in 0:num_clusters) {
#   SeuratOBJ@misc$temp <- FindConservedMarkers(SeuratOBJ, ident.1 = i, grouping.var = "orig.ident", min.cells.group = 0)
#   names(gene.conditions@misc)[names(gene.conditions@misc)=="temp"] <-
#     paste0(names(gene.conditions), ".cluster_", i, ".markers")
# }



######################. Plot conserved cell type markers with Doplot() ######################

unique(Idents(SeuratOBJ))
markers.to.plot <- c("MMRN1", "HTR2C", "EPHA5", "GPR151", "POU4F1", 
                     "AC109466.1", "AC008415.1", "GPR149", "GNG8", "LINC01876", "TLL1", "CD24", "AC004594.1")
# DotPlot(SeuratOBJ, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "orig.ident") +
#   RotatedAxis()

DotPlot(SeuratOBJ, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) +
  RotatedAxis()

######################. Identify differential expressed genes across conditions ######################

## We use AggregateExpression() to aggregate cells of a similar type and condition together to create “pseudobulk” profiles

colnames(SeuratOBJ@meta.data)

aggregate_ifnb <- AggregateExpression(SeuratOBJ, group.by = c("seurat_clusters.renamed", "orig.ident"), return.seurat = TRUE)
markers.to.plot <- c("MMRN1", "HTR2C", "EPHA5", "GPR151", "POU4F1")
markers.to.plot <- c("AC109466.1", "AC008415.1", "GPR149", "GNG8")
markers.to.plot <- c("LINC01876", "TLL1", "CD24", "AC004594.1")
markers.to.plot <- c("HTR2C")

unique(Cells(SeuratOBJ))
#p1 <- CellScatter(aggregate_ifnb, "S1_Hb_KDM_A", "Cell2", highlight = genes.to.label)
#p2 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)

SeuratOBJ@reductions

# FeaturePlot(SeuratOBJ, features = genes.to.label , split.by = "orig.ident", max.cutoff = 3, 
#             cols = c("grey","red"), reduction = "integrated.cca")

# Run umap 
SeuratOBJ <- RunUMAP(SeuratOBJ, dims = 1:30, reduction = "integrated.cca")
SeuratOBJ@reductions

# Plot in umap features for LHb/MHb marker genes 
FeaturePlot(SeuratOB, features = markers.to.plot , split.by = "orig.ident", max.cutoff = 3, 
            cols = c("grey","red"), reduction = "umap")

# Plot Violin plots for the same LHb/MHb marker genes 
plots <- VlnPlot(SeuratOBJ, features = markers.to.plot, split.by = "orig.ident", group.by = "seurat_clusters",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)


DoHeatmap(
  SeuratOBJ,
  features = NULL,
  cells = NULL,
  group.by = "orig.ident",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = NULL,
  label = TRUE,
  size = 5.5,
  hjust = 0,
  vjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
)



png_file <- paste0(s_sample, '_integrated.cca_pca_heatmap.png')
png_name <- here('plots/04_merge_seurats', png_file)
ggsave(p1, filename = png_name, height = 5, width = 10)


rds_name <- here('processed-data/04_merge_seurats', paste0(s_sample, '_PCA_CCA_Harmony.rds'))
# .../seurat.combined.data_counts_PCA_CCA_Harmony.rds
saveRDS(SeuratOBJ, file = rds_name)
message('Seurat combined saved in ', rds_name)   


# INTEGRATION methods for Seurat V5:  https://satijalab.org/seurat/articles/seurat5_integration (Oct 31, 2023)
# https://satijalab.org/seurat/articles/integration_introduction.html (Nov 16, 2023)



############ Reproducibility information ####################

library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
#"2023-04-04 12:42:26 EDT"
proc.time()
options(width = 120)
session_info()

