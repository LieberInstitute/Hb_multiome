
########################################################################
## Build tree clustering plot 
## INPUT: Seurat object with cluster information
## 
## OUPUT: Tree clustering plot
## 
## Date: March 29th, 2024 
########################################################################

# load libraries
library('Seurat')                             
library("here")
library("ape")
library("ggplot2")
library("bluster")
library("ComplexHeatmap")

library(here)

here::here()



############        Initials      ############

## read directory with Seurat objects
if (count_mtx_type=='data_counts') { 
  Seurat_base_name <- 'seurat.combined.data_counts_PCA'
} else { 
  Seurat_base_name <- 'seurat.combined.norm_counts_PCA' 
}
## Compose Seurat object name processed before
rds_name1 <- here('processed-data/02_merge_seurats', paste0(Seurat_base_name, '_CCA.rds'))
rds_name2 <- here('processed-data/02_merge_seurats', paste0(Seurat_base_name, '_Harmony.rds'))
rds_name1
rds_name2

# Check if plots directory exists, if not create it
if (!dir.exists(here("plots/04_DiffExpr_Clustering_seurat/"))) {
  dir.create(here("plots/04_DiffExpr_Clustering_seurat/"))
}
rds_path_out <- here('plots/04_DiffExpr_Clustering_seurat')
rds_path_out



############ Load Seurat to plot  ############


if (length(rds_name1) < 1 || length(rds_name2) < 1) { stop('\nNone seurats available.') }

message('Losing Seurats with clustering: \n', rds_name1, ' and \n', rds_name2)


SeuratOBJ_CCA <-readRDS(rds_name1) #CCA
SeuratOBJ_Harm <-readRDS(rds_name2) # Harmony
#str(SeuratOBJ_CCA)

head(SeuratOBJ_CCA@reductions$pca) 

head(SeuratOBJ_CCA@meta.data)

# # Set the top 6 marker genes from top50-ratio Habenula marker genes
# Mhb_Top50r <- c("CHAT", "LINC01307", "NEUROD1", "CHRNB4", "LINC02143", "AC114321.1")
# MHb <- FetchData(SeuratOBJ_CCA,c("ident","PC_1","nFeature_RNA","orig.ident", Mhb_Top50r))
# Lhb_Top50r <- c("HTR4", "BVES", "NRP1", "HTR2C", "CDH4", "COL25A1")
# LHb <- FetchData(SeuratOBJ_CCA,c("ident","PC_1","nFeature_RNA","orig.ident", Lhb_Top50r))
# 
# head(MHb,5)
# head(LHb,5)
# 
# VlnPlot(SeuratOBJ_CCA, Mhb_Top50r)
# VlnPlot(SeuratOBJ_CCA, Lhb_Top50r)


# ===============================================
# Tree Clustering for CCA


if (requireNamespace("ape", quietly = TRUE)) {
  SeuratOBJ_CCA_tree <- BuildClusterTree(SeuratOBJ_CCA)
  p1 <- PlotClusterTree(SeuratOBJ_CCA_tree, reduction = "pca") 
  p1
}


if (requireNamespace("ape", quietly = TRUE)) {
  SeuratOBJ_CCA_tree <- BuildClusterTree(SeuratOBJ_CCA, 
                                      reduction = "integrated.cca")
  p1 <- PlotClusterTree(SeuratOBJ_CCA_tree)
  p1
}

# set the marker genes identified as feature property to build the tree 
CCA_Hb = c("MMRN1", "GPR151", "POU4F1", "LINC01876", "CD24", "AC004594.1")
if (requireNamespace("ape", quietly = TRUE)) {
  SeuratOBJ_CCA_tree <- BuildClusterTree(SeuratOBJ_CCA, 
                                      reduction = "integrated.cca",
                                      features = CCA_Hb)
  p1 <- PlotClusterTree(SeuratOBJ_CCA_tree) 
  p1
}

# par(mfrow=c(1,2))
# VlnPlot(SeuratOBJ_CCA, CCA_Hb)



# ===============================================
# Tree Clustering for Harmony


if (requireNamespace("ape", quietly = TRUE)) {
  SeuratOBJ_Harm_tree <- BuildClusterTree(SeuratOBJ_Harm)
  p1 <- PlotClusterTree(SeuratOBJ_Harm_tree, 
                        reduction = "pca")
  p1
}

# SeuratOBJ_Harm@reductions$integrated.harmony
if (requireNamespace("ape", quietly = TRUE)) {
  SeuratOBJ_Harm_tree <- BuildClusterTree(SeuratOBJ_Harm, 
                                      reduction = "integrated.harmony")
  p1 <- PlotClusterTree(SeuratOBJ_Harm_tree)
  p1
}

par(mfrow=c(1,1))

# set the marker genes identified as feature property to build the tree 
Harm_Hb = c("MMRN1", "GPR151", "POU4F1", "EPHA5", "TLL1", "AC109466.1", "AC008415.1", "GPR149", "GNG8", "NEUROD1", "RASGPR1", "SLC5A7", "CHRNB3", "SCUBE1", "LINC02143", "CD24")
if (requireNamespace("ape", quietly = TRUE)) {
  SeuratOBJ_Harm_tree <- BuildClusterTree(SeuratOBJ_Harm, 
                                      reduction = "integrated.harmony",
                                      features = Harm_Hb)
  p1 <- PlotClusterTree(SeuratOBJ_Harm_tree) 
  p1
}


# par(mfrow=c(1,2))
# VlnPlot(SeuratOBJ_CCA, Harm_Hb)


# ===============================================

# we can look at the cells in the clusters and compare the average expression
# For example, I compare cluster g14 and cluster g15


CCA.avg <- AverageExpression(SeuratOBJ_CCA)
head(CCA.avg$RNA)
par(mfrow=c(2,3))

p1<- plot(CCA.avg$RNA[,"g18"], CCA.avg$RNA[,"g14"], pch=16,cex=0.8)
p2<- plot(CCA.avg$RNA[,"g18"], CCA.avg$RNA[,"g15"], pch=16,cex=0.8)
p3<- plot(CCA.avg$RNA[,"g18"], CCA.avg$RNA[,"g8"], pch=16,cex=0.8)
p4<- plot(CCA.avg$RNA[,"g14"], CCA.avg$RNA[,"g15"], pch=16,cex=0.8)
p5<- plot(CCA.avg$RNA[,"g14"], CCA.avg$RNA[,"g8"], pch=16,cex=0.8)
p6<- plot(CCA.avg$RNA[,"g15"], CCA.avg$RNA[,"g8"], pch=16,cex=0.8)


H.avg <- AverageExpression(SeuratOBJ_Harm)
head(H.avg$RNA)
par(mfrow=c(2,3))

p1<- plot(H.avg$RNA[,"g15"], H.avg$RNA[,"g14"], pch=16,cex=0.8)
p2<- plot(H.avg$RNA[,"g15"], H.avg$RNA[,"g12"], pch=16,cex=0.8)
p3<- plot(H.avg$RNA[,"g15"], H.avg$RNA[,"g6"], pch=16,cex=0.8)
p4<- plot(H.avg$RNA[,"g14"], H.avg$RNA[,"g12"], pch=16,cex=0.8)
p5<- plot(H.avg$RNA[,"g14"], H.avg$RNA[,"g6"], pch=16,cex=0.8)
p6<- plot(H.avg$RNA[,"g12"], H.avg$RNA[,"g6"], pch=16,cex=0.8)



message('Done!')


############ Reproducibility information ####################

library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
proc.time()
options(width = 120)
session_info()