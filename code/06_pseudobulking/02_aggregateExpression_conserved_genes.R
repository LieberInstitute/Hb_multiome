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

library(Seurat)                                 # 4.9.9.9045 2023-05-17 [1] Github (satijalab/seurat@7d1094c)
#options(Seurat.object.assay.version = 'v5')    # To use new Seurat v5: Please run: options(Seurat.object.assay.version = 'v5')
#library(Signac)                                 # 1.9.0.9000 2023-05-08 [1] Github (stuart-lab/signac@cf31022)
#library(EnsDb.Hsapiens.v86)
#library(BSgenome.Hsapiens.UCSC.hg38)
options(tidyverse.quiet = TRUE)
library(tidyverse)

library(ggplot2)
library(patchwork)
library(cowplot)
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

source(here("code/functions_custom", "remote_plot_functions.R"))    # Call to plot GEX assay
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
SeuratOBJ <- get_seurat(rds_name)


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
# [1] "0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14"
# [16] "15" "16" "17" "18"

newIdent <- paste("Clust", 0:18, sep = "_")
# [1] "Clust_0"  "Clust_1"  "Clust_2"  "Clust_3"  "Clust_4"  "Clust_5" 
# [7] "Clust_6"  "Clust_7"  "Clust_8"  "Clust_9"  "Clust_10" "Clust_11"
# [13] "Clust_12" "Clust_13" "Clust_14" "Clust_15" "Clust_16" "Clust_17"
# [19] "Clust_18"

# count cells by clusters
table(Idents(SeuratOBJ))
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 4122 1625 1615 1330 1311 1250 1100 1026  821  787  683  565  564  515  273  163 
# 16   17   18 
# 136   75   33 

# rename clusters to make them more readable
# require scCustomize/Wrapper funtion to rename clusters
if (FALSE) {
  SeuratOBJ <- Rename_Clusters(SeuratOBJ, new_idents = newIdent,
                         meta_col_name = "seurat_clusters.renamed")
}
colnames(SeuratOBJ@meta.data)
table(Idents(SeuratOBJ))
# Clust_0  Clust_1  Clust_2  Clust_3  Clust_4  Clust_5  Clust_6  Clust_7 
# 4122     1625     1615     1330     1311     1250     1100     1026 
# Clust_8  Clust_9 Clust_10 Clust_11 Clust_12 Clust_13 Clust_14 Clust_15 
# 821      787      683      565      564      515      273      163 
# Clust_16 Clust_17 Clust_18 
# 136       75       33 


# NEEDS TO BE FIXED AND SET ORDER CORRECTLY
# Idents(SeuratOBJ) <- factor(Idents(SeuratOBJ), levels = c("Oli-1", "Unknown-1", "Unknown-2", "Oli-Prec", "Unknown-3", 
#                                                           "Unknown-4", "LHb-1", "Ast-1", "Ast-2", "LHb-2",
#                                                           "Neu", "Medio_Thal", "Oli-2", "Unknown-5", "LHb-Mhb-1", 
#                                                           "MHb", "Endo", "Unknown-6", "LHb-Mhb-2"))
# 
# Idents(SeuratOBJ) <- factor(Idents(SeuratOBJ), levels = c("0", "1", "2", "3", "4", 
#                                                           "5", "6", "7", "8", "9",
#                                                           "10", "11", "12", "13", "14", 
#                                                           "15", "16", "17", "18"))




######################. Identify conserved cell type markers ######################



######################. Plot conserved cell type markers with Doplot() ######################


unique(Idents(SeuratOBJ))
markers.to.plot <- c("MMRN1", "HTR2C", "EPHA5", "GPR151", "POU4F1", 
                     "AC109466.1", "AC008415.1", "GPR149", "GNG8", "LINC01876", "TLL1", "CD24", "AC004594.1")
DotPlot(SeuratOBJ, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "orig.ident") +
  RotatedAxis()


######################. Identify differential expressed genes across conditions ######################

## We use AggregateExpression() to aggregate cells of a similar type and condition together to create “pseudobulk” profiles

colnames(SeuratOBJ@meta.data)

aggregate_ifnb <- AggregateExpression(SeuratOBJ, group.by = c("seurat_clusters", "orig.ident"), return.seurat = TRUE)
genes.to.label =  c("MMRN1", "HTR2C", "EPHA5", "GPR151", "POU4F1")

p1 <- CellScatter(aggregate_ifnb, "Cell1", "Cell2", highlight = genes.to.label)
p2 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)

SeuratOBJ@reductions

genes.to.label =  c("MMRN1", "EPHA5", "GPR151", "POU4F1")
FeaturePlot(SeuratOBJ, features = genes.to.label , split.by = "orig.ident", max.cutoff = 3, 
            cols = c("grey","red"), reduction = "integrated.cca")



plots <- VlnPlot(SeuratOBJ, features = genes.to.label, split.by = "orig.ident", group.by = "seurat_clusters",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)





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

