########################################################################
## Plot gene expression for RNA assay from CCA and Harmony clusters
## Authors. CSC
## Date. March 8th, 2024
## Last.Adaptation: xxx
##
## Input: Seurat integrated with correction
## Output: Seurat with pre-selected clusters and multiple plots
##
## NOTES:
## For slurm env: runsrun --x11 --pty --partition=interactive bash
########################################################################

library('Seurat')                                 # 4.9.9.9045 2023-05-17 [1] Github (satijalab/seurat@7d1094c)
## Additional required packages for rename clusters
library('scCustomize')
## Packages to plot
library('ggplot2')
library('patchwork')
library('cowplot')
theme_set(theme_cowplot())

library(here)

here::here()

if (!packageVersion("Seurat")=='4.9.9.9060') {
    stop
    message('This pipeline was implemented with Seurat v5 and Signac v1.11+ ')
    message('You need the laterst Seurat v5 (‘4.9.9.9060’)')
    message('Current available repository on: https://satijalab.org/seurat/articles/install.html  ') }

# Check if plot directory exists, if not create it
if (!dir.exists(here("plots/03_pseudobulking/"))) {
    dir.create(here("plots/03_pseudobulking/")) }



########################    Initials ########################  

## Select the count-mtx to merge (raw or normalized data)
count_mtx_type <- 'data_counts'
#count_mtx_type <- 'norm_counts' 
#Seurat_reduction <- 'CCA'
Seurat_reduction <- 'Harmony'

## load pre-existing Seurat
get_seurat <- function(name) {
    sobj <- readRDS(name)
    # verification of the integration
    print(table(sobj$orig.ident))
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


## Set pre-selected clusters and marker genes to plot
if (Seurat_reduction=='CCA') {
  
  clust_selected <- c(8, 15, 18)
  ## genes selected from Hb_pilot, defined as Hb general, LHb and MHb
  Hb <- c("MMRN1", "GPR151", "POU4F1") # Hb
  LHb <- c("LINC01876") # LHb
  MHb <- c("CD24", "AC004594.1") #MHb
  markers.to.plot <- c(Hb, LHb, MHb)
  
} else {
  
  clust_selected <- c(6, 12, 14, 15)
  ## genes selected from Hb_pilot, defined as Hb general, LHb and MHb
  ## marker genes pre-selected (from Hb_pilot)
  Hb <- c("MMRN1", "GPR151", "POU4F1") # Hb
  LHb <- c("EPHA5", "TLL1") # LHb
  MHb <- c("AC109466.1", "AC008415.1", "GPR149", "GNG8", "NEUROD1", "RASGRP1", "SLC5A7", "CHRNB3", "SCUBE1", "LINC02143", "CD24") #MHb
  markers.to.plot <- c(Hb, LHb, MHb)
  
}  


## Compose Seurat object name

## USE THIS OF YOU HAVE pre-existing Seurat with pre-selected clusters AND Jump to visualizations directly

if (count_mtx_type=='data_counts') { Seurat_base_name <- 'seurat.combined.data_counts_PCA' } else { Seurat_base_name <- 'seurat.combined.norm_counts_PCA' }
if (Seurat_reduction=='CCA') {
  rds_name <- here('processed-data/03_pseudobulking', paste0(Seurat_base_name, '_CCA__selected_clust.rds'))
} else {
  rds_name <- here('processed-data/03_pseudobulking', paste0(Seurat_base_name, '_Harmony_selected_clust.rds'))
}
rds_name
# EX. seurat.combined.data_counts_PCA_Harmony_selected_clust.rdS
# Load pre-existing seurat
SeuratOBJ_Hb  <- get_seurat(rds_name)

## NOTE: Jump to visualizations directly


# =============================================================================================================


## USE THIS OF YOU want to plot all the clusters

if (count_mtx_type=='data_counts') { Seurat_base_name <- 'seurat.combined.data_counts_PCA' } else { Seurat_base_name <- 'seurat.combined.norm_counts_PCA' }
if (Seurat_reduction=='CCA') {
  rds_name <- here('processed-data/02_merge_seurats', paste0(Seurat_base_name, '_CCA.rds'))
} else {
  rds_name <- here('processed-data/02_merge_seurats', paste0(Seurat_base_name, '_Harmony.rds'))
}
rds_name
# Load pre-existing seurat
SeuratOBJ <- get_seurat(rds_name)



# =============================================================================================================



# some verification
table(SeuratOBJ$orig.ident)
# S1_Hb_KDM S2_Hb_KDM 
# 8178      9816 

colnames(SeuratOBJ@meta.data)
SeuratOBJ@reductions

## Subset the clusters and marker genes selected
SeuratOBJ_Hb <- subset(SeuratOBJ, subset = seurat_clusters %in% clust_selected)
levels(Idents(SeuratOBJ))


####################  Rename clusters  ######################

oldIdent <- levels(Idents(SeuratOBJ_Hb))
# Ex. Harmony: [1] "6"  "12" "14" "15"
newIdent <- paste("C", 0:(length(oldIdent)-1), sep = "_")

newIdent <- lapply(oldIdent, function(x) { paste0("C_", x) })
unlist(newIdent)

## rename clusters to make them more readable
## require scCustomize/Wrapper funtion to rename clusters
SeuratOBJ_Hb <- Rename_Clusters(SeuratOBJ_Hb, new_idents = newIdent,
                         meta_col_name = "seurat_clusters.renamed")


rds_name <- paste0(Seurat_base_name,'_', Seurat_reduction, '_selected_clust.rds')
rds_name <- here('processed-data/03_pseudobulking', rds_name)
# file name: ~/.seurat.combined.data_counts_PCA_Harmony_pseudobulk.rds
saveRDS(SeuratOBJ_Hb, file = rds_name)
message(' Gene distribution plots completed! ')   




###################################.        PLOTS      #################################################

# check the count cells by clusters
table(Idents(SeuratOBJ_Hb))
# Ex. Harmony
# C_6 C_12 C_14 C_15 
# 830  180  102   99 

###################### DoPlot for gene expression across the selected clusters    ######################

## File name
pdf_file <- paste0(Seurat_base_name, '_', Seurat_reduction, '_GeneExpr_DoPlots_Hb.pdf')
pdf_name <- here('plots/03_pseudobulking', pdf_file)
pdf(file = pdf_name)

unique(Idents(SeuratOBJ_Hb))
# Ex. Harmony
# [1] C_6  C_12 C_14 C_15
# Levels: C_6 C_12 C_14 C_15

p1 <- DotPlot(SeuratOBJ_Hb, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) +
  RotatedAxis() 
# Define costume descriptions
L1 <- geom_vline(xintercept = length(Hb), linetype="dotted", color = "gray")
L2 <- geom_vline(xintercept = (length(Hb) + length(LHb)), linetype="dotted", color = "gray")
L3 <- geom_vline(xintercept = length(markers.to.plot), linetype="dotted", color = "gray")
main_title <- paste0("Gene expression across clusters from ", Seurat_reduction)
sub_title1 <-  paste0('Broad Hb | LHb | MHb')

p1 <- p1 + L1 + L2 + L3 
p1 + ggtitle(label = , main_title,
                            subtitle = sub_title1)

dev.off()




###################### Feature Plot for gene expression across the clusters selected in UMAP ###################### 


pdf_file <- paste0(Seurat_base_name, '_', Seurat_reduction, '_FeaturePlot_Hb.pdf')
pdf_name <- here('plots/03_pseudobulking', pdf_file)
pdf(file = pdf_name)
#par.save <- par(mfrow = c(4, 1))
# seurat.combined.data_counts_PCA_Harmony_FeaturePlot_Hb.pdf

Hb_type <- list(Hb)
Hb_type <- append(Hb_type, list(LHb))
Hb_type <- append(Hb_type, list(MHb))
Hb_title_lst <- c('Broad Habenula', 'Lateral Habenula', 'Median Habenula')


main_title <- paste0("Gene expression across clusters from ", Seurat_reduction)
#SeuratOBJ_Hb@reductions

# for loop to control numbe of plots per page (set to 3)
for (x in length(Hb_type)) {
  
  sub_title <- Hb_title_lst[x]
  genes.to.plot <- as.vector(unlist(Hb_type[x]))
  
  maxl <- length(genes.to.plot). #max number of plts per page 
  
  `%+=%` = function(e1,e2) eval.parent(substitute(e1 <- e1 + e2))
  x1 = 1
  x2 = 1
  x2 %+=% 2 ; x2
    
  while (x2 == maxl) {
    if (x2 > maxl) { 
        x2 <- maxl
        plot_feat <- genes.to.plot[x1:x2]
      } else {
        plot_feat <- genes.to.plot[x1:x2]
        x2 <- x2 + 1
        x1 <- x2
        x2 %+=% 2 ; x2
      }
    print(plot_feat)
    
    p2 <- FeaturePlot(SeuratOBJ_Hb, features = plot_feat, split.by = "orig.ident", 
                      reduction = "umap",
                      label=TRUE, label.size = 2, label.color = "black", repel = TRUE,
                      ncol = 3) + 
      plot_layout(ncol = 2, nrow = 3) + 
      plot_layout(axis_titles = "collect") +
      plot_annotation( title = main_title, subtitle = sub_title ) 
      # A patchworked ggplot object if combine = TRUE
      p2  
    
  }

}

#par(par.save)
dev.off()




###################### Violin plot for gene expression across the clusters selected  ###################### 


pdf_file <- paste0(Seurat_base_name, '_', Seurat_reduction, '_ViolinPlots_Hb.pdf')
pdf_name <- here('plots/03_pseudobulking', pdf_file)
main_title <- paste0("Gene expression across clusters from ", Seurat_reduction)

## Plot gene expression across the clusters selected in umap reduction

#Hb_type
plot_list = list()
for (x in 1:length(Hb_type)) {

  sub_title2 <- Hb_title_lst[x]
  genes.to.plot <- as.vector(unlist(Hb_type[x]))
  
  #plot_violinQC(SeuratOBJ_Hb, genes.to.plot, sub_title2)

  ## Plot Violin plots for the same LHb/MHb marker genes
  Vplots_Hb <- VlnPlot(SeuratOBJ_Hb, features = genes.to.plot, 
                       group.by = "seurat_clusters", pt.size = 0,
                       ncol = 3, combine = TRUE) + 
    plot_layout(ncol = 2, nrow = 3) + 
    plot_layout(axis_titles = "collect") +
    plot_annotation(title = main_title, subtitle = sub_title2, caption = 'Samples: S1 and S2')

    #plot_annotation( title = main_title, subtitle = sub_title2 ) 
  plot_list[[x]] <- Vplots_Hb  #wrap_plots(Vplots_Hb, ncol = 1)
  
}

pdf(file = pdf_name)
print(plot_list)
dev.off()



###################### Violin plot for gene expression across the clusters selected  ###################### 


pdf_file <- paste0(Seurat_base_name, '_', Seurat_reduction, '_GeneExpr_Heatmap_Hb.pdf')
pdf_name <- here('plots/03_pseudobulking', pdf_file)
pdf(file = pdf_name)

DoHeatmap(
  SeuratOBJ_Hb,
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

dev.off()


############ Reproducibility information ####################

library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
proc.time()
options(width = 120)
session_info()

