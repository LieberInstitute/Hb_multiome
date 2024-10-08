########################################################################
## Plot gene expression for RNA assay for integrated Seurats
## Authors. CSC
## Date. March 8th, 2024
## Last.Adaptation: xxx
##
## Input:  Seurat integrated with CCA / Harmony
## Output:  (1) DoPlot
##          (2) Feature plots in umap
##          (3) Violin plots 
##          (4) DoHeatmap
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
# library(RColorBrewer) #To use gradient color to umaps 
theme_set(theme_cowplot())

library('purrr')
library('tidyverse')
library('here')

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

# =============================================================================================================

## Compose Seurat object name

## USE THIS CHUNK IF YOU HAVE pre-existing Seurat with pre-selected clusters AND Jump to visualizations ...

if (count_mtx_type=='data_counts') { Seurat_base_name <- 'seurat.combined.data_counts_PCA' } else { Seurat_base_name <- 'seurat.combined.norm_counts_PCA' }
if (Seurat_reduction=='CCA') {
    rds_name <- here('processed-data/03_pseudobulking', paste0(Seurat_base_name, '_CCA_selected_clust.rds'))
} else {
    rds_name <- here('processed-data/03_pseudobulking', paste0(Seurat_base_name, '_Harmony_selected_clust.rds'))
}
rds_name
# EX. seurat.combined.data_counts_PCA_Harmony_selected_clust.rdS
# Load pre-existing seurat
SeuratOBJ_Hb  <- get_seurat(rds_name)

## NOTE: Jump to PLOTS label


# =============================================================================================================


## USE THIS CHUNK IF YOU HAVE pre-existing Seurat with none pre-selected clusters. ALL the clusters. 

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
table(SeuratOBJ$seurat_clusters)
SeuratOBJ@reductions

## Subset the clusters and marker genes selected
SeuratOBJ_Hb <- subset(SeuratOBJ, subset = seurat_clusters %in% clust_selected)
levels(Idents(SeuratOBJ_Hb))



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

message(' Seurat intgrated idents renamed! ')   




###################################.        PLOTS      #################################################

# check the count cells by clusters
table(Idents(SeuratOBJ_Hb))
# Ex. Harmony
# C_6 C_12 C_14 C_15 
# 830  180  102   99 



###################### DoPlot for gene expression across the selected clusters    ######################


pdf_file <- paste0(Seurat_base_name, '_', Seurat_reduction, '_GeneExpr_DoPlots_Hb.pdf')
pdf_name <- here('plots/03_pseudobulking', pdf_file)
pdf(file = pdf_name)

p1 <- DotPlot(SeuratOBJ_Hb, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) +
  RotatedAxis() 
# Define costume descriptions
L1 <- geom_vline(xintercept = length(Hb), linetype="dotted", color = "gray")
L2 <- geom_vline(xintercept = (length(Hb) + length(LHb)), linetype="dotted", color = "gray")
L3 <- geom_vline(xintercept = length(markers.to.plot), linetype="dotted", color = "gray")
main_title <- paste0("Gene expression across clusters from ", Seurat_reduction)
sub_title1 <-  paste0('Broad Hb | LHb | MHb')

p1 <- p1 + L1 + L2 + L3 
p1 + ggtitle(label = main_title,
                            subtitle = sub_title1)

dev.off()




###################### Feature Plot for gene expression across the clusters selected in UMAP ###################### 


Hb_type <- list(Hb)
Hb_type <- append(Hb_type, list(LHb))
Hb_type <- append(Hb_type, list(MHb))
Hb_title_lst <- c('Broad Habenula', 'Lateral Habenula', 'Median Habenula')
#SeuratOBJ_Hb@reductions

## Function to plot genes (features) in umap  

plot_Feature_UMAP <- function(g, st) {

    FeaturePlot(SeuratOBJ_Hb, features = g, split.by = "orig.ident",
                    reduction = "umap",
                    label=TRUE, label.size = 3, label.color = "black", repel = TRUE,
                    ncol = 3) +
    plot_layout(ncol = 2, nrow = 3, axis_titles = "collect") +
    plot_annotation( title = main_title, subtitle = st ) &
    theme(text = element_text(face = "bold"), 
          plot.title = element_text(size = 14),
          plot.subtitle = element_text(size = 10),
          axis.text.x=element_text(size=10),
          axis.text.y=element_text(size=10),
          axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10),
          axis.title.y.right = element_text(size=10,face="bold")) 
  
}

#plots2 <- lapply(X = plot_list, FUN = function(x) x + theme(plot.title = element_text(size = 1)))

## initialize variables 

main_title <- paste0("Gene expression across clusters from ", Seurat_reduction)
genes.to.plotALL = list()
plot_list = list()
st_lst = list()

## create the lists of genes and sub-titles to plot by page

for (x in 1:length(Hb_type)) {

  genes.tmp <- Hb_type[[x]]
  ## If sub-list of markers grater than 6 genes, split the list to set 6 plots per page in a grid of 2x3 
  if (length(genes.tmp) > 6) { genes.to.plot <- split(genes.tmp, ceiling(seq_along(genes.tmp) / 3)) } else { genes.to.plot <- list(Hb_type[[x]]) }
  genes.to.plotALL <- append(genes.to.plotALL, genes.to.plot)
  st_tmp <- rep(Hb_title_lst[x], length(genes.to.plot))
  st_lst <- append(st_lst, st_tmp)
  
}

#str(genes.to.plotALL)

## map the list(s) and build Feature plots in umaps for Hb, LHb and MHb

if (length(genes.to.plotALL) == length(st_lst)) {
  plot_list <- map2(genes.to.plotALL, st_lst, 
                    ~ plot_Feature_UMAP(g = .x, st = .y)) }


plot_list[[1]]

pdf_file <- paste0(Seurat_base_name, '_', Seurat_reduction, '_FeaturePlot_Hb.pdf')
pdf_name <- here('plots/03_pseudobulking', pdf_file)
pdf(file = pdf_name)
# Ex. file name: seurat.combined.data_counts_PCA_Harmony_FeaturePlot_Hb.pdf
print(plot_list)
dev.off()





###################### Violin plot for gene expression across the clusters selected  ###################### 


## create violin plots

plot_violin <- function(g, st) {
  
    VlnPlot(SeuratOBJ_Hb, features = g, group.by = "seurat_clusters", 
            pt.size = 0, combine = TRUE) + geom_boxplot() + 
    plot_layout(ncol = 3, nrow = 4, axis_titles = "collect", guides = "collect") +
    plot_annotation(title = main_title, subtitle = st)  & # , caption = 'Samples: S1 and S2'
      theme_bw()  &
      theme(text = element_text(face = "bold"),
            legend.position="none",
            plot.title = element_text(size = 14),
            plot.subtitle = element_text(size = 10),
            axis.text.x=element_text(size=8),
            axis.text.y=element_text(size=8),
            axis.title.x = element_blank(),
            axis.title.y = element_blank()
            )
}


## initialize variables

main_title <- paste0("Gene expression across clusters from ", Seurat_reduction)
plot_list = list()
genes.to.plotALL = list()
st_lst = list()

## create the lists of genes and sub-titles to plot by page

for (x in 1:length(Hb_type)) {
  
  genes.tmp <- Hb_type[[x]]
  ## Split list of markers grater than 6 genes, to set 6 plots per page in a grid of 1x6 
  if (length(genes.tmp) > 6) { genes.to.plot <- split(genes.tmp, ceiling(seq_along(genes.tmp) / 12)) } else { genes.to.plot <- list(Hb_type[[x]]) }
  genes.to.plotALL <- append(genes.to.plotALL, genes.to.plot)
  st_tmp <- rep(Hb_title_lst[x], length(genes.to.plot))
  st_lst <- append(st_lst, st_tmp)
  
}

## map the lists to build Violin plots for Hb, LHb and MHb

if (length(genes.to.plotALL) == length(st_lst)) {
  plot_list <- map2(genes.to.plotALL, st_lst, 
                    ~ plot_violin(g = .x, st = .y)) }

#plot_list[[3]]

pdf_file <- paste0(Seurat_base_name, '_', Seurat_reduction, '_ViolinPlots_Hb.pdf')
pdf_name <- here('plots/03_pseudobulking', pdf_file)
pdf(file = pdf_name)
# Ex. file name: seurat.combined.data_counts_PCA_Harmony_ViolinPlots_Hb.pdf
print(plot_list)
dev.off()



###################### Basic Heatmap for gene expression across the clusters selected (mpt sorted or balanced) ###################### 

## initialize variables

main_title <- paste0("Gene expression across clusters from ", Seurat_reduction)
sub_title <-  paste0('Broad Hb | LHb | MHb')

## build basic heatmap 

p1 <- DoHeatmap(
  SeuratOBJ_Hb,
  features = markers.to.plot,
  cells = NULL,
  group.by = "orig.ident",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = NULL,
  label = TRUE,
  size = 3,
  hjust = 0,
  vjust = 0,
  angle = 0,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
) &
theme(
      plot.title = element_text(size = 14),
      plot.subtitle = element_text(size = 10),
      axis.text.x=element_text(size=10),
      axis.text.y=element_text(size=10)
      ) 

p1 + ggtitle(label = main_title,
             subtitle = sub_title)
# L1 <- geom_hline(xintercept = length(Hb), linetype="dotted", color = "gray")
# L2 <- geom_hline(xintercept = (length(Hb) + length(LHb)), linetype="dotted", color = "gray")
# L3 <- geom_hline(xintercept = length(markers.to.plot), linetype="dotted", color = "gray")
# p1 <- p1 + L1 + L2 + L3  

pdf_file <- paste0(Seurat_base_name, '_', Seurat_reduction, '_GeneExpr_Heatmap_Hb.pdf')
pdf_name <- here('plots/03_pseudobulking', pdf_file)
# Ex. file name: seurat.combined.data_counts_PCA_Harmony_GeneExpr_Heatmap_Hb.pdf
pdf(file = pdf_name)
dev.off()




############ Reproducibility information ####################

library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
proc.time()
options(width = 120)
session_info()

