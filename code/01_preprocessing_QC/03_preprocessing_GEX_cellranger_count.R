########################################################################
## Build Seurat and calculate initial Quality Controls for CellRanger-count only GEX data
##
## Authors. CSC
## Date. Oct 04th, 2024
##
## NOTES: 
## For a ~10k cells cellranger dataset it is recommended ~40G free mem to process the TSS() ATAC score.  Without storing the base-resolution matrix of integration counts at each site you can use less memory, but does not allow plotting the accessibility profile at the TSS.
## For slurm env: $srun --pty --mem=40GB --x11 bash
########################################################################

library(Seurat)                                 # 4.9.9.9045 2023-05-17 [1] Github (satijalab/seurat@7d1094c)
#library(Signac)                                 # 1.9.0.9000 2023-05-08 [1] Github (stuart-lab/signac@cf31022)
#library(EnsDb.Hsapiens.v86)
#library(BSgenome.Hsapiens.UCSC.hg38)
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(here)

here::here()

cellrangerDir <- here("processed-data", "cellrangerGEX")
processedDir <- here("processed-data", "01_preprocessing_QC", "cellranger_count")
plotDir <- here("plots", "01_preprocessing_QC", "cellranger_count")
functionsDir <- here("code", "01_preprocessing_QC")

# Check processed_data and plot directories exists
#if (!dir.exists(processedDir_reanalyze)) { dir.create(processedDir_reanalyze) }
if (!dir.exists(processedDir)) { dir.create(processedDir) }
if (!dir.exists(plotDir)) { dir.create(plotDir) }

source(here(functionsDir, "remote_seurat_functions_v2.R"))  # Call functions to create and handle Seurat object
source(here(functionsDir, "remote_signac_functions_v2.R"))  # Call functions to create Signac object

# Function to plot GEX basic quality controls 
plot_GEX_QCs <- function(SeuratOBJ, sample_name, b_UMIscorr=FALSE) {    

  # @b_UMIscorr:  If TRUE only plots UMI/Counts by MT levels, otherwise plot all. 
  #               This is used to re-plot the correlation in QCed data
    
    # Assign the layer to the df
  genes_per_cell <- as.data.frame(SeuratOBJ[[]])
  if (b_UMIscorr) {
    
    p1 <- VlnPlot(object = SeuratOBJ, 
                  features = c('nCount_RNA','nFeature_RNA','percent.mt'), 
                  layer = "counts", 
                  group.by = "orig.ident", raster=FALSE) &
      theme(#legend.position = 'none',
        axis.text.x = element_text(angle=0, hjust=1, size=8),  #10
        axis.text.y = element_text(size=8),  #10
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) # &labs(title = "", x = 'Samples', y ="")
    png_name <- here(plotDir, paste0(sample_name,'_UMIs_Genes_MT.png'))  
    ggsave(p1, filename = png_name, height = 4, width = 7)

    ## Plot Genes Density per cell 
    p3 <- genes_per_cell %>%
      ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) +
      geom_density(alpha = 0.2) +
      scale_x_log10() +
      theme_classic() +
      theme(plot.title = element_text(hjust=0.5)) +
      geom_vline(xintercept = 300) +
      ylab("Log10(UMIs)") +
      xlab("Gene-counts") +
      ggtitle("Genes density by cell")   
    png_name <- here(plotDir, paste0(sample_name,'_Genes_Density.png'))  
    ggsave(p3, filename = png_name, height = 4, width = 4)

    # Plot Genes Distribution per cell 
    p4 <- genes_per_cell %>%
      ggplot(aes(x=orig.ident, y=log10(nFeature_RNA), fill=orig.ident)) +
      geom_boxplot(alpha = 0.7) +
      theme_classic() +
      theme(axis.text.x = element_text(vjust = 1, hjust=1)) +
      theme(plot.title = element_text(hjust=0.5)) +
      ylab("Log10(gene-counts)") +
      xlab("") +
      ggtitle("Genes distribution by cell")
    png_name <- here(plotDir, paste0(sample_name,'_Genes_Distribution.png'))  
    ggsave(p4, filename = png_name, height = 4, width = 4)

  }
    
  # Correlation btw genes and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
  #p5 <- get_plt_UMIS_genes_MT_geomlm(genes_per_cell) # (df_genes_per_cell, 500, 500) 
  p5 <- genes_per_cell %>%
    ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt, group.by = 'orig.ident')) + # MTRatio
    #    ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=MTRatio)) + # MTRatio
    geom_point() +
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method=lm) +
    scale_x_log10() +
    scale_y_log10() +
    theme_classic() +
    ylab("log10(genes-counts)") +
    xlab("log10(UMI-counts)") +
    ggtitle('UMIs/Genes by MT levels')
  png_name <- here(plotDir, paste0(sample_name,'_UMIS_per_MT.png'))  
  ggsave(p5, filename = png_name, height = 4, width = 4)
    
  # pALL <- p1 + p3 + p4 + p5 +
  #     plot_annotation(paste0(s_sample,' Quality Scores Before Quality Controls')) &
  #     theme(plot.tag = element_text(size = 10)) 
  # png_file <- paste0(sample_name,'_ALL.pdf')
  # png_name <- here(plotDir, png_file)  
  # ggsave(pALL, filename = png_name, height = 12, width = 7)
  
  message('Plots for RNA quality controls saved!')
    
}



########################    Initials ########################  

## commandArgs scans the arguments which have been supplied when the current R script was invoked (from shell sh)
sample_args <- commandArgs(trailingOnly = TRUE)
sample_name <- sample_args[1]
# For testing: sample_name <- "8C_Hb_KDM"
#              sample_name <- "12C_Hb_KDM"
s_tissue <- "human"
message('Processing `Cell Ranger count (GEX)` for sample `', sample_name, "`")


# Create preliminary plots
b_get_GEX_plots <- TRUE         



####  1.  Create a Seurat object containing the RNA     


filtered_barcode_path <- here(cellrangerDir, sample_name, "outs", "filtered_feature_bc_matrix.h5")
#barcode_csv_path <- here(cellrangerDir, sample_name, "outs", "per_barcode_metrics.csv")

rna_counts <- Read10X_h5(filtered_barcode_path)
message("Processing GEX sample (gex x cells):")
#dim(rna_counts)

SeuratOBJ <- CreateSeuratObject(
  counts = rna_counts,
  assay = "RNA",
  project = sample_name,
  #meta.data = meta
)
SeuratOBJ

message('Seurat object created successfully!')

## Add additional meta-data: Chr-Mitochondrial levels

SeuratOBJ$log10GenesPerUMI <- log10(SeuratOBJ$nFeature_RNA) / log10(SeuratOBJ$nCount_RNA)
## For humans or mouse. GRCh38 and mm10, respectively
if (s_tissue=='human') {
  SeuratOBJ[["percent.mt"]] <- PercentageFeatureSet(SeuratOBJ, pattern = "^MT-")
} else { # it is mouse
  SeuratOBJ[["percent.mt"]] <- PercentageFeatureSet(SeuratOBJ, pattern = "^Mt")
}
if (s_tissue=='human') {
  SeuratOBJ[["percent.ribo"]] <- PercentageFeatureSet(SeuratOBJ, pattern = "^RP[LS]")
} else { # it is mouse
  SeuratOBJ[["percent.ribo"]] <- PercentageFeatureSet(SeuratOBJ, pattern = "^Rp[ls]")
}
SeuratOBJ[["MTRatio"]] <- SeuratOBJ$percent.mt / 100
x <- as.data.frame.character(summary(SeuratOBJ$MTRatio))

message('Mitochondrial and Ribosomal percentage levels added')
colnames(SeuratOBJ@meta.data)

## Optional meta-data for alculate percentages of largest genes by single cell
# if (b_additional_feat==TRUE) { SeuratOBJ <- l_get_perc_largest_genes(SeuratOBJ) }


########  ################################################# ########
########        2. Get Visualizations for the GEX           ########
########  ################################################# ########

# Build UMI, Genes, MITO and RIBO violin plots, number of cells per sample, UMI/transcripts per cell plots, and
#   Distribution of genes per cell histogram

## Before QCed data, plot GEX basic quality controls
if (b_get_GEX_plots) {
    base_name <- paste0(sample_name)
    plot_GEX_QCs(SeuratOBJ, base_name, TRUE)
    ## Calculate basic interquartile range for basic GEX stats
    source(here(functionsDir, "remote_seurat_functions_v2.R"))  # Call functions to create and handle Seurat object
    x <- get_basic_stats_GEX(SeuratOBJ, s_tissue)
    #write_csv(as.data.frame(x), here(processedDir_reanalyze, paste0(crARC_Sample_r, "_metrics_summary.csv")))
    write_csv(as.data.frame(x), here(processedDir, paste0(sample_name, "_CRcount_metrics_summary.csv")))
}


# Only one Seurat object available
rds_name <- here(processedDir, paste0(sample_name,'_CR.rds'))
saveRDS(SeuratOBJ, file = rds_name)

message('Saving Seurat none filtered.')
message('Preprocess completed!')


############ Reproducibility information ####################

library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
#"2023-04-04 12:42:26 EDT"
proc.time()
options(width = 120)
session_info()

