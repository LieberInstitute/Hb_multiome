########################################################################
## Measure Quality Controls for GEX and ATAC assays from CellRanger-ARC data using Seurat & Signac packages
##
## Authors. CSC / HT
## Date. April 21st, 2023
## Last Update: August 2024
##
## Input: Truncated H5, meta-data.cvs and fragments.tvs
## Output: rds Seurat objects and plots
##
## NOTES: 
## For a ~10k cells cellranger dataset it is recommended ~40G free mem to process the TSS() ATAC score.  Without storing the base-resolution matrix of integration counts at each site you can use less memory, but does not allow plotting the accessibility profile at the TSS.
## For slurm env: $srun --pty --mem=40GB --x11 bash
########################################################################

library(Seurat)                                 # 4.9.9.9045 2023-05-17 [1] Github (satijalab/seurat@7d1094c)
#options(Seurat.object.assay.version = 'v5')    # To use new Seurat v5: Please run: options(Seurat.object.assay.version = 'v5')
library(Signac)                                 # 1.9.0.9000 2023-05-08 [1] Github (stuart-lab/signac@cf31022)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
options(tidyverse.quiet = TRUE)
library(tidyverse)
#library(SeuratDisk)                             
library(here)

here::here()

processedDir <- here("processed-data", "01_preprocessing_QC")
plotDir <- here("plots", "01_preprocessing_QC")
functionsDir <- here("code", "functions_custom")

# Check processed_data and plot directories
if (!dir.exists(processedDir)) { dir.create(processedDir) }
if (!dir.exists(plotDir)) { dir.create(plotDir) }

source(here(functionsDir, "remote_file_caller.R"))       # Call functions to read paths
source(here(functionsDir, "remote_seurat_functions.R"))  # Call functions to create and handle Seurat object
source(here(functionsDir, "remote_signac_functions.R"))  # Call functions to create Signac object
source(here(functionsDir, "remote_plot_functions.R"))    # Call to plot GEX assay
source(here(functionsDir, "remote_plot_functions_ATAC.R"))       # Call to plot ATAC assay
source(here(functionsDir, "remote_filtering_functions.R"))       # Call functions to subset the Seurat object

# Function to plot GEX QCs 
plot_GEX_QCs <- function(SeuratO, sample_name, b_UMIscorr=FALSE) {    

    #@b_UMIscorr: If TRUE only plots UMI/Counts by MT levels, otherwise plot all
    
    # Assign the layer to the df
    df_genes_per_cell <- as.data.frame(SeuratO[[]])
    # Bar plot, number of cells per sample 
    #p2 <- get_plt_cells_by_sample(df_genes_per_cell)
    
    if (b_UMIscorr) {     
        # Violin plot with UMIs, Genes, ^MT and RIBO levels
        p1 <- get_Vplots_main_GEX(SeuratO)  
        png_name <- here(plotDir, paste0(sample_name,'_UMIs_Genes_MT.png'))  
        ggsave(p1, filename = png_name, height = 4, width = 7)
        message('Violin plots for UMIs, Genes and MT levels saved!')
        
        # Plot Genes Density per cell 
        p3 <- get_plt_genes_per_cell_density(df_genes_per_cell)
        png_name <- here(plotDir, paste0(sample_name,'_Genes_Density.png'))  
        ggsave(p3, filename = png_name, height = 4, width = 4)
        message('Genes density plot saved!')
        
        # Plot Genes Distribution per cell 
        p4 <- get_plt_genes_per_cell_boxplot(df_genes_per_cell)
        png_name <- here(plotDir, paste0(sample_name,'_Genes_Distribution.png'))  
        ggsave(p4, filename = png_name, height = 4, width = 4)
        message('Genes distribution plot saved!')
    }
    
    # Correlation btw genes and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
    p5 <- get_plt_UMIS_genes_MT_geomlm(df_genes_per_cell) # (df_genes_per_cell, 500, 500) 
    png_name <- here(plotDir, paste0(sample_name,'_UMIS_per_MT.png'))  
    ggsave(p5, filename = png_name, height = 4, width = 4)
    message('UMI/Genes by MT plot saved!')
    
    # pALL <- p1 + p3 + p4 + p5 +
    #     plot_annotation(paste0(s_sample,' Quality Scores Before Quality Controls')) &
    #     theme(plot.tag = element_text(size = 10)) 
    # png_file <- paste0(sample_name,'_ALL.pdf')
    # png_name <- here(plotDir, png_file)  
    # ggsave(pALL, filename = png_name, height = 12, width = 7)
    
    message('QCs reference saved')
    
}

# Some descriptive stats for further analysis
table_descriptive_stats_GEX <- function(SeuratO, sample_name, sample_tissue) {    
    
    tab_stats <- get_basic_stats_GEX(SeuratO, sample_tissue)
    message('Exporting table with QC quantiles for sample ', sample_name)
    s_file_name <- here(processedDir, paste0(sample_name,'_GEX_MITO_stats.csv'))
    write.csv(tab_stats, file=s_file_name, quote=TRUE, row.names=FALSE)

}
    
# p5 plot warning ----- CSC
# Warning message:
#     The following aesthetics were dropped during statistical transformation: colour
# ℹ This can happen when ggplot fails to infer the correct grouping structure in the data.
# ℹ Did you forget to specify a `group` aesthetic or to convert a numerical variable into a factor? 

########################    Initials ########################  

## commandArgs scans the arguments which have been supplied when the current R script was invoked (from shell sh)
sample_tmp <- commandArgs(trailingOnly = TRUE)
#sample_tmp <- args[1]
# testing
#sample_tmp <- 'S1_Hb_KDM,human'  # testing HUMAN tissue
#sample_tmp <- 'S2_Hb_KDM,human'  # testing HUMAN tissue
#sample_tmp <- '2_HPC_KDM,human'  # testing HUMAN tissue
#sample_tmp <- '3_HPC_KDM,mouse'  # testing MOUSE tissue 
#sample_tmp <- 'hippo42_1,human'  # testing HUMAN tissue 
sample_data = unlist(strsplit(sample_tmp,","))

s_sample <- sample_data[[1]]
s_tissue <- sample_data[[2]]
message('Processing sample: ', s_sample, ' from ', s_tissue, ' tissue.')

# Create preliminary plots
b_get_GEX_plots <- TRUE         
#if (b_get_GEX_plots) { source(here("code/functions_custom", "remote_plot_functions.R")) }

# Remove mitochondrial levels (by sample dynamically)
b_get_filtered_GEX <- TRUE      
if (b_get_filtered_GEX) {
    source(here(functionsDir, "remote_filtering_functions.R"))    # Filter Seurat assay by MT levels (GEX assay)
    i_filtering_method <- 1         # Pick up probabilities method to remove mito levels   
}

# Create attac and get quality control measures
b_get_ATAC_QC <- TRUE


########  ################################################# ######## 
########  1.  Create a Seurat object containing the RNA     
########      Recommended 40G of free_mem to 3k-10k cells 
########  ################################################# ######## 

# Read filtered barcode matrix, meta-data and fragment file names from cellranger-ARC
s_bc_mtx <- get_filtered_barcode_mtx(s_sample)      # Read H5 file
meta_path <- get_metadata_path(s_sample)            # Read csv file (meta-data path)
s_frag_namefile <- get_ATAC_barcode_tsv(s_sample)   # Read tvs file (fragments)

# Create Seurat Object 

SeuratOBJ <- get_seurat_obj(s_sample, s_bc_mtx, s_tissue, meta_path, FALSE)
# Syntax: function(seuratName, s_bc_mtx, s_tissue, s_meta, b_additional_feat = FALSE)
#SeuratOBJ@meta.data
print(SeuratOBJ)

message('Seurat object created successfully!')
message('Saving Seurat ...')
rds_name <- here(processedDir, paste0(s_sample,'.rds'))
saveRDS(SeuratOBJ, file = rds_name)


########  ################################################# ######## 
########        2. Get Visualizations for the GEX           ######## 
########  ################################################# ######## 

# Build and plot:
#       UMI, Genes, MITO and RIBO violin plots
#       Number of cells per sample
#       UMI/transcripts per cell
#       Distribution of genes per cell (histogram)

# Before QC any assay, plot the data ang get some descriptive stats
if (b_get_GEX_plots) { 
 
    base_name <- paste0(s_sample, '_None_QC', TRUE)
    plot_GEX_QCs(SeuratOBJ, base_name)
    # Calculate and save some descriptive stats for further analysis
    table_descriptive_stats_GEX(SeuratOBJ, base_name, s_tissue)
    
}
    
########  ################################################# ######## 
##        3. Preprocess Seurat Object based on the Mitochondrial percentage in the GEX assay
##              METHOD 1: (M1p)  calculate cut-off based on probabilities
##              METHOD 2: (M2sd) calculate cut-off at +/-2 Standard Deviations (95%)
########  ################################################# ######## 


# Create a second Seurat filtered by custom method of 1) probabilities or 2) SD
if (b_get_filtered_GEX) {
    # Preferred method 1. filter by probabilities
    
    # Validate or Load Seurat object from disk if available
    # Note. GEX assay need to be active

    if (i_filtering_method==1) {
        # filter by probabilities: M1
        SeuratOBJ.filtered <- get_seurat_GEX_filteringM1p(SeuratOBJ, FALSE)
    } else {    
        # filter by SD: M2
        SeuratOBJ.filtered <- get_seurat_GEX_filteringM2sd(SeuratOBJ, FALSE)
    }
    # Rename origin ident to identify the filtered Seurat object
    base_name <- paste0(levels(SeuratOBJ$`orig.ident`[1]), '_M', i_filtering_method)
    SeuratOBJ.filtered$orig.ident <- base_name
    #SeuratOBJ.filtered$`orig.ident`[1]

    plot_GEX_QCs(SeuratOBJ.filtered, base_name, TRUE)
    # Calculate and save some descriptive stats for further analysis
    table_descriptive_stats_GEX(SeuratOBJ.filtered, base_name, s_tissue)
    
    # Create a list of Seurat Objects to attach ATAC assay and process by QC
    lst_seurats <- list(SeuratOBJ, SeuratOBJ.filtered)
    
    message('Seurat filtered by MT level successfully!')
    
} else {
    
    # Only one Seurat object available
    lst_seurats <- list(SeuratOBJ)    
    
}


# message(paste0('Saving new Seurat filtered.'))
# rds_name <- here(processedDir, paste0(s_sample,'filteredM',i_filtering_method,'.rds'))
# saveRDS(SeuratOBJ.filtered, file = rds_name)


########  ################################################# ############ 
########  4. Create ATAC assay and attach it to Seurat object     ##### 
########  ################################################# ############ 

# Validate or Load Seurat object from disk if available
#SeuratOBJ <- load_seurat_obj(SeuratOBJ, s_seurat_name) 
#head(SeuratOBJ, n = 3)

# General use: Create gene annotations for hg38 and extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# show(annotations) / # View(head(annotations,n=5)) /# names(genomeStyles('Homo_sapiens'))
seqlevelsStyle(annotations) <- "UCSC"
length(extractSeqlevelsByGroup(species = 'Homo_sapiens', style = 'UCSC', group = 'all')) #group: sex, auto, circular
genome(annotations) <- "hg38"
# names(mcols(annotations))
#[1] "tx_id"        "gene_name"    "gene_id"      "gene_biotype" "type"   


####### Create and calculate chromatin QC metrics #######

# Parse a list of Seurat objects  (lst_seurats)

message(length(lst_seurats), ' Seurat objects to process...')

for (S in lst_seurats) {

    rm('SeuratOBJ')
    # pull base name to label objects and plots 
    # testing: S <- lst_seurats[[1]]
    message('Processing ATAC for sample ', levels(S$`orig.ident`[1]))
    base_name <- paste0(levels(S$`orig.ident`[1]),'_QC')
    
    # Create the chromatin assay with annotations and attach it to the Seurat object 
    start.time = Sys.time()
    SeuratOBJ <- get_create_atac_objs(S, s_bc_mtx, s_frag_namefile, annotations, TRUE) 
    end.time = Sys.time()
    print(end.time - start.time)        # 15s / 40G free_mem / 3k Cells / Object.size 806.2 Mb
    #f_InspectSeurat(SeuratOBJ)
    #tbs_atac <- get_basic_stats_ATAC(SeuratOBJ)

    if (b_get_ATAC_QC) {
        # Calculate Nucleosome Signal
        
        DefaultAssay(SeuratOBJ) <- "ATAC"

        # Calculate the strength of the nucleosome signal per cell
        SeuratOBJ <- NucleosomeSignal(SeuratOBJ)
        SeuratOBJ$nucleosome_group <- ifelse(SeuratOBJ$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
        p1_NS <- FragmentHistogram(object = SeuratOBJ, group.by = 'nucleosome_group')
        #head(SeuratOBJ, n = 3) 

        # Calculate the "Transcription Start Site (TSS)" enrichment score
        tryCatch( {
            # Important NOTES: https://github.com/stuart-lab/signac/issues/374 
            #           1. Some times jumps an issue of ** memory allocation ** when "fast=FALSE", FALSE argument is mandatory to
            #              compute the TSS enrichment scores and visualize with TSSPlot(). You need more memory
            #           2. Error in `colnames<-`(`*tmp*`, value = seq_len(length.out = region.width) -  : attempt to set 'colnames' ...
            #              This is a vague message that would happen if no fragments are found in the set of TSS regions. 
            #              You could double-checking that the correct gene annotations is being used or you have a low ATAC quality. 
            
            SeuratOBJ <- TSSEnrichment(SeuratOBJ, fast = FALSE) 
            # Group by cells with TSS enrichment scores in two groups.
            SeuratOBJ$high.tss <- ifelse(SeuratOBJ$TSS.enrichment > 2, 'High', 'Low')
            #colnames(SeuratOBJ@meta.data)
            p1_TSS <- TSSPlot(SeuratOBJ, group.by = 'high.tss') + NoLegend()
            png_file_TSS <- paste0(base_name,'_TSS.png')
            png_name <- here(plotDir, png_file_TSS)
            ggsave(p1_TSS, filename = png_name, height = 4, width = 4)
            
        }
        , error = function(e) { print('An error occurred. Check the annotation or you probably have ATAC quality loss.') } )
        
        # Add blacklist ratio and fraction of reads in peaks
        SeuratOBJ$blacklist_fraction <- FractionCountsInRegion(
            object = SeuratOBJ,
            assay = 'ATAC',
            regions = blacklist_hg38
        )
        # Add blacklist ratio and fraction of reads in peaks
        SeuratOBJ$pct_reads_in_peaks <- SeuratOBJ$atac_peak_region_fragments / SeuratOBJ$atac_fragments * 100
        SeuratOBJ$blacklist_ratio <- SeuratOBJ$blacklist_fraction / SeuratOBJ$atac_peak_region_fragments
        #        Error in `x[[i, drop = TRUE]]`:
        #       ! 'blacklist_fraction' not found in this Seurat object
        # Plot Peaks in black ratio and ATAC main feature scoreds
        p1_BlackR <- get_Vplots_blackR_ATAC(SeuratOBJ)
        p1_ATAC <- get_Vplots_main_ATAC(SeuratOBJ)
        
        png_file_NS <- paste0(base_name,'_Fragment_Distribution_grp.png')
        png_file_BlackR <- paste0(base_name,'_reads_in_peaks.png')
        png_file_ATAC <- paste0(base_name,'_ATAC_QCs.png')
        
        png_name <- here(plotDir, png_file_NS)
        ggsave(p1_NS, filename = png_name, height = 4, width = 4)
        png_name <- here(plotDir, png_file_BlackR)
        ggsave(p1_BlackR, filename = png_name, height = 4, width = 5)
        png_name <- here(plotDir, png_file_ATAC)
        ggsave(p1_ATAC, filename = png_name, height = 4, width = 7)
  
        message('ATAC QCs plots saved!')  
    
    }
    # Save RDS Object
    rds_name <- here(processedDir, paste0(base_name,'_ATAC.rds'))
    saveRDS(SeuratOBJ, file = rds_name)
    message('ATAC RDS object saved!')   
}



# Chunk of code to filter by conventional or custom profile moved to: 01b_preprocessing_GEX_ATAC.R script


############ Reproducibility information ####################

library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
#"2023-04-04 12:42:26 EDT"
proc.time()
options(width = 120)
session_info()


# > session_info()
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.1 Patched (2023-07-19 r84711)
# os       Rocky Linux 9.2 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2024-02-06
# pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3/bin/ (via rmarkdown)
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package                     * version    date (UTC) lib source
# abind                         1.4-5      2016-07-21 [2] CRAN (R 4.3.1)
# AnnotationDbi               * 1.62.2     2023-07-02 [2] Bioconductor
# AnnotationFilter            * 1.24.0     2023-04-25 [2] Bioconductor
# backports                     1.4.1      2021-12-13 [2] CRAN (R 4.3.1)
# base64enc                     0.1-3      2015-07-28 [2] CRAN (R 4.3.1)
# beeswarm                      0.4.0      2021-06-01 [2] CRAN (R 4.3.1)
# Biobase                     * 2.60.0     2023-04-25 [2] Bioconductor
# BiocFileCache                 2.8.0      2023-04-25 [2] Bioconductor
# BiocGenerics                * 0.46.0     2023-04-25 [2] Bioconductor
# BiocIO                        1.10.0     2023-04-25 [2] Bioconductor
# BiocParallel                  1.34.2     2023-05-22 [2] Bioconductor
# biomaRt                       2.56.1     2023-06-09 [2] Bioconductor
# Biostrings                  * 2.68.1     2023-05-16 [2] Bioconductor
# biovizBase                    1.48.0     2023-04-25 [2] Bioconductor
# bit                           4.0.5      2022-11-15 [2] CRAN (R 4.3.1)
# bit64                         4.0.5      2020-08-30 [2] CRAN (R 4.3.1)
# bitops                        1.0-7      2021-04-24 [2] CRAN (R 4.3.1)
# blob                          1.2.4      2023-03-17 [2] CRAN (R 4.3.1)
# BPCells                       0.1.0      2023-10-02 [1] Github (bnprks/BPCells@ac4376d)
# BSgenome                    * 1.68.0     2023-04-25 [2] Bioconductor
# BSgenome.Hsapiens.UCSC.hg38 * 1.4.5      2023-09-25 [1] Bioconductor
# cachem                        1.0.8      2023-05-01 [2] CRAN (R 4.3.1)
# checkmate                     2.3.0      2023-10-25 [1] CRAN (R 4.3.1)
# cli                           3.6.1      2023-03-23 [2] CRAN (R 4.3.1)
# cluster                       2.1.4      2022-08-22 [3] CRAN (R 4.3.1)
# codetools                     0.2-19     2023-02-01 [3] CRAN (R 4.3.1)
# colorspace                    2.1-0      2023-01-23 [2] CRAN (R 4.3.1)
# cowplot                       1.1.1      2020-12-30 [2] CRAN (R 4.3.1)
# crayon                        1.5.2      2022-09-29 [2] CRAN (R 4.3.1)
# curl                          5.1.0      2023-10-02 [1] CRAN (R 4.3.1)
# data.table                    1.14.8     2023-02-17 [2] CRAN (R 4.3.1)
# DBI                           1.1.3      2022-06-18 [2] CRAN (R 4.3.1)
# dbplyr                        2.3.3      2023-07-07 [2] CRAN (R 4.3.1)
# DelayedArray                  0.26.7     2023-07-28 [2] Bioconductor
# deldir                        1.0-9      2023-05-17 [2] CRAN (R 4.3.1)
# dichromat                     2.0-0.1    2022-05-02 [2] CRAN (R 4.3.1)
# digest                        0.6.33     2023-07-07 [2] CRAN (R 4.3.1)
# dotCall64                     1.1-0      2023-10-17 [1] CRAN (R 4.3.1)
# dplyr                       * 1.1.4      2023-11-17 [1] CRAN (R 4.3.1)
# ellipsis                      0.3.2      2021-04-29 [2] CRAN (R 4.3.1)
# EnsDb.Hsapiens.v86          * 2.99.0     2023-09-25 [1] Bioconductor
# ensembldb                   * 2.24.0     2023-04-25 [2] Bioconductor
# evaluate                      0.23       2023-11-01 [1] CRAN (R 4.3.1)
# fansi                         1.0.5      2023-10-08 [1] CRAN (R 4.3.1)
# farver                        2.1.1      2022-07-06 [2] CRAN (R 4.3.1)
# fastDummies                   1.7.3      2023-07-06 [1] CRAN (R 4.3.1)
# fastmap                       1.1.1      2023-02-24 [2] CRAN (R 4.3.1)
# fastmatch                     1.1-4      2023-08-18 [2] CRAN (R 4.3.1)
# filelock                      1.0.2      2018-10-05 [2] CRAN (R 4.3.1)
# fitdistrplus                  1.1-11     2023-04-25 [1] CRAN (R 4.3.1)
# forcats                     * 1.0.0      2023-01-29 [2] CRAN (R 4.3.1)
# foreign                       0.8-85     2023-09-09 [3] CRAN (R 4.3.1)
# Formula                       1.2-5      2023-02-24 [2] CRAN (R 4.3.1)
# future                        1.33.0     2023-07-01 [2] CRAN (R 4.3.1)
# future.apply                  1.11.0     2023-05-21 [1] CRAN (R 4.3.1)
# generics                      0.1.3      2022-07-05 [2] CRAN (R 4.3.1)
# GenomeInfoDb                * 1.36.3     2023-09-07 [2] Bioconductor
# GenomeInfoDbData              1.2.10     2023-07-20 [2] Bioconductor
# GenomicAlignments             1.36.0     2023-04-25 [2] Bioconductor
# GenomicFeatures             * 1.52.2     2023-08-25 [2] Bioconductor
# GenomicRanges               * 1.52.0     2023-04-25 [2] Bioconductor
# ggbeeswarm                    0.7.2      2023-04-29 [2] CRAN (R 4.3.1)
# ggplot2                     * 3.4.4      2023-10-12 [1] CRAN (R 4.3.1)
# ggrastr                       1.0.2      2023-06-01 [2] CRAN (R 4.3.1)
# ggrepel                       0.9.4      2023-10-13 [1] CRAN (R 4.3.1)
# ggridges                      0.5.4      2022-09-26 [2] CRAN (R 4.3.1)
# globals                       0.16.2     2022-11-21 [2] CRAN (R 4.3.1)
# glue                          1.6.2      2022-02-24 [2] CRAN (R 4.3.1)
# goftest                       1.2-3      2021-10-07 [1] CRAN (R 4.3.1)
# gridExtra                     2.3        2017-09-09 [1] CRAN (R 4.3.1)
# gtable                        0.3.4      2023-08-21 [2] CRAN (R 4.3.1)
# hdf5r                         1.3.8      2023-01-21 [2] CRAN (R 4.3.1)
# here                        * 1.0.1      2020-12-13 [2] CRAN (R 4.3.1)
# Hmisc                         5.1-1      2023-09-12 [2] CRAN (R 4.3.1)
# hms                           1.1.3      2023-03-21 [2] CRAN (R 4.3.1)
# htmlTable                     2.4.2      2023-10-29 [1] CRAN (R 4.3.1)
# htmltools                     0.5.7      2023-11-03 [1] CRAN (R 4.3.1)
# htmlwidgets                   1.6.2      2023-03-17 [2] CRAN (R 4.3.1)
# httpuv                        1.6.12     2023-10-23 [1] CRAN (R 4.3.1)
# httr                          1.4.7      2023-08-15 [2] CRAN (R 4.3.1)
# ica                           1.0-3      2022-07-08 [1] CRAN (R 4.3.1)
# igraph                        1.5.1      2023-08-10 [2] CRAN (R 4.3.1)
# IRanges                     * 2.34.1     2023-06-22 [2] Bioconductor
# irlba                         2.3.5.1    2022-10-03 [2] CRAN (R 4.3.1)
# jsonlite                      1.8.7      2023-06-29 [2] CRAN (R 4.3.1)
# KEGGREST                      1.40.0     2023-04-25 [2] Bioconductor
# KernSmooth                    2.23-22    2023-07-10 [3] CRAN (R 4.3.1)
# knitr                         1.45       2023-10-30 [1] CRAN (R 4.3.1)
# labeling                      0.4.3      2023-08-29 [2] CRAN (R 4.3.1)
# later                         1.3.1      2023-05-02 [2] CRAN (R 4.3.1)
# lattice                       0.21-8     2023-04-05 [3] CRAN (R 4.3.1)
# lazyeval                      0.2.2      2019-03-15 [2] CRAN (R 4.3.1)
# leiden                        0.4.3      2022-09-10 [1] CRAN (R 4.3.1)
# lifecycle                     1.0.4      2023-11-07 [1] CRAN (R 4.3.1)
# listenv                       0.9.0      2022-12-16 [2] CRAN (R 4.3.1)
# lmtest                        0.9-40     2022-03-21 [2] CRAN (R 4.3.1)
# lubridate                   * 1.9.3      2023-09-27 [1] CRAN (R 4.3.1)
# magrittr                      2.0.3      2022-03-30 [2] CRAN (R 4.3.1)
# MASS                          7.3-60     2023-05-04 [3] CRAN (R 4.3.1)
# Matrix                        1.6-1.1    2023-09-18 [3] CRAN (R 4.3.1)
# MatrixGenerics                1.12.3     2023-07-30 [2] Bioconductor
# matrixStats                   1.1.0      2023-11-07 [1] CRAN (R 4.3.1)
# memoise                       2.0.1      2021-11-26 [2] CRAN (R 4.3.1)
# mgcv                          1.9-0      2023-07-11 [3] CRAN (R 4.3.1)
# mime                          0.12       2021-09-28 [2] CRAN (R 4.3.1)
# miniUI                        0.1.1.1    2018-05-18 [2] CRAN (R 4.3.1)
# munsell                       0.5.0      2018-06-12 [2] CRAN (R 4.3.1)
# nlme                          3.1-163    2023-08-09 [3] CRAN (R 4.3.1)
# nnet                          7.3-19     2023-05-03 [3] CRAN (R 4.3.1)
# parallelly                    1.36.0     2023-05-26 [2] CRAN (R 4.3.1)
# patchwork                   * 1.1.3      2023-08-14 [2] CRAN (R 4.3.1)
# pbapply                       1.7-2      2023-06-27 [2] CRAN (R 4.3.1)
# pillar                        1.9.0      2023-03-22 [2] CRAN (R 4.3.1)
# pkgconfig                     2.0.3      2019-09-22 [2] CRAN (R 4.3.1)
# plotly                        4.10.3     2023-10-21 [1] CRAN (R 4.3.1)
# plyr                          1.8.9      2023-10-02 [1] CRAN (R 4.3.1)
# png                           0.1-8      2022-11-29 [1] CRAN (R 4.3.1)
# polyclip                      1.10-6     2023-09-27 [1] CRAN (R 4.3.1)
# prettyunits                   1.1.1      2020-01-24 [2] CRAN (R 4.3.1)
# progress                      1.2.2      2019-05-16 [2] CRAN (R 4.3.1)
# progressr                     0.14.0     2023-08-10 [1] CRAN (R 4.3.1)
# promises                      1.2.1      2023-08-10 [2] CRAN (R 4.3.1)
# ProtGenerics                  1.32.0     2023-04-25 [2] Bioconductor
# purrr                       * 1.0.2      2023-08-10 [2] CRAN (R 4.3.1)
# R6                            2.5.1      2021-08-19 [2] CRAN (R 4.3.1)
# ragg                          1.2.5      2023-01-12 [2] CRAN (R 4.3.1)
# RANN                          2.6.1      2019-01-08 [2] CRAN (R 4.3.1)
# rappdirs                      0.3.3      2021-01-31 [2] CRAN (R 4.3.1)
# RColorBrewer                  1.1-3      2022-04-03 [2] CRAN (R 4.3.1)
# Rcpp                          1.0.11     2023-07-06 [2] CRAN (R 4.3.1)
# RcppAnnoy                     0.0.21     2023-07-02 [2] CRAN (R 4.3.1)
# RcppHNSW                      0.5.0      2023-09-19 [2] CRAN (R 4.3.1)
# RcppRoll                      0.3.0      2018-06-05 [1] CRAN (R 4.3.1)
# RCurl                         1.98-1.12  2023-03-27 [2] CRAN (R 4.3.1)
# readr                       * 2.1.4      2023-02-10 [2] CRAN (R 4.3.1)
# reshape2                      1.4.4      2020-04-09 [2] CRAN (R 4.3.1)
# restfulr                      0.0.15     2022-06-16 [2] CRAN (R 4.3.1)
# reticulate                    1.34.0     2023-10-12 [1] CRAN (R 4.3.1)
# rjson                         0.2.21     2022-01-09 [2] CRAN (R 4.3.1)
# rlang                         1.1.2      2023-11-04 [1] CRAN (R 4.3.1)
# rmarkdown                     2.25       2023-09-18 [2] CRAN (R 4.3.1)
# ROCR                          1.0-11     2020-05-02 [2] CRAN (R 4.3.1)
# rpart                         4.1.19     2022-10-21 [3] CRAN (R 4.3.1)
# rprojroot                     2.0.4      2023-11-05 [1] CRAN (R 4.3.1)
# Rsamtools                     2.16.0     2023-04-25 [2] Bioconductor
# RSpectra                      0.16-1     2022-04-24 [2] CRAN (R 4.3.1)
# RSQLite                       2.3.1      2023-04-03 [2] CRAN (R 4.3.1)
# rstudioapi                    0.15.0     2023-07-07 [2] CRAN (R 4.3.1)
# rtracklayer                 * 1.60.1     2023-08-15 [2] Bioconductor
# Rtsne                         0.16       2022-04-17 [2] CRAN (R 4.3.1)
# S4Arrays                      1.0.6      2023-08-30 [2] Bioconductor
# S4Vectors                   * 0.38.2     2023-09-22 [1] Bioconductor
# scales                        1.2.1      2022-08-20 [2] CRAN (R 4.3.1)
# scattermore                   1.2        2023-06-12 [1] CRAN (R 4.3.1)
# sctransform                   0.4.1      2023-10-19 [1] CRAN (R 4.3.1)
# sessioninfo                 * 1.2.2      2021-12-06 [2] CRAN (R 4.3.1)
# Seurat                      * 4.9.9.9067 2023-10-02 [1] Github (satijalab/seurat@99b9ded)
# SeuratObject                * 4.9.9.9091 2023-09-25 [1] Github (mojaveazure/seurat-object@c51dd86)
# shiny                         1.7.5.1    2023-10-14 [1] CRAN (R 4.3.1)
# Signac                      * 1.11.9000  2023-09-25 [1] Github (stuart-lab/signac@4f60de7)
# sp                          * 2.1-1      2023-10-16 [1] CRAN (R 4.3.1)
# spam                          2.10-0     2023-10-23 [1] CRAN (R 4.3.1)
# spatstat.data                 3.0-3      2023-10-24 [1] CRAN (R 4.3.1)
# spatstat.explore              3.2-5      2023-10-22 [1] CRAN (R 4.3.1)
# spatstat.geom                 3.2-7      2023-10-20 [1] CRAN (R 4.3.1)
# spatstat.random               3.2-1      2023-10-21 [1] CRAN (R 4.3.1)
# spatstat.sparse               3.0-3      2023-10-24 [1] CRAN (R 4.3.1)
# spatstat.utils                3.0-4      2023-10-24 [1] CRAN (R 4.3.1)
# stringi                       1.8.1      2023-11-13 [1] CRAN (R 4.3.1)
# stringr                     * 1.5.1      2023-11-14 [1] CRAN (R 4.3.1)
# SummarizedExperiment          1.30.2     2023-06-06 [2] Bioconductor
# survival                      3.5-7      2023-08-14 [3] CRAN (R 4.3.1)
# systemfonts                   1.0.4      2022-02-11 [2] CRAN (R 4.3.1)
# tensor                        1.5        2012-05-05 [1] CRAN (R 4.3.1)
# textshaping                   0.3.6      2021-10-13 [2] CRAN (R 4.3.1)
# tibble                      * 3.2.1      2023-03-20 [2] CRAN (R 4.3.1)
# tidyr                       * 1.3.0      2023-01-24 [2] CRAN (R 4.3.1)
# tidyselect                    1.2.0      2022-10-10 [2] CRAN (R 4.3.1)
# tidyverse                   * 2.0.0      2023-02-22 [2] CRAN (R 4.3.1)
# timechange                    0.2.0      2023-01-11 [2] CRAN (R 4.3.1)
# tzdb                          0.4.0      2023-05-12 [2] CRAN (R 4.3.1)
# utf8                          1.2.4      2023-10-22 [1] CRAN (R 4.3.1)
# uwot                          0.1.16     2023-06-29 [2] CRAN (R 4.3.1)
# VariantAnnotation             1.46.0     2023-04-25 [2] Bioconductor
# vctrs                         0.6.4      2023-10-12 [1] CRAN (R 4.3.1)
# vipor                         0.4.5      2017-03-22 [2] CRAN (R 4.3.1)
# viridisLite                   0.4.2      2023-05-02 [2] CRAN (R 4.3.1)
# withr                         2.5.2      2023-10-30 [1] CRAN (R 4.3.1)
# xfun                          0.41       2023-11-01 [1] CRAN (R 4.3.1)
# XML                           3.99-0.14  2023-03-19 [2] CRAN (R 4.3.1)
# xml2                          1.3.5      2023-07-06 [2] CRAN (R 4.3.1)
# xtable                        1.8-4      2019-04-21 [2] CRAN (R 4.3.1)
# XVector                     * 0.40.0     2023-04-25 [2] Bioconductor
# yaml                          2.3.7      2023-01-23 [2] CRAN (R 4.3.1)
# zlibbioc                      1.46.0     2023-04-25 [2] Bioconductor
# zoo                           1.8-12     2023-04-13 [2] CRAN (R 4.3.1)
# 
# [1] /users/csoto/R/4.3
# [2] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/library

