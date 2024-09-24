########################################################################
## Measure Quality Controls for GEX and ATAC assays from CellRanger-ARC data using Seurat & Signac packages
##
## Authors. CSC
## Date. Sep 23rd, 2024
## Last Update: XXX
##
## NOTES: 
## For a ~10k cells cellranger dataset it is recommended ~40G free mem to process the TSS() ATAC score.  Without storing the base-resolution matrix of integration counts at each site you can use less memory, but does not allow plotting the accessibility profile at the TSS.
## For slurm env: $srun --pty --mem=40GB --x11 bash
########################################################################

library(Seurat)                                 # 4.9.9.9045 2023-05-17 [1] Github (satijalab/seurat@7d1094c)
library(Signac)                                 # 1.9.0.9000 2023-05-08 [1] Github (stuart-lab/signac@cf31022)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(here)

here::here()

cellrangerDir <- here("processed-data", "cellrangerARC")
cellrangerDir_reanalyze <- here("processed-data", "cellrangerARC_reanalyze")

processedDir <- here("processed-data", "01_preprocessing_QC", "cellrangerARC")
processedDir_reanalyze <- here("processed-data", "01_preprocessing_QC", "cellrangerARC_reanalyze")

plotDir_reanalyze <- here("plots", "01_preprocessing_QC", "cellrangerARC_reanalyze")
functionsDir <- here("code", "01_preprocessing_QC")

# Check processed_data and plot directories exists
if (!dir.exists(processedDir_reanalyze)) { dir.create(processedDir_reanalyze) }
if (!dir.exists(plotDir_reanalyze)) { dir.create(plotDir_reanalyze) }

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
    png_name <- here(plotDir_reanalyze, paste0(sample_name,'_UMIs_Genes_MT.png'))  
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
    png_name <- here(plotDir_reanalyze, paste0(sample_name,'_Genes_Density.png'))  
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
    png_name <- here(plotDir_reanalyze, paste0(sample_name,'_Genes_Distribution.png'))  
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
    # if (i_vline > 0) {
    #     geom_vline(xintercept = i_vline, linetype=2)} +
    # if (i_hline > 0) {
    #     geom_hline(yintercept = i_hline, linetype=2)} +
    # facet_wrap(~orig.ident) +
    ylab("log10(genes-counts)") +
    xlab("log10(UMI-counts)") +
    ggtitle('UMIs/Genes by MT levels')
  png_name <- here(plotDir_reanalyze, paste0(sample_name,'_UMIS_per_MT.png'))  
  ggsave(p5, filename = png_name, height = 4, width = 4)
    
  # pALL <- p1 + p3 + p4 + p5 +
  #     plot_annotation(paste0(s_sample,' Quality Scores Before Quality Controls')) &
  #     theme(plot.tag = element_text(size = 10)) 
  # png_file <- paste0(sample_name,'_ALL.pdf')
  # png_name <- here(plotDir, png_file)  
  # ggsave(pALL, filename = png_name, height = 12, width = 7)
  
  message('Plots for RNA quality controls saved!')
    
}

get_Vplots_main_ATAC  <- function(seuratOBJ) {
  # TSS.enrichment will fail if you do not have enough memory
  if (!is.null(seuratOBJ@meta.data$TSS.enrichment)) { 
    p1 <- VlnPlot(object = seuratOBJ, features = c("nCount_ATAC", "nFeature_ATAC", "TSS.enrichment"), 
                  group.by = "orig.ident") # , ncol = 4 
  } else {
    p1 <- VlnPlot(object = seuratOBJ, features = c("nCount_ATAC", "nFeature_ATAC", "nucleosome_signal"), 
                  group.by = "orig.ident") # , ncol = 3 
  }
  p1 <- p1 &
    theme(#legend.position = 'none',
      axis.text.x = element_text(angle=0, hjust=1, size=8),  #10
      axis.text.y = element_text(size=8),  #10
      axis.title.x = element_blank(),
      axis.title.y = element_blank()) # &labs(title = "", x = 'Samples', y ="")
  
  return(p1)
  
} 

get_Vplots_blackR_ATAC  <- function(seuratOBJ) {
  
  if ((!is.null(seuratOBJ@meta.data$pct_reads_in_peaks)) && (!is.null(seuratOBJ@meta.data$blacklist_ratio))) { 
    p1 <- VlnPlot(object = seuratOBJ, features = c("pct_reads_in_peaks","blacklist_ratio"), group.by = "orig.ident", ncol = 2)  
  } else {
    p1 <- VlnPlot(object = seuratOBJ, features = c("pct_reads_in_peaks","blacklist_ratio"), group.by = "orig.ident", ncol = 2) 
  }
  p1 <- p1 &
    theme(#legend.position = 'none',
      axis.text.x = element_text(angle=0, hjust=1, size=8),  #10
      axis.text.y = element_text(size=8),  #10
      axis.title.x = element_blank(),
      axis.title.y = element_blank()) # &labs(title = "", x = 'Samples', y ="")
  
  return(p1)
  
} 



########################    Initials ########################  

## commandArgs scans the arguments which have been supplied when the current R script was invoked (from shell sh)
sample_tmp <- commandArgs(trailingOnly = TRUE)
# sample_tmp <- args[1]
# For testing: sample_tmp <- "S3_Hb_KDM_reanalysis, S3_Hb_KDM"
#              sample_tmp <- "4S_Hb_KDM_reanalysis, 4S_Hb_KDM"

sample_data = unlist(strsplit(sample_tmp,","))
crARC_Sample_r <- trimws(sample_data[[1]])  # Cell Ranger ARC reanalyze sample sub-directory name
crARC_Sample <- trimws(sample_data[[2]])    # Cell Ranger ARC sample sub-dir name
message('Processing `Cell Ranger ARC reanalyze` sample `', crARC_Sample_r, "` corresponding to `Cell Ranger ARC` sample `", crARC_Sample, "`")


# Create preliminary plots
b_get_GEX_plots <- TRUE         

# Remove mitochondrial levels (by sample dynamically)
b_get_filtered_GEX <- FALSE      
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

filtered_barcode_path <- here(cellrangerDir_reanalyze, crARC_Sample_r, "outs", "filtered_feature_bc_matrix.h5")
barcode_csv_path <- here(cellrangerDir, crARC_Sample, "outs", "per_barcode_metrics.csv")
fragments_tsv_path <- here(cellrangerDir, crARC_Sample, "outs", "atac_fragments.tsv.gz")

mtx <- Read10X_h5(filtered_barcode_path)             
rna_counts <- mtx$`Gene Expression`
metadata <- read.csv(file = barcode_csv_path, header = TRUE, row.names = 1)
meta_tmp = c('atac_peak_region_fragments','atac_fragments')
meta = metadata[meta_tmp]
message("Processing GEX sample (gex x cells):")
dim(rna_counts)


SeuratOBJ <- CreateSeuratObject(
  counts = rna_counts,
  assay = "RNA",
  project = crARC_Sample_r,
  meta.data = meta
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

message('Mitochondrial and Ribosomal percentage levels added')

## Optional meta-data for alculate percentages of largest genes by single cell
# if (b_additional_feat==TRUE) { SeuratOBJ <- l_get_perc_largest_genes(SeuratOBJ) }  


########  ################################################# ######## 
########        2. Get Visualizations for the GEX           ######## 
########  ################################################# ######## 

# Build UMI, Genes, MITO and RIBO violin plots, number of cells per sample, UMI/transcripts per cell plots, and 
#   Distribution of genes per cell histogram

## Before QCed data, plot GEX basic quality controls
if (b_get_GEX_plots) {
    base_name <- paste0(crARC_Sample_r, '_None_QC')
    plot_GEX_QCs(SeuratOBJ, base_name, TRUE)
    ## Calculate basic interquartile range for basic GEX stats
    # table_descriptive_stats_GEX(SeuratOBJ, base_name, s_tissue)
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

    rds_name <- here(processedDir, paste0(s_sample,'filteredM',i_filtering_method,'.rds'))
    saveRDS(SeuratOBJ.filtered, file = rds_name)
    message('Saving Seurat filtered.')
    
} else {
    
    # Only one Seurat object available
    lst_seurats <- list(SeuratOBJ)    
    rds_name <- here(processedDir_reanalyze, paste0(crARC_Sample_r,'.rds'))
    saveRDS(SeuratOBJ, file = rds_name)
    message('Saving Seurat none filtered.')
    
}


########  ################################################# ############ 
########  4. Create ATAC assay and attach it to Seurat object     ##### 
########  ################################################# ############ 

# Validate or Load Seurat object from disk if available
#SeuratOBJ <- load_seurat_obj(SeuratOBJ, s_seurat_name) 
#head(SeuratOBJ, n = 3)

## General use: Create gene annotations for hg38 and extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# show(annotations) / # View(head(annotations,n=5)) /# names(genomeStyles('Homo_sapiens'))
seqlevelsStyle(annotations) <- "UCSC"
length(extractSeqlevelsByGroup(species = 'Homo_sapiens', style = 'UCSC', group = 'all')) #group: sex, auto, circular
genome(annotations) <- "hg38"
# names(mcols(annotations))
#[1] "tx_id"        "gene_name"    "gene_id"      "gene_biotype" "type"   


##### Create and calculate chromatin QC metrics #######

## Parse a list of Seurat objects (lst_seurats). lst_seurats may containt:
##      (1) A new Seurat Obj from a specific experiment, and/or
##      (2) A Seurat Obj filtered by some GEX quality threshold (m1 or m2) 

message(length(lst_seurats), ' Seurat objects to process...')

for (S in lst_seurats) {
    ## S <- lst_seurats[[1]]
    rm('SeuratOBJ')
    # pull base name to label objects and plots 
    # testing: S <- lst_seurats[[1]]
    message('Processing ATAC for sample ', levels(S$`orig.ident`[1]))
    base_name <- paste0(levels(S$`orig.ident`[1]),'_QC')
    
    # Create the chromatin assay with annotations and attach it to the Seurat object 
    start.time = Sys.time()
    
    atac_counts <- mtx$Peaks
    SeuratOBJ <- get_create_atac_objs(S, atac_counts, fragments_tsv_path, annotations, TRUE) 
    print(start.time - Sys.time())        # 15s / 40G free_mem / 3k Cells / Object.size 806.2 Mb
    # f_InspectSeurat(SeuratOBJ)
    # tbs_atac <- get_basic_stats_ATAC(SeuratOBJ)

    if (b_get_ATAC_QC) {
        # Calculate Nucleosome Signal
        
        DefaultAssay(SeuratOBJ) <- "ATAC"

        # Calculate the strength of the nucleosome signal per cell
        SeuratOBJ <- NucleosomeSignal(SeuratOBJ)
        SeuratOBJ$nucleosome_signal
        plt_NS <- ggplot(SeuratOBJ@meta.data, aes(x=nucleosome_signal)) + 
          geom_histogram(aes(y=..density..), colour="black", fill="white")+
          geom_density(alpha=.5, color="darkblue", fill="lightblue") 
        #SeuratOBJ$nucleosome_group <- ifelse(SeuratOBJ$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
        SeuratOBJ$nucleosome_group <- ifelse(SeuratOBJ$nucleosome_signal > 2, 'NS_FAIL', 'NS_PASS')
        plt_NSgrp <- FragmentHistogram(object = SeuratOBJ, group.by = 'nucleosome_group')
        summary(SeuratOBJ$nucleosome_signal)
        #head(SeuratOBJ, n = 3) 
        message("NS score: ")
        addmargins(table(SeuratOBJ$nucleosome_group))

        # Calculate the "Transcription Start Site (TSS)" enrichment score
        tryCatch( {
            # Important NOTES: https://github.com/stuart-lab/signac/issues/374 
            #           1. Some times jumps an issue of ** memory allocation ** when "fast=FALSE", FALSE argument is mandatory to
            #              compute the TSS enrichment scores and visualize with TSSPlot(). You need more memory
            #           2. Error in `colnames<-`(`*tmp*`, value = seq_len(length.out = region.width) -  : attempt to set 'colnames' ...
            #              This is a vague message that would happen if no fragments are found in the set of TSS regions. 
            #              You could double-checking that the correct gene annotations is being used or you have a low ATAC quality
            22
            ## Filtering by ATAC basic quality controls 
          
            set.seed(24092024)  
            SeuratOBJ <- TSSEnrichment(SeuratOBJ, fast = FALSE) 
            summary(SeuratOBJ$TSS.enrichment)
            plt_TSS <- ggplot(SeuratOBJ@meta.data, aes(x=TSS.enrichment)) + 
              geom_histogram(aes(y=..density..), colour="black", fill="white")+
              geom_density(alpha=.5, color="darkblue", fill="lightblue")
            # Group by cells with TSS enrichment scores in two groups.
            #SeuratOBJ$high.tss <- ifelse(SeuratOBJ$TSS.enrichment > 1, 'High', 'Low')
            SeuratOBJ$high.tss <- ifelse(SeuratOBJ$TSS.enrichment > 1, 'TSS_PASS', 'TSS_FAIL')
            message("TSS scores: ")
            addmargins(table(SeuratOBJ$high.tss))
            pass_TSS <- sum(SeuratOBJ$high.tss == "TSS_PASS")
            fail_TSS <- (length(Cells(SeuratOBJ)) - pass_TSS)
            TSS_total <- length(Cells(SeuratOBJ))
            TSS_cap <- paste("TSS total:", TSS_total)
            TSS_caption1 <- paste("TSS PASS:", pass_TSS,  "(", round(pass_TSS*100 / TSS_cap, digits = 2), "%)")
            TSS_caption2 <- paste("TSS FAIL:", fail_TSS,  "(", round(fail_TSS*100 / TSS_cap, digits = 2), "%)")
            
            #colnames(SeuratOBJ@meta.data)  
            subtitle = "This is the subtitile."
            plt_TSSgrp <- TSSPlot(SeuratOBJ, group.by = 'high.tss') + NoLegend() 
            p1_TSS <- (plt_TSS + labs(title = "TSS distribution and TSS Scores")) / plt_TSSgrp + 
              theme(plot.caption = element_text(hjust = 0)) +
              labs(caption = paste(TSS_cap, "\n", TSS_caption1, "\n", TSS_caption2))
            png_file_TSS <- paste0(base_name,'_TSS.png')
            png_name <- here(plotDir_reanalyze, png_file_TSS)
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
        # Plot Peaks in black ratio and ATAC main feature scoreds
        p1_BlackR <- get_Vplots_blackR_ATAC(SeuratOBJ)
        p1_ATAC <- get_Vplots_main_ATAC(SeuratOBJ)
        
        png_file_NS <- paste0(base_name,'_Fragment_Distribution_grp.png')
        png_file_BlackR <- paste0(base_name,'_reads_in_peaks.png')
        png_file_ATAC <- paste0(base_name,'_ATAC_QCs.png')
        
        png_name <- here(plotDir_reanalyze, png_file_NS)
        p1_NS <- plt_NS + plt_NSgrp
        ggsave(p1_NS, filename = png_name, height = 4, width = 4)
        png_name <- here(plotDir_reanalyze, png_file_BlackR)
        ggsave(p1_BlackR, filename = png_name, height = 4, width = 5)
        png_name <- here(plotDir_reanalyze, png_file_ATAC)
        ggsave(p1_ATAC, filename = png_name, height = 4, width = 7)
  
        message('ATAC QCs plots saved!')  
    
    }
    # Save RDS Object
    rds_name <- here(processedDir_reanalyze, paste0(base_name,'_ATAC.rds'))
    saveRDS(SeuratOBJ, file = rds_name)
    message('ATAC RDS object saved!')   

}



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

