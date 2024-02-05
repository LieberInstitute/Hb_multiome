########################################################################
## Measure Quality Controls for GEX and ATAC assays from CellRanger-ARC data using Seurat & Signac packages
##
## Authors. CSC / HT
## Date. April 21st, 2023
## Last.Adaptation: Jan.22, 2024
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

if (!packageVersion("Seurat")=='4.9.9.9060') {
    stop
    message('This pipeline was implemented with Seurat v5 and Signac v1.11+ ')
    message('You need the laterst Seurat v5 (‘4.9.9.9060’)')
    message('Current available repository on: https://satijalab.org/seurat/articles/install.html  ') }

# Check if processed_data directory exists, if not create it
if (!dir.exists(here("processed-data/01_preprocessing_QC/"))) {
    dir.create(here("processed-data/01_preprocessing_QC/"))
}
# Check if plot directory exists, if not create it
if (!dir.exists(here("plots/01_preprocessing_QC/"))) {
    dir.create(here("plots/01_preprocessing_QC/"))
}

source(here("code/functions_custom", "remote_file_caller.R"))       # Call functions to read paths
source(here("code/functions_custom", "remote_seurat_functions.R"))  # Call functions to create and handle Seurat object
source(here("code/functions_custom", "remote_signac_functions.R"))  # Call functions to create Signac object
source(here("code/functions_custom", "remote_plot_functions.R"))    # Call to plot GEX assay
source(here("code/functions_custom", "remote_filtering_functions.R"))   # Call functions to subset the Seurat object

# Function to plot GEX QCs 
plot_GEX_QCs <- function(SeuratO, sample_name) {    
    
    # Violin plot with UMIs, Genes, ^MT and RIBO levels
    p1 <- get_Vplots_main_GEX(SeuratO)  
    png_file <- paste0(sample_name,'_UMIs_Genes_MT.png')
    png_name <- here('plots/01_preprocessing_QC', png_file)  
    ggsave(p1, filename = png_name, height = 4, width = 7)
    message('Violin plots for UMIs, Genes and MT levels saved!')
    
    # Assign the layer to the df
    df_genes_per_cell <- as.data.frame(SeuratO[[]])
    # Bar plot, number of cells per sample 
    #p2 <- get_plt_cells_by_sample(df_genes_per_cell)
    
    # Plot Genes Density per cell 
    p3 <- get_plt_genes_per_cell_density(df_genes_per_cell)
    png_file <- paste0(sample_name,'_Genes_Density.png')
    png_name <- here('plots/01_preprocessing_QC', png_file)  
    ggsave(p3, filename = png_name, height = 4, width = 4)
    message('Genes density plot saved!')
    
    # Plot Genes Distribution per cell 
    p4 <- get_plt_genes_per_cell_boxplot(df_genes_per_cell)
    png_file <- paste0(sample_name,'_Genes_Distribution.png')
    png_name <- here('plots/01_preprocessing_QC', png_file)  
    ggsave(p4, filename = png_name, height = 4, width = 4)
    message('Genes distribution plot saved!')
    
    # Correlation btw genes and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
    p5 <- get_plt_UMIS_genes_MT_geomlm(df_genes_per_cell) # (df_genes_per_cell, 500, 500) 
    png_file <- paste0(sample_name,'_UMIS_per_MT.png')
    png_name <- here('plots/01_preprocessing_QC', png_file)  
    ggsave(p5, filename = png_name, height = 4, width = 4)
    message('UMI/Genes by MT plot saved!')
    
    # pALL <- p1 + p3 + p4 + p5 +
    #     plot_annotation(paste0(s_sample,' Quality Scores Before Quality Controls')) &
    #     theme(plot.tag = element_text(size = 10)) 
    # png_file <- paste0(sample_name,'_ALL.pdf')
    # png_name <- here('plots/01_preprocessing_QC', png_file)  
    # ggsave(pALL, filename = png_name, height = 12, width = 7)
    
    message('QCs reference saved')
    
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
message('Processing sample: ',s_sample, ' from ', s_tissue, ' tissue.')

# Create preliminary plots
b_get_GEX_plots <- TRUE         
#if (b_get_GEX_plots) { source(here("code/functions_custom", "remote_plot_functions.R")) }

# Remove mitochondrial levels (by sample dynamically)
b_get_filtered_GEX <- TRUE      
if (b_get_filtered_GEX) {
    source(here("code/functions_custom", "remote_filtering_functions.R"))    # Filter Seurat assay by GEX
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
str(SeuratOBJ)
print(SeuratOBJ)

message('Seurat object created successfully!')
message('Saving Seurat ...')
rds_name <- here('processed-data/01_preprocessing_QC', paste0(s_sample,'.rds'))
saveRDS(SeuratOBJ, file = rds_name)


########  ################################################# ######## 
########        2. Get Visualizations for the GEX           ######## 
########  ################################################# ######## 

# Build and plot:
#       UMI, Genes, MITO and RIBO violin plots
#       Number of cells per sample
#       UMI/transcripts per cell
#       Distribution of genes per cell (histogram)

# Before QC any assay
if (b_get_GEX_plots) { 
    base_name <- paste0(s_sample, '_None_QC')
    plot_GEX_QCs(SeuratOBJ, base_name)
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
    
    # Some stats and plot correlation 
    get_basic_stats_GEX(SeuratOBJ.filtered, s_tissue)
    df_genes_per_cell <- as.data.frame(SeuratOBJ.filtered[[]])
    # Plot correlation between genes and number of UMIs to visualize strong presence of cells with low numbers of genes/UMIs
    p1 <- get_plt_UMIS_genes_MT_geomlm(df_genes_per_cell, 500, 500)
    png_file <- paste0(s_sample,'_corr_filteredM',i_filtering_method,'.png')
    png_name <- here('plots/01_preprocessing_QC', png_file)  
    ggsave(p1, filename = png_name, height = 4, width = 4)
    message('UMI/Counts by MT plot saved!')
    
    # Create a list of Seurat Objects to attach ATAC assay and process by QC
    lst_seurats <- c(SeuratOBJ, SeuratOBJ.filtered)

} else {
    # Only one Seurat Onb available
    lst_seurats <- c(SeuratOBJ)    
}

message('Seurat filtered successfully!')
# message(paste0('Saving new Seurat filtered.'))
# rds_name <- here('processed-data/01_preprocessing_QC', paste0(s_sample,'filteredM',i_filtering_method,'.rds'))
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

# Parse a list of Seurat objects  
for (S in lst_seurats) {

    message('Processing ATAC for sample ', s_sample)
    
    # Create the chromatin assay with annotations and attach it to the Seurat object 
    start.time = Sys.time()
    SeuratOBJ <- get_create_atac_objs(S, s_bc_mtx, s_frag_namefile, annotations, TRUE) 
    end.time = Sys.time()
    print(end.time - start.time)        # 15s / 40G free_mem / 3k Cells / Object.size 806.2 Mb
    # 38s / 40G free_mem / 10k cells / 
    SeuratOBJ@meta.data$`orig.ident`[1]
    #f_InspectSeurat(SeuratOBJ)
    # get_basic_stats_ATAC(SeuratOBJ)
    
    ####### Evaluate chromatin assay #######
    # Build base-name for plots and rds objects
    if (b_get_filtered_GEX) {
        base_name <- paste0(s_sample, '_M', i_filtering_method,'_QC')
    } else {            
        base_name <- paste0(s_sample, '_None_QC')
    } 
    
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
            
            SeuratOBJ <- TSSEnrichment(SeuratOBJ, fast = TRUE) 
            # Group by cells with TSS enrichment scores in two groups.
            SeuratOBJ$high.tss <- ifelse(SeuratOBJ$TSS.enrichment > 2, 'High', 'Low')
            #colnames(SeuratOBJ@meta.data)
            p1_TSS <- TSSPlot(SeuratOBJ, group.by = 'high.tss') + NoLegend()
            png_file_TSS <- paste0(base_name,'_TSS.png')
            png_name <- here('plots/01_preprocessing_QC', png_file_TSS)
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
        p1_BlackR <- VlnPlot(SeuratOBJ, features = c("pct_reads_in_peaks","blacklist_ratio"), ncol = 2)
        p1_ATAC <- VlnPlot(SeuratOBJ, features = c("nCount_ATAC", "nFeature_ATAC", "nucleosome_signal", "TSS.enrichment"), ncol = 4)
        
        png_file_NS <- paste0(base_name,'_Fragment_Distribution_grp.png')
        png_file_BlackR <- paste0(base_name,'_reads_in_peaks.png')
        png_file_ATAC <- paste0(base_name,'_ATAC_QCs.png')
        
        png_name <- here('plots/01_preprocessing_QC', png_file_NS)
        ggsave(p1_NS, filename = png_name, height = 4, width = 4)
        png_name <- here('plots/01_preprocessing_QC', png_file_BlackR)
        ggsave(p1_BlackR, filename = png_name, height = 4, width = 4)
        png_name <- here('plots/01_preprocessing_QC', png_file_ATAC)
        ggsave(p1_ATAC, filename = png_name, height = 4, width = 7)
  
        message('ATAC QCs plots saved!')  
    
    }
    # Save RDS Object
    rds_name <- here('processed-data/01_preprocessing_QC', paste0(base_name,'_ATAC.rds'))
    saveRDS(SeuratOBJ, file = rds_name)
    message('ATAC RDS object saved!')     
}



# if (b_get_ATAC) {     
# 
#     # Validate or Load Seurat object from disk if available
#     #SeuratOBJ <- load_seurat_obj(SeuratOBJ, s_seurat_name) 
#     #head(SeuratOBJ, n = 3)
#     
#     # General use: Create gene annotations for hg38 and extract gene annotations from EnsDb
#     annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
#     # show(annotations) / # View(head(annotations,n=5)) /# names(genomeStyles('Homo_sapiens'))
#     seqlevelsStyle(annotations) <- "UCSC"
#     length(extractSeqlevelsByGroup(species = 'Homo_sapiens', style = 'UCSC', group = 'all')) #group: sex, auto, circular
#     genome(annotations) <- "hg38"
#     
#     # names(mcols(annotations))
#     #[1] "tx_id"        "gene_name"    "gene_id"      "gene_biotype" "type"   
#     
#     # Create the chromatin assay with annotations and attach it to the seurat object 
#     start.time = Sys.time()
#     SeuratOBJ <- get_create_atac_objs(SeuratOBJ, 
#                                       s_bc_mtx, s_frag_namefile, 
#                                       annotations, TRUE) 
#     end.time = Sys.time()
#     print(end.time - start.time)        # 15s / 40G free_mem / 3k Cells / Object.size 806.2 Mb
#                                         # 38s / 40G free_mem / 10k cells / 
#     SeuratOBJ@meta.data$`orig.ident`[1]
#     #f_InspectSeurat(SeuratOBJ)
#     # get_basic_stats_ATAC(SeuratOBJ)
#     # ATAC counts loaded successfully
#     # ATAC counts annotation attached successfully
#     # Computing hash
#     # Checking for 9816 cell barcodes
#     # Chromatin assay completed successfully
#     # ChromatinAssay data with 64288 features for 9816 cells
#     # Variable features: 0 
#     # Genome: 
#     #     Annotation present: TRUE 
#     # Motifs present: FALSE 
#     # Fragment files: 1 
#     
#     if (b_get_filtered_GEX) {
#         # Preferred method 1. filter by probabilities
#         SeuratOBJ.filtered <- get_create_atac_objs(SeuratOBJ.filtered, 
#                                           s_bc_mtx, s_frag_namefile, 
#                                           annotations, TRUE) 
#         # ATAC counts loaded successfully
#         # ATAC counts annotation attached successfully
#         # Computing hash
#         # Checking for 9816 cell barcodes
#         # Chromatin assay completed successfully
#         # ChromatinAssay data with 64288 features for 9816 cells
#         # Variable features: 0 
#         # Genome: 
#         #     Annotation present: TRUE 
#         # Motifs present: FALSE 
#         # Fragment files: 1 
#     }    
# 
# }


    
# message('ATAC attached successfully!')



########  ################################################# ######## 
########  5. QC ATAC metrics: Signac package                ######## 
########  ################################################# ######## 

# Get next ATAC measures:
#           1. Fragment size distribution
#           2. Nucleosome signal
#           3. TSS profile
#           4. TSS enrichment at nucleosome signal
#           5. Linear model in the #UMIs against the #Genes by ^MT levels 


# Evaluate chromatin assay

# if (b_get_ATAC_QC) {
# 
#     if (b_get_filtered_GEX) {
#         
#         # Calculate Nucleosome Signal and TSS Scores 
#         
#         DefaultAssay(SeuratOBJ.filtered) <- "ATAC"
#         
#         # Preferred method 1. filter by probabilities
#         SeuratOBJ.filtered <- NucleosomeSignal(SeuratOBJ.filtered)
#         SeuratOBJ.filtered$nucleosome_group <- ifelse(SeuratOBJ.filtered$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
#         p1 <- FragmentHistogram(object = SeuratOBJ.filtered, group.by = 'nucleosome_group')
#         base_name <- paste0(s_sample, '_M1_QC')
#         png_file <- paste0(base_name,'_Fragment_Distribution_grp.png')
#         png_name <- here('plots/01_preprocessing_QC', png_file)  
#         ggsave(p1, filename = png_name, height = 4, width = 4)
#         message('Fragments Distribution plot saved!')    
#         head(SeuratOBJ.filtered, n = 3)    
#         
#     } else {
#         
#         DefaultAssay(SeuratOBJ) <- "ATAC"
#         
#         # Calculate the strength of the nucleosome signal per cell
#         SeuratOBJ <- NucleosomeSignal(SeuratOBJ)
#         SeuratOBJ$nucleosome_group <- ifelse(SeuratOBJ$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
#         p1 <- FragmentHistogram(object = SeuratOBJ, group.by = 'nucleosome_group')
#         base_name <- paste0(s_sample, '_None_QC')
#         png_file <- paste0(base_name,'_Fragment_Distribution_grp.png')
#         png_name <- here('plots/01_preprocessing_QC', png_file)  
#         ggsave(p1, filename = png_name, height = 4, width = 4)
#         message('Fragments Distribution plot saved!')        
#         head(SeuratOBJ, n = 3) 
#         
#     }
# 
# }
    
# if (b_get_ATAC_QC) {
# 
# 
#     
#     head(SeuratOBJ, n = 3)
#     
#     # Calculate the "Transcription Start Site (TSS)" enrichment score for each cell, as defined by ENCODE.
#     tryCatch( {
#         
#         # Important NOTES: https://github.com/stuart-lab/signac/issues/374 
#         #           1. Some times jumps an issue of ** memory allocation ** when "fast=FALSE", FALSE argument is mandatory to
#         #              compute the TSS enrichment scores and visualize with TSSPlot(). You need more memory
#         #           2. Error in `colnames<-`(`*tmp*`, value = seq_len(length.out = region.width) -  : attempt to set 'colnames' ...
#         #              This is a vague message that would happen if no fragments are found in the set of TSS regions. 
#         #              You could double-checking that the correct gene annotations is being used or you have a low ATAC quality. 
#         
#         SeuratOBJ <- TSSEnrichment(SeuratOBJ, fast = TRUE)
#         SeuratOBJ <- TSSEnrichment(SeuratOBJ, fast = FALSE, verbose = TRUE)
#         # Note
#         # Group by cells with TSS enrichment scores in two groups.
#         SeuratOBJ$high.tss <- ifelse(SeuratOBJ$TSS.enrichment > 2, 'High', 'Low')
#         #colnames(SeuratOBJ@meta.data)
#         TSSPlot(SeuratOBJ, group.by = 'high.tss') + NoLegend()
#     }
#     , error = function(e) {print('An error occurred. Check the annotation or you probably have ATAC quality loss.') })
#     
#     # Add blacklist ratio and fraction of reads in peaks
#     SeuratOBJ$blacklist_fraction <- FractionCountsInRegion(
#         object = SeuratOBJ,
#         assay = 'ATAC',
#         regions = blacklist_hg38
#     )
#     # Add blacklist ratio and fraction of reads in peaks
#     SeuratOBJ$pct_reads_in_peaks <- SeuratOBJ$atac_peak_region_fragments / SeuratOBJ$atac_fragments * 100
#     SeuratOBJ$blacklist_ratio <- SeuratOBJ$blacklist_fraction / SeuratOBJ$atac_peak_region_fragments
#     # Plot reads in black ratio
#     
#     p1 <- VlnPlot(SeuratOBJ, features = c("pct_reads_in_peaks","blacklist_ratio"), ncol = 2)
#     base_name <- paste0(s_sample, '_None_QC')
#     png_file <- paste0(base_name,'_reads_in_peaks.pdf')
#     png_name <- here('plots/01_preprocessing_QC', png_file)  
#     ggsave(p1, filename = png_name, height = 4, width = 4)
#     message('Black ratio reads in peaks plot saved!')  
#     
#     p1 <- VlnPlot(SeuratOBJ, features = c("nCount_ATAC", "nFeature_ATAC", "nucleosome_signal", "TSS.enrichment"), ncol = 4)
#     base_name <- paste0(s_sample, '_None_QC')
#     png_file <- paste0(base_name,'_ATAC.png')
#     png_name <- here('plots/01_preprocessing_QC', png_file)  
#     ggsave(p1, filename = png_name, height = 4, width = 7)
#     message('ATAC quality general plots saved!')  
#     
# }

# message('ATAC quality processed successfully!')
# message(paste0('Saving new Seurat with ATAC attached'))
# rds_name <- here('processed-data/01_preprocessing_QC', paste0(s_sample,'filteredM',i_filtering_method,'_ATAC.rds'))
# saveRDS(SeuratOBJ, file = rds_name)


# Chunk of code to filter by conventional standar and custom profile moved to: 01b_preprocessing_GEX_ATAC.R script


############ Reproducibility information ####################

library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
#"2023-04-04 12:42:26 EDT"
proc.time()
options(width = 120)
session_info()

