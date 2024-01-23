########################################################################
## Measure Quality Controls for GEX and ATAC assays from CellRanger-ARC data using Seurat & Signac packages
##
## Authors. CSC / HT
## Date. April 21st, 2023
## Last.Adaptation: Jan.22, 2024
##
## Input: Truncated H5, meta-data.cvs and fragments.tvs
## Output: rds Seurat objects and plots
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
source(here("code/functions_custom", "remote_plot_functions.R"))


########################    Initials ########################  

## commandArgs scans the arguments which have been supplied when the current R script was invoked (from shell sh)
sample_tmp <- commandArgs(trailingOnly = TRUE)
#sample_tmp <- args[1]
# testing
#sample_tmp <- 'I_am_an_error_file,dog'  # testing error file and tissue
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
b_get_GEX_plots <- FALSE         
#if (b_get_GEX_plots) { source(here("code/functions_custom", "remote_plot_functions.R")) }

# Remove mitochondrial levels (by sample dynamically)
b_get_filtered_GEX <- TRUE      
if (b_get_filtered_GEX) {
    source(here("code/functions_custom", "remote_filtering_functions.R"))    # Filter Seurat assay by GEX
    i_filtering_method <- 1         # Pick up probabilities method to remove mito levels   
}

# Create attac and get quality control measures
b_get_ATAC <- TRUE
b_get_ATAC_QC <- TRUE


########  ################################################# ######## 
########  1.  Create a Seurat object containing the RNA     
########      Recommended 30G of free_mem to 3k-10k cells 
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
########        2. Get Visualizations for the GEX             ######## 
########  ################################################# ######## 

# Build and plot:
#       UMI, Genes, MITO and RIBO violin plots
#       Number of cells per sample
#       UMI/transcripts per cell
#       Distribution of genes per cell (histogram)

# Before QC any assay
if (b_get_GEX_plots) { 

    # Load a valid Seurat object if available
    #SeuratOBJ <- load_seurat_obj(SeuratOBJ, s_seurat_name)      

    # To check what are the combined datasets run next command
    #table(pbmc10k.combined.GEX$orig.ident)
    
    # Violin plot with UMIs, Genes, ^MT and RIBO levels
    p1 <- get_Vplots_main_GEX(SeuratOBJ)  
    
    # Assign the layer to the df
    df_genes_per_cell <- as.data.frame(SeuratOBJ[[]])
    head(df_genes_per_cell, n = 3)
    
    # Plot number of cells per sample (Geom_bar)
    p2 <- get_plt_cells_by_sample(df_genes_per_cell)
    
    # Plot the number UMIs/transcripts per cell
    p3 < get_plt_UMIs_per_cell(df_genes_per_cell)
    
    # Plot the distribution of genes detected per cell (density or quantiles)
    p4 <- get_plt_genes_per_cell_density(df_genes_per_cell)
    p1 <- get_plt_genes_per_cell_boxplot(df_genes_per_cell)
    # Correlation btw genes and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
    p1 <- get_plt_UMIS_genes_MT_geomlm(df_genes_per_cell, 500, 500)

}

    
########  ################################################# ######## 
##        3. Preprocess Seurat Object based on the Mitochondrial percentage in the GEX assay
##              METHOD 1: (M1p)  calculate cut-off based on probabilities
##              METHOD 2: (M2sd) calculate cut-off at +/-2 Standard Deviations (95%)
########  ################################################# ######## 


if (b_get_filtered_GEX) {
    # Preferred method 1. filter by probabilities
    
    # Validate or Load Seurat object from disk if available
    # Note. GEX assay need to be active

    if (i_filtering_method==1) {
        # filter by probabilities: M1
        SeuratOBJ.filtered.M1 <- get_seurat_GEX_filteringM1p(SeuratOBJ, FALSE)
        # Some stats and plot correlation 
        get_basic_stats_GEX(SeuratOBJ.filtered.M1, s_tissue)
        df_genes_per_cell <- as.data.frame(SeuratOBJ.filtered.M1[[]])
        
    } else {    

        # filter by SD: M2
        SeuratOBJ.filtered.M2 <- get_seurat_GEX_filteringM2sd(SeuratOBJ, FALSE)
        # Some stats and plot correlation 
        get_basic_stats_GEX(SeuratOBJ.filtered.M2, s_tissue)
        df_genes_per_cell <- as.data.frame(SeuratOBJ.filtered.M2[[]])
    }

    # Plot correlation between genes and number of UMIs to visualize strong presence of cells with low numbers of genes/UMIs
    p1 <- get_plt_UMIS_genes_MT_geomlm(df_genes_per_cell, 500, 500) #+
    # plot_annotation(paste0(s_sample,': Genes/UMIs correlation based on method1')) &
    #     plot_annotation(tag_levels = '1') &
    #     theme(plot.tag = element_text(color = "blue", size = 10)) 
    p1
    stitle <- paste0(s_sample,'corr_filteredM',i_filtering_method,'.png')
    plt_name <- here('plot/01_preprocessing_QC', stitle)  
    ggsave(p1, filename = plt_name)
    
}

message('Seurat filtered successfully!')
message(paste0('Saving new Seurat filtered.'))
rds_name <- here('processed-data/01_preprocessing_QC', paste0(s_sample,'filteredM',i_filtering_method,'.rds'))
saveRDS(SeuratOBJ, file = rds_name)


########  ################################################# ############ 
########  4. Create ATAC assay and attache it to Seurat object     ##### 
########  ################################################# ############ 

if (b_get_ATAC) {     

    # Validate or Load Seurat object from disk if available
    #SeuratOBJ <- load_seurat_obj(SeuratOBJ, s_seurat_name) 
    #head(SeuratOBJ, n = 3)
    
    # General use: Create gene annotations for hg38 and extract gene annotations from EnsDb
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
    # show(annotations) / # View(head(annotations,n=5)) /# names(genomeStyles('Homo_sapiens'))
    seqlevelsStyle(annotations) <- "UCSC"
    length(extractSeqlevelsByGroup(species = 'Homo_sapiens', style = 'UCSC', group = 'all')) #group: sex, auto, circular
    genome(annotations) <- "hg38"
    
    # Create the chromatin assay with annotations and attach it to the seurat object 
    start.time = Sys.time()
    SeuratOBJ <- get_create_atac_objs(SeuratOBJ, 
                                      s_bc_mtx, s_frag_namefile, 
                                      annotations, TRUE) 
    end.time = Sys.time()
    print(end.time - start.time)        # 15s / 40G free_mem / 3k Cells / Object.size 806.2 Mb
                                        # 38s / 40G free_mem / 10k cells / 
    SeuratOBJ@meta.data$`orig.ident`[1]
    #f_InspectSeurat(SeuratOBJ)
    get_basic_stats_ATAC(SeuratOBJ)

}
    
message('ATAC attached successfully!')



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
if (b_get_ATAC_QC) {

    # Validate or Load Seurat object from disk if available
    #SeuratOBJ <- load_ATAC_obj(SeuratOBJ, s_seurat_name)
    
    DefaultAssay(SeuratOBJ) <- "ATAC"
    # Calculate fragment size distribution
    FragmentHistogram(object = SeuratOBJ)
    # Calculate the strength of the nucleosome signal per cell
    SeuratOBJ <- NucleosomeSignal(SeuratOBJ)
    head(SeuratOBJ, n = 3)
    # Group by cells with high or low nucleosome signal strength. 
    SeuratOBJ$nucleosome_group <- ifelse(SeuratOBJ$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
    FragmentHistogram(object = SeuratOBJ, group.by = 'nucleosome_group')
    head(SeuratOBJ, n = 3)
    
    # Calculate the "Transcription Start Site (TSS)" enrichment score for each cell, as defined by ENCODE.
    tryCatch( {
        
        # Important NOTES: https://github.com/stuart-lab/signac/issues/374 
        #           1. Some times jumps an issue of ** memory allocation ** when "fast=FALSE", FALSE argument is mandatory to
        #              compute the TSS enrichment scores and visualize with TSSPlot(). You need more memory
        #           2. Error in `colnames<-`(`*tmp*`, value = seq_len(length.out = region.width) -  : attempt to set 'colnames' ...
        #              This is a vague message that would happen if no fragments are found in the set of TSS regions. 
        #              You could double-checking that the correct gene annotations is being used or you have a low ATAC quality. 
        
        SeuratOBJ <- TSSEnrichment(SeuratOBJ, fast = FALSE)
        # Note
        # Group by cells with TSS enrichment scores in two groups.
        SeuratOBJ$high.tss <- ifelse(SeuratOBJ$TSS.enrichment > 2, 'High', 'Low')
        #colnames(SeuratOBJ@meta.data)
        TSSPlot(SeuratOBJ, group.by = 'high.tss') + NoLegend()
    }
    , error = function(e) {print('An error occurred. Check the annotation or you probably have ATAC quality loss.') })
    
    # Add blacklist ratio and fraction of reads in peaks
    SeuratOBJ$blacklist_fraction <- FractionCountsInRegion(
        object = SeuratOBJ,
        assay = 'ATAC',
        regions = blacklist_hg38
    )
    # Add blacklist ratio and fraction of reads in peaks
    SeuratOBJ$pct_reads_in_peaks <- SeuratOBJ$atac_peak_region_fragments / SeuratOBJ$atac_fragments * 100
    SeuratOBJ$blacklist_ratio <- SeuratOBJ$blacklist_fraction / SeuratOBJ$atac_peak_region_fragments
    # Plot reads in black ratio
    
    VlnPlot(SeuratOBJ, features = c("pct_reads_in_peaks","blacklist_ratio"), ncol = 2)
    
    VlnPlot(SeuratOBJ, features = c("nCount_ATAC", "nFeature_ATAC", "nucleosome_signal", "TSS.enrichment"), ncol = 4)

}

message('ATAC quality processed successfully!')
message(paste0('Saving new Seurat with ATAC attached'))
rds_name <- here('processed-data/01_preprocessing_QC', paste0(s_sample,'filteredM',i_filtering_method,'_ATAC.rds'))
saveRDS(SeuratOBJ, file = rds_name)

# Clean dataset

#get_preprocessed_subset()          # functions_custom/remote_filtering_functions.R    

############ Reproducibility information ####################

library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
#"2023-04-04 12:42:26 EDT"
proc.time()
options(width = 120)
session_info()

