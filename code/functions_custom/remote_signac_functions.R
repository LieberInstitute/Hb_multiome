########################################################################
##
## FUNCTIONS TO HANDLE SIGNAC OBJECTS 
##
## Authors. CSC
## Date. May 18th, 2023
##
########################################################################

## not used longer (csc. oct.1.2023), integrated in 04_preprocessing_ATAC.R
# get_ATAC_obj  <- function(seuratOBJ, s_bc_mtx, s_frag, ann) { 
# 
#     # Pull the full paths to build the chromatic assay.
#     # NOTE. Since v2.0.0 to v2.0.2 is used the same name standard. Be careful in future cellranger-arc releases. (Sep, 2023)
# 
#     # Pull ATAC data from barcode mtx
#     message(s_bc_mtx)
#     mtx <- Read10X_h5(s_bc_mtx)
#     atac_counts <- mtx$Peaks
#     message('ATAC counts loaded successfully')
# 
#     # Convert a genomic coordinate string to a GRanges object ands assign it to the ATAC counts mtx
#     if(!("GRanges" %in% class(ann))) stop("Annotation class not found")
#     grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))      # length(grange.counts)
#     grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
#     atac_counts <- atac_counts[as.vector(grange.use), ]
#     message('ATAC annotation attached successfully')
# 
#     # Ensure seurat_obj is a Seurat object
#     if (!("Seurat" %in% class(seuratOBJ))) { stop("Seurat object not found") }
# 
#     # Create chromatin assay with annotations
#     chrom_assay <- CreateChromatinAssay(counts = atac_counts,
#                                         sep = c(":", "-"),
#                                         fragments = s_frag,
#                                         annotation = ann)
#     if(!("ChromatinAssay" %in% class(chrom_assay))) stop("Chromatin assay not found")
#     message('Chromatin assay completed successfully')
# 
#     seuratOBJ[["ATAC"]] <- chrom_assay
# 
#     # Annotations of the object are set
#     Annotation(seuratOBJ[["ATAC"]]) <- ann
#     message('Chromatin assay attached to seurat object successfully')
# 
#     return(seuratOBJ)
# }


# This is an old version (v1) of the function get_create_atac_objs()
get_create_atac_objs  <- function(seuratOBJ, s_bc_mtx, s_frag, ann, b_save_signal_obj = FALSE) {

    # Pull ATAC data from barcode mtx
    mtx <- Read10X_h5(s_bc_mtx)
    atac_counts <- mtx$Peaks
    message('ATAC counts loaded successfully')
    # Convert a genomic coordinate string to a GRanges object ands assign it to the ATAC counts mtx
    if(!("GRanges" %in% class(ann))) stop("Annotation class not found")
    grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))      # length(grange.counts)
    grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
    atac_counts <- atac_counts[as.vector(grange.use), ]
    message('ATAC counts annotation attached successfully')
    
    # Ensure seurat_obj is a Seurat object
    if (!("Seurat" %in% class(seuratOBJ))) { stop("Seurat object not found") } 
#    UpdateSeuratObject(seuratOBJ)
        
    # Create chromatin assay
    chrom_assay <- CreateChromatinAssay(counts = atac_counts, 
                                        sep = c(":", "-"), 
                                        fragments = s_frag, 
                                        annotation = ann)
    if(!("ChromatinAssay" %in% class(chrom_assay))) stop("Chromatin assay not found")
    message('Chromatin assay completed successfully')
    print(chrom_assay)

    seuratOBJ[["ATAC"]] <- chrom_assay
    # Annotations of the object are set
    Annotation(seuratOBJ[["ATAC"]]) <- ann
    message('Chromatin assay attached to seurat object successfully')
    
    # Saving the object 
    if (b_save_signal_obj) {
        # Check if processed_data directory exists, if not create it
        if (!dir.exists(here("processed-data"))) {
            dir.create(here("processed-data"))
        }
        # Save the Seurat object
        # New Seurat functions fails ... ----- need to be check why SaveH5Seurat library is not working well
        # seuratName <- paste0(seuratName, '.h5Seurat')
        # UpdateSeuratObject(seur_obj)
        # SaveH5Seurat(seur_obj, filename = here("processed-data", seuratName), overwrite = TRUE)
        seuratName = toString(unique(seuratOBJ@meta.data$orig.ident))
        seuratName <- paste0(seuratName, '.ATAC.RDS')
        write_rds(seuratOBJ, here("processed-data", seuratName), compress = ('gz'))
        # save(seuratOBJ, file = here("processed-data", seuratName))
        message('Process completed successfully. Seurat object with ATAC saved as ', seuratName)
    }
    print(seuratOBJ)
    #print(head(seuratOBJ, n = 3))
    return(seuratOBJ) 
}


f_TSS_enrichment_NS  <- function(seuratOBJ) {
    seuratOBJ$blacklist_fraction <- FractionCountsInRegion(
        object = seuratOBJ,
        assay = 'ATAC',
        regions = blacklist_hg38
    )
    # add blacklist ratio and fraction of reads in peaks
    seuratOBJ$pct_reads_in_peaks <- seuratOBJ$atac_peak_region_fragments / seuratOBJ$atac_fragments * 100
    seuratOBJ$blacklist_ratio <- seuratOBJ$blacklist_fraction / seuratOBJ$atac_peak_region_fragments
    ##### VlnPlot(seuratOBJ, features = c("pct_reads_in_peaks","blacklist_ratio"), ncol = 2)
    # Classified the TSS enrichment scores in two groups
    seuratOBJ$high.tss <- ifelse(seuratOBJ$TSS.enrichment > 2, 'High', 'Low')
    head(seuratOBJ, n =3)
    colnames(seuratOBJ@meta.data)
    
    return(seur_obj) 
}

get_basic_stats_ATAC <- function(seuratOBJ) {
    ## Print basic stats from a Seurat ATAC assay
    ## INPUT: 
    ## @seuratOBJ: Seurat obj 
    ## @sname: string describing the object (optional / provide short-name)
    ## OUTPUT:
    ## @s: long string with basic statistics about the ATAC Seurat Assay
    
    # Get the name project (ex. pbmc3k)
    tryCatch({
        # Get the name project (ex. pbmc3k)
        s_ProjName = toString(unique(seuratOBJ@meta.data$orig.ident))
        # Initialize metrics
        s = paste('\nATAC ASSAY #################')        
        s = paste(s, '\n\nObject name:  ', s_ProjName, sep = '') 
        s = paste(s, '\nCells:  ', toString(table(seuratOBJ$orig.ident)))        
        # Discontinuous sample quantile: inverse of empirical distribution function
        # tmp = quantile(seuratOBJ$nCount_ATAC, na.rm = TRUE, type = 1)       #  if na.rm = true, any NA and NaN's are removed
        # Continuous sample quantile types 4: linear interpolation of the Empirical Cumulative Distribution Function 
        tmp = quantile(seuratOBJ$nCount_ATAC, na.rm = TRUE, type = 4)       #  if na.rm = true, any NA and NaN's are removed
        s = paste(s, '\nQuantiles in nCount_ATAC:\n', toString(tmp))
        tmp = quantile(seuratOBJ$nFeature_ATAC, na.rm = TRUE, type = 4) 
        s = paste(s, '\nQuantiles in nFeature_ATAC:\n', toString(tmp))
        s = cat(paste(s, '\n', sep = ''))           
        return(s) }
        , error = function(e) {print('An error ocurred. Verify you have an ATAC object active') })
}
