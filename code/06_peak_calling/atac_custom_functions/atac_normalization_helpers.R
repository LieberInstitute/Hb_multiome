########################################################################
## Helper for atac analysis 
##
## Authors. CSC
## Date. August 29, 2025
## Note. function based on:
# Seurat                      * 5.0.1      2023-11-17 [1] CRAN (R 4.3.2)
# SeuratObject                * 5.1.0      2025-04-22 [1] CRAN (R 4.3.2)
# Signac                      * 1.14.0     2024-08-21 [1] CRAN (R 4.3.2)
# module load conda_R/4.3.x
# Note: incompatibility between version could broke the chromatin assays
########################################################################

## Requires:
# library("Seurat")
# library("Signac")

## fast verification
if (!requireNamespace("Signac", quietly = TRUE)) {
    stop("Package 'Signac' is required but not installed.")
}
if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required but not installed.")
}

# Rebuild TF-IDF normalization (with previous thresholds), runSVD & make GC content correction
# based on:
# https://github.com/LieberInstitute/Hb_multiome/blob/a87c3512395b0488af9413174796787412e99567/code/03_pseudobulking/08_harmony_CR_ARCr.R#L99-L105 

global_rebuild_atac_normalization <- function(
        SeuratOBJ,
        assay_name
) {
    
    message("Starting normalization ... ")
    
    SeuratOBJ <- RunTFIDF(
        SeuratOBJ,
        assay = assay_name,
        method = 1,
        scale.factor = 10000
    )
    
    SeuratOBJ <- FindTopFeatures(
        SeuratOBJ,
        assay = assay_name,
        min.cutoff = 'q5',
        verbose = TRUE
    )
    
    SeuratOBJ <- RunSVD(
        SeuratOBJ,
        assay = assay_name
    )
    
    # message("Starting GC content correction ... ")
    # 
    # genome <- BSgenome.Hsapiens.UCSC.hg38
    # # Normalize naming style before running RegionStats
    # seqlevelsStyle(SeuratOBJ) <- "UCSC"
    # # Keep only standard chromosomes
    # SeuratOBJ[[assay_name]] <- SeuratOBJ[[assay_name]][
    #     seqnames(granges(SeuratOBJ)) %in% standardChromosomes(granges(SeuratOBJ)), ]
    # # Suppress scaffolds
    # # suppressWarnings(
    # #     SeuratOBJ <- RegionStats(SeuratOBJ, genome = genome, assay = assay_name)
    # # )
    # 
    # SeuratOBJ <- RegionStats(
    #     object = SeuratOBJ,
    #     assay = assay_name,
    #     genome = genome
    # )
    
    message("Normalization done!")
    
    return(SeuratOBJ)
    
}


# Set Seurat identities from a metadata column

global_set_idents_from_meta <- function(seurat_obj, meta_col, level_order = NULL, na_fill = "Unknown") {
    
    message("Set new Seurat idents ... ")
    
    ## double check level exist on meta-data
    if (!meta_col %in% colnames(seurat_obj@meta.data)) {
        stop("Meta column '", meta_col, "' not found in SeuratOBJ@meta.data")
    }
    # extract target vector
    target_vec <- as.character(seurat_obj[[meta_col]][, 1])
    
    # decide levels and keep appearance order
    if (is.null(level_order)) {
        level_order <- sort(unique(target_vec))
    }
    
    # only update if different from current Idents
    current_idents <- as.character(Idents(seurat_obj))
    if (!identical(current_idents, target_vec)) {
        seurat_obj <- SetIdent(seurat_obj, value = factor(target_vec, levels = level_order))
        message("Idents set from meta column '", meta_col, "'.")
    } else {
        message("Idents already match '", meta_col, "', nothing to do.")
    }
    
    message("Set new idents done!")
    
    return(seurat_obj)
    
}

