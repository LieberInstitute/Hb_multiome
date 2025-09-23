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

## code assumes peak_called_in column exists and contains a comma-separated list of cell types
global_peaks_ranges_cellType <- function(
        pb_obj = SeuratOBJ_pb,
        atac_assay_name = PSEUDO_ATAC_ASSAY,
        cluster_name
) {
    
    # Use the pseudobulk ATAC assay that exists
    #stopifnot(atac_assay_name %in% Assays(pb_obj))
    DefaultAssay(pb_obj) <- atac_assay_name
    
    # subset to cluster
    seurat_subset <- subset(pb_obj, idents = cluster_name)
    n_samples <- ncol(seurat_subset)
    message("Processing ", cluster_name, " (", n_samples, " pseudobulk samples)")
    if (n_samples < 3) stop("Too few samples for cluster ", cluster_name, "")
    
    # peaks for assay
    peaks_gr <- granges(seurat_subset[[atac_assay_name]])
    peaks_gr$peak_id <- Signac::GRangesToString(peaks_gr)
    
    if (!"peak_called_in" %in% colnames(mcols(peaks_gr))) {
        warning("No peak_called_in metadata; returning all peaks in subset")
        return(peaks_gr)
    }
    
    peaks_map <- as.data.frame(peaks_gr) |>
        transmute(peak_id = peak_id,
                  peak_called_in = as.character(peak_called_in)) |>
        mutate(peak_called_in = str_split(peak_called_in, "\\s*,\\s*")) |>
        tidyr::unnest(peak_called_in) |>
        mutate(peak_called_in = trimws(peak_called_in)) |>
        filter(!is.na(peak_called_in), peak_called_in != "") |>
        distinct()
    #head(peaks_map)
    
    cluster_peaks <- peaks_map |> 
        filter(peak_called_in == cluster_name) |> 
        pull(peak_id) |> 
        intersect(rownames(seurat_subset[[atac_assay_name]]))
    
    if (length(cluster_peaks) == 0) stop("No peaks for ", cluster_name)
    
    # drop ultra-sparse peaks in this cluster (improves stability)
    counts_mat <- GetAssayData(
        seurat_subset, 
        assay = atac_assay_name, 
        layer = "counts")[cluster_peaks, , drop = FALSE]
    
    min_cells_sub <- max(3, floor(0.05 * n_samples))   # ≥5% or at least 3
    keep_peaks_sub <- rownames(counts_mat)[Matrix::rowSums(counts_mat > 0) >= min_cells_sub]
    if (length(keep_peaks_sub) == 0) stop("No peaks pass support filter in ", cluster_name, ".")
    
    # preserves annotation
    peak_ranges <- Signac::StringToGRanges(keep_peaks_sub)
    
    return(peak_ranges)
    
}


## conversion from Seurat to SingleCellExperiment (SCE) is a standard and 
## necessary step for using Bioconductor packages like edgeR and limma
global_convert_atac_subset_to_sce <- function(
        SeuratOBJ_pb, 
        PSEUDO_ATAC_ASSAY, 
        peaks_to_keep, 
        meta
) {
    
    # subset by peaks
    Seurat_subset <- subset(
        SeuratOBJ_pb[[PSEUDO_ATAC_ASSAY]],
        features = peaks_to_keep
    )
    Seurat_subset
    
    ## convert object into sce with meta.data
    counts_mat <- GetAssayData(Seurat_subset, layer = "counts")
    
    Seurat_subset2 <- CreateSeuratObject(
        counts = counts_mat,
        assay = PSEUDO_ATAC_ASSAY,
        meta.data = meta
    )
    Seurat_subset2
    # > Seurat_subset2
    # An object of class Seurat 
    # 5224 features across 169 samples within 1 assay 
    # Active assay: ATAC_macs2_merged_pseudo (5224 features, 0 variable features)
    # 1 layer present: counts
    
    sce_pb <- as.SingleCellExperiment(
        Seurat_subset2,
        assay = PSEUDO_ATAC_ASSAY
    )
    
    # copy Seurat's "logcounts" layer (if it existed) into SCE
    # if ("logcounts" %in% Layers(Seurat_subset2[[PSEUDO_ATAC_ASSAY]])) {
    #     logcounts(sce_pb) <- GetAssayData(Seurat_subset2, assay = PSEUDO_ATAC_ASSAY, layer = "logcounts")
    # }
    
    # Note. actual voomLmFit still use the raw counts - no need to compute logcounts / (for QC or visualization only)
    if (!"logcounts" %in% assayNames(sce_pb)) {
        # create logcounts from raw counts
        dge <- DGEList(counts = assay(sce_pb, "counts"))
        dge <- calcNormFactors(dge, method = "TMM") # gives scaling factors for library sizes
        lcpm <- edgeR::cpm(dge, log = FALSE, prior.count = 0, normalized.lib.sizes = TRUE)
        logcounts(sce_pb) <- log2(lcpm + 1)
    }
    
    #dim(sce_pb)
    #colData(sce_pb) #Endo test: [1] 5224  169
    #rowData(sce_pb)
    #table(sce_pb$cellType)
    #table(sce_pb$donor)
    
    return(sce_pb)
}


global_rebuild_rna_normalization <- function(
        SeuratOBJ,
        assay_name = "RNA"
) {
    
    message("Starting normalization ... ")
    
    SeuratOBJ <- CreateSeuratObject(counts = pb_rna_counts, assay = assay_name)
    SeuratOBJ <- NormalizeData(SeuratOBJ) # RNA log-normalize per cluster
    SeuratOBJ <- FindVariableFeatures(SeuratOBJ, selection.method = "vst") 
    SeuratOBJ <- ScaleData(SeuratOBJ, features = rownames(SeuratOBJ))

    message("RNA normalization done!")
    
    return(SeuratOBJ)
    
}

    
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
    
    message("ATAC normalization done!")
    
    return(SeuratOBJ)
    
}


## Set Seurat identities from a metadata column
global_set_idents_from_meta <- function(
        seurat_obj, 
        meta_col, 
        level_order = NULL, 
        na_fill = "Unknown"
) {
    
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

