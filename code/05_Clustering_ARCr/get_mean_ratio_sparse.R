## Adaptation from Deconvobuddies
# minimal changes to avoid dense conversion, with sparse-safe mean computation
# https://github.com/LieberInstitute/DeconvoBuddies/blob/1d16d86a52adbba0e5492625a0c1c3d376d9c6bd/R/get_mean_ratio.R

library("purrr")
library("dplyr")
library("tidyverse")
library("stringr")
library("matrixStats")

get_mean_ratio_sparse <- function(sce, cellType_col, assay_name = "logcounts", gene_ensembl = NULL, 
                                  gene_name = NULL) {
    cellType.target <- NULL
    cellType <- NULL
    ratio <- NULL
    rank_ratio <- NULL
    anno_ratio <- NULL
    
    stopifnot(cellType_col %in% colnames(colData(sce)))
    stopifnot(assay_name %in% names(SummarizedExperiment::assays(sce)))
    
    cell_types <- unique(sce[[cellType_col]])
    names(cell_types) <- cell_types
    ct_table <- table(sce[[cellType_col]])
    if (any(ct_table < 10)) 
        warning("One or more cell types has < 10 cells, this may result in unstable marker genes results. 
                Check details of get_mean_ratio() for more info")
    
    sce_assay <- assay(sce, assay_name)  # Here, NO dense conversion
    cell_means <- map(cell_types, ~as.data.frame(Matrix::rowMeans(sce_assay[, 
                                                                            sce[[cellType_col]] == .x])))
    cell_means <- do.call("rbind", cell_means)
    colnames(cell_means) <- "mean"
    cell_means$cellType <- rep(cell_types, each = nrow(sce))
    cell_means$gene <- rep(rownames(sce), length(cell_types))
    
    ratio_tables <- map(cell_types, ~.get_ratio_table(.x, sce, 
                                                      sce_assay, cellType_col, cell_means))
    ratio_tables <- rename(mutate(do.call("rbind", ratio_tables), 
                                  anno_ratio = paste0(cellType.target, "/", cellType, ": ", 
                                                      base::round(ratio, 3))), cellType.2nd = cellType, 
                           mean.2nd = mean, MeanRatio = ratio, MeanRatio.rank = rank_ratio, 
                           MeanRatio.anno = anno_ratio)
    
    if (!is.null(gene_ensembl)) {
        if (gene_ensembl %in% colnames(SummarizedExperiment::rowData(sce))) {
            ratio_tables$gene_ensembl <- SummarizedExperiment::rowData(sce)[ratio_tables$gene, 
            ][[gene_ensembl]]
        } else {
            warning("'", gene_ensembl, "' not in col rowData, gene_ensembl not included in output")
        }
    }
    if (!is.null(gene_name)) {
        if (gene_name %in% colnames(SummarizedExperiment::rowData(sce))) {
            ratio_tables$gene_name <- SummarizedExperiment::rowData(sce)[ratio_tables$gene, 
            ][[gene_name]]
        } else {
            warning("'", gene_name, "' not in col rowData, gene_name not included in output")
        }
    }
    return(ratio_tables)
}
.get_ratio_table <- function(x, sce, sce_assay, cellType_col, cell_means) {
    # RCMD Fix
    mean.target <- NULL
    gene <- NULL
    ratio <- NULL
    cellType.target <- NULL
    cellType <- NULL
    
    # filter target median != 0
    #median_index <- matrixStats::rowMedians(sce_assay[, sce[[cellType_col]] == x]) != 0
    median_index <- matrixStats::rowMedians(as.matrix(sce_assay[, sce[[cellType_col]] == x])) != 0
    
    # message("Median == 0: ", sum(!median_index))
    # filter for target means
    target_mean <- cell_means[cell_means$cellType == x, ]
    target_mean <- target_mean[median_index, ]
    colnames(target_mean) <- c("mean.target", "cellType.target", "gene")
    
    nontarget_mean <- cell_means[cell_means$cellType != x, ]
    
    ratio_table <- dplyr::left_join(target_mean, nontarget_mean, by = "gene") |>
        mutate(ratio = mean.target / mean) |>
        dplyr::group_by(gene) |>
        arrange(ratio) |>
        dplyr::slice(1) |>
        dplyr::select(gene, cellType.target, mean.target, cellType, mean, ratio) |>
        arrange(-ratio) |>
        dplyr::ungroup() |>
        mutate(rank_ratio = dplyr::row_number())
    
    return(ratio_table)
}