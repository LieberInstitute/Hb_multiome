########################################################################
## Differential Accessibility (DA) using voomLmFit / pseudobulk multiome assays with merged peaks
##
## Authors. CSC
## Date. Sep 16, 2025
## Recommended resources mem=30GB
########################################################################

suppressPackageStartupMessages({
    library("edgeR")     # for DGEList, filterByExpr, voomLmFit
    library("limma")     # for contrasts.fit, eBayes, topTable, plot helpers
    library("Matrix")    # if counts are sparse; otherwise base matrix is fine
    library("data.table")
})
library("Seurat")
library("spatialLIBD")
library("GenomicRanges")
library("ggplot2")
library("ggrepel") 
library("Seurat")
library("Signac")
library("dplyr")
library("stringr")
library("here")


# Testing spearman at 5e4 on macs2 peaks 
resolution_level = "Mid"
# p_met = "spearman"
# w_size = "5e5"
sig_thresh <- 0.1    # FDR cutoff
lfc_thresh <- 0.25   # logFC cutoff / log2

## ATAC function's helper used globally
source(here("code", "06_peak_calling", "multiome_custom_functions", "multiome_idents_normalization_helper.R"))

# Inputs:
# counts:   matrix of raw counts, rows = peaks, cols = pseudobulk samples
#           (e.g., aggregated by cluster x donor)
# coldata:  data.frame with sample metadata, nrow = ncol(counts)
#           must include at least: sample_id, cluster_id, group, donor
# peaks_df: data.frame/GRanges of peaks for annotation (optional)

# Check/create directories
inputRDS_Dir <- here(
    "processed-data",
    "06_peak_calling",
    "12_pseudobulk_MACS2"
)
# these are mac2 peaks merged and normalized 
Seurat_base_name <- "Mid_pseudobulk.spearman.5e5_merged_peaks.rds"

output_Dir <- here(
    "processed-data",
    "06_peak_calling",
    "17_pseudobulk_DARs_MACS2_reduced_voomLmFit"
)
plotDir <- here(
    "plots",
    "06_peak_calling",
    "17_pseudobulk_DARs_MACS2_reduced_voomLmFit"
)
PSEUDO_ATAC_ASSAY <- "ATAC_macs2_merged_pseudo"

## Check directories
if (!dir.exists(plotDir)) {
    dir.create(plotDir)
}
if (!dir.exists(output_Dir)) {
    dir.create(output_Dir)
}


##==============================================================================

## extract peaks cell-type specific
## code assumes peak_called_in column exists and contains a comma-separated list of cell types
peaks_ranges_cellType <- function(
        pb_obj = SeuratOBJ_pb,
        atac_assay_name = PSEUDO_ATAC_ASSAY,
        cluster_name #  = "Endo"
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
convert_atac_subset_to_sce <- function(
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


registration_stats_enrichment_voomLmFit <-
    function(
        sce_pseudo,
        block_cor,
        covars = NULL,                              # c("age", "sex", "ethnicity")
        var_registration = "registration_variable", # cellType = 18 ct
        var_sample_id = "registration_sample_id",   # donor = 169 samples
        gene_ensembl = NULL,
        gene_name = NULL
        # sce_pseudo = sce_pb
        # covars = covars_vec
        # var_registration = "registration_variable"
        # var_sample_id = "registration_sample_id"
        # block_cor = block_cor
    ) {
        ## For each cluster, test it against the rest
        cluster_idx <- split(
            seq(along = sce_pseudo[[var_registration]]),
            sce_pseudo[[var_registration]]
        )
        
        message(Sys.time(), " computing enrichment statistics")
        eb0_list_cluster <- lapply(cluster_idx, function(x) {
            res <- rep(0, ncol(sce_pseudo))
            res[x] <- 1
            if (!is.null(covars)) {
                res_formula <-
                    eval(str2expression(paste(
                        "~",
                        "res",
                        "+",
                        paste(covars, collapse = " + ")
                    )))
            } else {
                res_formula <- eval(str2expression(paste("~", "res")))
            }
            # res_formula = ~res + age + sex + ethnicity
            m <- model.matrix(res_formula, data = colData(sce_pseudo))
            
            ## voomLmFit, is a shorthand that includes: log2-counts per million (logCPM) -> Weighting -> Linear model fitting
            ## run voomLmFit for the pseudobulked data, referring donor to duplicateCorrelation; 
            ## using an adaptive span (number of genes, based on the number of genes in the dge) for smoothing the mean-variance trend
            res <- limma::eBayes(edgeR::voomLmFit(
                assay(sce_pseudo, "counts"),
                design = m,
                block = sce_pseudo[[var_sample_id]], # ge. donor
                correlation = block_cor, # Use the estimated correlation
                sample.weights = TRUE
            ))
            # if (is.finite(block_cor)) {
            #     res <- limma::eBayes(edgeR::voomLmFit(
            #         assay(sce_pseudo, "counts"),
            #         design = m,
            #         block = sce_pseudo[[var_sample_id]],
            #         correlation = block_cor
            #     ))
            # } else {
            #     res <- limma::eBayes(edgeR::voomLmFit(
            #         assay(sce_pseudo, "counts"),
            #         design = m
            #     ))
            # }
            return(res)
        })
        
        message(Sys.time(), " extract and reformat enrichment results")
        
        ## Extract the p-values
        pvals0_contrasts_cluster <-
            sapply(eb0_list_cluster, function(x) {
                x$p.value[, 2, drop = FALSE]
            })
        rownames(pvals0_contrasts_cluster) <- rownames(sce_pseudo)
        
        ## Extract t-statistics
        t0_contrasts_cluster <- sapply(eb0_list_cluster, function(x) {
            x$t[, 2, drop = FALSE]
        })
        rownames(t0_contrasts_cluster) <- rownames(sce_pseudo)
        
        ## Extract logFC
        logFC_contrasts_cluster <- sapply(eb0_list_cluster, function(x) {
            x$coefficients[, 2, drop = FALSE]
        })
        rownames(logFC_contrasts_cluster) <- rownames(sce_pseudo)
        
        ## Compute FDRs
        fdrs0_contrasts_cluster <-
            apply(pvals0_contrasts_cluster, 2, p.adjust, "fdr")
        
        ## Merge into one data.frame
        results_specificity <-
            spatialLIBD:::f_merge(
                p = pvals0_contrasts_cluster,
                fdr = fdrs0_contrasts_cluster,
                t = t0_contrasts_cluster,
                logFC = logFC_contrasts_cluster
            )
        
        ## Add gene info
        results_specificity$ensembl <-
            rowData(sce_pseudo)[[gene_ensembl]]
        results_specificity$gene <- rowData(sce_pseudo)[[gene_name]]
        
        ## Done!
        return(results_specificity)
    }


########################################################################
## Differential Accessibility Analysis using voomLmFit
## Input: Seurat object with pseudobulk RNA+ATAC assay with merged peaks
########################################################################

## Load Seurat / macs peaks / filtered genes. And make verification

seurat_name <- here(inputRDS_Dir, Seurat_base_name)
SeuratOBJ_pb <- readRDS(here(inputRDS_Dir, Seurat_base_name))

# check-ins
SeuratOBJ_pb
# An object of class Seurat 
# 366933 features across 169 samples within 2 assays 
# Active assay: ATAC_macs2_merged_pseudo (351037 features, 333643 variable features)
# 2 layers present: counts, data
# 1 other assay present: RNA
# 1 dimensional reduction calculated: lsi

# ## get rna counts
# DefaultAssay(SeuratOBJ_pb) <- "RNA"
# rna_counts <- GetAssayData(SeuratOBJ_pb, assay="RNA", layer="data")
# length(rownames(rna_counts)) # [1] 36601


## =============================================================================
## add additional meta-data

# get matrix of peaks: rows = peaks, cols = pseudobulk samples
atac_counts <- GetAssayData(SeuratOBJ_pb, assay = PSEUDO_ATAC_ASSAY, layer ="counts")
length(rownames(atac_counts)) # [1] 351037
head(atac_counts)

## get metadata, include at least: sample_id, cluster_id, group, donor
meta <- SeuratOBJ_pb@meta.data
colnames(meta)
head(meta)

## add meta-data in the expected format
meta$sample_id <- row.names(meta) # ge. LHb.1_S05-Hb-r
meta$donor <- sapply(strsplit(row.names(meta), "_", fixed = TRUE),  "[[", 2) # ge. S05-Hb-r
table(meta$donor)
# S03-Hb-r S04-Hb-r S05-Hb-r S06-Hb-r S07-Hb-r S08-Hb-r S09-Hb-r S10-Hb-r 
# 17       15       16       18       18       17       17       18 
# S11-Hb-r S12-Hb-r 
# 17       16

## add additional meta-data
meta$sex <- "M"
meta$ethnicity <- "EA.CAUC" # waiting confirmation KDM
meta$sex[grepl("S03|S08|S11", meta$donor)] <- "F"
meta$ethnicity[grepl("S04|S05", meta$donor)] <- "AA"
meta$age[grepl("S03", meta$donor)] <- 41.3
meta$age[grepl("S04", meta$donor)] <- 40.57
meta$age[grepl("S05", meta$donor)] <- 46.72
meta$age[grepl("S06", meta$donor)] <- 33.39
meta$age[grepl("S07", meta$donor)] <- 48.88
meta$age[grepl("S08", meta$donor)] <- 37.33
meta$age[grepl("S09", meta$donor)] <- 39.98
meta$age[grepl("S10", meta$donor)] <- 47.95
meta$age[grepl("S11", meta$donor)] <- 65
meta$age[grepl("S12", meta$donor)] <- 57.7

SeuratOBJ_pb@meta.data <- meta
head(meta)


## Ensure we are contrasting by cellType
if (!"cellType" %in% colnames(SeuratOBJ_pb@meta.data)) {
    SeuratOBJ_pb$cellType <- sub("^[^_]+_", "", SeuratOBJ_pb$orig.ident)  # keep part after first underscore
}
length(SeuratOBJ_pb$cellType)
# [1] 169

Idents(SeuratOBJ_pb) <- "cellType"
cluster_ids <- levels(SeuratOBJ_pb)
unique(Idents(SeuratOBJ_pb))

message("Pseudobulk groups (cell-types):")
cluster_ids

stopifnot(all(colnames(atac_counts) == meta$sample_id))

# ====/


message("Starting registration_stats_enrichment_voomLmFit ... ")

## refs:
## current version https://github.com/LieberInstitute/spatialLIBD/blob/40da043d0235e01a12a7f52a0b367d3850bad9e8/R/registration_stats_enrichment.R#L40
## My previous manual implementation: https://github.com/LieberInstitute/Habenula_Visium/blob/2d21e39f51c9e46ebddbcf57f959f70c27d78678/code/05_brain_area_differential_expression/05_pseudobulk_DEG_contrast.R#L176-L255 

# Inputs for atac
# logcounts(sce_pseudo) peaks (it fits limma on log-scale data)
# var_registration: in my case cellType
# var_sample_id: blocking by donor factor
# covars:age, sex and ethnicity

DefaultAssay(SeuratOBJ_pb) <- PSEUDO_ATAC_ASSAY
# Get all peaks in the assay
all_peaks <- granges(SeuratOBJ_pb)
length(all_peaks) # [1] 351037
FDR_thr = 0.20

for (clus in cluster_ids) {
    # test: 
    # clus = "Endo"
    message("Subset peaks for [", clus, "] cellType")
    
    # Get peaks cellType specific
    peak_ranges_ct <- peaks_ranges_cellType(
        SeuratOBJ_pb,
        PSEUDO_ATAC_ASSAY,
        cluster_name = clus
    ) 
    length(peak_ranges_ct) # 5224
    
    # Get the assay’s ranges and rownames
    ap_gr   <- granges(SeuratOBJ_pb[[PSEUDO_ATAC_ASSAY]])
    # GRanges object with 351037 ranges and 2 metadata columns:
    # seqnames          ranges strand |        revmap
    # <Rle>       <IRanges>  <Rle> | <IntegerList>
    # [1]       chr1   181329-181534      * |             1

    ap_rows <- rownames(SeuratOBJ_pb[[PSEUDO_ATAC_ASSAY]])
    # [1] "chr1-181329-181534"   "chr1-191217-191619"   "chr1-629146-629354"  
    # [4] "chr1-629811-630032"   "chr1-630189-630389"   "chr1-632189-632410"  

    # Find exact matches by genomic coordinates (use type="any" for overlaps)
    hits <- findOverlaps(ap_gr, peak_ranges_ct, type = "equal")
    idx  <- queryHits(hits)
    # [1]      4      5      7     53     60     73    116    168    189    226
    
    # Map hits to peaks IDs actually present in the assay
    peaks_to_keep <- ap_rows[idx]
    peaks_to_keep <- unique(as.character(peaks_to_keep))
    # [1495] "chr5-72307684-72308762"   "chr5-72319550-72320868"  
    # [1497] "chr5-72815810-72817418"   "chr5-73119548-73120878"  
    
    # checks before subsetting seurat
    stopifnot(length(peaks_to_keep) > 0)
    stopifnot(all(peaks_to_keep %in% ap_rows))
    stopifnot(sum(peaks_to_keep %in% rownames(
        GetAssayData(SeuratOBJ_pb, assay = PSEUDO_ATAC_ASSAY, layer = "counts")
    )) == length(peaks_to_keep))
    
    # cat("Assay: ", DefaultAssay(SeuratOBJ_pb), "\n")
    # cat("Duplicates in assay rownames: ", anyDuplicated(rownames(SeuratOBJ_pb[[PSEUDO_ATAC_ASSAY]])) > 0, "\n")
    # cat("Overlaps found: ", length(idx), "\n")
    # cat("Example peaks_to_keep:\n"); print(utils::head(peaks_to_keep))
    # cat("Present in counts rows: ",
    #     sum(peaks_to_keep %in% rownames(GetAssayData(SeuratOBJ_pb, assay = PSEUDO_ATAC_ASSAY, layer = "counts"))),
    #     "/", length(peaks_to_keep), "\n")
    
    ## subset seurat and convert to sce
    sce_pb <- convert_atac_subset_to_sce(
            SeuratOBJ_pb, 
            PSEUDO_ATAC_ASSAY, 
            peaks_to_keep, 
            SeuratOBJ_pb@meta.data
            )
    assayNames(sce_pb)
    
    # added this chunk to patch: registration_stats_enrichment() which tries to pull row annotations via 
    # - rowData(sce_pseudo)[[gene_ensembl]] and/or [[gene_name]]. 
    # - If none of those arguments is provided, it resolves to NULL in rowData, and gets error
    # Give rowData (peak ranges) to explicit ID columns
    if (!"peak_id" %in% colnames(rowData(sce_pb))) {
        rowData(sce_pb)$peak_id <- rownames(sce_pb)
    }
    # Provide a second column for 'gene_ensembl' argument (duplicate peak_id)
    if (!"peak_ensembl" %in% colnames(rowData(sce_pb))) {
        rowData(sce_pb)$peak_ensembl <- rownames(sce_pb)
    }
    
    ## required columns
    sce_pb$registration_variable <- factor(sce_pb$cellType)          # group to test
    sce_pb$registration_sample_id <- factor(sce_pb$donor)            # block by donor
    sce_pb$ethnicity <- factor(sce_pb$ethnicity, levels = c("AA", "EA.CAUC"))
    sce_pb$sex <- factor(sce_pb$sex, levels = c("F", "M"))
    sce_pb$age <- as.numeric(sce_pb$age)
    
    # Inspect the registration model
    covars_vec <- c("age", "sex", "ethnicity")
    reg_mod <- registration_model(
        sce_pseudo = sce_pb,
        covars = covars_vec,
        var_registration = "registration_variable" #cellType
    )
    #head(reg_mod)  # inspect column names / coding
    
    ## estimate donor-level block correlation
    nrow(assay(sce_pb, "counts")) # [1] 5224
    
    block_cor <- registration_block_cor(
        sce_pseudo = sce_pb,
        registration_model = reg_mod,
        var_sample_id = "registration_sample_id"  ##cellType
    )
    # Endo ct test:
    # 2025-09-19 13:06:16.032861 run duplicateCorrelation()
    # 2025-09-19 13:06:22.540479 The estimated correlation is: 0.00931803725048961
    
    ## Run enrichment t-stats (1-vs-all for each cell type)
    rowData(sce_pb)
    head(rowData(sce_pb)) # have peak_id and peak_ensembl

    res_enrich <- registration_stats_enrichment_voomLmFit(
        sce_pseudo = sce_pb,
        block_cor = block_cor,
        covars = covars_vec,
        var_registration = "registration_variable",
        var_sample_id = "registration_sample_id",
        gene_ensembl = "peak_ensembl",     # must exist in rowData(sce_pb)
        gene_name = "peak_id"              # carry peak IDs into the output
    )
    
    head(res_enrich)
    #                       t_stat_Astrocyte t_stat_Endo t_stat_Excit.Thal
    # chr1-629811-630032          -2.233975   -1.980414          3.315414
    # chr1-630189-630389          -2.287579   -1.959382          3.321134
    # chr1-633694-634122          -2.257702   -2.094713          3.321954

    ## rename columns 
    res_enrich <- res_enrich |>
        rename(
            peak_macs2 = ensembl,
            peak_id    = gene
        )
    
    # save enrichment stats
    f_name <- here(output_Dir, paste0("voomlmFit_DAR_peaks_ALL_in_", clus, ".csv"))
    write.csv(res_enrich, f_name, row.names = FALSE)
    
    # Identify top up/down peaks based on FDR
    fdr_cols <- grep("^fdr_", colnames(res_enrich), value = TRUE)
    
    # cluster-specific significant peaks
    res_sig_ct <- res_enrich |>
        filter(.data[[paste0("fdr_", clus)]] < FDR_thr) |>
        mutate(
            logFC_ct = .data[[paste0("logFC_", clus)]],
            direction = case_when(
                logFC_ct >  0 ~ "Up",    # opening
                logFC_ct <  0 ~ "Down",  # closing
                TRUE ~ "NS"              # should not occur if you filtered
            )
        )
    
    dim(res_sig_ct)
    head(res_sig_ct[, c("peak_id", paste0("fdr_", clus), "logFC_ct", "direction")])
    
    if (nrow(res_sig_ct) > 0) { 
        f_name <- here(output_Dir, paste0("voomlmFit_DAR_peaks_", clus, ".csv"))
        write.csv(res_sig_ct, f_name, row.names = FALSE)
        message("Enrichment statistics saved [", clus, "]")    
    } else { 
        next 
    }

}

message("Enrichment statistics done!")


# library("slurmjobs")
# job_single(
#   "17_pseudobulk_DARs_MACS2_reduced_voomLmFit",
#   create_shell = TRUE,
#   partition = "katun",
#   memory = "30G",
#   cores = 2,
#   logdir = "logs",
#   command = "Rscript 17_pseudobulk_DARs_MACS2_reduced_voomLmFit.R",
#   create_logdir = FALSE
# )


## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
