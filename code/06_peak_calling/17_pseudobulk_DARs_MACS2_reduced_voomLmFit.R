########################################################################
## Differential Accessibility (DA) using voomLmFit / pseudobulk multiome assays with merged peaks
## ## refs:
## adapted from https://github.com/LieberInstitute/spatialLIBD/blob/40da043d0235e01a12a7f52a0b367d3850bad9e8/R/registration_stats_enrichment.R#L40
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
library("Signac")
library("spatialLIBD")
library("GenomicRanges")
library("ggplot2")
library("ggrepel") 
library("dplyr")
library("stringr")
library("here")


# Testing spearman at 5e4 on macs2 peaks 
resolution_level = "Mid"
# lfc_thresh <- 0.25 # 2^0.25 ≈1.189
FDR_thr = 0.20

## read input arguments
args = commandArgs(trailingOnly = TRUE)
clus <- args[1]
if (is.na(clus) || !nzchar(clus)) stop("Missing cluster_name argument")

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

registration_stats_enrichment_voomLmFit <- function(
        sce_pseudo,
        covars = NULL,                              # c("age", "sex", "ethnicity")
        var_registration = "registration_variable", # cellType = 18 ct
        var_sample_id = "registration_sample_id",   # donor = 169 samples
        gene_ensembl = NULL,
        gene_name = NULL
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
            # res <- limma::eBayes(edgeR::voomLmFit(
            #     assay(sce_pseudo, "counts"),
            #     design = m,
            #     block = sce_pseudo[[var_sample_id]], # ge. donor
            #     # correlation = block_cor, #recent version of the edgeR package does not accept the correlation argument 
            #     #  - it computes the inter-block correlation internally 
            #     sample.weights = TRUE
            # ))
            # add trend = TRUE to manage small size-samples (<10)
            vfit <- edgeR::voomLmFit(
                assay(sce_pseudo, "counts"),
                design = m,
                block = sce_pseudo[[var_sample_id]],
                sample.weights = TRUE
            )
            res <- limma::eBayes(vfit, trend = TRUE)
            
            #plotSA(res_enrich, main = "Mean–variance trend (voom+trend)")
            message("Residual DF summary:")
            print(summary(res$df.residual))
            
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


## =============================================================================
## add meta-data and setup ATAC assay

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
table(SeuratOBJ_pb$cellType)
# Astrocyte       Endo Excit.Thal Inhib.Thal      LHb.1    LHb.1.3  LHb.1.3.4 
# 10         10         10         10         10          7         10 
# LHb.2.7      LHb.4      LHb.7      MHb.1    MHb.1.2      MHb.2      MHb.3 
# 10         10          6         10         10         10         10 
# Microglia      Oligo        OPC       Thal 
# 10         10         10          6 

Idents(SeuratOBJ_pb) <- "cellType"
cluster_ids <- levels(SeuratOBJ_pb)
unique(Idents(SeuratOBJ_pb))

message("Pseudobulk groups (cell-types):")
cluster_ids
# [1] "Astrocyte"  "Endo"       "Excit.Thal" "Inhib.Thal" "LHb.1"     
# [6] "LHb.1.3"    "LHb.1.3.4"  "LHb.2.7"    "LHb.4"      "LHb.7"     
# [11] "MHb.1"      "MHb.1.2"    "MHb.2"      "MHb.3"      "Microglia" 
# [16] "Oligo"      "OPC"        "Thal"  

stopifnot(all(colnames(atac_counts) == meta$sample_id))



message("Starting registration_stats_enrichment_voomLmFit ... ")

# Inputs for atac:
#        var_registration: in my case cellType
#        var_sample_id: blocking by donor factor
#        covars:age, sex and ethnicity

DefaultAssay(SeuratOBJ_pb) <- PSEUDO_ATAC_ASSAY
# Get all peaks in the assay
all_peaks <- granges(SeuratOBJ_pb)
length(all_peaks) # [1] 351037


##==============================================================================

message("Starting stats enrichment for cluster", clus)
    
# Get peaks cellType specific

peak_ranges_ct <- global_peaks_ranges_cellType(
    SeuratOBJ_pb,
    PSEUDO_ATAC_ASSAY,
    cluster_name = clus
) 
message("Subset ", length(peak_ranges_ct), " peaks for [", clus, "] cellType") # ge. ENDO: 5224

# Get the assay’s ranges and rownames
ap_gr   <- granges(SeuratOBJ_pb[[PSEUDO_ATAC_ASSAY]])
ap_rows <- rownames(SeuratOBJ_pb[[PSEUDO_ATAC_ASSAY]])
# Find exact matches by genomic coordinates (use type="any" for overlaps)
hits <- findOverlaps(ap_gr, peak_ranges_ct, type = "equal")
idx  <- queryHits(hits)
# Map hits to peaks IDs actually present in the assay
peaks_to_keep <- ap_rows[idx]
peaks_to_keep <- unique(as.character(peaks_to_keep))
message( length(peak_ranges_ct), " / ", length(peaks_to_keep))

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

## convert seurat to sce
sce_pb <- global_convert_atac_subset_to_sce(
        SeuratOBJ_pb, 
        PSEUDO_ATAC_ASSAY, 
        peaks_to_keep, 
        SeuratOBJ_pb@meta.data
        )
assayNames(sce_pb)

# added chunk to patch registration_stats_enrichment() which tries to pull row annotations via 
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
sce_pb$registration_variable <- factor(sce_pb$cellType)          # group to test / 18 ct
sce_pb$registration_sample_id <- factor(sce_pb$donor)            # block by donor / 10 donors
sce_pb$ethnicity <- factor(sce_pb$ethnicity, levels = c("AA", "EA.CAUC"))
sce_pb$sex <- factor(sce_pb$sex, levels = c("F", "M"))
sce_pb$age <- as.numeric(sce_pb$age)

# Inspect the registration model
covars_vec <- c("age", "sex", "ethnicity")
reg_mod <- registration_model(
    sce_pseudo = sce_pb,
    covars = covars_vec,
    var_registration = "registration_variable"
)
#head(reg_mod)  # inspect column names / coding

## Run enrichment t-stats (1-vs-all for each cell type)

head(rowData(sce_pb)) # have peak_id and peak_ensembl

res_enrich <- registration_stats_enrichment_voomLmFit(
    sce_pseudo = sce_pb,
    covars = covars_vec,
    var_registration = "registration_variable",
    var_sample_id = "registration_sample_id",
    gene_ensembl = "peak_ensembl",     # must exist in rowData(sce_pb)
    gene_name = "peak_id"              # carry peak IDs into the output
)
# head(res_enrich)

## rename columns 
res_enrich <- res_enrich |>
    rename(
        peak_macs2 = ensembl,
        peak_id    = gene
    )

# save enrichment stats
f_name <- here(output_Dir, paste0("voomlmFit_DAR_peaks_ALL_in_", clus, ".csv"))
write.csv(res_enrich, f_name, row.names = FALSE)

# # Identify top up/down peaks based on FDR
# fdr_cols <- grep("^fdr_", colnames(res_enrich), value = TRUE)
# 
# # cluster-specific significant peaks
# res_sig_ct <- res_enrich |>
#     filter(.data[[paste0("fdr_", clus)]] < FDR_thr) |>
#            #abs(.data[[paste0("logFC_", clus)]]) > lfc_thresh) |>
#     mutate(
#         logFC_ct = .data[[paste0("logFC_", clus)]],
#         direction = case_when(
#             logFC_ct >  0 ~ "Up",    # opening
#             logFC_ct <  0 ~ "Down",  # closing
#             TRUE ~ "NS"              # should not occur if logFC = 0
#         )
#     )
# 
# if (nrow(res_sig_ct) > 0) { 
#     f_name <- here(output_Dir, paste0("voomlmFit_DAR_peaks_", clus, ".csv"))
#     write.csv(res_sig_ct, f_name, row.names = FALSE)
#     message("Enrichment statistics saved [", clus, "]")    
# }


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
