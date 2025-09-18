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
library("spatialLIBD")
library("ggplot2")
library("ggrepel") 
library("future")
library("future.apply")
library("Seurat")
library("Signac")
library("dplyr")
library("here")


# Testing spearman at 5e4 on macs2 peaks 
resolution_level = "Mid"
p_met = "spearman"
w_size = "5e5"
sig_thresh <- 0.1    # FDR cutoff
lfc_thresh <- 0.25   # logFC cutoff / log2


if (length(resolution_level)) {
    message("Processing job for peak-method:\n",
            p_met,
            "\nWindow-size\n",
            w_size)
} else {
    message("Input arguments missed")
    stop()
}

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


## Check directories
if (!dir.exists(plotDir)) {
    dir.create(plotDir)
}
if (!dir.exists(output_Dir)) {
    dir.create(output_Dir)
}

# List all files matching the specific clustering resolution level
lst_peak_files <- list.files(
    path = inputRDS_Dir,
    pattern = "*.rds",
)
#lst_peak_files = list.files(path = input_cvsDir)
message("Link peak-genes files found:")
lst_peak_files
# [1] "mtx_merged_peaks_cell_level_Mid_resolution.rds"   
# [2] "Seurat_peaks_merged_cell_level_Mid_resolution.rds"


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

## (1) get matrix of raw peaks: rows = peaks, cols = pseudobulk samples
PSEUDO_ATAC_ASSAY <- "ATAC_macs2_merged_pseudo"
atac_counts <- GetAssayData(SeuratOBJ_pb, assay = PSEUDO_ATAC_ASSAY, layer ="counts")
length(rownames(atac_counts)) # [1] 351037
head(atac_counts)

## (2) get data.frame with sample metadata, nrow = ncol(counts)
#  must include at least: sample_id, cluster_id, group, donor
meta <- SeuratOBJ_pb@meta.data
colnames(meta)
# [1] "orig.ident"                        "nCount_RNA"                       
# [3] "nFeature_RNA"                      "nCount_ATAC_macs2_merged_pseudo"  
# [5] "nFeature_ATAC_macs2_merged_pseudo"
head(meta)

## add meta-data in the expected format
meta$sample_id <- row.names(meta)
meta$donor <- sapply(strsplit(row.names(meta), "_", fixed = TRUE), 
                "[[", 2)
tail(meta)

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


## Ensure we are contrasting by cell-type
if (!"cell_type" %in% colnames(SeuratOBJ_pb@meta.data)) {
    SeuratOBJ_pb$cell_type <- sub("^[^_]+_", "", SeuratOBJ_pb$orig.ident)  # keep part after first underscore
}
length(SeuratOBJ_pb$cell_type)
# [1] 169

Idents(SeuratOBJ_pb) <- "cell_type"
cluster_ids <- levels(SeuratOBJ_pb)
unique(Idents(SeuratOBJ_pb))

message("Pseudobulk groups (cell-types):")
cluster_ids


# ====/

stopifnot(all(colnames(atac_counts) == meta$sample_id))

# # Tweakables
# min_cpm        <- 1            # expression filter threshold (logical CPM rule)
# min_samples    <- 3            # require in >= this many samples overall or per-group
# use_block      <- TRUE         # TRUE = account for repeated measures (e.g., donor)
# # use_samp_wts   <- TRUE         # TRUE = estimate sample quality weights
# # trend_ebayes   <- TRUE         # TRUE = eBayes(trend=TRUE) often good for counts
# # robust_ebayes  <- TRUE         # TRUE = robust empirical Bayes
# fdr_cutoff     <- 0.10


## ============================================================================/


## convert Seurat object into sce
sce_pb <- as.SingleCellExperiment(SeuratOBJ_pb,
                                  layer = PSEUDO_ATAC_ASSAY)
sce_pb
dim(sce_pb)
colData(sce_pb)
## SingleCellExperiment() attempts to transfer data from the Seurat object's feature metadata to the SCE object's rowData. 
## - But Seurat meta.data are saved in Seurat@meta.data, thus feature-level metadata need to be added directly to Seurat
rowData(sce_pb)
# DataFrame with 351037 rows and 0 columns
table(sce_pb$cell_type)
# Astrocyte       Endo Excit.Thal Inhib.Thal      LHb.1    LHb.1.3  LHb.1.3.4 
#   10         10         10         10         10          7         10 
# LHb.2.7      LHb.4      LHb.7      MHb.1    MHb.1.2      MHb.2      MHb.3 
#   10         10          6         10         10         10         10 
# Microglia      Oligo        OPC       Thal 
#   10          10         10          6 
table(sce_pb$donor)
# S03-Hb-r S04-Hb-r S05-Hb-r S06-Hb-r S07-Hb-r S08-Hb-r S09-Hb-r S10-Hb-r 
# 17       15       16       18       18       17       17       18 
# S11-Hb-r S12-Hb-r 
# 17       16 

## required columns
sce_pb$registration_variable <- factor(sce_pb$cell_type)          # group to test
sce_pb$registration_sample_id <- factor(sce_pb$donor)             # block by donor
sce_pb$ethnicity <- factor(sce_pb$ethnicity, levels = c("AA", "EA.CAUC"))
sce_pb$sex <- factor(sce_pb$sex, levels = c("F", "M"))
sce_pb$age <- as.numeric(sce_pb$age)

## Provide logcounts for ATAC pseudobulk
dge <- DGEList(counts = assay(sce_pb, "counts"))
dge <- calcNormFactors(dge, method = "TMM")
sizeFactors(sce_pb) <- dge$samples$norm.factors * dge$samples$lib.size / exp(mean(log(dge$samples$lib.size)))

# Compute CPM normalized counts, then log-transform
lcpm <- scuttle::calculateCPM(sce_pb, size.factors = sizeFactors(sce_pb))
logcounts(sce_pb) <- log2(lcpm + 1)

## annotate rows
rowData(sce_pb)$peak_id <- rownames(sce_pb)

# Inspect the registration model
covars_vec <- c("age", "sex", "ethnicity")
reg_mod <- registration_model(
    sce_pseudo = sce_pb,
    covars = covars_vec,
    var_registration = "registration_variable"
)
head(reg_mod)  # inspect column names / coding





message(Sys.time(), " - Loop voomlmFit by cluster")

## run voomLmFit: Transform count data to log2-counts per million (logCPM), estimate voom precision weights and fit limma linear models while allowing for loss of residual degrees of freedom due to exact zeros

## ============================================================================/
## version adapter from https://github.com/LieberInstitute/DeconvoBuddies/blob/d128d498c18318d05528bf75c6fa8436f1bab8c6/R/findMarkers_1vAll.R
## And from my previous implementation: https://github.com/LieberInstitute/Habenula_Visium/blob/2d21e39f51c9e46ebddbcf57f959f70c27d78678/code/05_brain_area_differential_expression/05_pseudobulk_DEG_contrast.R#L176-L255 

# Input
# logcounts(sce_pseudo) peaks (it fits limma on log-scale data)
# var_registration: in my case cellType
# var_sample_id: blocking by donor factor
# covars:age, sex and ethnicity













## ============================================================================/
## version from LFF_sspatial_ERC project (Louise)
## https://github.com/LieberInstitute/LFF_spatial_ERC/blob/b0b174759f0f539d7d6e7362df9205d7daebbf22/code/12_voomLmFit/01_Clusterwise_voomLmFit.R
# 
# lmf_summary <- map_dfr(clusters, function(clus){
#     # clus = "Endo"
#     
#     dge <- sce_pb[,sce_pb$registration_variable == clus]
#     # class: SingleCellExperiment 
#     # dim: 351037 10 
#     # metadata(0):
#     #     assays(2): counts logcounts
#     # rownames(351037): chr1-181329-181534 chr1-191217-191619 ...
#     # KI270728.1-1791302-1791701 KI270728.1-1792066-1792310
#     # rowData names(0):
#     #     colnames(10): Endo_S03-Hb-r Endo_S04-Hb-r ... Endo_S11-Hb-r
#     # Endo_S12-Hb-r
#     # colData names(13): orig.ident nCount_RNA ... ident
#     # registration_variable
#     # reducedDimNames(1): LSI
#     # mainExpName: ATAC_macs2_merged_pseudo
#     # altExpNames(1): RNA
#     
#     ## set up the variables
#     colData(dge)$sex <- as.factor(colData(dge)$sex)
#     levels(dge$sex) # [1] "F" "M"
#     colData(dge)$ethnicity <- as.factor(colData(dge)$ethnicity)
#     levels(dge$ethnicity) # [1] "AA"      "EA.CAUC"
#     colData(dge)$age <- as.numeric(colData(dge)$age)
#     
#     #des <- model.matrix(~0 + APOE_syn + Sex + Age + Anc_Afr + pseudo_expr_chrM_ratio, data = colData(dge))
#     des <- model.matrix(~0 + sex + ethnicity + age , data = colData(dge))
#     # sexF, sexM, ethnicityAA, ethnicityEA.CAUC, age
#     colnames(des) <- make.names(colnames(des))
#     # ~0 + → no intercept
#     # sex (2 levels: F, M) → 2 columns: sexF, sexM
#     # ethnicity (2 levels: "AA", "EA/CAUC") → 2 columns: ethnicityAA, ethnicityEA.CAUC
#     # age (numeric) → 1 column: age
#     
#     des_df <- as.data.frame(des)
#     des_df
#     # sexF sexM ethnicityEA.CAUC   age
#     # Endo_S03-Hb-r    1    0                1 41.30
#     # Endo_S04-Hb-r    0    1                0 40.57
#     # Endo_S05-Hb-r    0    1                0 46.72
#     # Endo_S06-Hb-r    0    1                1 33.39
#     # Endo_S07-Hb-r    0    1                1 48.88
#     # Endo_S08-Hb-r    1    0                1 37.33
#     # Endo_S09-Hb-r    0    1                1 39.98
#     # Endo_S10-Hb-r    0    1                1 47.95
#     # Endo_S11-Hb-r    1    0                1 65.00
#     # Endo_S12-Hb-r    0    1                1 57.70
#     
#     # filter low expression genes
#     dge <- edgeR::calcNormFactors(dge)
#     keep <- edgeR::filterByExpr.DGEList(dge,design=des)
#     dge <- dge[keep,,keep.lib.sizes=FALSE]
#     dge <- edgeR::calcNormFactors(dge)
#     head(dge, n=2)
#     
#     # >     head(dge, n=2)
#     # An object of class "DGEList"
#     # $counts
#     # Endo_S03-Hb-r Endo_S04-Hb-r Endo_S05-Hb-r Endo_S06-Hb-r
#     # chr1-629811-630032             9             9            27            28
#     # chr1-633694-634122            31            68           135           157
#     # Endo_S07-Hb-r Endo_S08-Hb-r Endo_S09-Hb-r Endo_S10-Hb-r
#     # chr1-629811-630032            17            29            22            23
#     # chr1-633694-634122           105           129           105           134
#     # Endo_S11-Hb-r Endo_S12-Hb-r
#     # chr1-629811-630032             7             4
#     # chr1-633694-634122            15             8
#     # 
#     # $samples
#     # group lib.size norm.factors orig.ident nCount_RNA nFeature_RNA
#     # Endo_S03-Hb-r     1      910    0.9578798       Endo     104052        10740
#     # Endo_S04-Hb-r     1     2478    1.3085681       Endo     208249        14020
#     # Endo_S05-Hb-r     1     5381    0.3267022       Endo     102209        12480
#     # Endo_S06-Hb-r     1     4926    1.2079846       Endo     275224        14319
#     # Endo_S07-Hb-r     1     1619    0.9246938       Endo     180751        13337
#     # Endo_S08-Hb-r     1     1239    0.9702263       Endo      97452        11324
#     # Endo_S09-Hb-r     1     4095    1.3104865       Endo     222285        13848
#     # Endo_S10-Hb-r     1     2434    1.1383460       Endo     135145        12464
#     # Endo_S11-Hb-r     1      365    1.3083337       Endo      48230         9096
#     # Endo_S12-Hb-r     1      404    1.1544705       Endo      20285         6345
#     # nCount_ATAC_macs2_merged_pseudo nFeature_ATAC_macs2_merged_pseudo
#     # Endo_S03-Hb-r                           71264                             56601
#     # Endo_S04-Hb-r                          201904                            121829
#     # Endo_S05-Hb-r                          186304                            123898
#     # Endo_S06-Hb-r                          347620                            171788
#     # Endo_S07-Hb-r                           83383                             62606
#     # Endo_S08-Hb-r                           63774                             47724
#     # Endo_S09-Hb-r                          227054                            123037
#     # Endo_S10-Hb-r                          121119                             75853
#     # Endo_S11-Hb-r                           27142                             24644
#     # Endo_S12-Hb-r                           10583                              9726
#     # sample_id    donor sex ethnicity   age cell_type ident
#     # Endo_S03-Hb-r Endo_S03-Hb-r S03-Hb-r   F   EA/CAUC 41.30      Endo  Endo
#     # Endo_S04-Hb-r Endo_S04-Hb-r S04-Hb-r   M        AA 40.57      Endo  Endo
#     # Endo_S05-Hb-r Endo_S05-Hb-r S05-Hb-r   M        AA 46.72      Endo  Endo
#     # Endo_S06-Hb-r Endo_S06-Hb-r S06-Hb-r   M   EA/CAUC 33.39      Endo  Endo
#     # Endo_S07-Hb-r Endo_S07-Hb-r S07-Hb-r   M   EA/CAUC 48.88      Endo  Endo
#     # Endo_S08-Hb-r Endo_S08-Hb-r S08-Hb-r   F   EA/CAUC 37.33      Endo  Endo
#     # Endo_S09-Hb-r Endo_S09-Hb-r S09-Hb-r   M   EA/CAUC 39.98      Endo  Endo
#     # Endo_S10-Hb-r Endo_S10-Hb-r S10-Hb-r   M   EA/CAUC 47.95      Endo  Endo
#     # Endo_S11-Hb-r Endo_S11-Hb-r S11-Hb-r   F   EA/CAUC 65.00      Endo  Endo
#     # Endo_S12-Hb-r Endo_S12-Hb-r S12-Hb-r   M   EA/CAUC 57.70      Endo  Endo
#     # registration_variable
#     # Endo_S03-Hb-r                  Endo
#     # Endo_S04-Hb-r                  Endo
#     # Endo_S05-Hb-r                  Endo
#     # Endo_S06-Hb-r                  Endo
#     # Endo_S07-Hb-r                  Endo
#     # Endo_S08-Hb-r                  Endo
#     # Endo_S09-Hb-r                  Endo
#     # Endo_S10-Hb-r                  Endo
#     # Endo_S11-Hb-r                  Endo
#     # Endo_S12-Hb-r                  Endo
#     
#     #message(Sys.time(), sprintf(" - voomLmFit - cluster: %s, block= '%s', ncol: %s, ngene: %i", clus, batch, ncol(dge), nrow(dge$genes)))
#     message(Sys.time(), sprintf(" - voomLmFit - cluster: %s, block= '%s', ncol: %s, nATAC: %i", clus, batch, ncol(dge), nrow(dge$counts)))
#     
#     
#     ## run voomLmFit for the pseudobulked data, referring donor to duplicateCorrelation; 
#     ## using an adaptive span (number of genes, based on the number of genes in the dge) for smoothing the mean-variance trend
#     #v.swt <- voomLmFit(dge,design = des,block = as.factor(dge$samples[[batch]]),adaptive.span = T,sample.weights = T)
#     
#     v.swt <- voomLmFit(dge$counts,                               # numeric matrix containing raw counts
#                        design = des,                             # rows corresponding to samples and columns to coefficients to be estimated
#                        block = as.factor(dge$samples[[batch]]),  # blocking variable on the samples / donor
#                        #adaptive.span = T,                        # width of smoothing window for the lowess mean-variance trend. Proportion between 0 and 1
#                        sample.weights = T)                       # if TRUE then empirical sample quality weights will be estimated
#     
#     # First sample weights (min/max) 0.2176917/1.8742648
#     # First intra-block correlation  0
#     # Final sample weights (min/max) 0.2533203/1.9176781
#     # Final intra-block correlation  0
#     
#     # This function returns df.residual values that are less than or equal to those from lmFit 
#     # - and sigma values that are greater than or equal to those from lmFit
#     v.swt
#     
#     # make these more readable
#     cont <- makeContrasts(
#         MaleVsFemale = sexM - sexF,
#         #AfrVsCauc    = ethnicityAA - ethnicityEA.CAUC,
#         levels = des
#     )
#     
#     # cont <- makeContrasts(
#     #     # main
#     #     carrier = "-0.5*(APOE_E2.E2 + APOE_E2.E3) + 0.5*(APOE_E3.E4 + APOE_E4.E4)",
#     #     E4E4 = "-APOE_E4.E4 + (APOE_E2.E2 + APOE_E2.E3 + APOE_E3.E4)/3",
#     #     ## apoe pairwise
#     #     apoe_E2E2_E4E4 = "-APOE_E2.E2 + APOE_E4.E4",
#     #     apoe_E3E4_E4E4 = "-APOE_E3.E4 + APOE_E4.E4",
#     #     apoe_E2E3_E4E4 = "-APOE_E2.E3 + APOE_E4.E4",
#     #     apoe_E2E2_E3E4 = "-APOE_E2.E2 + APOE_E3.E4",
#     #     apoe_E2E2_E2E3 = "-APOE_E2.E2 + APOE_E2.E3",
#     #     apoe_E2E3_E3E4 = "-APOE_E2.E3 + APOE_E3.E4",
#     #     # heterozygous vs. homozygous
#     #     anyE2_E4E4 = "- 0.5*(APOE_E2.E3 + APOE_E2.E2) + APOE_E4.E4",
#     #     E2E2_anyE4 = "-APOE_E2.E2 + 0.5*(APOE_E3.E4 + APOE_E4.E4)",
#     #     E2E3_anyE4 = "-APOE_E2.E3 + 0.5*(APOE_E3.E4 + APOE_E4.E4)",
#     #     # other
#     #     Sex="SexM",
#     #     Anc="Anc_Afr",
#     #     levels=des
#     # )
# 
#     v.swt.fit <- contrasts.fit(v.swt,contrasts=cont)
#     v.swt.fit.e <- eBayes(v.swt.fit)
#     
#     ## run top table over contrasts
#     v.swt.e.tt <- purrr::map(colnames(cont), ~topTable(v.swt.fit.e,coef = .x, number=Inf, adjust.method = "BH") |>
#                                  mutate(data_type = opt$datatype, 
#                                         cluster = clus,
#                                         contrast = .x, 
#                                         .before = 1) |>
#                                  arrange(adj.P.Val)) 
#     
#     names(v.swt.e.tt) <- colnames(cont)
#     
#     message("Done - Save data")
#     saveRDS(v.swt.e.tt, file = here(data_dir, sprintf("voomLmFit_%s_%s.rds", opt$datatype, clus)))
#     return(purrr::map_int(v.swt.e.tt, ~sum(.x$adj.P.Val < 0.05)))
# })
#
# lmf_summary <- lmf_summary |>
#     add_column(cluster = clusters, .before=1)
#
#write.csv(lmf_summary, file = here(data_dir, sprintf("vlmf_FDR05_summary-%s.csv", opt$datatype)), row.names = FALSE)
#
## ============================================================================/



# library("slurmjobs")
# job_single(
#   "16_pseudobulk_DARs_MACS2_reduced",
#   create_shell = TRUE,
#   partition = "katun",
#   memory = "30G",
#   cores = 2,
#   logdir = "logs",
#   command = "Rscript 16_pseudobulk_DARs_MACS2_reduced.R",
#   create_logdir = FALSE
# )


## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
