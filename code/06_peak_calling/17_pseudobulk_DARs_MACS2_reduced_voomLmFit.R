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
p_met = "spearman"
w_size = "5e5"
sig_thresh <- 0.1    # FDR cutoff
lfc_thresh <- 0.25   # logFC cutoff / log2

## ATAC function's helper used globally
source(here("code", "06_peak_calling", "multiome_custom_functions", "multiome_idents_normalization_helper.R"))

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
PSEUDO_ATAC_ASSAY <- "ATAC_macs2_merged_pseudo"

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



##==============================================================================

## extract peaks cell-type specific
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

## get matrix of peaks: rows = peaks, cols = pseudobulk samples

atac_counts <- GetAssayData(SeuratOBJ_pb, assay = PSEUDO_ATAC_ASSAY, layer ="counts")
length(rownames(atac_counts)) # [1] 351037
head(atac_counts)

## get data.frame with sample metadata, nrow = ncol(counts)
#  must include at least: sample_id, cluster_id, group, donor
meta <- SeuratOBJ_pb@meta.data
colnames(meta)
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


# ====/

message("Subset peaks for specific cellType ... ")

stopifnot(all(colnames(atac_counts) == meta$sample_id))

DefaultAssay(SeuratOBJ_pb) <- PSEUDO_ATAC_ASSAY

# Get all peaks in the assay
all_peaks <- granges(SeuratOBJ_pb)
length(all_peaks) # [1] 351037

# Get peaks cellType specific
peak_ranges_ct <- peaks_ranges_cellType(
    SeuratOBJ_pb,
    PSEUDO_ATAC_ASSAY,
    cluster_name = "Endo"
) 
length(peak_ranges_ct) # 5224

# Get the assay’s ranges and rownames
ap_gr   <- granges(SeuratOBJ_pb[[PSEUDO_ATAC_ASSAY]])
# GRanges object with 351037 ranges and 2 metadata columns:
# seqnames          ranges strand |        revmap
# <Rle>       <IRanges>  <Rle> | <IntegerList>
# [1]       chr1   181329-181534      * |             1
# [2]       chr1   191217-191619      * |             2
# [3]       chr1   629146-629354      * |             3

ap_rows <- rownames(SeuratOBJ_pb[[PSEUDO_ATAC_ASSAY]])
# [1] "chr1-181329-181534"   "chr1-191217-191619"   "chr1-629146-629354"  
# [4] "chr1-629811-630032"   "chr1-630189-630389"   "chr1-632189-632410"  
# [7] "chr1-633694-634122"   "chr1-777714-777962"   "chr1-778314-779308"  

# Find exact matches by genomic coordinates
#    (use type="any" for overlaps instead of exact start/end equality)
hits <- findOverlaps(ap_gr, peak_ranges_ct, type = "equal")
idx  <- queryHits(hits)
# [1]      4      5      7     53     60     73    116    168    189    226
# [11]    288    290    293    324    340    353    361    503    535    557
# [21]    919    989    994   1018   1222   1234   1258   1259   1267   1428

# Map hits to peaks IDs actually present in the assay
peaks_to_keep <- ap_rows[idx]
peaks_to_keep <- unique(as.character(peaks_to_keep))
# [1495] "chr5-72307684-72308762"   "chr5-72319550-72320868"  
# [1497] "chr5-72815810-72817418"   "chr5-73119548-73120878"  
# [1499] "chr5-73564688-73566247"   "chr5-75236269-75237381"  

# checks before subsetting
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


# subset by peaks
Seurat_subset <- subset(
    SeuratOBJ_pb[[PSEUDO_ATAC_ASSAY]],
    features = peaks_to_keep
)
Seurat_subset
# ChromatinAssay data with 5224 features for 169 cells
# Variable features: 5216 
# ...

## convert object into sce with meta.data
class(Seurat_subset)
counts_mat <- GetAssayData(Seurat_subset, layer = "counts")
logcounts_mat <- GetAssayData(Seurat_subset, layer = "data")  # “data” in Seurat usually holds lognorm values

Seurat_subset2 <- CreateSeuratObject(
    counts = counts_mat,
    assay = PSEUDO_ATAC_ASSAY,
    meta.data = SeuratOBJ_pb@meta.data
)
Seurat_subset2
# > Seurat_subset2
# An object of class Seurat 
# 5224 features across 169 samples within 1 assay 
# Active assay: ATAC_macs2_merged_pseudo (5224 features, 0 variable features)
# 1 layer present: counts

# Add logcounts back
Seurat_subset2[[PSEUDO_ATAC_ASSAY]] <- SetAssayData(
    Seurat_subset2[[PSEUDO_ATAC_ASSAY]],
    layer = "logcounts",
    new.data = logcounts_mat
)
Seurat_subset2
# An object of class Seurat 
# 5224 features across 169 samples within 1 assay 
# Active assay: ATAC_macs2_merged_pseudo (5224 features, 0 variable features)
# 2 layers present: counts, logcounts

sce_pb <- as.SingleCellExperiment(
    Seurat_subset2,
    assay = PSEUDO_ATAC_ASSAY
)

sce_pb
dim(sce_pb)
colData(sce_pb)
## SingleCellExperiment() attempts to transfer data from the Seurat object's feature metadata to the SCE object's rowData. 
## - But Seurat meta.data are saved in Seurat@meta.data, thus feature-level metadata need to be added directly to Seurat
rowData(sce_pb)
# DataFrame with 351037 rows and 0 columns
table(sce_pb$cellType)
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


##=======

## refs:
## current version https://github.com/LieberInstitute/spatialLIBD/blob/40da043d0235e01a12a7f52a0b367d3850bad9e8/R/registration_stats_enrichment.R#L40
## My previous manual implementation: https://github.com/LieberInstitute/Habenula_Visium/blob/2d21e39f51c9e46ebddbcf57f959f70c27d78678/code/05_brain_area_differential_expression/05_pseudobulk_DEG_contrast.R#L176-L255 

# Inputs for atac
# logcounts(sce_pseudo) peaks (it fits limma on log-scale data)
# var_registration: in my case cellType
# var_sample_id: blocking by donor factor
# covars:age, sex and ethnicity

# added this chunck to patch: registration_stats_enrichment() tries to pull row annotations via rowData(sce_pseudo)[[gene_ensembl]] and/or [[gene_name]]. 
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

## extract/transform peak logcounts 

if (!"logcounts" %in% assayNames(sce_pb)) {
    # create logcounts from raw counts
    dge <- DGEList(counts = assay(sce_pb, "counts"))
    dge <- calcNormFactors(dge, method = "TMM")
    sizeFactors(sce_pb) <- dge$samples$norm.factors *
        dge$samples$lib.size /
        exp(mean(log(dge$samples$lib.size)))
    
    # CPM normalize then log-transform
    lcpm <- edgeR::cpm(dge, log = FALSE, prior.count = 0, normalized.lib.sizes = TRUE)
    logcounts(sce_pb) <- log2(lcpm + 1)
    
} else {

    message("logcounts already present, skipping normalization")
    
}

    
# Inspect the registration model
covars_vec <- c("age", "sex", "ethnicity")
reg_mod <- registration_model(
    sce_pseudo = sce_pb,
    covars = covars_vec,
    var_registration = "registration_variable" #cellType
)
head(reg_mod)  # inspect column names / coding


## estimate donor-level block correlation
nrow(assay(sce_pb, "counts"))
# 351037

block_cor <- registration_block_cor(
    sce_pseudo = sce_pb,
    registration_model = reg_mod,
    var_sample_id = "registration_sample_id"  ##cellType
)
# Endo ct test:
# 2025-09-19 13:06:16.032861 run duplicateCorrelation()
# 2025-09-19 13:06:22.540479 The estimated correlation is: 0.00931803725048961

## Run enrichment t-stats (1-vs-all for each cell type)

# Add column with peak IDs
rowData(sce_pb)
#rowData(sce_pb)$peak_id <- rownames(sce_pb)
head(rowData(sce_pb))
# DataFrame with 6 rows and 1 column
# peak_id
# <character>
#     chr1-181329-181534 chr1-181329-181534
# chr1-191217-191619 chr1-191217-191619
# chr1-629146-629354 chr1-629146-629354
# chr1-629811-630032 chr1-629811-630032
# chr1-630189-630389 chr1-630189-630389
# chr1-632189-632410 chr1-632189-632410

res_enrich <- registration_stats_enrichment(
    sce_pseudo = sce_pb,
    block_cor = block_cor,
    covars = covars_vec,
    var_registration = "registration_variable",
    var_sample_id = "registration_sample_id",
    #gene_ensembl = NULL,              # not genes here
    gene_ensembl = "peak_ensembl",     # must exist in rowData(sce_pb)
    gene_name = "peak_id"             # carry peak IDs into the output
)

head(res_enrich)
# save summary
f_name <- here(output_Dir, "stats_enrichment_lmFit_atac_peaks_by_cluster.csv")
write.csv(res_enrich, f_name, row.names = FALSE)

message("Enrichment statistics done!")


# message(Sys.time(), " - Loop voomlmFit by cluster")





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
#     # sample_id    donor sex ethnicity   age cellType ident
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
