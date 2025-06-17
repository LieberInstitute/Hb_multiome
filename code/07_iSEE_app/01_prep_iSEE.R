library("SingleCellExperiment")
library("Seurat")
library("here")
library("lobstr")
library("Matrix")
library("tidyverse")
library("sessioninfo")

## Load RDS multiome
message(Sys.time(), "- load Harmony corrected Seurat")

Seurat_base_name <- here("processed-data", "05_Clustering_ARCr", "05_rename_idents",
                         "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium.rds")

Seurat_base_name
SeuratOBJ <- readRDS(Seurat_base_name)

## verification
SeuratOBJ
levels(SeuratOBJ)
# [1] "C.01.undeterminated" "C.02.DD_Oligo"       "C.03.undeterminated"
# [4] "C.04.undeterminated" "C.05.DD_LHb"         "C.06.DD_Exit.Thal"  
# [7] "C.07.DD_MHb"         "C.08.undeterminated" "C.09.undeterminated"
# [10] "C.10.DD_MHb"         "C.11.DD_MHb"         "C.12.undeterminated"
# [13] "C.13.no-match"       "C.14.DD_MHb"         "C.15.DD_Exit.Thal"  
# [16] "C.16.DD_MHb"         "C.17.DD_Exit.Thal"   "C.18.DD_LHb"        
# [19] "C.19.DD_Inhib.Thal"  "C.20.DD_Astrocyte"   "C.21.DD_Astrocyte"  
# [22] "C.22.undeterminated" "C.23.DD_LHb"         "C.24.DD_LHb"        
# [25] "C.25.undeterminated" "C.26.DD_OPC"         "C.27.DD_Microglia"  
# [28] "C.28.DD_Inhib.Thal"  "C.29.DD_Endo"        "C.30.DD_LHb"        
# [31] "C.31.DD_Exit.Thal"   "C.32.undeterminated" "C.33.DD_LHb"        
# [34] "C.34.DD_Oligo"       "C.35.undeterminated" "C.36.DD_MHb"        
# [37] "C.37.undeterminated" "C.38.DD_Inhib.Thal"  "C.39.DD_Inhib.Thal" 
# [40] "C.40.DD_LHb"         "C.41.DD_Microglia"   "C.42.no-match"    

## Add active cluster identity (levels(SeuratOBJ)) as a metadata to easy cell-type identification
SeuratOBJ$cluster_annotation <- Idents(SeuratOBJ)
head(SeuratOBJ@meta.data[, "cluster_annotation"])

## verify that RNA assay has gene symbols
head(rownames(SeuratOBJ[["RNA"]]))
# [1] "MIR1302-2HG" "FAM138A"     "OR4F5"       "AL627309.1"  "AL627309.3" 
# [6] "AL627309.2" 

head(SeuratOBJ@meta.data)

## make slim Seurat with ony RNA modality
SeuratOBJ <- DietSeurat(SeuratOBJ, 
                        assays = "RNA")
SeuratOBJ
# An object of class Seurat 
# 36601 features across 55702 samples within 1 assay 
# Active assay: RNA (36601 features, 2000 variable features)
# 3 layers present: data, counts, scale.data
head(rownames(SeuratOBJ[["RNA"]]))
# [1] "MIR1302-2HG" "FAM138A"     "OR4F5"       "AL627309.1"  "AL627309.3" 
# [6] "AL627309.2" 

# Remove all reductions from SeuratOBJ to avoid mismatch issues with sce convertion
for (red in Reductions(SeuratOBJ)) {
    SeuratOBJ[[red]] <- NULL
}
Reductions(SeuratOBJ)
# NULL

## convert Seurat object into sce
sce <- as.SingleCellExperiment(SeuratOBJ)
## verification
sce
# class: SingleCellExperiment 
# dim: 36601 55702 
# metadata(0):
#     assays(3): counts logcounts scaledata
# rownames(36601): MIR1302-2HG FAM138A ... AC007325.4 AC007325.2
# rowData names(0):
#     colnames(55702): S04_AAACAGCCAGAATGAC-1 S04_AAACAGCCAGCAAGGC-1 ...
# S09_TTTGTTGGTCATGCAA-1 S09_TTTGTTGGTTGTTCAC-1
# colData names(29): orig.ident nCount_RNA ... C.leiden_wnn ident
# reducedDimNames(0):
#     mainExpName: RNA
# rownames(sce) 

head(rownames(sce))
# [1] "MIR1302-2HG" "FAM138A"     "OR4F5"       "AL627309.1"  "AL627309.3" 
# [6] "AL627309.2" 

## Drop data we don't need for iSEE. We keep only logcount
assayNames(sce)
# [1] "counts"    "logcounts" "scaledata"
assays(sce)$counts <- NULL
assays(sce)$scaledata <- NULL
assayNames(sce)
# [1] "logcounts"

class(logcounts(sce))
# [1] "DelayedMatrix"
# ======
# [1] "dgCMatrix"
# attr(,"package")
# [1] "Matrix"

message(Sys.time(), "- Convert logcounts to sparse Matrix")
logcounts(sce) <- as(logcounts(sce), "sparseMatrix")  

## Drop metadata we don't need
metadata(sce) <- list()
# sce$path <- NULL
# sce$total <- NULL

## sourcing official color palette
# load(here("processed-data", "04_snRNA-seq", "cell_type_colors.Rdata"))
# sn_colors<- cell_type_colors

## Check final size
lobstr::obj_size(sce)
# 2.36 GB

#rownames(sce) <- rowData(sce)$Symbol
head(rownames(sce))
# [1] "MIR1302-2HG" "FAM138A"     "OR4F5"       "AL627309.1"  "AL627309.3" 
# [6] "AL627309.2" 

# # ===== I do not have mean-ratio genes, trying to adpapt to work with the DGE from Seurat
# #### Add MeanRatio Marker Gene Details ####
# load(here("processed-data", "04_snRNA-seq", "16_sn_MeanRatio", "MarkerStats_cell_type_fine.Rdata"))
# marker_stats |> dplyr::count(cellType.target)
# 
# marker_anno <- marker_stats |>
#     filter(MeanRatio.rank <= 50 & MeanRatio > 1) |>
#     select(gene,
#            cellType.target,
#            MeanRatio.rank,
#            MeanRatio,
#            MeanRatio.anno) |>
#     column_to_rownames("gene")
# 
# rowData(sce) <- cbind(rowData(sce), marker_anno[rownames(sce),])
# 
# rowData(sce)[which(rowData(sce)$MeanRatio.rank ==1),]
# # ============================================================================

# saveRDS(sce, file = here("code", "06_iSEE_app", "sce_ERC_iSEE.rds"))
saveRDS(sce, here("code", "07_iSEE_app", "sce_Habenula_iSEE.rds"))

# slurmjobs::job_single('01_prep_iSEE', create_shell = TRUE, memory = '25G', command = "Rscript 01_prep_iSEE.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
