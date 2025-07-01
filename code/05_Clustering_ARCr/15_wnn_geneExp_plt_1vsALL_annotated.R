########################################################################
## Compute and Plot gene expression plots for top marker genes for one cell type 
## I used Deconvobuddies::findMarkers_1vAll(), a convenient wrapped (Scran/Deconvobuddies) to compute test.type="binom" (1vsALL)
## - Used: Default direction = "up"
##
## Authors. CSC
## Date. Jun 30, 2025
##
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
## DeconvoBuddies 1.1+ is not available for Bioconductor 3.18 (conda_R/4.3.x), so I am using conda_R/4.4.x instead
########################################################################

library("SingleCellExperiment")
library("DeconvoBuddies")
library("purrr")
library("dplyr")
library("stringr")
library("ggplot2")
library("here")


## directories

#inputSCE_Dir <- "~/Habenula_Visium/processed-data/05_snRNA-seq_model_stats/"
inputSCE_Dir <- here(
    "processed-data", 
    "08_spatial_registration_vs_multiome_snRNA-seq"
)
processedDir <- here(
    "processed-data",
    "05_Clustering_ARCr",
    "15_wnn_geneExp_plt_1vsALL_annotated"
)
plotDir <- here(
    "plots",
    "05_Clustering_ARCr",
    "15_wnn_geneExp_plt_1vsALL_annotated"
)

## Check directories
if (!dir.exists(plotDir)) {
    dir.create(plotDir)
}
if (!dir.exists(processedDir)) {
    dir.create(processedDir)
}

# no longer required - fixed on recent deconvobuddies release
#source(here("code", "05_Clustering_ARCr", "get_mean_ratio_sparse.R"))
#get_mean_ratio_sparse

#===============================================================================

## SingleCellExperiment object derived from Seurat rna modality
sce_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_v4.rds"
sce_name <- here(inputSCE_Dir, sce_base_name)
sce_name
title_name <- str_extract(sce_base_name, regex("C\\.\\w*\\_r2"))
title_name

# Load Seurat
sce <- readRDS(sce_name)
## verification
sce
# class: SingleCellExperiment 
# dim: 29690 55702 
# metadata(0):
#     assays(3): counts logcounts scaledata
# rownames(29690): MIR1302-2HG FAM138A ... AC007325.4 AC007325.2
# rowData names(2): gene_symbol gene_id
# colnames(55702): S04_AAACAGCCAGAATGAC-1 S04_AAACAGCCAGCAAGGC-1 ...
# S09_TTTGTTGGTCATGCAA-1 S09_TTTGTTGGTTGTTCAC-1
# colData names(31): orig.ident nCount_RNA ... cluster_ann ident
# reducedDimNames(0):
#     mainExpName: RNA
# altExpNames(0):

head(rownames(sce))

colnames(colData(sce))

## briefly check cell types and counts
cluster_counts_df <- as.data.frame(table(sce$cluster_ann))
head(cluster_counts_df)

# check rowData ensembl names and gene id(s)
head(rowData(sce)$gene_symbol)
head(rowData(sce)$gene_id)

#===============================================================================
## Prepare data: remove cell types with fewer than 10 cells

celltypes <- colData(sce)$cluster_ann
celltype_counts <- table(celltypes)
low_ct <- names(celltype_counts[celltype_counts <= 10])

# filter into the get_mean_ratio() function
if (length(low_ct)==TRUE) {
    message("Removing cell types with less<10 cells: ", low_ct)
    # Remove cell types with <= 10 cells
    valid_types <- names(celltype_counts[celltype_counts > 10])
    sce <- sce[, colData(sce)$cluster_ann %in% valid_types]
}

## =============================================================================
## Function to compute 1vsALL across several sizes of datasets
##==============================================================================

## set some settings
genes_of_interest <- c("GPR151", "TAC3", "POU4F1")
all_genes <- nrow(sce)
# sizes <- c(2000, 4000, 6000, 8000, 10000, 12000, all_genes)
## for speeding the process, I will only tests
sizes <- c(all_genes)

message("Subset sizes: ")
sizes


# Wrapper function to compute and extract ranks
get_marker_ranks_1vsALL <- function(sce, 
                                    gene_subset, 
                                    n, 
                                    celltype_regex = NULL,
                                    genes_of_interest = NULL) {
    # celltype_regex = NULL: No filtering on cellType.target
    # celltype_regex = "LHb|MHb": Filters clusters with names matching that regex
    # Optional gene filtering via genes_of_interest
    
    sce_sub <- sce[gene_subset, ]
    
    message("Processing mean-ratio for subset:", n)
    
    # Apply test.type="binom" (Scran/Deconvobuddies)
    One_vsALL_df <-  findMarkers_1vAll(
        sce_sub,
        assay_name = "logcounts",
        cellType_col = "cluster_ann",
        mod = NULL, 
        add_symbol = FALSE,
        verbose = TRUE,
        direction = "up"
    ) 
    
    message("1vsALL done!")
    
    return(One_vsALL_df)
}


rank_results <- lapply(sizes, function(n) {
    if (as.integer(n) < as.integer(all_genes)) {
        gene_subset <- head(order(Matrix::rowMeans(assay(sce, "logcounts")), decreasing = TRUE), n)
    } else {
        message("Processing full dataset!")
        gene_subset <- head(order(Matrix::rowMeans(assay(sce, "logcounts")), decreasing = TRUE))
    }
    get_marker_ranks_1vsALL(
        sce,
        gene_subset,
        n
    )
})

## verification

length(rank_results)
names(rank_results) <- paste0("size_", sizes)
names(rank_results)
# [1] "size_29690"

rank_summary <- bind_rows(rank_results)
summary(rank_summary)
# gene               logFC           log.p.value          log.FDR       
# Length:246         Min.   :-2.66888   Min.   :-2532.12   Min.   :-2530.3  
# Class :character   1st Qu.:-0.23933   1st Qu.: -195.72   1st Qu.: -195.1  
# Mode  :character   Median : 0.13408   Median :  -22.79   Median :  -22.5  
# Mean   :-0.02792   Mean   : -224.27   Mean   : -223.7  
# 3rd Qu.: 0.47191   3rd Qu.:    0.00   3rd Qu.:    0.0  
# Max.   : 1.56867   Max.   :    0.00   Max.   :    0.0  
# std.logFC        cellType.target    std.logFC.rank std.logFC.anno    
# Min.   :-2.53057   Length:246         Min.   :1.0    Length:246        
# 1st Qu.:-0.22443   Class :character   1st Qu.:2.0    Class :character  
# Median : 0.16215   Mode  :character   Median :3.5    Mode  :character  
# Mean   :-0.03141                      Mean   :3.5                      
# 3rd Qu.: 0.47591                      3rd Qu.:5.0                      
# Max.   : 1.27905                      Max.   :6.0    


message("Check clusters present:")

print(sort(unique(rank_summary$cellType.target)))
length(unique(rank_summary$cellType.target))
# 41

##==============================================================================

message("save marker stats / rank_results")

# marker_ranks_6000 <- rank_results[["size_6000"]]
# marker_ranks_8000 <- rank_results[["size_8000"]]
# marker_ranks_10000 <- rank_results[["size_10000"]]
# marker_ranks_12000 <- rank_results[["size_12000"]]
marker_ranks_all <- rank_results[["size_29690"]]
f_name <- here(processedDir, "marker_1vsALL_29k.RData")
save(marker_ranks_all, file = f_name)

message("Mean ratio results saved !!!")


##==============================================================================
## plots the top n marker genes for a specified cell type based off of the stats table from get_mean_ratio()

message("Prepare data to plot top 'n' marker genes")

sorted_levels <- sort(unique(colData(sce)$cluster_ann))
levels(sce$cluster_ann)
head(sce[["cluster_ann"]])

marker_stats <- rank_results[["size_29690"]]
print(marker_stats, n=50)
# # A tibble: 6 × 8
# # Groups:   cellType.target [1]
# gene      logFC log.p.value   log.FDR std.logFC cellType.target std.logFC.rank
# <chr>     <dbl>       <dbl>     <dbl>     <dbl> <chr>                    <int>
#     1 SNHG14  3.92e-1   -1.41e+ 2 -1.39e+ 2   3.71e-1 C.25.Excit.Thal          1
# 2 MALAT1  1.18e-1   -1.33e+ 1 -1.22e+ 1   2.20e-1 C.25.Excit.Thal              2
# 3 CNTNA… -1.12e-4   -6.91e- 1 -3.27e-11  -8.63e-5 C.25.Excit.Thal              3
# 4 NRXN1  -1.27e-1   -3.33e- 4 -3.27e-11  -1.17e-1 C.25.Excit.Thal              4
# 5 LRP1B  -2.24e-1   -4.10e- 6 -3.27e-11  -1.79e-1 C.25.Excit.Thal              5
# 6 CADM2  -2.71e-1   -3.27e-11 -3.27e-11  -2.64e-1 C.25.Excit.Thal              6
# # ℹ 1 more variable: std.logFC.anno <chr>

valid_clusters <- unique(marker_stats$cellType.target)
sorted_ct_valid <- sorted_levels[sorted_levels %in% valid_clusters]

message("Plotting 1vsALL by cell type")

f_name <- here(plotDir ,"VPlot_1vsALL_wnn_cluster_genes29k.pdf")
pdf(file = f_name, width = 8.5, height = 11)  # standard letter size

for (ct in sorted_ct_valid) {
    message("Plotting wnn cluster: ", ct)
    
    # Filter top 10 genes for this cluster with log.FDR < -0.05
    top_genes <- marker_stats |>
        filter(cellType.target == ct, log.FDR < 0.05) |>
        arrange(std.logFC.rank) |>
        slice_head(n = 10) |>
        pull(gene)
    
    if (length(top_genes) == 0) {
        message("No significant genes found for ", ct, ". Skipping.")
        next
    } else { message(length(top_genes), " passed the filter.") }
    
    ## plot expression of top 10 genes
    p1 <- plot_gene_express(
        sce = sce,
        category = "cluster_ann",
        genes = top_genes
    ) + 
        ggtitle(paste0("Cluster: ", ct, " — Top 10 marker genes")) +
        labs(caption= paste0("Top10 on log.FDR < 0.05 / Binom-Test, Direction UP"))
    print(p1)
}

dev.off()

message("Top 10 mean-ratio plots done!")

