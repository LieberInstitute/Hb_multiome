########################################################################
## Plot Coverage Plots on 3 groups of clusters: MHb, LHb and None Hb
##
## Authors. CSC
## Date. May09, 2025
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
## Note. Seurat objects were created with module load conda_R/4.3.x
## - Processing the Seurat objects with newer versions could have unexpected results on the chromatin object
########################################################################

library("Seurat")
library("Signac")
library("ggplot2")
library("purrr")
library("patchwork")
library("tidyverse")
library("stringr")
library("here")

## input directories

here()

# Check/create directories

## clusters renamed for Spatial-Registration on Visium project
inputRDS_Dir <- here(
    "processed-data",
    "05_Clustering_ARCr",
    "08_wnn_gene_expression_plts_renamed_idents"
)
inputCVS_Dir <- here(
    "processed-data",
    "05_Clustering_ARCr",
    "02_Hb_celltypes_from_seurat_reanalyze_v3",
    "cvs_files_markers"
)
plotDir <- here(
    "plots",
    "06_peak_calling"
)

## Check directories
if (!dir.exists(plotDir)) {
    dir.create(plotDir)
}

## Load Seurat
# Use Seurat with clusters renamed for Spatial-Registration on Visium project
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium.rds"
seurat_name <- here(inputRDS_Dir, Seurat_base_name)
SeuratOBJ <- readRDS(here(inputRDS_Dir, Seurat_base_name))
#levels(SeuratOBJ)

DefaultAssay(SeuratOBJ) <- "ATAC"
class(SeuratOBJ[["ATAC"]])
# ensure the identities are set
Idents(SeuratOBJ) <- "seurat_clusters"


## Read DEG to select genes for the coverage. Top 5 genes highly expressed

DEG_file_name <- "WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_cellTypes_integrated_top50.csv"
DEG_file_name <- here(inputCVS_Dir, DEG_file_name)
dge_cluster_names <- read.csv(DEG_file_name)
dge_cluster_names <- dge_cluster_names |> drop_na(cell_type)
head(dge_cluster_names)
# p_val avg_log2FC pct.1 pct.2 p_val_adj cluster     gene            cell_type
# 1     0   4.336764 0.938 0.100         0       1 OTX2-AS1        DD_Inhib.Thal
# 2     0   3.849508 0.900 0.083         0       1      KIT        DD_Inhib.Thal
# 3     0   3.682782 0.927 0.122         0       1    MEIS2      LB_Thalamus/MDm
unique(dge_cluster_names$cluster)
# [1]  1  2  3  4  5  6  7  8  9 10 11 12 14 15 16 17 18 19 20 21 22 23 24 25 26
# [26] 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41


## Identified and subset clusters annotated in e groups

message("Cluster-IDs from `WNN`")

all_clusters <- unlist(levels(SeuratOBJ))
#unique(SeuratOBJ$seurat_clusters)

# Create 3 list with clusters to plot: LHb, MHb, None

lhb_clusters <- all_clusters[grepl("LHb", all_clusters)]
lhb_cluster_ids <- str_extract(lhb_clusters, "(?<=C\\.)\\d+")

mhb_clusters <- all_clusters[grepl("MHb", all_clusters)]
mhb_cluster_ids <- str_extract(mhb_clusters, "(?<=C\\.)\\d+")

noneHb_clusters <- all_clusters[!grepl("MHb|LHb", all_clusters)]
noneHb_clusters_ids <- str_extract(noneHb_clusters, "(?<=C\\.)\\d+")



################################################################################
# Create ATAC Coverage Plot with 3 columns side to side
# - compares canonical Hb gene markers
################################################################################

## Create 3 subsets of Seurat objects to plot: LHb, MHb, None

SeuratOBJ_LHb <- subset(SeuratOBJ, idents = lhb_clusters)
# some verification: fragments assigned after subset  and check gene annotations loaded 
#Fragments(SeuratOBJ_LHb)
#SeuratOBJ_LHb[["ATAC"]]@annotation
SeuratOBJ_MHb <- subset(SeuratOBJ, idents = mhb_clusters)
SeuratOBJ_None <- subset(SeuratOBJ, idents = noneHb_clusters)
rm("SeuratOBJ")


## function to build coverage plots from Seurat counts

make_coverage_plot <- function(seurat_subset, 
                                gene_list, 
                                region_size,
                                group_by = "seurat_clusters") {
    # Check inputs
    if (length(gene_list)==0) {
        stop(" gene_list empty")
    }
    
    # Generate coverage plots
    plt1 <- CoveragePlot(
        object = seurat_subset,
        region = gene_list,
        features = gene_list,
        group_by = group_by,
        extend.upstream = region_size,
        extend.downstream = region_size,
        peaks = TRUE,
        links = TRUE
    )
    plt1 <- plt1 +
        labs(title = paste0("Muti-coverage plot for gene: ", gene_list)) +
        theme(
            text = element_text(size = 6),
            axis.text.x = element_text(size = 6),
            axis.text.y = element_text(size = 5),
            plot.title = element_text(hjust = 0.5)
        )
    return(plt1)
}



## function to build coverage plots from Seurat pseudobulk data

make_coverage_plot_pseudobulk_by_cluster <- function(seurat_subset, 
                               gene_list, 
                               region_size,
                               group_by = "seurat_clusters") {
    # Check inputs
    if (length(gene_list)==0) {
        stop(" gene_list empty")
    }
    
    # Generate coverage plots
    plt1 <- CoveragePlot(
        object = seurat_subset,
        region = gene_list,
        features = gene_list,
        group.by = "seurat_clusters",  # pseudobulk by cluster
        extend.upstream = region_size,
        extend.downstream = region_size,
        peaks = TRUE,
        links = TRUE
    )
    plt1 <- plt1 +
        labs(title = paste0("Muti-coverage plot for gene: ", gene_list)) +
        theme(
            text = element_text(size = 6),
            axis.text.x = element_text(size = 6),
            axis.text.y = element_text(size = 5),
            plot.title = element_text(hjust = 0.5)
        )
    return(plt1)
}

## prepara data to plot canonical genes to plot

# set width and heigh to plot
hg = 6
wd = 14

# window size in the coverage graph track to loop
open_window_sizes = c(500,1000,2000)

# List of Seurat objects
seurat_objs <- list(
    MHb = SeuratOBJ_MHb,
    LHb = SeuratOBJ_LHb,
    None = SeuratOBJ_None
)
# Create all combinations
param_grid <- cross2(seurat_objs, open_window_sizes)

# Generate individual coverage plot for each gene in the list of 'features'
features <- c("POU4F1", "GPR151", "TAC3")

# ## testing raw vs pseudobulk
# 
# plt <- make_coverage_plot(SeuratOBJ, "POU4F1", 2000)
# plt1r <- make_coverage_plot(subset(SeuratOBJ, idents = lhb_clusters[1]), "POU4F1", 4000)
# 
# # Balancing number of cells per cluster. Downsampling for visualization purposes
# balanced_obj <- subset(SeuratOBJ, downsample = 200)  # max 200 cells/cluster
# plt1 <- make_coverage_plot_pseudobulk_by_cluster(balanced_obj, "POU4F1", 2000)
# 
# # manual comparison - subset lost cell identities
# plt1 <- make_coverage_plot_pseudobulk_by_cluster(subset(balanced_obj, idents = lhb_clusters[1]), "POU4F1", 2000)
# plt2 <- make_coverage_plot_pseudobulk_by_cluster(subset(balanced_obj, idents = mhb_clusters[1]), "POU4F1", 2000)


################################################################################
# Create 3 ATAC Coverage Plots for Hb genes ("POU4F1", "GPR151", "TAC3")
# - Plot both, raw accessibility signals, and pseudobulk signals for comparison purposes 
# - Plots are arranged by columns
################################################################################

for (gene in features) {
    
    # gene="POU4F1"
    feature <- gene
    message("Gene to track: ", feature)
    
    ## Plot raw (individual-cell) accessibility signals
    ## The plot has still pseudobulk data, but globally, not by cluster, condition, or any group.
    
    coverage_plots <- map(param_grid, function(params) {
        seurat_obj <- params[[1]]
        window_size <- params[[2]]
        make_coverage_plot(seurat_obj, feature, window_size)
    })

    # Save plots to PDF
    if (length(coverage_plots)>1) {

        f_name <- paste0("Multiome_peaks_gene_", feature,".pdf")
        pdf(file = here(plotDir, f_name), width = wd, height = hg)

        # Loop through in chunks of 3
        for (i in seq(1, length(coverage_plots), by = 3)) {
            plots_chunk <- coverage_plots[i:min(i+2, length(coverage_plots))]
            combined_plot <- wrap_plots(plotlist = plots_chunk, ncol = 3)
            print(combined_plot)
        }

        dev.off()
    }
    

    ## Plot pseudobulk (individual-cell) accessibility signals by cluster
    
    ## Note. I can see peaks in the raw (all-cells) CoveragePlot(), but not after pseudobulking by cluster
    # The issue is likely related to cell mismatches. e.g., cells not linked properly to fragment (need to reattach the  fragment file)
    
    coverage_plots <- map(param_grid, function(params) {
        seurat_obj <- params[[1]]
        window_size <- params[[2]]
        make_coverage_plot_pseudobulk_by_cluster(seurat_obj, feature, window_size)
    })
    
    
    # Save plots to PDF
    if (length(coverage_plots)>1) {
        
        f_name <- paste0("Multiome_peaks_gene_", feature,"_pseudobulk_by_cluster.pdf")
        pdf(file = here(plotDir, f_name), width = 14, height = 6)
        
        # Loop through in chunks of 3
        for (i in seq(1, length(coverage_plots), by = 3)) {
            plots_chunk <- coverage_plots[i:min(i+2, length(coverage_plots))]
            combined_plot <- wrap_plots(plotlist = plots_chunk, ncol = 3)
            print(combined_plot)
        }
        
        dev.off()
    }
    
    
}


################################################################################
# Create One ATAC Coverage Plot for all clusters
################################################################################

# Generate plot for all canonical in 1 plot

features <- c("POU4F1", "GPR151")
# window size in the coverage graph track to loop
window_size = 500
# List of Seurat objects
seurat_objs <- list(
    MHb = SeuratOBJ_MHb,
    LHb = SeuratOBJ_LHb
)

# size window to track in the plots
open_window_sizes = c(500)
# Create all combinations
param_grid <- cross2(seurat_objs, open_window_sizes)
#param_grid[[1]]

# # Generate names for each plot
# plot_names <- map_chr(param_grid, function(params) {
#     obj_name <- names(seurat_objs)[sapply(seurat_objs, identical, params[[1]])]
#     paste0(obj_name, "_win", params[[2]])
# })
# plot_names

message("Genes to track: ", paste(features, collapse = ","))

coverage_plots <- map(param_grid, function(params) {
    seurat_obj <- params[[1]]
    window_size <- params[[2]]
    make_coverage_plot(seurat_objs, features, window_size)
})
length(coverage_plots)
# # Add names
# names(coverage_plots) <- paste0(param_grid$obj, "_win", param_grid$window)
# names(coverage_plots) 

# Save plots
if (length(coverage_plots)>1) {
    
    f_name <- paste0("CoveragePlot_all_clusters_canonical_", paste(features, collapse = "_"), ".pdf")
    pdf(file = here(plotDir, f_name), width = wd, height = hg)
    
    # Loop through each plot individually
    for (i in seq_along(coverage_plots)) {
        # Add title
        plot_title <- if (!is.null(names(coverage_plots))) {
            names(coverage_plots)[i]
        } else {
            paste("Plot", i)
        }
        combined_plot <- coverage_plots[[i]] + plot_annotation(title = plot_title)
        print(combined_plot)
    }
    dev.off()
    
}



################################################################################
# Create 3 ATAC Coverage Plots for Hb genes ("POU4F1", "GPR151", "TAC3")
# - Plot raw accessibility signals
# - Plots are arranged in one column
################################################################################

## load seurat with merged clusters for visualization purposes 
rm("SeuratOBJ")
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_MERGED.rds"
seurat_name <- here(inputRDS_Dir, Seurat_base_name)
SeuratOBJ_merged <- readRDS(here(inputRDS_Dir, Seurat_base_name))

DefaultAssay(SeuratOBJ_merged) <- "ATAC"
class(SeuratOBJ_merged[["ATAC"]])

# ensure the identities are set
Idents(SeuratOBJ_merged) <- "merged_cluster"
colnames(SeuratOBJ_merged@meta.data)
table(SeuratOBJ_merged$merged_cluster)
# LHb_merged   MHb_merged No-Hb_merged 
# 6883        10944        37875 

# set width and heigh to plot
hg = 6
wd = 6

# Define the genes
genes_to_plot <- c("TAC3", "GPR151", "POU4F1")

"TAC3" %in% rownames(SeuratOBJ_merged[["RNA"]])  # Should be TRUE

# Create one plot per gene category (LHb, MHb and No-Hb)

# plt <- make_coverage_plot(SeuratOBJ_merged, "GPR151", 500)
# plt <- make_coverage_plot(SeuratOBJ_merged, "GPR151", 500, "merged_cluster")

coverage_plots <- lapply(genes_to_plot, function(gene) {
    make_coverage_plot(SeuratOBJ_merged, gene, 500, "merged_cluster") + patchwork::plot_annotation(title = gene)
})

# Combine plots vertically
combined_plot <- wrap_plots(coverage_plots, ncol = 1)

# Save to PDF
pdf(here(plotDir, "multi_peaks_by_MEGED_Hb_category_vertical.pdf"), width = 10, height = 12)
print(combined_plot)
dev.off()


message("Coverage plots completed")

# library("slurmjobs")
# job_single(
#   "03_coverage_3columns.R",
#   create_shell = TRUE,
#   partition = "katun",
#   memory = "60G",
#   cores = 2,
#   logdir = "logs",
#   command = "Rscript -e \"options(width = 120); sessioninfo::session_info()\"",
#   create_logdir = TRUE
# )

## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
