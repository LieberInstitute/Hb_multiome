########################################################################
## Plot Coverage Plots on WNN clusters
##
## Authors. CSC
## Date. March 24, 2025
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
## Note. Seurat objects were created with module load conda_R/4.3.x
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
levels(SeuratOBJ)

DefaultAssay(SeuratOBJ) <- "ATAC"
class(SeuratOBJ[["ATAC"]])
Seurat_base_name <- str_extract(seurat_name, regex("C\\.\\w+"))
# C.leiden_lsi_r2_renamed_visium


## Read DEG to prepare coverage plots of the top 5 genes highly expressed

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


# Create 3 subsets of Seurat objects to plot: LHb, MHb, None

SeuratOBJ_LHb <- subset(SeuratOBJ, idents = lhb_clusters)
SeuratOBJ_MHb <- subset(SeuratOBJ, idents = mhb_clusters)
SeuratOBJ_None <- subset(SeuratOBJ, idents = noneHb_clusters)
rm("SeuratOBJ")


##### function to build coverage plots

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


# Make and combine plots in a single row

## canonical genes to plot

hg = 6
wd = 14

# size window to track in the plots
open_window_sizes = c(500,1000,2000)
# List of Seurat objects
seurat_objs <- list(
    MHb = SeuratOBJ_MHb,
    LHb = SeuratOBJ_LHb,
    None = SeuratOBJ_None
)
# Create all combinations
param_grid <- cross2(seurat_objs, open_window_sizes)

# Generate plot for POU4F1

features <- c("POU4F1", "GPR151", "TAC3")

for (gene in features) {
    
    feature <- gene
    coverage_plots <- map(param_grid, function(params) {
        seurat_obj <- params[[1]]
        window_size <- params[[2]]
        make_coverage_plot(seurat_obj, feature, window_size)
    })
    
    # Save plots to PDF
    if (length(coverage_plots)>1) {
        
        f_name <- paste0("CoveragePlot_all_clusters_canonical_", feature,".pdf")
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


# =======

# Generate plot for all canonical in 1 plot

# size window to track in the plots
open_window_sizes = c(500,1000)
# Create all combinations
param_grid <- cross2(seurat_objs, open_window_sizes)

coverage_plots <- map(param_grid, function(params) {
    seurat_obj <- params[[1]]
    window_size <- params[[2]]
    make_coverage_plot(seurat_obj, feature, window_size)
})

# Save plots to PDF
if (length(coverage_plots)>1) {
    
    f_name <- paste0("CoveragePlot_all_clusters_canonical_", paste(features, collapse = "_"), ".pdf")
    pdf(file = here(plotDir, f_name), width = 14, height = 6)
    # Loop through in chunks of 3
    for (i in seq(1, length(coverage_plots), by = 3)) {
        plots_chunk <- coverage_plots[i:min(i+2, length(coverage_plots))]
        combined_plot <- wrap_plots(plotlist = plots_chunk, ncol = 3)
        print(combined_plot)
    }
    
    dev.off()
}


message("Coverage plots completed")

# library("slurmjobs")
# job_single(
#   "01_coverage_basic",
#   create_shell = TRUE,
#   partition = "katun",
#   memory = "30G",
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
