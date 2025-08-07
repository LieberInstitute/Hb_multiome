########################################################################
## Baseline reference for region-specific accessibility
## Use pseudobulk-style visualization at ±500 bp
##     for Hb cannonical genes and top5 Hb DGE genes by cluster
##
## Output: coverage plots from WNN
##
## Authors. CSC
## Date. March 24, 2025
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
## Note. Seurat objects were created with module load conda_R/4.3.x
########################################################################

library("Seurat")
library("Signac")
library("ggplot2")
library("BSgenome.Hsapiens.UCSC.hg38")
library("tidyverse")
library("purrr")
library("here")


# Check/create directories

## clusters renamed for Spatial-Registration on Visium project
inputRDS_Dir <- here(
  "processed-data",
  "05_Clustering_ARCr",
  "17_wnn_clustering_final_ct"
)
inputCVS_Dir <- here(
  "processed-data",
  "05_Clustering_ARCr",
  "02_Hb_celltypes_from_seurat_reanalyze_v3",
  "cvs_files_markers"
)
plotDir <- here(
  "plots",
  "06_peak_calling",
  "01_coverage_basic"
)

## Check directories
if (!dir.exists(plotDir)) {
  dir.create(plotDir)
}

## Load Seurat
# Use Seurat with clusters renamed for Spatial-Registration on Visium project
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_HD.rds"
seurat_name <- here(inputRDS_Dir, Seurat_base_name)
SeuratOBJ <- readRDS(here(inputRDS_Dir, Seurat_base_name))
levels(SeuratOBJ)

DefaultAssay(SeuratOBJ) <- "ATAC"
class(SeuratOBJ[["ATAC"]])
# make a readable base-name for plots
Seurat_base_name <- str_extract(seurat_name, regex("C\\.\\w+"))
Seurat_base_name <- sub("_renamed_visium$", "", Seurat_base_name)
Seurat_base_name
# C.leiden_lsi_r2

##==============================================================================

## Identifies cis-regulatory elements by linking chromatin-accessible peaks to gene expression using correlation (and optionally accounting for covariates).

## GC content correction
genome <- BSgenome.Hsapiens.UCSC.hg38

SeuratOBJ <- RegionStats(
    object = SeuratOBJ,
    genome = genome,
    assay = "ATAC"  
)

## see all chromosomes
table(seqnames(granges(SeuratOBJ)))
## see what are considered standard chromosomes
standardChromosomes(granges(SeuratOBJ))

## Even though this filtering doesn’t change anything in this dataset, it ensures reproducibility
## remove the features that correspond to chromosome scaffolds or other sequences instead of the (22+2) standard chromosomes
peaks.keep <- seqnames(granges(SeuratOBJ)) %in% standardChromosomes(granges(SeuratOBJ))
tryCatch(
    {
        SeuratOBJ <- SeuratOBJ[as.vector(peaks.keep), ]
    }, error = function(e) {
        message(e)
    })


## find peaks that are correlated with the expression of nearby genes
atac <- LinkPeaks(
    object = SeuratOBJ,
    peak.assay = "ATAC",
    expression.assay = "RNA",
    # genes.use = top5_genes_habenula,      # supply vector of gene names
    method = "pearson",        # Default settings
    distance = 1e5             # Only consider peaks within ±100 kb of gene TSS (cis-window)
)


## =============================================================================
## Coverage plots with Habenula canonical genes and top5 DEG

message("Processing coverage plots for `POU4F1`, `GPR151` ... ")

## setup Extend Region
# based on distance for promoter or small gene
upstream = 1000
downstream = 1000

## make plot
get_coveragePlot <- function(
        SeuratOBJ,
        features, 
        upstream = 500, 
        downstream = 500,
        seurat_name
){
    plt1 <- CoveragePlot(
        object = SeuratOBJ, # Seurat_subset
        region = features,
        features = features,
        extend.upstream = upstream,
        extend.downstream = downstream,
        peaks = TRUE,
        links = TRUE
    ) +
        labs(title = paste0("Clusters from WNN: ", seurat_name)) +
        theme(
            text = element_text(size = 8),
            axis.text.x = element_text(size = 7),
            axis.text.y = element_text(size = 7),
            plot.title = element_text(hjust = 0.5)
        )
    return(plt1)
}


hb_cannonical_genes <- c("GPR151",  "POU4F1", "TAC3")

walk(hb_cannonical_genes, function(gene) {
    
    message("Plotting gene: ", gene)
    
    # Generate the plot
    p1 <- get_coveragePlot(
        SeuratOBJ,
        gene, upstream, downstream,
        paste0("Clusters from WNN: ", Seurat_base_name)
    )
    
    # save the plot
    f_name <- paste0(Seurat_base_name, "_coverage_", upstream, "bp_", gene, ".pdf")
    pdf(file = here(plotDir, f_name))
    print(p1)
    dev.off()
    
})



## =============================================================================
## Coverage plots of the top 5 genes highly expressed by cluster

# Read All DEG
DEG_file_name <- "WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_cellTypes_integrated_top50.csv"
DEG_file_name <- here(inputCVS_Dir, DEG_file_name)
df_cluster_names <- read.csv(DEG_file_name)
df_cluster_names <- df_cluster_names |> drop_na(cell_type)
head(df_cluster_names)

## Subset clusters annotated for `habenula`.

message("Cluster-IDs from `WNN`")

## set TRUE if you wish to plot only Hb clusters
hb_clusters_only = FALSE

if (hb_clusters_only) {
    Seurat_subset <- subset(SeuratOBJ, idents = grep("MHb|LHb", Idents(SeuratOBJ), value = TRUE))
    hb_clusters <- levels(Seurat_subset)
    hb_clusters
    cluster_numbers <- as.numeric(sub("^C\\.(\\d+)\\..*$", "\\1", hb_clusters))
} else {
    all_clusters <- levels(SeuratOBJ)
    all_clusters
    cluster_numbers <- as.numeric(sub("^C\\.(\\d+)\\..*$", "\\1", all_clusters))   
}
    

##  Use length of clusters to extract clusters IDs
cluster_numbers

## filter the top 5
#unique(df_cluster_names$cluster)

top5 <- df_cluster_names |>
  filter(cluster %in% cluster_numbers) |>
  group_by(cluster) |>
  top_n(n = 5, wt = avg_log2FC)
head(top5)


##==============================================================================

## plot top5 DEG on WNN

for (clus in unique(top5$cluster)) {
  # testing: clus = 5
  
  PDF_name <- paste0(Seurat_base_name, "_coverage", upstream, "bp_", clus, ".pdf") #"_coverage_hb_cluster_",
  pdf(file = here(plotDir, PDF_name))
  
  message("Processing cluster: ", clus, "; Save as: ", PDF_name)

  top5_cluster <- top5 |>
    filter(cluster == clus)

  walk(
    seq_along(top5_cluster$gene),
    ~ {
      tryCatch(
        {
        message(paste0("Processing gene ", top5_cluster$gene[.x]))
        ## make plot gene .x
        p1 <- make_coveragePlot(SeuratOBJ,
                                top5_cluster$gene[.x], upstream, downstream,
                                "")
        print(p1)
        },
        error = function(e) {
          message(paste0(
            "Error occurred while processing gene ",
            top5_cluster$gene[.x],
            ": ",
            e$message
          ))
        }
      )
    }
  )

  dev.off()
  
}


# Error occurred while processing gene AC109466.1: Gene not found
# Error occurred while processing gene LINC02143: Gene not found

message("Coverage plots for top 5 genes completed!")





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
