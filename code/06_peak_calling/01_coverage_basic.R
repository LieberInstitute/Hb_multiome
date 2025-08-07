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
## I use BSgenome.Hsapiens.UCSC.hg38 for extracting DNA motifs, k-mers, sequence-based features and compute Tn5 bias correction
library("BSgenome.Hsapiens.UCSC.hg38")  # full reference genome sequence / actual DNA bases (A/T/C/G) for each chromosome
## I use EnsDb.Hsapiens.v86 for extracting gene names, positions, TSSs, chr locations, etc.
library("EnsDb.Hsapiens.v86")           # Gene annotation (GTF-style)
library("tidyverse")
library("dplyr")
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

## preprocessed peaks

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


## set a subset of genes to test
hb_cannonical_genes <- c("GPR151",  "POU4F1", "TAC3")

## find peaks that are correlated with the expression of nearby genes
atac <- LinkPeaks(
    object = SeuratOBJ,
    peak.assay = "ATAC",
    expression.assay = "RNA",
    genes.use = hb_cannonical_genes,
    # genes.use = top5_genes_habenula,      # supply vector of gene names
    method = "pearson",        # Default settings
    distance = 1e5             # Only consider peaks within ±100 kb of gene TSS (cis-window)
)


##==============================================================================
## explore and filter strong peak-gene links

## inspect data
head(Links(atac))
# GRanges object with 5 ranges and 5 metadata columns:
#     seqnames              ranges strand |     score        gene
# <Rle>           <IRanges>  <Rle> | <numeric> <character>
# [1]     chr5 146497025-146516190      * | 0.0515208      GPR151
# [2]     chr5 146516043-146516190      * | 0.0979604      GPR151
# [3]    chr13   78596294-78603560      * | 0.0569620      POU4F1
# [4]    chr13   78597457-78603560      * | 0.0579676      POU4F1
# [5]    chr13   78603450-78603560      * | 0.0594102      POU4F1
# peak    zscore      pvalue
# <character> <numeric>   <numeric>
#     [1] chr5-146496517-14649..   2.71746 3.28927e-03
# [2] chr5-146515546-14651..   5.29336 6.00436e-08
# [3] chr13-78595767-78596..   4.35176 6.75248e-06
# [4] chr13-78597000-78597..   4.33686 7.22667e-06
# [5] chr13-78602870-78604..   5.15894 1.24175e-07

link_df <- as.data.frame(Links(atac))
summary(link_df$score)

write.csv(
    link_df,
    file = here(plotDir, "all_peak_gene_links.csv"),
    row.names = FALSE
)

## compute the distance between each peak and its linked gene's TSS and add it to the link_df
# Get TSS per gene
gene_coords <- genes(EnsDb.Hsapiens.v86)
#head(gene_coords)
tss_coords <- resize(gene_coords, width = 1, fix = "start")
tss_coords <- keepStandardChromosomes(tss_coords, pruning.mode = "coarse")
#head(tss_coords)
# ensure chromosome names in TSS match UCSC-style peaks (chr1, chr2, etc.)
seqlevelsStyle(tss_coords) <- "UCSC"
head(tss_coords)

## Match genes in link_df to their TSS
# Make sure gene names match
link_df$gene <- as.character(link_df$gene)
tss_df <- as.data.frame(tss_coords)
tss_df <- tss_df[, c("seqnames", "start", "gene_id", "gene_name")]
head(tss_df)
#                   seqnames start         gene_id gene_name
# ENSG00000223972     chr1 11869 ENSG00000223972   DDX11L1
# ENSG00000227232     chr1 29570 ENSG00000227232    WASH7P
# ENSG00000278267     chr1 17436 ENSG00000278267 MIR6859-1

# Join link_df with TSS info by gene name
link_df <- left_join(link_df, tss_df, by = c("gene" = "gene_name"))
head(link_df, n = 3)
# remove rows where chromosomes don’t match
link_df <- link_df[as.character(link_df$seqnames.x) == as.character(link_df$seqnames.y), ]

## Compute the distance to TSS
# Compute center of each peak
link_df$peak_center <- (link_df$start.x + link_df$end) / 2

# Compute absolute distance to TSS
link_df$distance <- abs(link_df$peak_center - link_df$start.y)  # `start.y` is the TSS
# gives directionality (upstream = negative, downstream = positive).
link_df$signed_distance <- link_df$peak_center - link_df$start.y

h1 <- hist(link_df$distance / 1000, breaks = 100,
     main = "Distance from Peaks to TSS",
     xlab = "Distance (kb)",
     col = "lightblue")

f_name <- "histogram_scores_pearson_1000bp.pdf"
print(h1)
dev.off()
    
g1 <- ggplot(link_df, aes(x = distance / 1000, y = score)) +
    geom_point(alpha = 0.3, color = "steelblue") +
    geom_hline(yintercept = 0.3, linetype = "dashed", color = "red") +
    labs(
        x = "Distance from TSS (kb)",
        y = "Correlation Score",
        title = "Peak-Gene Correlation vs. Distance"
    ) +
    theme_minimal()

f_name <- here(plotDir, "distribution_scores_pearson_1000bp.pdf")
print(g1)
dev.off()



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
