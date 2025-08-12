########################################################################
## Pseudobulk ATAC peaks from your SeuratOBJ and then explore LinkPeaks() correlation scores
##
## Tables to asses links at different thresholds 
##
## Authors. CSC
## Date. March 24, 2025
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
## Note. Seurat objects were created with module load conda_R/4.3.x
########################################################################

library("Seurat")
library("Signac")
library("Matrix")
library("ComplexHeatmap")
## I use BSgenome.Hsapiens.UCSC.hg38 for extracting DNA motifs, k-mers, sequence-based features and compute Tn5 bias correction
# library("BSgenome.Hsapiens.UCSC.hg38")  # full reference genome sequence / actual DNA bases (A/T/C/G) for each chromosome
library("stringr")
library("tidyverse")
library("dplyr")
library("purrr")
library("here")


## read input arguments
args = commandArgs(trailingOnly = TRUE)
p_met <- args[2]
w_size <- args[4]
# 0: pearson, 1e5
# 1: pearson, 5e4
# 2: spearman, 1e5
# 3: spearman, 5e4

# for testing
# p_met = "spearman"
# w_size = 5e4

## p_met:
# pearson -> peak-scores<0.2: likely due scATAC counts are ultra‑sparse; scRNA is zero‑inflated. Pearson r’s of 0.05–0.2 are common even for real links
# spearman -> as enhancer → gene relationships aren’t strictly linear; Pearson seems to underestimates. I will try spearman, more robust to nonlinearity/zeros

if (length(p_met) && length(w_size)) {
    message(
        "Processing job for peak-method:\n",
        p_met,
        "\nWindow-size\n",
        w_size
    )
    ## Use numeric comparison first, then assign string labels
    w_size_label <- case_when(
        isTRUE(all.equal(w_size, 25000))  ~ "2.5e4",
        isTRUE(all.equal(w_size, 50000))  ~ "5e4",
        isTRUE(all.equal(w_size, 100000)) ~ "1e5",
        TRUE                              ~ "00"
    )
    f_sufix <- paste0(".", p_met, ".", w_size_label, ".cells_filtered_2perc")
    message("Processing: ", f_sufix)
} else {
    message("Input arguments missed")
    stop()
}


# Check/create directories

## clusters renamed for Spatial-Registration on Visium project
inputRDS_Dir <- here(
    "processed-data",
    "05_Clustering_ARCr",
    "17_wnn_clustering_final_ct"
)
input_cvsDir <- here(
    "processed-data",
    "06_peak_calling",
    #"00_link_peaks"
    "06_exploratory_peak_scores"
)
plotDir <- here(
    "plots",
    "06_peak_calling",
    "02_link_peaks_pseudobulk"
)
cvsDir <- here(
    "processed-data",
    "06_peak_calling",
    "02_link_peaks_pseudobulk"
)

if (!dir.exists(plotDir)) {
    dir.create(plotDir)
}
if (!dir.exists(cvsDir)) {
    dir.create(cvsDir)
}


## Load Seurat
# Use Seurat with clusters renamed for Spatial-Registration on Visium project
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_HD.rds"
seurat_name <- here(inputRDS_Dir, Seurat_base_name)
SeuratOBJ <- readRDS(here(inputRDS_Dir, Seurat_base_name))
levels(SeuratOBJ)

message("Seurat loaded ...")

## filter genes to those expressed in 2% of cells
DefaultAssay(SeuratOBJ) <- "RNA"
rna_counts <- GetAssayData(SeuratOBJ, assay="RNA", layer="data")
length(rownames(rna_counts)) # [1] 36601
#length(rownames(rna_counts)[Matrix::rowSums(rna_counts > 0)]) # 34738
#length(rownames(rna_counts)[Matrix::rowSums(rna_counts > 0) > 0.02]) # 34738
keep_genes <- rownames(rna_counts)[Matrix::rowSums(rna_counts > 0) > 0.02 * ncol(rna_counts)]
length(keep_genes) # in count: [1] 14526

message("filter genes to those expressed in 2% of cells ...")

message("Extracting ATAC counts for pseudobulk ...")

DefaultAssay(SeuratOBJ) <- "ATAC"

atac_counts <- GetAssayData(SeuratOBJ, assay = "ATAC", layer = "counts")
head(atac_counts)
# get cell grouping info to pseudobulk by cluster
groups <- SeuratOBJ$cluster_ann
# Split cells by group and sum counts within each group

set.seed(11082025)

pseudobulk_counts <- map(unique(groups), function(g) {
    cells_in_group <- WhichCells(SeuratOBJ, idents = g)
    if (length(cells_in_group) == 1) {
        atac_counts[, cells_in_group]
    } else {
        Matrix::rowSums(atac_counts[, cells_in_group, drop = FALSE])
    }
}) |> 
    set_names(unique(groups)) |> 
    reduce(cbind)  # purrr's way to combine into a matrix

message("Pseudobulk completed!")

class(pseudobulk_counts)
# name columns by cluster
pseudobulk_counts <- as.matrix(pseudobulk_counts)
colnames(pseudobulk_counts) <- unique(groups)
message("Dimension on atac pseudobulk mtx is:")
dim(pseudobulk_counts)

# one column per cluster, representing aggregated ATAC peak counts
pseudobulk_counts[1:5, 1:5]

message("pseudobulk done!!")


message("Explore score distribution from LinkPeaks cvs file ...")

# load previous link peak-gene csv
#gene_peaks_csv <- here(input_cvsDir, paste0("all_peak_gene_links", f_sufix, ".csv"))
# peak_gene_links_with_TSS_and_CC.pearson.2.5e4.cells_filtered_2perc.csv
gene_peaks_csv <- here(input_cvsDir, paste0("peak_gene_links_with_TSS_and_CC", f_sufix, ".csv"))

if (file.exists(gene_peaks_csv)) {
    link_df <- read.csv(gene_peaks_csv)
    message("File loaded!")
} else {
    stop(paste("File not found:", gene_peaks_csv))
}
colnames(link_df)
head(link_df)

ggplot(link_df, aes(x = score)) +
    geom_histogram(binwidth = 0.02, fill = "steelblue", color = "white") +
    labs(
        title = "Distribution of Peak–Gene Correlation Scores",
        x = "Correlation score",
        y = "Count"
    ) +
    theme_minimal()

# Find top correlated peaks
high_conf <- link_df %>% filter(score > 0.1, FDR < 0.1)
high_conf
nrow(high_conf)

# Subset pseudobulk counts to these peaks
high_conf_counts <- pseudobulk_counts[rownames(pseudobulk_counts) %in% high_conf$peak, ]

# Heatmap of these peaks across clusters
Heatmap(
    high_conf_counts,
    name = "Counts",
    cluster_rows = TRUE,
    cluster_columns = TRUE
)


# --- Packages ---
## build a pseudobulk score–distance heatmap
# Assign each peak to a wnn cluster where it’s most accessible
# For each cluster, bin links by distance from TSS (e.g., every 10 kb up to 100 kb)
# Compute a summary per bin (median score or fraction above a threshold)
# Plot a ComplexHeatmap: rows = clusters, columns = distance bins

library(circlize)

DefaultAssay(SeuratOBJ) <- "ATAC"
colnames(link_df)
# Pseudobulk: assign each peak to the cluster with max accessibility ---

avg_atac <- AggregateExpression(
    SeuratOBJ, assays = "ATAC", group.by = "cluster_ann", layers = "counts"
)$ATAC

# Map peak -> cluster where it is most accessible
peak_cluster <- data.frame(
    peak = rownames(avg_atac),
    cluster_ann = colnames(avg_atac)[max.col(avg_atac, ties.method = "first")]
)

# Attach cluster to links (keep only links whose peaks are in the assay)
links_cl <- link_df %>%
    inner_join(peak_cluster, by = "peak") %>%
    filter(!is.na(cluster_ann), !is.na(distance_kb), !is.na(score))
head(links_cl)

# Keep a fixed cluster order (all clusters)
all_clusters <- sort(unique(SeuratOBJ$cluster_ann))
links_cl$cluster_ann <- factor(links_cl$cluster_ann, levels = all_clusters)

# Bin links by distance. Breaks: 0-10,10-20,...,90-100, >100
bin_breaks  <- c(seq(0, 100, by = 10), Inf)   # length = 12
edges       <- seq(0, 100, by = 10)           # 0,10,...,100
bin_labels  <- c(paste(head(edges, -1), tail(edges, -1), sep = "-"), ">100")

length(bin_labels) # 11 -> OK
stopifnot(length(bin_labels) == length(bin_breaks) - 1)

links_cl <- links_cl %>%
    mutate(dist_bin = cut(distance_kb, breaks = bin_breaks, labels = bin_labels, right = TRUE))

# Summaries per (cluster, dist_bin) ---
# A) Median score
summ_median <- links_cl %>%
    group_by(cluster_ann, dist_bin) %>%
    summarize(median_score = median(score, na.rm = TRUE), n = dplyr::n(), .groups = "drop") %>%
    complete(cluster_ann = all_clusters, dist_bin = factor(bin_labels, levels = bin_labels),
             fill = list(median_score = NA_real_, n = 0))

mat_median <- summ_median %>%
    select(cluster_ann, dist_bin, median_score) %>%
    pivot_wider(names_from = dist_bin, values_from = median_score) %>%
    as.data.frame()
rownames(mat_median) <- as.character(mat_median$cluster_ann)
mat_median$cluster_ann <- NULL
mat_median <- as.matrix(mat_median)

# B) Fraction with score ≥ 0.2 (change threshold if you like)
thr <- 0.2
summ_frac <- links_cl %>%
    group_by(cluster_ann, dist_bin) %>%
    summarize(frac_ge_thr = mean(score >= thr, na.rm = TRUE), n = dplyr::n(), .groups = "drop") %>%
    complete(cluster_ann = all_clusters, dist_bin = factor(bin_labels, levels = bin_labels),
             fill = list(frac_ge_thr = 0, n = 0))

mat_frac <- summ_frac %>%
    select(cluster_ann, dist_bin, frac_ge_thr) %>%
    pivot_wider(names_from = dist_bin, values_from = frac_ge_thr) %>%
    as.data.frame()
rownames(mat_frac) <- as.character(mat_frac$cluster_ann)
mat_frac$cluster_ann <- NULL
mat_frac <- as.matrix(mat_frac)

# Heatmaps (ComplexHeatmap)
col_fun_med  <- circlize::colorRamp2(
    c(min(mat_median, na.rm = TRUE), 0, max(mat_median, na.rm = TRUE)),
    c("#2166AC", "white", "#B2182B")
)
col_fun_frac <- circlize::colorRamp2(c(0, 0.5, 1), c("#F7FBFF", "#6BAED6", "#08306B"))

mat_median[is.na(mat_median)] <- 0 
ht_median <- Heatmap(
    mat_median,
    name = "Median score",
    col = col_fun_med,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    row_title = "Clusters",
    column_title = "Distance from TSS (kb, binned)",
    na_col = "grey95",
    heatmap_legend_param = list(direction = "horizontal"),
    column_names_rot = 0
)

mat_frac[is.na(mat_frac)] <- 0
ht_frac <- Heatmap(
    mat_frac,
    name = paste0("Frac score ≥ ", thr),
    col = col_fun_frac,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    row_title = "Clusters",
    column_title = "Distance from TSS (kb, binned)",
    heatmap_legend_param = list(direction = "horizontal"),
    column_names_rot = 0
)
mat_median[is.na(mat_median)] <- 0 

pdf(file.path(plotDir, "pseudobulk_score_distance_heatmaps.pdf"), width = 10, height = 8)
draw(ht_median %v% ht_frac, heatmap_legend_side = "bottom", merge_legend = FALSE)
dev.off()

# (optional) PNG export with better font handling
png(file.path(plotDir, "pseudobulk_score_distance_heatmaps.png"), width = 2200, height = 1600, res = 200)
draw(ht_median %v% ht_frac, heatmap_legend_side = "bottom", merge_legend = FALSE)
dev.off()





## Set atac
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

message("Starting GC content correction ... ")

## GC content correction
genome <- BSgenome.Hsapiens.UCSC.hg38

SeuratOBJ <- RegionStats(
    object = SeuratOBJ,
    genome = genome,
    assay = "ATAC"
)
SeuratOBJ

message("GC content correction done!")

message("Searching LinkPeaks ... ")

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
#hb_cannonical_genes <- c("GPR151",  "POU4F1", "TAC3")

message("Computing link-peaks correlations ...")

## find peaks that are correlated with the expression of nearby genes
atac <- LinkPeaks(
    object = SeuratOBJ,
    peak.assay = "ATAC",
    expression.assay = "RNA",
    genes.use = keep_genes,
    method = p_met,
    distance = as.numeric(w_size)             # Only consider peaks within ±100 kb of gene TSS (cis-window)
)

message("Computing link-peaks correlations completed!")

##==============================================================================
## explore and filter strong peak-gene links

## inspect data
head(Links(atac))
# GRanges object with 5 ranges and 5 metadata columns:
#     seqnames              ranges strand |     score        gene
#        <Rle>           <IRanges>  <Rle> | <numeric> <character>
# [1]     chr5 146497025-146516190      * | 0.0515208      GPR151
# [2]     chr5 146516043-146516190      * | 0.0979604      GPR151
# [3]    chr13   78596294-78603560      * | 0.0569620      POU4F1
# [4]    chr13   78597457-78603560      * | 0.0579676      POU4F1
# [5]    chr13   78603450-78603560      * | 0.0594102      POU4F1
# peak    zscore      pvalue
# <character> <numeric>   <numeric>
# [1] chr5-146496517-14649..   2.71746 3.28927e-03
# [2] chr5-146515546-14651..   5.29336 6.00436e-08
# [3] chr13-78595767-78596..   4.35176 6.75248e-06
# [4] chr13-78597000-78597..   4.33686 7.22667e-06
# [5] chr13-78602870-78604..   5.15894 1.24175e-07

link_df <- as.data.frame(Links(atac))
summary(link_df$score)

write.csv(
    link_df,
    file = here(cvsDir, paste0("all_peak_gene_links", f_sufix, ".csv")),
    row.names = FALSE
)

message("Peak-genes table saved!")

message("All done!!!")


# library("slurmjobs")
# job_single(
#   "00_link_peaks",
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
