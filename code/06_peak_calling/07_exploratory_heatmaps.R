########################################################################
## explore and filter strong peak-gene links for evaluation 
## - Make several visualization to evaluate Peak scores
## - Make table with several confidence Peak scores
##
## Authors. CSC
## Date. Aug 11, 2025
## Recommended resources on interactive mode: srun --pty --mem=30GB --x11 bash
## Note. Seurat objects were created with module load conda_R/4.4.x
########################################################################

library("Seurat")
library("Signac")
library("ggplot2")
library("patchwork")
library("tidyverse")
library("dplyr")
library("here")


## read input arguments
args = commandArgs(trailingOnly = TRUE)
p_met <- args[2]
w_size <- args[4]
# 0: pearson, 1e5
# 1: pearson, 5e4
# 2: spearman, 1e5
# 3: spearman, 2.5e4

# for testing
# p_met = "spearman"
# w_size = 2.5e4

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
    f_sufix <- paste0(".", p_met, ".", format(w_size_label, scientific = TRUE), ".cells_filtered_2perc")
    message("Processing: ", f_sufix)
} else {
    message("Input arguments missed")
    stop()
}


# Check/create directories

## clusters renamed for Spatial-Registration on Visium project
input_cvsDir <- here(
    "processed-data",
    "06_peak_calling",
    "06_exploratory_peak_scores"
)
plotDir <- here(
  "plots",
  "06_peak_calling",
  "07_exploratory_heatmaps"
)
csvDir <- here(
    "processed-data",
    "06_peak_calling",
    "07_exploratory_heatmaps"
)
## Seurat with wnn final ct
inputRDS_Dir <- here(
    "processed-data",
    "05_Clustering_ARCr",
    "17_wnn_clustering_final_ct"
)
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_HD.rds"


## Check directories
if (!dir.exists(plotDir)) {
  dir.create(plotDir)
}
if (!dir.exists(csvDir)) {
    dir.create(csvDir)
}

## for testing: ================================================================
#link_df2 <- read.csv(file = here(input_cvsDir, "peak_gene_links_with_TSS_and_CC.spearman.2.5e4.cells_filtered_2perc.csv"))
## for testing: ================================================================

# load link peak-gene with TSS anc CC
gene_peaks_csv <- here(input_cvsDir, paste0("peak_gene_links_with_TSS_and_CC", f_sufix, ".csv"))

if (file.exists(gene_peaks_csv)) {
    link_df2 <- read.csv(gene_peaks_csv)
    message("File loaded!")
} else {
    stop(paste("File not found:", gene_peaks_csv))
}

message("Link gene-peak scores with TSS loaded ...")

## inspect data
head(link_df2)
table(link_df2$tier)
summary(link_df2)
table(link_df2$gene_strand, useNA = "ifany")
    # -    + 
    # 2714 2786 

#===============================================================================

message("Building plots ...")

#===============================================================================


## Prepare Seurat for Heatmap
# Use Seurat with clusters renamed for Spatial-Registration on Visium project
seurat_name <- here(inputRDS_Dir, Seurat_base_name)
SeuratOBJ <- readRDS(here(inputRDS_Dir, Seurat_base_name))

message("WNN clusters:")
levels(SeuratOBJ)

## filter genes to those expressed in 2% of cells
DefaultAssay(SeuratOBJ) <- "RNA"
rna_counts <- GetAssayData(SeuratOBJ, assay="RNA", layer="data")
length(rownames(rna_counts)) # [1] 36601
keep_genes <- rownames(rna_counts)[Matrix::rowSums(rna_counts > 0) > 0.02 * ncol(rna_counts)]

message("Genes kept after filtering those expressed in 2% of cells ")
length(keep_genes) # in count: [1] 14526

## Set atac as defaul assay
DefaultAssay(SeuratOBJ) <- "ATAC"
class(SeuratOBJ[["ATAC"]])

set.seed(12082025)

## assign each peak to a cluster (by accessibility)
## average peaks × clusters
avg_atac <- AggregateExpression(
    SeuratOBJ,
    assays = "ATAC",
    group.by = "cluster_ann",
    layers = "counts" # fragment counts per peak per cell
)$ATAC  # matrix: peaks x clusters

head(avg_atac)

##==============================================================================
## Make complex Heatmap

# for each peak, which cluster has highest average accessibility
peak_cluster <- data.frame(
    peak = rownames(avg_atac),
    cluster_ann = colnames(avg_atac)[max.col(avg_atac, ties.method = "first")]
)
head(peak_cluster)
nrow(peak_cluster) # [1] 262951

# keep only links whose peaks are known in the assay
links_cl <- link_df %>%
    inner_join(peak_cluster, by = "peak") %>%
    filter(!is.na(cluster_ann))

head(links_cl)
nrow(links_cl)
table(links_cl$tier)

## filter links of desired
links_cl <- links_cl %>% filter(tier=="Exploratory (0.10–0.20, FDR<0.10)")
nrow(links_cl) # ge. 468
table(links_cl$tier)
#links_cl <- links_cl %>% filter(score > 0.3, pvalue < 0.05)

# bin definition
bin_breaks <- c(0,1,2,3,5,10,Inf)
bin_labels <- c("1","2","3","4–5","6–10",">10")


##==============================================================================
## Distribution of number of linked peaks per gene, by cluster
# Count linked peaks per (gene, cluster)
peaks_per_gene <- links_cl %>%
    group_by(cluster_ann, gene) %>%
    summarise(n_peaks = n_distinct(peak), .groups = "drop") %>%
    mutate(bin = cut(
        n_peaks, breaks = bin_breaks,
        labels = bin_labels, right = TRUE
    ))
# plot summarize how many peaks each gene has in each cluster
raw_links <- nrow(links_cl)

## verification
peaks_per_gene[peaks_per_gene$cluster_ann=="C.19.Inhib.Thal", ]
# cluster_ann     gene      n_peaks bin  
# <chr>           <chr>       <int> <fct>
# 1 C.19.Inhib.Thal GAD1            1 1    
# 2 C.19.Inhib.Thal GAD2            1 1    
# 3 C.19.Inhib.Thal KIT             1 1    
# 4 C.19.Inhib.Thal LINC01210       1 1    
# 5 C.19.Inhib.Thal MEIS2           2 2    
# 6 C.19.Inhib.Thal OTX2-AS1        1 1    
# 7 C.19.Inhib.Thal SOX14           1 1   
nrow(peaks_per_gene[peaks_per_gene$cluster_ann=="C.19.Inhib.Thal", ])


# Make a cluster × bin table
dist_pg <- peaks_per_gene %>%
    count(cluster_ann, bin, name = "n_genes") %>%
    complete(cluster_ann, bin, fill = list(n_genes = 0)) %>%
    group_by(cluster_ann) %>%
    mutate(freq = n_genes / sum(n_genes)) %>%   # normalize by row
    ungroup()
head(dist_pg)

x_label <- paste0("Peaks per gene (bins)\n raw-links(", nrow(links_cl), ")")

# Make heatmap
g_dp <- ggplot(dist_pg, aes(x = bin, y = cluster_ann, fill = freq)) +
    geom_tile(color = "grey85") +
    scale_fill_viridis_c(name = "Fraction of genes", option = "C") +
    labs(
        title = "Peaks per gene",
        subtitle = "Exploratory (0.10–0.20, FDR<0.10)",
        x = x_label, y = "cluster_ann"
    ) +
    theme_minimal()


##==============================================================================

## Distribution of number of linked genes per peak, by cluster
# Count linked genes per (peak, cluster)
genes_per_peak <- links_cl %>%
    group_by(cluster_ann, peak) %>%
    summarise(n_genes = n_distinct(gene), .groups = "drop") %>%
    mutate(bin = cut(
        n_genes, breaks = c(0,1,2,3,5,10,Inf),
        labels = c("1","2","3","4-5","6-10",">10"), right = TRUE
    ))

dist_gp <- genes_per_peak %>%
    count(cluster_ann, bin, name = "n_peaks") %>%
    complete(cluster_ann, bin, fill = list(n_peaks = 0)) %>%
    group_by(cluster_ann) %>%
    mutate(freq = n_peaks / sum(n_peaks)) %>%
    ungroup()

x_label <- paste0("Genes per peak (bins)\n raw-links(", nrow(links_cl), ")")

g_dg <- ggplot(dist_gp, aes(x = bin, y = cluster_ann, fill = freq)) +
    geom_tile(color = "grey85") +
    scale_fill_viridis_c(name = "Fraction of peaks", option = "C") +
    labs(
        title = "Genes per peak",
        subtitle = "Exploratory (0.10–0.20, FDR<0.10)",
        x = x_label, y = "cluster_ann"
    ) +
    theme_minimal()

# both plots set to use same color mapping & labels
g_dp <- g_dp + labs(color = "Tier") + theme(legend.position = "bottom")
g_dg <- g_dg + labs(color = "Tier") + theme(legend.position = "bottom")

combined_plot <- g_dp + g_dg + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

pdf(here(plotDir, paste0("heatmaps_link_peak_genes_peaks", f_sufix, ".pdf")), width = 10, height = 6)
combined_plot +
    plot_annotation(
        title = "Peak-Gene Link Distributions",
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12)),
        caption = f_sufix
    ) 
dev.off()

##==============================================================================

##  Top genes by total linked peaks, across clusters
# Pick top N genes by total linked peaks (across all clusters)

N <- 50

top_genes <- links_cl %>%
    group_by(gene) %>%
    summarise(total_peaks = n_distinct(peak), .groups = "drop") %>%
    slice_max(total_peaks, n = N) %>%
    pull(gene)

mat_pg <- links_cl %>%
    filter(gene %in% top_genes) %>%
    group_by(cluster_ann, gene) %>%
    summarise(n_peaks = n_distinct(peak), .groups = "drop") %>%
    mutate(gene = factor(gene, levels = top_genes)) %>%
    complete(cluster_ann, gene, fill = list(n_peaks = 0))

# Heatmap: genes (rows) × clusters (cols)
g1 <- ggplot(mat_pg, aes(x = cluster_ann, y = gene, fill = n_peaks)) +
    geom_tile() +
    scale_fill_viridis_c(name = "# linked peaks") +
    labs(
        title = paste("Top", N, "genes: linked peaks per gene by cluster"),
        x = "WNN cluster", y = "Gene"
    ) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 7))
g1

ggsave(here(plotDir,
            paste0("heatmap_Top50_genes_peaks", f_sufix, ".pdf")),
       g1, width = 8, height = 5)



# ## ====
# 
# library(ComplexHeatmap)
# library(circlize)
# 
# bin_breaks  <- c(0,1,2,3,5,10,Inf)
# bin_labels  <- c("1","2","3","4-5","6-10",">10")   # <-- ASCII hyphen
# 
# #1) Heatmap: Distribution of peaks per gene (binned), by cluster
# #Counts how many distinct peaks are linked to each gene within each cluster, bins those counts, and shows the fraction per cluster.
# 
# peaks_per_gene <- links_cl %>%
#     group_by(cluster_ann, gene) %>%
#     summarise(n_peaks = n_distinct(peak), .groups = "drop") %>%
#     mutate(bin = cut(n_peaks, breaks = bin_breaks, labels = bin_labels, right = TRUE))
# 
# dist_pg <- peaks_per_gene %>%
#     count(cluster_ann, bin, name = "n_genes") %>%
#     complete(cluster_ann, bin = factor(bin_labels, levels = bin_labels), fill = list(n_genes = 0)) %>%
#     group_by(cluster_ann) %>%
#     mutate(freq = n_genes / sum(n_genes)) %>%
#     ungroup()
# 
# # matrix: rows = clusters, cols = bins
# mat_pg <- dist_pg %>%
#     select(cluster_ann, bin, freq) %>%
#     pivot_wider(names_from = bin, values_from = freq) %>%
#     as.data.frame()
# 
# rownames(mat_pg) <- mat_pg$cluster_ann
# mat_pg$cluster_ann <- NULL
# mat_pg <- as.matrix(mat_pg)  # numeric matrix
# 
# # color scale
# col_fun_pg <- circlize::colorRamp2(c(0, max(mat_pg, na.rm = TRUE)), c("#F7FBFF", "#08306B"))
# 
# ht_pg <- Heatmap(
#     mat_pg,
#     name = "Frac genes",
#     col = col_fun_pg,
#     cluster_rows = TRUE,
#     cluster_columns = FALSE,
#     row_title = "WNN clusters",
#     column_title = "Peaks per gene (binned)",
#     heatmap_legend_param = list(direction = "horizontal"),
#     column_names_rot = 0
# )
# 
# # Save as PDF
# pdf(here(plotDir, "heatmap_frac_genes", f_sufix, ".pdf"),
#     width = 8, height = 6)
# draw(ht_pg)
# dev.off()
# 
# 
# ##2) Heatmap: Distribution of genes per peak (binned), by cluster
# ##Counts how many distinct genes each peak links to within each cluster, bins those counts, shows the fraction per cluster.
# 
# genes_per_peak <- links_cl %>%
#     group_by(cluster_ann, peak) %>%
#     summarise(n_genes = n_distinct(gene), .groups = "drop") %>%
#     mutate(bin = cut(n_genes, breaks = bin_breaks, labels = bin_labels, right = TRUE))
# 
# dist_gp <- genes_per_peak %>%
#     count(cluster_ann, bin, name = "n_peaks") %>%
#     complete(cluster_ann, bin = factor(bin_labels, levels = bin_labels), fill = list(n_peaks = 0)) %>%
#     group_by(cluster_ann) %>%
#     mutate(freq = n_peaks / sum(n_peaks)) %>%
#     ungroup()
# 
# mat_gp <- dist_gp %>%
#     select(cluster_ann, bin, freq) %>%
#     pivot_wider(names_from = bin, values_from = freq) %>%
#     as.data.frame()
# 
# rownames(mat_gp) <- mat_gp$cluster_ann
# mat_gp$cluster_ann <- NULL
# mat_gp <- as.matrix(mat_gp)
# 
# col_fun_gp <- circlize::colorRamp2(c(0, max(mat_gp, na.rm = TRUE)), c("#FFF5F0", "#7F0000"))
# 
# ht_gp <- Heatmap(
#     mat_gp,
#     name = "Frac peaks",
#     col = col_fun_gp,
#     cluster_rows = TRUE,
#     cluster_columns = FALSE,
#     row_title = "WNN clusters",
#     column_title = "Genes per peak (binned)",
#     heatmap_legend_param = list(direction = "horizontal"),
#     column_names_rot = 0
# )
# 
# # Save as PDF
# pdf(here(plotDir, "heatmap_frac_peaks", f_sufix, ".pdf"),
#     width = 8, height = 6)
# draw(ht_gp)
# dev.off()
# 
# #3) Heatmap: Top N genes by linked peaks across clusters
# #Rows = genes (top N by total linked peaks), columns = clusters, values = # linked peaks.
# 
# top_genes <- links_cl %>%
#     group_by(gene) %>%
#     summarise(total_peaks = n_distinct(peak), .groups = "drop") %>%
#     slice_max(total_peaks, n = N) %>%
#     pull(gene)
# 
# mat_top <- links_cl %>%
#     filter(gene %in% top_genes) %>%
#     group_by(cluster_ann, gene) %>%
#     summarise(n_peaks = n_distinct(peak), .groups = "drop") %>%
#     complete(cluster_ann, gene, fill = list(n_peaks = 0)) %>%
#     pivot_wider(names_from = cluster_ann, values_from = n_peaks) %>%
#     as.data.frame()
# 
# rownames(mat_top) <- mat_top$gene
# mat_top$gene <- NULL
# mat_top <- as.matrix(mat_top)
# 
# # Optional: row-wise z-score to show relative enrichment per gene
# mat_top_z <- t(scale(t(mat_top)))  # center/scale per gene
# 
# col_fun_top <- circlize::colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
# 
# ## rows = the top N genes (most variable or most linked, depending on your earlier filtering)
# ht_top <- Heatmap(
#     mat_top_z,
#     name = "Z (gene)",
#     col = col_fun_top,
#     cluster_rows = TRUE,
#     cluster_columns = TRUE,
#     show_row_dend = TRUE,
#     show_column_dend = TRUE,
#     column_title = "Cluster annotation",           # X-axis label
#     column_title_side = "bottom",                  # put at bottom if desired
#     row_title = paste0("Top ", N, " genes"),       # Y-axis label
#     show_row_names = TRUE,
#     row_names_gp = grid::gpar(fontsize = 7)
# )
# 
# ht_distributions <- ht_pg %v% ht_gp
# ht_distributions
# 
# # draw in one page
# pdf(file.path(plotDir, paste0("heatmaps_distributions_complexheatmap", f_sufix, ".pdf")),
#               width = 8, height = 10)
# draw(ht_distributions, heatmap_legend_side = "right", merge_legend = TRUE)
# dev.off()
# 
# # top-genes heatmap on its own
# pdf(file.path(plotDir, paste0("heatmap_top_genes_complexheatmap", f_sufix, ".pdf")),
#     width = 8, height = 10)
# draw(ht_top, heatmap_legend_side = "right")
# dev.off()

message("Plots done!!!")


# library("slurmjobs")
# job_single(
#   "06_exploratory_peak_scores",
#   create_shell = TRUE,
#   partition = "katun",
#   memory = "30G",
#   cores = 2,
#   logdir = "logs",
#   command = "Rscript 06_exploratory_peak_scores.R",
#   create_logdir = TRUE
# )

## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
