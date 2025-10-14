########################################################################
## Explore FDR-Scores and Top overlapping (Heatmap): 
##
## Authors. CSC
## Date. Oct 14, 2025
## Recommended resources on interactive mode: srun --pty --mem=20GB --x11 bash
########################################################################

library("dplyr")
library("purrr")
library("ggplot2")
library("patchwork")
library("ggrepel")
library("tidyverse")
library("tidyr")
library("stringr")
library("here")

#===============================================================================
# resolution_level = "Mid"      # 18 cell-types
#===============================================================================

## setup variable names

# resolution_level = "Mid" 
FDR = 0.2 # actual value to run the DAR-Links
# we do not use log FC for this exploratory analysis

## Set directory names
inputCSV_Overlaps_Dir <- here(
    "processed-data",
    "06_peak_calling"
    #"20_explore_overlapping_linkpeaks_DARs"
    #"19_Linkage_DARs_analysis"
)
processedDir <- here(
    "processed-data",
    "06_peak_calling",
    "21_overlaping_FDRscores_TopHeatmap"
)
plotDir <- here(
    "plots",
    "06_peak_calling",
    "21_overlaping_FDRscores_TopHeatmap"
)

if (!dir.exists(processedDir)) { dir.create(processedDir) }
if (!dir.exists(plotDir)) { dir.create(plotDir) }

##==============================================================================

message("Loading Unique-Overlap with directionality ...")

cluster_specific = "LHb.4"
top_g = 20

## load raw overlaps df
overlap_df <- read.csv(here(inputCSV_Overlaps_Dir, "19_Linkage_DARs_analysis", 
                             paste0("Overlaps_LinkPeak_DARs_FDR", FDR, ".csv")))
shared_2ct_df <- read.csv(here(inputCSV_Overlaps_Dir, "20_explore_overlapping_linkpeaks_DARs", 
                             paste0("overlaps_summary_linkPeak_DARs_unique_ct_FDR", FDR, ".csv")))
                            
## 2-shared cell type-specific Overlaps
# overlaps_shared_2ct_df <- read.csv(here(inputCSV_Overlaps_Dir, 
#                              paste0("overlaps_summary_linkPeak_DARs_shared_2ct_FDR", FDR, ".csv")))

all_unique_df <- overlap_df |>
    distinct(peak_id_links, .keep_all = TRUE)
nrow(all_unique_df) # 5603
head(all_unique_df)

peaks_2n_cell_types <- shared_2ct_df |>
    filter(n_cell_types == 1 | n_cell_types == 2) |>
    pull(peak_id)
length(peaks_2n_cell_types) # 1967 / 3205

# Keep only those 1,967 unique LinkPeaks
unique_df <- all_unique_df |>
    filter(peak_id_links %in% peaks_2n_cell_types)

length(unique_df$peak_id_links) # 1967 / 3205
table(unique_df$cell_type)


## =============================================================================
## subset by cell_type and cluster, focus only on cis-links where both accessibility and expression specificity occur in the same cluster
# Conceptually: 
# cell_type → From the ATAC side: where the chromatin accessibility change occurs (DARs)
# cluster → From the RNA side: where the correlated gene was expressed or correlated in the LinkPeaks model


all_clusters <- unique(uniques_df$cell_type)
# sort by MHh and LHb first 
clusters_sorted <- c(
    sort(grep("MHb|LHb", all_clusters, value = TRUE)), 
    sort(grep("MHb|LHb", all_clusters, value = TRUE, invert = TRUE))
)
clusters_sorted

## subset df by cell-Type
subset_cell_type <- function(unique_df, cluster_specific) {
    message("Subsetting [", cluster_specific, "] cell_type")
    subset_uniques <- unique_df |> filter(cell_type == cluster_specific)
    message("Rows found: ", nrow(subset_uniques))
    message("Unique peaks: ", length(unique(subset_uniques$peak_id_links)))
    return(subset_uniques)
}

subsetted_df <- purrr::map_df(
    clusters_sorted,
    ~ subset_cell_type(uniques_df, .x)
)

length(subsetted_cellTypes_df)


## =============================================================================
## Generate scatterplot


make_scattered_plot_dars_cc_real <- function(
        subset_uniques,
        fdr_cutoff = 0.05
) {
    # FDR thresholds for LinkPeaks (CC) and DARs
    thr_CC <- 0.1 # magnitude threshold for correlation strength (|CCscore|)
    thr_DAR <- fdr_cutoff # threshold for accessibility significance
    clus_name <- unique(subset_uniques[["cell_type"]])
    
    message("Building plot for [", clus_name, "] (real FDR values)")
    message("Using FDR threshold =", thr_DAR)
    
    # Categorize points
    plot_data <- subset_uniques |>
        mutate(
            sig_CC  = abs(CCscore) > thr_CC,     # strong correlation in either direction
            sig_DAR = fdr_dars < thr_DAR,        # significant accessibility
            category = case_when(
                sig_CC & sig_DAR & logFC > 0  ~ "Active CRE (+)",
                sig_CC & sig_DAR & logFC < 0  ~ "Repressive CRE (-)",
                sig_CC & !sig_DAR              ~ "Shared CRE",
                !sig_CC & sig_DAR              ~ "Unlinked OCR",
                TRUE                           ~ "Non-significant"
            )
        )
    
    top_hits <- plot_data |>
        filter(category %in% c("Active CRE (+)", "Repressive CRE (-)")) |>
        arrange(FDR_CC) |>
        head(15)
    
    x_max <- min(1, max(plot_data$FDR_CC, na.rm = TRUE) * 1.05)
    y_max <- min(1, max(plot_data$fdr_dars, na.rm = TRUE) * 1.05)
    # define tick units
    axis_breaks <- seq(0, max(x_max, y_max), by = 0.05)
    axis_labels <- sprintf("%.2f", axis_breaks)
    
    g1 <- ggplot(plot_data, aes(x = CCscore, y = logFC, color = category)) +
        geom_hline(yintercept = 0, linetype = "solid", color = "grey70") +
        geom_vline(xintercept = 0, linetype = "solid", color = "grey70") +
        geom_vline(xintercept = c(-thr_CC, thr_CC), linetype = "dashed", color = "darkgrey") +
        geom_point(alpha = 0.7, size = 1.6) +
        geom_text_repel(
            data = top_hits,
            aes(label = gene_name),
            size = 3,
            max.overlaps = 15
        ) +
        scale_color_manual(values = c(
            "Active CRE (+)"  = "#E64B35FF",
            "Repressive CRE (-)" = "#4DBBD5FF",
            "Shared CRE"      = "#00A087FF",
            "Unlinked OCR"    = "#3C5488FF",
            "Non-significant" = "lightgrey"
        )) +
        labs(
            title = paste0(clus_name, " | CCscore vs logFC"),
            subtitle = paste0("Dashed lines: |CCscore| > ", thr_CC, 
                              " | FDR < ", thr_DAR),
            x = "LinkPeaks correlation (CCscore)",
            y = "log2 Fold Change (Accessibility)"
        ) +
        theme_minimal(base_size = 12) +
        theme(
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold"),
            legend.position = "bottom",
            legend.title = element_blank()
        )
    
    return(g1)
}

make_scattered_plot_dars_cc <- function(subset_uniques) {
    thr_CC <- 0.1
    thr_DAR <- 0.1
    thr_val <- thr_DAR
    clus_name <- unique(subset_uniques[["cell_type"]])
    
    message("Building plot for [", clus_name, "]")
    message("Using FDR threshold =", thr_val)
    
    plot_data <- subset_uniques |>
        mutate(
            neglog_FDR_CC = -log10(FDR_CC),
            neglog_fdr_dars = -log10(fdr_dars),
            sig_CC = FDR_CC < thr_CC,
            sig_DAR = fdr_dars < thr_DAR,
            category = case_when(
                sig_CC & sig_DAR & logFC > 0  ~ "Active CRE (+)",
                sig_CC & sig_DAR & logFC < 0  ~ "Repressive CRE (-)",
                sig_CC & !sig_DAR              ~ "Shared CRE",
                !sig_CC & sig_DAR              ~ "Unlinked OCR",
                TRUE                           ~ "Non-significant"
            )
        )
    
    top_hits <- plot_data |>
        filter(category %in% c("Active CRE (+)", "Repressive CRE (-)")) |>
        arrange(FDR_CC) |>
        head(20)
    
    g1 <- ggplot(plot_data, aes(x = neglog_FDR_CC, y = neglog_fdr_dars, color = category)) +
        geom_point(alpha = 0.7, size = 1.6) +
        geom_text_repel(
            data = top_hits,
            aes(label = gene_name),
            size = 3,
            max.overlaps = 15
        ) +
        geom_vline(xintercept = -log10(thr_CC), linetype = "dashed", color = "darkgrey") +
        geom_hline(yintercept = -log10(thr_DAR), linetype = "dashed", color = "darkgrey") +
        scale_color_manual(values = c(
            "Active CRE (+)" = "#E64B35FF",
            "Repressive CRE (-)" = "#4DBBD5FF",
            "Shared CRE" = "#00A087FF",
            "Unlinked OCR" = "#3C5488FF",
            "Non-significant" = "lightgrey"
        )) +
        +
        scale_x_continuous(limits = c(0, 0.5)) +
        scale_y_continuous(limits = c(0, 0.5)) +
        labs(
            title = "Peak–Gene Correlation vs. Differential Accessibility",
            subtitle = paste0(clus_name, " | FDR = ", thr_val),
            x = expression(-log[10](FDR[CC])),
            y = expression(-log[10](FDR[DARs]))
        ) +
        theme_minimal(base_size = 12)
    
    return(g1)
}

length(subsetted_df) # 18
table(subsetted_df$cell_type)
# Astrocyte Excit.Thal Inhib.Thal      LHb.1  LHb.1.3.4    LHb.2.7      LHb.4 
# 11        160        501         13         20         91        163 
# MHb.1    MHb.1.2      MHb.2  Microglia      Oligo 
# 7          2        258          1         13

head(subsetted_df)
subsetted_list_df <- split(subsetted_df, subsetted_df$cell_type)
length(subsetted_list_df)

# log10
scattered_plt_cell_type <- purrr::map(
    subsetted_list_df,
    ~ make_scattered_plot_dars_cc(.x)
)

scattered_plt_cell_type[1] 

pdf("ScatteredPlots_byCellType.pdf", width = 6, height = 5)
walk(scattered_plt_cell_type, print)
dev.off()

# real values
scattered_plt_cell_type_real_values <- purrr::map(
    subsetted_list_df,
    ~ make_scattered_plot_dars_cc_real(.x, FDR)
)
scattered_plt_cell_type_real_values[1]


## =============================================================================
## Classify links as within or cross: stacked bar (proportions and counts)

# Each row in links_summary is a peak–gene pair labeled according to whether
# the peak accessibility (cell_type) and correlated gene expression (cluster) belong to the same cluster.
links_summary <- uniques_df |>
    mutate(link_type = if_else(cell_type == cluster, "within_cluster", "cross_cluster"))
head(links_summary)

## Summarize counts per cell type
summary_counts <- links_summary |>
    count(cell_type, link_type) |>
    group_by(cell_type) |>
    mutate(
        total = sum(n),
        proportion = n / total
    )
summary_counts


## stacked bar (proportions) show, for each cell_type what fraction of its significant LinkPeaks 
## - connect to genes in the same cluster (green) vs other clusters

g1 <- ggplot(summary_counts, aes(x = cell_type, y = proportion, fill = link_type)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("within_cluster" = "#1b9e77", "cross_cluster" = "#d95f02")) +
    labs(
        title = "Within vs Cross-Cluster LinkPeaks-DARs overlaps by Cell Type",
        x = "Cell Type",
        y = "Proportion of Links",
        fill = "Link Type"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

f_name <- paste0("Cross-Cluster_proportions_FDR", FDR, ".pdf")
ggsave(here(plotDir, f_name),
       g1, width = 8, height = 8)


## stacked bar (absolute counts)
g1 <- ggplot(summary_counts, aes(x = cell_type, y = n, fill = link_type)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(y = "Number of Links") +
labs(
    title = "Within- vs Cross-Cluster LinkPeaks-DARs overlaps by Cell Type",
    x = "Cell Type",
    y = "Number of Links",
    fill = "Link Type"
) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

f_name <- paste0("Cross-Cluster_abs_counts_FDR", FDR, ".pdf")
ggsave(here(plotDir, f_name),
       g1, width = 8, height = 8)


## =============================================================================

## Highlighting reciprocal relationships
# matrix heatmap of regulatory interactions between cell-type pairs
# directionality and strength of inter-cluster links.

# heatmap h1:
# How many cross-cluster connections exist? 
# h1 - Excludes diagonal (cell_type != cluster), showing only cross-cluster links
# Shows how much cross-regulation occurs between clusters

pair_counts <- links_summary |>
    count(cell_type, cluster) |>
    filter(cell_type != cluster)
head(pair_counts)

f_name = here(plotDir, paste0("heatmap_CrossCluster_LinkCounts_h1_FDR", FDR,".pdf"))
pdf(f_name, width = 7, height = 6)
pheatmap::pheatmap(
    pivot_wider(pair_counts, names_from = cluster, values_from = n, values_fill = 0)[,-1],
    cluster_rows = TRUE, cluster_cols = TRUE,
    main = "Cross-Cluster Link Counts (peaks→genes)"
)
dev.off()

# ======

# heatmap h2:
# How strongly each cell type favors certain gene clusters?
# h2 Keeps all pairs, including within-cluster (diagonal)
# Highlights preferential link patterns (who regulates whom more strongly) across all clusters

pair_counts <- links_summary |>
    count(cell_type, cluster)
mat <- pair_counts |>
    pivot_wider(names_from = cluster, values_from = n, values_fill = 0)
mat_norm <- mat |>
    column_to_rownames("cell_type") |>
    as.matrix()
# normalize by row totals
mat_prop <- mat_norm / rowSums(mat_norm)
# # center/scale each row
mat_z <- t(scale(t(mat_norm)))  

# Y-axis (rows)	ATAC-defined cell_type	Where the open chromatin peak is located (e.g., LHb.4).
# X-axis (columns)	RNA-defined cluster	Where the correlated gene expression occurs (e.g., Excit.Thal).
# Color	Strength or proportion of links	Fraction (or Z-score) of LHb.4 peaks linked to genes in each expression cluster.
f_name = here(plotDir, paste0("heatmap_within_crossCluster_LinkCounts_h2_FDR", FDR,".pdf"))
pdf(f_name, width = 7, height = 6)
pheatmap::pheatmap(
    mat_z,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    main = "Z-scored Link Density (peaks→genes)"
)
dev.off()

    


message("All plots done!!!")





# library("slurmjobs")
# job_single(
#   "19_Linkage_DARs_analysis",
#   create_shell = TRUE,
#   partition = "katun",
#   memory = "30G",
#   cores = 2,
#   logdir = "logs",
#   command = "Rscript 19_Linkage_DARs_analysis.R",
#   create_logdir = FALSE
# )

## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()



