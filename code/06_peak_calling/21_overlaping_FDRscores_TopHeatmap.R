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
    "06_peak_calling",
    "20_explore_overlapping_linkpeaks_DARs"
)
processedDir <- here(
    "processed-data",
    "06_peak_calling",
    "21_explore_top_overlapping_linkpeaks_DARs"
)
plotDir <- here(
    "plots",
    "06_peak_calling",
    "21_explore_top_overlapping_linkpeaks_DARs"
)

if (!dir.exists(processedDir)) { dir.create(processedDir) }
if (!dir.exists(plotDir)) { dir.create(plotDir) }

##==============================================================================

message("Loading Unique-Overlap with directionality ...")

ct = "LHb.4"
top_g = 20

## load unique and shared df(s)
uniques_df <- read.csv(here(inputCSV_Overlaps_Dir, 
                            paste0("overlaps_unique_linkPeak_DARs_detail_FDR", FDR, ".csv")))
## 2-shared cell type-specific Overlaps
# overlaps_shared_2ct_df <- read.csv(here(inputCSV_Overlaps_Dir, 
#                              paste0("overlaps_summary_linkPeak_DARs_shared_2ct_FDR", FDR, ".csv")))

colnames(uniques_df)
head(uniques_df)
nrow(uniques_df) # 2904
length(unique(uniques_df$peak_id_links)) # [1] 1967

## verification, unique DARs
unique_summary_collapsed <- uniques_df |>
    group_by(cell_type, peak_id_links, accessibility) |> 
    summarise(.groups = "drop") # assign accessibility (More, Less, Mixed) based on all links

## verification
uniques_df |>
    count(peak_id_links, sort = TRUE) |>
    filter(n > 1)
# peak_id_links   n
# 1 chr1-2683852-2684422 108
# 2 chr1-161612222-161613071  96
# 3 chr11-65497297-65497968  91


## =============================================================================
## subset by cell_type and cluster, focus only on cis-links where both accessibility and expression specificity occur in the same cluster
# Conceptually: 
# cell_type → From the ATAC side: where the chromatin accessibility change occurs (DARs)
# cluster → From the RNA side: where the correlated gene was expressed or correlated in the LinkPeaks model

cluster_specific = "LHb.4"

subset_uniques <- uniques_df |> filter(cell_type == cluster_specific)
nrow(subset_uniques)
length(unique(subset_uniques$peak_id_links)) # 194
head(subset_uniques)
table(subset_uniques$cluster)

## focus only on cis-links where both accessibility and expression specificity occur in the same cluster
subset_uniques <- uniques_df |> 
    filter(cell_type == cluster)
head(subset_uniques)
subset_uniques <- subset_uniques |> filter(cell_type == cluster_specific)
head(subset_uniques)
nrow(subset_uniques)


## =============================================================================
## Generate scatterplot

# Define thresholds for visual reference
thr_CC <- 0.1      # significance for LinkPeaks correlation
thr_DAR <- 0.1     # significance for differential accessibility

# # Extract numeric value from FDR_threshold column, e.g. "FDR0.2" → 0.2
# thr_val <- uniques_df$FDR_threshold[1] |>
#     gsub("FDR", "", x = _) |>
#     as.numeric()
thr_val <- thr_DAR

thr_CC <- thr_val
thr_DAR <- thr_val

cat("Using FDR threshold =", thr_val, "\n")

plot_data <- uniques_df %>%
    mutate(
        neglog_FDR_CC = -log10(FDR_CC),
        neglog_fdr_dars = -log10(fdr_dars),
        sig_CC = FDR_CC < thr_CC,
        sig_DAR = fdr_dars < thr_DAR,
        sig_both = sig_CC & sig_DAR, 
        category = case_when(
            sig_CC & sig_DAR & logFC > 0  ~ "Active CRE (+)",
            sig_CC & sig_DAR & logFC < 0  ~ "Repressive CRE (−)",
            sig_CC & !sig_DAR              ~ "Shared CRE",
            !sig_CC & sig_DAR              ~ "Unlinked OCR",
            TRUE                           ~ "Non-significant"
        )
    )

top_hits <- plot_data %>%
    filter(category %in% c("Active CRE (+)", "Repressive CRE (−)")) %>%
    arrange(FDR_CC) %>%
    head(20)

ggplot(plot_data, aes(x = neglog_FDR_CC, y = neglog_fdr_dars, color = category)) +
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
        "Repressive CRE (−)" = "#4DBBD5FF",
        "Shared CRE" = "#00A087FF",
        "Unlinked OCR" = "#3C5488FF",
        "Non-significant" = "lightgrey"
    )) +
    labs(
        title = "Peak–Gene Correlation vs. Differential Accessibility",
        subtitle = paste0("Dashed lines: FDR = ", thr_val),
        x = expression(-log[10](FDR[CC])),
        y = expression(-log[10](FDR[DARs]))
    ) +
    theme_minimal(base_size = 12)


# ggplot(plot_data, aes(x = neglog_FDR_CC, y = neglog_fdr_dars)) +
#     geom_point(aes(color = sig_both), alpha = 0.7, size = 1.5) +
#     geom_vline(xintercept = -log10(thr_CC), linetype = "dashed", color = "darkgrey") +
#     geom_hline(yintercept = -log10(thr_DAR), linetype = "dashed", color = "darkgrey") +
#     facet_wrap(~ cell_type, scales = "free_y") +
#     scale_color_manual(values = c("TRUE" = "#E64B35FF", "FALSE" = "lightgrey"),
#                        name = paste0("FDR < ", thr_val, " in both")) +
#     labs(
#         title = "LinkPeaks vs DARs Significance per Cell Type",
#         x = expression(-log[10](FDR[CC])),
#         y = expression(-log[10](FDR[DARs]))
#     ) +
#     theme_minimal(base_size = 11) +
#     theme(strip.text = element_text(face = "bold"))



## Classify links as within or cross
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

ggplot(summary_counts, aes(x = cell_type, y = proportion, fill = link_type)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("within_cluster" = "#1b9e77", "cross_cluster" = "#d95f02")) +
    labs(
        title = "Within- vs Cross-Cluster LinkPeaks-DARs overlaps by Cell Type",
        x = "Cell Type",
        y = "Proportion of Links",
        fill = "Link Type"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

## stacked bar (absolute counts)
ggplot(summary_counts, aes(x = cell_type, y = n, fill = link_type)) +
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



