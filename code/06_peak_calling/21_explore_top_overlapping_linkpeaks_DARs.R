########################################################################
## Overlapping Scores: 
## Distribution of differential accessibility (logFC) values for genes whose linked peaks overlap with DARs
##
## Authors. CSC
## Date. Oct 06, 2025
## Recommended resources on interactive mode: srun --pty --mem=15GB --x11 bash
########################################################################

library("dplyr")
library("purrr")
library("ggplot2")
library("patchwork")
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

## load full overlaps df
overlaps_df <- read.csv(here(inputCSV_Overlaps_Dir, 
                             paste0("overlaps_unique_with_direction_ct_FDR", FDR, ".csv")))
## test
colnames(overlaps_df)
head(overlaps_df)
nrow(overlaps_df) # 3201
ct = "LHb.4|MHb.2"
top_g = 20


make_correlations_logFC_CCScore <- function(
        ct,
        top_g,
        metric # "logFC" or "FDR_CC"
) {
    
    # Identify all Habenula (Hb) cell types
    hb_ct <- unique_overlaps |>
        distinct(cell_type) |>
        pull(cell_type) |>
        grep(pattern = ct, value = TRUE)
    
    # Filter data for Hb cell types and calculate magnitude
    mhb_lhb_data <- unique_overlaps |>
        filter(cell_type %in% hb_ct) |>
        mutate(
            abs_logFC = abs(logFC),
            xmetric_val = if (metric == "logFC") abs_logFC else FDR_CC
        )
    
    # Select top genes per cell type
    if (metric == "logFC") {
        top_n_per_ct <- mhb_lhb_data |>
            group_by(cell_type, gene_name) |>
            slice_max(order_by = xmetric_val, n = 1, with_ties = FALSE) |>
            ungroup() |>
            group_by(cell_type) |>
            slice_max(order_by = xmetric_val, n = top_g) |>
            ungroup()
    } else {
        top_n_per_ct <- mhb_lhb_data |>
            group_by(cell_type, gene_name) |>
            slice_min(order_by = xmetric_val, n = 1, with_ties = FALSE) |>
            ungroup() |>
            group_by(cell_type) |>
            slice_min(order_by = xmetric_val, n = top_g) |>
            ungroup()
    }
    
    # Prepare plot data
    plot_data <- mhb_lhb_data |>
        inner_join(top_n_per_ct |> select(cell_type, gene_name),
                   by = c("cell_type", "gene_name")) |>
        group_by(cell_type) |>
        mutate(gene_name = fct_reorder(gene_name, xmetric_val, .fun = max, .desc = FALSE)) |>
        ungroup()
    
    # Transform FDR_CC scale for clarity
    if (metric == "FDR_CC") {
        plot_data <- plot_data |> 
            mutate(FDR_CC = -log10(pmax(FDR_CC, 1e-300)))  # avoid Inf for FDR=0
        y_label <- expression(-log[10]~"(FDR_CC)")
    } else {
        y_label <- "Log Fold Change (logFC)"
    }
    
    sub_label <- if (metric == "logFC") 
        "ranked by maximum |logFC|" 
    else 
        "ranked by lowest FDR_CC (most significant links)"
    
    # Generate violin plot
    g_violin_faceted <- ggplot(
        plot_data,
        aes(x = gene_name, y = !!sym(metric), fill = direction)
    ) +
        geom_violin(trim = FALSE, scale = "width", alpha = 0.7) +
        geom_jitter(aes(color = direction),
                    height = 0, width = 0.1, size = 1.5, alpha = 0.6) +
        coord_flip() +
        facet_wrap(~cell_type, scales = "free_y", ncol = 1) +
        scale_fill_manual(values = c("Down" = "#1F77B4", "Up" = "#D62728")) +
        scale_color_manual(values = c("Down" = "black", "Up" = "black"), guide = "none") +
        geom_hline(
            yintercept = if (metric == "logFC") 0 else -log10(0.05),
            linetype = "dashed",
            color = "grey50"
        ) +
        labs(
            title = "Top unique Link–DARs",
            subtitle = sub_label,
            x = "Gene Name",
            y = y_label,
            fill = "Direction"
        ) +
        theme_minimal() +
        theme(
            legend.position = "bottom",
            strip.background = element_rect(fill = "grey90", color = "grey50"),
            strip.text = element_text(face = "bold"),
            axis.title.y = element_text(size = 10),
            axis.title.x = element_text(size = 10)
        )
    
    return(g_violin_faceted)
}

# Generate plots
g1 <- make_correlations_logFC_CCScore(ct, top_g, "logFC")
g2 <- make_correlations_logFC_CCScore(ct, top_g, "FDR_CC")

g_combined <- (g1 + g2) + plot_annotation(tag_levels = 'A')

ggsave(here(plotDir, "top_unique_overlaps_ranked_by_logFC_FDR_CC.pdf"),
       g_combined, width = 8, height = 8)


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



