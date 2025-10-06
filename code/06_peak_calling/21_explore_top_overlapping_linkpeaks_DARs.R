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

# Identify all Habenula (Hb) cell types
hb_ct <- unique_overlaps |>
    distinct(cell_type) |>
    pull(cell_type) |>
    grep(pattern = ct, value = TRUE)
    #grep(pattern = "MHb|LHb", value = TRUE)

# Filter data for Hb cell types and calculate magnitude
mhb_lhb_data <- unique_overlaps |>
    filter(cell_type %in% hb_ct) |>
    mutate(abs_logFC = abs(logFC))

# Identify the Top 10 LinkPeaks/Genes PER CELL TYPE by Magnitude
top_n_per_ct <- mhb_lhb_data |>
    # Group by cell type AND gene to get the single max magnitude for each gene/cell-type combo
    group_by(cell_type, gene_name) |>
    slice_max(order_by = abs_logFC, n = 1, with_ties = FALSE) |>
    ungroup() |>
    # Now group only by cell_type to select the top 10 genes within each cell type
    group_by(cell_type) |>
    slice_max(order_by = abs_logFC, n = top_g) |>
    ungroup()

# Filter the original data using the list of (cell_type, gene_name) pairs and prepare for plotting
plot_data <- mhb_lhb_data |>
    # Inner join to filter to only the top 10 genes in each cell type
    inner_join(
        top_n_per_ct |> select(cell_type, gene_name),
        by = c("cell_type", "gene_name")
    ) |>
    # Order genes within each cell type by magnitude for cleaner visualization
    group_by(cell_type) |>
    mutate(gene_name = fct_reorder(gene_name, abs_logFC, .fun = max, .desc = FALSE)) |>
    ungroup()


# Generate the Faceted Violin Plot
g_violin_faceted <- ggplot(
    plot_data,
    aes(x = gene_name, y = logFC, fill = direction)
) +
    geom_violin(
        trim = FALSE,
        scale = "width",
        alpha = 0.7
    ) +
    geom_jitter(
        aes(color = direction),
        height = 0,
        width = 0.1,
        size = 1.5,
        alpha = 0.6
    ) +
    coord_flip() +
    facet_wrap(~cell_type, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c("Down" = "#1F77B4", "Up" = "#D62728")) +
    scale_color_manual(values = c("Down" = "black", "Up" = "black"), guide = "none") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    labs(
        #title = "Distribution of LogFC for Top 10 Habenula Unique Links",
        subtitle = "Top Unique Link-DARs per Hb Sub-type\nranked by maximum |LogFC|",
        x = "Gene Name",
        y = "Log Fold Change (LogFC)",
        fill = "Direction"
    ) +
    theme_minimal() +
    theme(
        legend.position = "bottom",
        #plot.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "grey50"),
        strip.text = element_text(face = "bold"),
        axis.title.y = element_text(size = 10), 
        axis.title.x = element_text(size = 10)
    )

print(g_violin_faceted)

ggsave(here(plotDir, paste0("overlaps_distribution_top", top_g, "_logFC.pdf")),
       g_violin_faceted, width = 4, height = 8)



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



