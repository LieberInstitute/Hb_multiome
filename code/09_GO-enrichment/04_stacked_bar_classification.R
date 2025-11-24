
library(tidyverse)
library(dplyr)
library(ggplot2)
library(forcats)
library(here)

## =============================================================================
## setup variables and file names
FDR = 0.2
inputCSV_Overlaps_Dir <- here(
    "processed-data",
    "06_peak_calling"
)

message("Loading Overlap ...")

## load classified overlaps with 2-shared ct
overlap_df <- read.csv(here(inputCSV_Overlaps_Dir, "21_overlaping_FDRscores_TopHeatmap", 
                            "overlaps_linkPeak_DARs_classified_thr_CC0.3_thr_DAR0.1.csv"))

nrow(overlap_df) #3205
colnames(overlap_df)

## =============================================================================
## Classify Regulatory Element Types
## Uses strand-aware distance (in kb) to assign RE categories

overlap_df <- overlap_df |>
    mutate(
        RE_type = case_when(
            abs(distance_kb) <= 2 ~ "Promoter (±2 kb)",
            abs(distance_kb) <= 10 ~ "Near-Promoter (2–10 kb)",
            abs(distance_kb) <= 100 ~ "Enhancer (10–100 kb)",
            abs(distance_kb) <= 500 ~ "Long-range (100–500 kb)",
            TRUE ~ "Distal (>500 kb)"
        ),
        RE_type = factor(
            RE_type,
            levels = c(
                "Promoter (±2 kb)",
                "Near-Promoter (2–10 kb)",
                "Enhancer (10–100 kb)",
                "Long-range (100–500 kb)",
                "Distal (>500 kb)"
            )
        )
    )
head(overlap_df)

re_colors <- c(
    "Promoter (±2 kb)"          = "#E69F00",
    "Near-Promoter (2–10 kb)"   = "#56B4E9",
    "Enhancer (10–100 kb)"      = "#009E73",
    "Long-range (100–500 kb)"   = "#CC79A7",
    "Distal (>500 kb)"          = "grey70"
)


## =============================================================================

## summarize counts per category × RE_type
count_plot_df_category <- function(
        overlap_df,
        RE_type
    ) {
    plot_df_category <- overlap_df |>
        count(category, RE_type) |>
        group_by(category) |>
        mutate(percentage = 100 * n / sum(n)) |>
        ungroup()
    
    ## order categories by total number of peaks
    plot_df_category$category <- fct_reorder(plot_df_category$category, plot_df_category$n, .fun = sum)
    return(plot_df_category)
}


## =============================================================================

## stacked bar plot for category
plot_df_category <- count_plot_df_category(overlap_df, RE_type)

p1 <- ggplot(plot_df_category, aes(x = category, y = percentage, fill = RE_type)) +
    geom_bar(stat = "identity", width = 0.75, color = "black", linewidth = 0.2) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
    scale_fill_manual(values = re_colors) +
    labs(
        title = "Regulatory Element Composition per Category",
        x = "Peak Category",
        y = "Percentage (%)",
        fill = "Regulatory Element Type"
    ) +
    theme_bw(base_size = 12) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold")
    )
p1


## =============================================================================

# Recompute stacked cumulative percentages to make a pecentiage composition plot, due some RE categories have percentage = 0
plot_wide <- plot_df_category |>
    arrange(category, RE_type) |>
    group_by(category) |>
    mutate(cum_low  = cumsum(percentage) - percentage,
           cum_high = cumsum(percentage)) |>
    ungroup() |> 
    mutate(idx = as.numeric(category))

## Ribbon plot across categories 
p2 <- ggplot(plot_wide, aes(x = idx)) +
    geom_ribbon(aes(ymin = cum_low, ymax = cum_high, fill = RE_type),
                color = "black", linewidth = 0.15, alpha = 0.95) +
    scale_x_continuous(
        breaks = unique(plot_wide$idx),
        labels = levels(plot_df_category$category)   # display category names
    ) +
    scale_y_continuous(
        labels = scales::percent_format(scale = 1),
        expand = c(0, 0), limits = c(0, 100)
    ) +
    scale_fill_manual(values = re_colors) +
    theme_minimal(base_size = 13) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()
    ) +
    labs(
        y = "% Composition",
        title = "Regulatory Element Composition Across Categories"
    )
p2


library(sessioninfo)
session_info()