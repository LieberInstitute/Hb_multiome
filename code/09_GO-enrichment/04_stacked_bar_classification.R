
library(tidyverse)
library(dplyr)
library(ggplot2)
library(forcats)
library(here)

## setup variables and file names
FDR = 0.2
inputCSV_Overlaps_Dir <- here(
    "processed-data",
    "06_peak_calling"
)

message("Loading Overlap ...")

## load classified overlaps with 2-shared ct
overlap_df2 <- read.csv(here(inputCSV_Overlaps_Dir, "21_overlaping_FDRscores_TopHeatmap", 
                            "overlaps_linkPeak_DARs_classified_thr_CC0.3_thr_DAR0.1.csv"))

nrow(overlap_df2) #3205
colnames(overlap_df2)
# [1] "source_df"           "peak_id_links"       "CCscore"            
# [4] "gene_name"           "gene_id"             "FDR_CC"             
# [7] "cluster"             "tss"                 "gene_strand"        
# [10] "distance"            "distance_kb"         "signed_distance"    
# [13] "signed_by_strand"    "peak_id"             "cell_type"          
# [16] "FDR_threshold"       "logFC"               "fdr_dars"           
# [19] "overlap_type"        "n_cell_types"        "sig_CC"             
# [22] "sig_DAR"             "category"            "type_classification"

overlap_df2 <- overlap_df2 |>
    mutate(RE_type = case_when(
        abs(distance_kb) <= 2 ~ "Promoter",
        abs(distance_kb) <= 10 ~ "Near-Promoter",
        abs(distance_kb) <= 100 ~ "Enhancer",
        abs(distance_kb) <= 500 ~ "Long-range",
        TRUE ~ "Distal >500kb"
    ))
head(overlap_df2)


plot_df <- overlap_df2 |>
    count(category, RE_type) |>
    group_by(category) |>
    mutate(pct = n / sum(n) * 100)

head(plot_df)
# category                RE_type           n   pct
# <chr>                   <chr>         <int> <dbl>
#     1 Discordant Linked DAR   Enhancer    221 35.8 
# 2 Discordant Linked DAR   Long-range      365 59.2 
# 3 Discordant Linked DAR   Near-Promoter    20  3.24
# 4 Discordant Linked DAR   Promoter         11  1.78
# 5 Linked OCR (+) enriched Enhancer        105 37.8 
# 6 Linked OCR (+) enriched Long-range      155 55.8 


## Classify Regulatory Element Types
## Uses strand-aware distance (in kb) to assign RE categories

overlap_df2 <- overlap_df2 |>
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

## summarize counts per category × RE_type
plot_df <- overlap_df2 |>
    count(category, RE_type) |>
    group_by(category) |>
    mutate(percentage = 100 * n / sum(n)) |>
    ungroup()

## order categories by total number of peaks
plot_df$category <- fct_reorder(plot_df$category, plot_df$n, .fun = sum)

## likely we woild like to order categories by % Enhancer
# plot_df <- plot_df |>
#     group_by(category) |>
#     mutate(order_metric = ifelse(RE_type == "Enhancer", percentage, NA)) |>
#     fill(order_metric, .direction = "downup") |>  # fill NA with available value
#     ungroup() |>
#     arrange(order_metric)

plot_df$category <- factor(plot_df$category,
                           levels = unique(plot_df$category))

## convert to wide cumulative for area stacking
plot_wide <- plot_df |>
    arrange(category, RE_type) |>
    group_by(category) |>
    mutate(cum_pct = cumsum(percentage)) |>
    ungroup()
dim(plot_wide)
head(plot_wide)
table(plot_wide$RE_type)
table(plot_wide$category)
summary(plot_wide$cum_pct)


# stacked bar plot for category
ggplot(plot_df, aes(x = category, y = percentage, fill = RE_type)) +
    geom_bar(stat = "identity", width = 0.75, color = "black", linewidth = 0.2) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
        scale_fill_brewer(palette = "RdYlBu") +
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

## Maybe is better a variance components style plot 



library(sessioninfo)
session_info()