########################################################################
## Summaries/Prepare plots for interpret Peaks Overlapping Links-DARs
##
## Authors. CSC
## Date. Oct 02, 2025
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
    "19_Linkage_DARs_analysis"
)
processedDir <- here(
    "processed-data",
    "06_peak_calling",
    "20_explore_overlapping_linkpeaks_DARs"
)
plotDir <- here(
    "plots",
    "06_peak_calling",
    "20_explore_overlapping_linkpeaks_DARs"
)

if (!dir.exists(processedDir)) { dir.create(processedDir) }
if (!dir.exists(plotDir)) { dir.create(plotDir) }

##==============================================================================

message("Loading Unique-Overlap Hits ...")

# ## load summary overlaps df
# link_summary <- read.csv(here(inputCSV_Overlaps_Dir, 
#                               paste0("summary_LinkPeak_overlap_stats_FDR", FDR, ".csv")))

## load full overlaps df
overlaps_df <- read.csv(here(inputCSV_Overlaps_Dir, 
                             paste0("Overlaps_LinkPeak_DARs_FDR", FDR, ".csv")))
head(overlaps_df)
# colnames(overlaps_df)
# overlaps_df |>
#     filter(peak_id == "chr1-1745993-1746550") |>
#     select(peak_id, cell_type, FDR_CC, fdr_dars) |>
#     distinct()
# length(unique(overlaps_df$peak_id)) # 5535


##==============================================================================
## stacked bar plot showing Unique vs. Shared overlaps per cell type

# Prepare classification
link_summary_classification <- overlaps_df |>
    group_by(peak_id) |>
    summarise(
        n_cell_types = n_distinct(cell_type),
        .groups = "drop"
    ) |>
    mutate(
        overlap_type = ifelse(n_cell_types == 1, "Unique", "Shared")
    ) |>
    select(peak_id, overlap_type) # Keep only peak_id and classification
head(link_summary_classification)


## Join classification back to cell-type level and count distinct LinkPeaks --> save result

unique_shared_overlaps <- overlaps_df |>
    # Ensure each LinkPeak is counted once per cell type
    select(peak_id, cell_type) |>
    distinct() |>
    # Merge with the classification calculated above
    left_join(link_summary_classification, by = "peak_id") 
nrow(unique_shared_overlaps)
# 14611

# save detailed unique and shared overlaps
f_name <- here(processedDir, paste0("overlaps_detailed_unique_shared_FDR", FDR, ".csv"))
write.csv(unique_shared_overlaps, f_name, row.names = FALSE)


# Join classification back to cell-type level and count distinct LinkPeaks
plot_df_ct <- unique_shared_overlaps |>
    group_by(cell_type, overlap_type) |>
    summarise(
        n_links = n(),
        .groups = "drop"
    )
head(plot_df_ct)

# Calculate total counts for each overlap type
total_counts <- plot_df_ct |>
    group_by(overlap_type) |>
    summarise(total_n = sum(n_links))
# Create named vectors for colors and labels
overlap_colors <- c("Shared" = "grey", "Unique" = "#D62728")
# Define the custom labels using the calculated totals
overlap_labels <- c(
    "Shared" = paste0("Shared (N=", total_counts$total_n[total_counts$overlap_type == "Shared"], ")"),
    "Unique" = paste0("Unique (N=", total_counts$total_n[total_counts$overlap_type == "Unique"], ")")
)

g_cell_type_overlap <- plot_df_ct |>
    ggplot(aes(
        # Order cell types by the total number of links (sum of Unique + Shared)
        x = fct_reorder(cell_type, n_links, .fun = sum), 
        y = n_links, 
        fill = overlap_type
    )) +
    # Use geom_col (or geom_bar(stat="identity")) for counts
    geom_col(position = position_stack(reverse = TRUE), color = "black") +
    geom_text(
        aes(label = n_links), 
        position = position_stack(vjust = 0.5, reverse = TRUE), # Center the labels
        size = 4
    ) +
    labs(
        title = "Unique vs. Shared Overlaps by Cell Type",
        subtitle = paste0("FDR thr = ", FDR),
        x = NULL,
        y = "Number of Overlaps",
        fill = "Overlap Type"
    ) +
    scale_fill_manual(
        values = overlap_colors, 
        labels = overlap_labels
    ) +
    coord_flip() + # Flip for readability
    theme_minimal() +
    theme( # legend right / bottom 
        legend.position.inside = c(0.95, 0.05), 
        legend.justification = c("right", "bottom"), 
        legend.background = element_rect(colour = "gray80", fill = "white") 
    )

f_name = paste0("overlaps_unique_replicated_cellType_FDR", FDR, ".pdf")
ggsave(here(plotDir, f_name),
       g_cell_type_overlap, width = 8, height = 8)


##==============================================================================
## Summarize unique LinkPeak with n_DARs and n_cell_types

link_summary <- overlaps_df |>
    group_by(peak_id) |>
    summarise(
        n_DARs       = n(),                         # total DAR rows linked to this LinkPeak
        n_cell_types = n_distinct(cell_type),       # distinct od distinc cell types represented
        cell_types   = paste(unique(cell_type), collapse = "; "), # list cell types
        .groups = "drop"
    ) |>
    # 2138 = LinkPeaks overlapping DARs from exactly one cell type (new definition, consistent with multiplicity)
    mutate(overlap_type = ifelse(n_cell_types == 1, "Unique", "Shared"))
    # 1572 = LinkPeaks overlapping exactly one DAR (old definition)

table(link_summary$overlap_type)
# Shared Unique 
# 3636   1967
#table(link_summary$n_cell_types==2)

f_name <- here(processedDir, paste0("overlaps_summary_LinkPeak-DARs_FDR", FDR, ".csv"))
write.csv(link_summary, f_name, row.names = FALSE)

message("Summary overlaps done!")


##==============================================================================
## make combined plot with: "LinkPeak Overlap: Unique vs. Multi-Cell Type DARs" and "Multiplicity of Cell Type Overlap"

message("Plotting exploratory LinkPeaks overlap  ... ")

## Plot "Count of Unique vs. Replicated Overlaps"
# Counts from the overlap_type column
g_overlap_type <- link_summary |>
    ggplot(aes(x = overlap_type, fill = overlap_type)) +
    geom_bar(color = "black") +
    geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 3) +
    scale_fill_manual(
        values = c("Unique" = "#D62728", "Shared" = "grey") 
    ) +
    labs(
        #title = "LinkPeak Overlap: Unique vs. Multi-Cell Type DARs",
        x = "Classification",
        y = NULL,
        fill = "Overlap Type"
    ) +
    theme_minimal() +
    theme(legend.position = "none")


## Distribution of Cell Type Multiplicity

g_multiplicity <- link_summary |>
    mutate(n_cell_types_factor = factor(n_cell_types)) |>
    # create the color classification column, I want first red
    mutate(
        bar_color = ifelse(n_cell_types == 1, "#D62728", "grey")
    ) |>
    ggplot(aes(x = n_cell_types_factor, fill = bar_color)) +
    geom_bar(color = "black") + 
    # use the actual values in the 'bar_color' column as colors
    scale_fill_identity() + 
    geom_text(
        stat = "count", 
        aes(label = after_stat(count)), 
        vjust = -0.5, 
        size = 4
    ) +
    labs(
        #title = "Multiplicity of Cell Type Overlap (DARs per LinkPeak)",
        x = "Number of Distinct Cell-Types Overlapping",
        y = NULL
    ) +
    scale_y_continuous(limits = c(0, 4000)) +
    theme_minimal() +
    theme(legend.position = "none") +
    # textual annotation for simple display only
    plot_annotation(
        caption = paste(
            "Note: The n=1 bar is defined by a single DAR row match (n_DARs=1) for consistency with other analyses.",
            "Key insight: A large portion of overlapping LinkPeaks are differentially accessible in multiple cell types (n > 1).",
            sep = "\n"
        )
    )

# Adjust widths (1, 2) and keep a single y-axis title
combined_plot <- (g_overlap_type + g_multiplicity) + plot_layout(widths = c(1, 3)) + 
    plot_layout(axes = "collect_y") & 
    plot_annotation(
        title = "LinkPeak-DARs Overlaps",
        subtitle = NULL,
        caption = NULL,
        theme = theme(
            plot.margin = margin(5, 5, 5, 5)
        )
    ) &
    labs(y = "Number of LinkPeaks-DARs")

f_name <- paste0("overlaps_combined_unique_distinct_shared_ct_FDR", FDR, ".pdf")
ggsave(here(plotDir, f_name),
       combined_plot, width = 7, height = 7)



##==============================================================================
## make combined plot with: g_multiplicity+g_cell_type_overlap

combined_unique_dup_plots <- (g_cell_type_overlap + g_multiplicity) + plot_layout(widths = c(1, 2)) + 
    plot_layout(axes = "collect_y") & 
    plot_annotation(
        #title = "Duplicted Overlaps",
        subtitle = NULL,
        caption = NULL,
        theme = theme(
            plot.margin = margin(5, 5, 5, 5)
        )
    ) &
    labs(y = "Number of Overlaps")

f_name <- paste0("overlaps_combined_unique_shared_ct_detail_FDR", FDR, ".pdf")

ggsave(here(plotDir, f_name),
       combined_unique_dup_plots, width = 10, height = 8)



##==============================================================================
## Make "Distribution of Cell Type-Specific LinkPeaks (Unique Overlaps)" plot

message("Creating barplot showing the frequency of each unique cell type .." )

# Filter the link_summary dataframe to include only "Unique" overlaps
uniques_df <- link_summary |>
    filter(overlap_type == "Unique") |>
    arrange(desc(n_DARs))
nrow(uniques_df) # 2138

g_uniques <- uniques_df |>
    #ggplot(aes(x = fct_infreq(cell_types), fill = cell_types)) + # fct_infreq orders bars by count
    ggplot(aes(
        # Use .desc = FALSE (Ascending count) to achieve largest bar at the top of the flipped plot.
        x = fct_reorder(cell_types, cell_types, .fun = length, .desc = FALSE), 
        fill = cell_types
    )) + 
    geom_bar(color = "black",  fill = "#D62728") +
    geom_text(
        stat = "count", 
        aes(label = after_stat(count)), 
        hjust = -0.5, 
        size = 4
    ) +
    coord_flip() + # Flip coordinates for readable cell type labels
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + # add 10% space on right
    labs(
        title = "Distribution of Cell Type-Specific Unique LinkPeaks-DARs",
        subtitle = paste0("Total Overlaps: ", nrow(uniques_df), " | FDR thr = ", FDR),
        x = NULL,
        y = "Number of Overlaps"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "none" # Remove legend since fill is redundant with the y-axis
    )

f_name <- paste0("overlaps_linkPeak_unique_ct_specific_", FDR, ".pdf")
ggsave(here(plotDir, f_name),
       g_uniques, width = 7, height = 8)


##==============================================================================
# Calculate and plot top 3 most frequent cell type combinations for each multiplicity level. I have 10 levels 
# plot the results for the first 6 multiplicity levels (e.g., n=1 through n=6).

message("Calculate the top 3 most frequent cell type combinations for each multiplicity level")

# aggregate your link_summary data to find the top 3 most common cell type combinations 
# for each multiplicity level (n_cell_types)

# Calculate up to 3 rows per multiplicity level ( I have set 10 level )
head(link_summary)
top_combinations_key <- link_summary |>
    # Group by the multiplicity level (the bar) and the specific cell type combination
    group_by(n_cell_types, cell_types) |>
    # Count how many unique LinkPeaks share that exact combination
    summarise(n_peaks = n(), .groups = 'drop_last') |>
    # Get the top 3 most frequent combinations
    slice_max(order_by = n_peaks, n = 3) |>
    mutate(
        display_text = paste0(cell_types, " (n=", n_peaks, ")"),
        `Multiplicity (n)` = n_cell_types
    ) |>
    ungroup() |>
    arrange(`Multiplicity (n)`, desc(n_peaks))

nrow(top_combinations_key)
head(top_combinations_key)

f_name <- here(processedDir, paste0("overlaps_summary_linkPeak_multiplicity_top3_freq_ct_", FDR, ".csv"))
write.csv(link_summary, f_name, row.names = FALSE)

message("Summary overlaps done!")

# make combined plot
g_top_combinations <- ggplot(
    #top_combinations_key,
    # Filter to include only the first 6 multiplicity levels (n=1 to n=5)
    top_combinations_key |> filter(`Multiplicity (n)` <= 5),
    aes(x = reorder(display_text, n_peaks), y = n_peaks, fill = factor(`Multiplicity (n)`))
) +
    geom_col(show.legend = FALSE) +
    coord_flip() +
    facet_wrap(~`Multiplicity (n)`, scales = "free_y", ncol = 1) +
    labs(
        title = "Top 3 Cell-Type Combinations\nper Multiplicity Level",
        subtitle = "first 6 multiplicity levels",
        x = "Cell-Type Combination",
        y = "Number of LinkPeaks Overlaps"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0, face = "bold"), # left-aligned title
        strip.text = element_text(face = "bold"),
        axis.text.y = element_text(size = 9)                 # smaller y-axis labels
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

f_name <- paste0("overlaps_linkPeak_5multiplicity_top3_frequent_ct_", FDR, ".pdf")
ggsave(here(plotDir, f_name),
       g_top_combinations, width = 8, height = 8)



#===========================================================================
## Process and Plot Up/Down LinkPekas-DARs overlapping per cell type

message("Processing Up/Down LinkPekas-DARs overlapping per cell type ... ")

colnames(overlaps_df)
head(overlaps_df)

# Define direction for all the overlaps: Up / Down
# peak_id_links → the LinkPeak ID (from gr_links)
# peak_id → the DAR peak ID (from gr_dars)

overlaps_df <- overlaps_df |>
    mutate(direction = case_when(
        logFC > 0 ~ "Up",
        logFC < 0 ~ "Down",
        TRUE ~ "Neutral"
    )) |>
    filter(direction %in% c("Up", "Down"))

nrow(overlaps_df) # [1] 20650
n_uniques <- nrow(uniques_df) # 2138
n_uniques

# Restrict overlaps to Unique LinkPeaks only
unique_overlaps <- overlaps_df |>
    filter(peak_id_links %in% (uniques_df |> pull(peak_id)))
head(unique_overlaps)
nrow(unique_overlaps)

f_name <- here(processedDir, paste0("overlaps_unique_with_direction_ct_FDR", FDR, ".csv"))
write.csv(unique_overlaps, f_name, row.names = FALSE)

message("Saved unique overlaps!")



## ===== Plot 1: Peak-level (collapse to one direction per LinkPeak) =====

unique_summary_collapsed <- unique_overlaps |>
    group_by(cell_type, peak_id_links) |>
    summarise(
        # direction = ifelse(mean(logFC) > 0, "Up", "Down"),
        direction = case_when(
            all(logFC > 0) ~ "Up",
            all(logFC < 0) ~ "Down",
            TRUE ~ "Mixed"   # at least one Up and one Down DAR
        ),
        .groups = "drop"
    ) |>
    count(cell_type, direction) |>
    group_by(cell_type) |>
    mutate(total = sum(n)) |>
    ungroup() |>
    mutate(
        #n_signed = ifelse(direction == "Down", -n, n)
        n_signed = case_when(
            direction == "Down"  ~ -n,
            direction == "Up"    ~  n,
            direction == "Mixed" ~  n   # show Mixed as positive
        )
    )
head(unique_summary_collapsed)
# cell_type  direction     n total n_signed
# <chr>      <chr>     <int> <int>    <int>
# 1 Astrocyte  Down         45    66      -45
# 2 Astrocyte  Up           21    66       21
# 3 Endo       Up            1     1        1
# 4 Excit.Thal Down        127   301     -127
table(unique_summary_collapsed$direction)
# Down   Up 
# 13   15 

cell_order1 <- unique_summary_collapsed |>
    group_by(cell_type) |>
    summarise(total = sum(abs(n_signed)), .groups = "drop") |>
    arrange(total) |>
    pull(cell_type)

unique_summary_collapsed$cell_type <- factor(unique_summary_collapsed$cell_type, levels = cell_order1)

## Plot with Mixed in the center
g_collapsed <- ggplot(unique_summary_collapsed,
                      aes(x = cell_type, y = n_signed, fill = direction)) +
    geom_col() +
    geom_text(
        aes(label = n,
            hjust = ifelse(direction == "Down", 1.1, 
                           ifelse(direction == "Up", -0.1, 0.5))),
        size = 3, color = "black"
    ) +
    scale_y_continuous(labels = abs, expand = expansion(mult = c(0.15, 0.15))) +
    labs(
        title = "Up/Down/Mixed Link-DARs",
        subtitle = paste0("Overlaps (Unique LinkPeaks - ", n_uniques ,")"),
        x = "Cell Type",
        y = NULL,
        fill = "Direction"
    ) +
    theme_minimal() +
    coord_flip()


## ===== Plot 2: DAR-level  =====

unique_summary_div <- unique_overlaps |>
    group_by(cell_type, direction) |>
    summarise(n = n(), .groups = "drop") |>
    group_by(cell_type) |>
    mutate(total = sum(n)) |>
    ungroup() |>
    mutate(n_signed = ifelse(direction == "Down", -n, n))

n_DARS_peaks_total_DARS <- sum(unique_summary_div$n) 

cell_order2 <- unique_summary_div |>
    group_by(cell_type) |>
    summarise(total = sum(abs(n_signed)), .groups = "drop") |>
    arrange(total) |>
    pull(cell_type)

unique_summary_div$cell_type <- factor(unique_summary_div$cell_type, levels = cell_order2)

g_darlevel <- ggplot(unique_summary_div,
                     aes(x = cell_type, y = n_signed, fill = direction)) +
    geom_col() +
    geom_text(
        aes(label = abs(n_signed),
            hjust = ifelse(direction == "Down", 1.1, -0.1)),
        size = 3, color = "black"
    ) +
    scale_y_continuous(labels = abs, expand = expansion(mult = c(0.15, 0.15))) +
    labs(
        title = "Up/Down DAR-level",
        subtitle = paste0("Multiple DARs per LinkPeak - ", n_DARS_peaks_total_DARS ,")"),
        x = "Cell Type",
        y = NULL,
        fill = "Direction"
    ) +
    theme_minimal() +
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank() 
    ) +
    coord_flip()

combined_plot <- (g_collapsed | g_darlevel) +
    plot_layout(
        guides = "collect")

ggsave(here(plotDir, paste0("overlaps_diverg_linkPeaks_vs_DARlevel_FDR", FDR, ".pdf")),
       combined_plot, width = 8, height = 8)


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


#===========================================================================
## Other strategies considered: 

# # count Up/Down per cell type
# unique_summary <- unique_overlaps |>
#     group_by(cell_type, direction) |>
#     summarise(n = n(), .groups = "drop")
# 
# unique_summary
# 
# # prepare data for diverging plot
# unique_summary_div <- unique_summary |>
#     group_by(cell_type) |>
#     mutate(total = sum(n)) |>
#     ungroup() |>
#     mutate(n_signed = ifelse(direction == "Down", -n, n))
# 
# cell_order <- unique_summary_div |>
#     group_by(cell_type) |>
#     summarise(total = sum(abs(n_signed)), .groups = "drop") |>
#     arrange(total) |>
#     pull(cell_type)
# 
# unique_summary_div$cell_type <- factor(unique_summary_div$cell_type, levels = cell_order)
# 
# # diverging barplot
# g_unique_div <- ggplot(unique_summary_div, aes(x = cell_type, y = n_signed, fill = direction)) +
#     geom_col() +
#     geom_text(
#         aes(label = abs(n_signed),
#             hjust = ifelse(direction == "Down", 1.1, -0.1)),
#         size = 3, color = "black"
#     ) +
#     scale_y_continuous(labels = abs, expand = expansion(mult = c(0.15, 0.15))) +
#     labs(
#         title = "Up/Down DAR-Link Overlaps (Unique LinkPeaks, n=2138)",
#         x = "Cell Type",
#         y = "Number of Overlaps",
#         fill = "Direction"
#     ) +
#     theme_minimal() +
#     coord_flip()
# 
# g_unique_div
# 
# f_name <- paste0("overlaps_diverg_plot_unique_ct_", FDR, ".pdf")
# ggsave(here(plotDir, f_name),
#        g_unique_div, width = 7, height = 7)

# prepare_data_to_plot <- function(
#         df_summary,
#         count_col = "n"
# ) {
#     # make a df for diverging plots (Up → right, Down → left), flip the sign of n for Down peaks
#     message("Preparing data ... ")
#     cell_totals <- df_summary |>
#         group_by(cell_type) |>
#         summarise(total = sum(n), .groups = "drop")
#     df_summary_div <- df_summary |>
#         mutate(n_signed = ifelse(direction == "Down", -n, n)) |>
#         left_join(cell_totals, by = "cell_type")
# 
#     return(df_summary_div)
# 
# }
# 
# # Unique Link peaks (per Link cluster & direction)
# # method = "majority" → assigns direction by majority of overlaps
# # method = "exclusive" → keeps only Links with exclusively Up or exclusively Down DARs
# # method = "any" → counts a Link in both categories (original behaviour, usually symmetric)
# summarise_unique_links <- function(
#         overlaps_df,
#         method = c("majority", "exclusive", "any")
# ) {
#     message("Summarizing data ...")
#     method <- match.arg(method)
#     df <- overlaps_df|>
#         group_by(cluster, peak_id) |>
#         summarise(
#             up_count   = sum(direction == "Up"),
#             down_count = sum(direction == "Down"),
#             .groups = "drop"
#         )
#     if (method == "majority") {
#         df <- df |>
#             mutate(direction = ifelse(up_count >= down_count, "Up", "Down")) |>
#             count(cluster, direction, name = "n")
#     }
#     if (method == "exclusive") {
#         df <- df |>
#             filter(xor(up_count > 0, down_count > 0)) %>%   # keep Links with only one direction
#             mutate(direction = ifelse(up_count > 0, "Up", "Down")) |>
#             count(cluster, direction, name = "n")
#     }
#     if (method == "any") {
#         df <- overlaps_df |>
#             group_by(cluster, direction) |>
#             summarise(n = n_distinct(peak_id), .groups = "drop")
#     }
# 
#     df |> rename(cell_type = cluster)
# }
# 
# 
# make_div_prop_barplots <- function(
#         summary_df,
#         summary_div_df,
#         FDR,
#         suffix
# ) {
#     # # testing:
#     # summary_df = DAR_link_summary
#     # summary_div_df = DAR_link_summary_div
#     # FDR = FDR
#     # suffix = "DAR–Link_pairwise"
#     message("Ploting ...")
# 
#     g4_overlap <- ggplot(summary_div_df,
#                          aes(x = reorder(cell_type, total), y = n_signed, fill = direction)) +
#         geom_col() +
#         # add labels outside bars
#         geom_text(aes(label = abs(n_signed),
#                       hjust = ifelse(direction == "Down", 1.1, -0.1)), # left for Down, right for Up
#                   position = position_identity(),
#                   size = 3,
#                   color = "#3D3936") +
#         scale_y_continuous(
#             labels = abs,
#             expand = expansion(mult = c(0.15, 0.15))) + # add padding
#         labs(
#             title = paste("Overlaps by cell type -", suffix),
#             subtitle = paste0("FDR = ", FDR),
#             x = "Cell Type",
#             y = "# overlapping",
#             fill = "Direction"
#         ) +
#         theme_minimal() +
#         coord_flip()
# 
#     # define factor levels by total counts
#     cell_totals <- summary_df |>
#         group_by(cell_type) |>
#         summarise(total = sum(n), .groups = "drop")
#     cell_order <- cell_totals |>
#         arrange(total) |>
#         pull(cell_type)
#     summary_df$cell_type <- factor(summary_df$cell_type, levels = cell_order)
#     summary_div_df$cell_type <- factor(summary_div_df$cell_type, levels = cell_order)
# 
#     # Proportion consistency checks
#     stopifnot(all.equal(
#         summary_df |> group_by(cell_type) |> summarise(total = sum(n)),
#         summary_div_df |> group_by(cell_type) |> summarise(total = sum(abs(n_signed)))
#     ))
# 
#     g4b_overlap <- summary_df |>
#         # Calculate the necessary proportions and cumulative sums
#         group_by(cell_type) |>
#         mutate(
#             total = sum(n),
#             prop = n / total
#         ) |>
#         # Add a cumulative proportion for position calculation *before* plotting
#         arrange(cell_type, direction) |> # Ensure consistent stacking order
#         mutate(
#             c_prop = cumsum(prop),
#             mid_prop = c_prop - (prop / 2) # Calculate the center position
#         ) |>
#         ungroup() |>
#         ggplot(aes(x = cell_type, y = prop, fill = direction)) +
#         geom_col(position = "fill") + # Draw the 100% stacked bars
#         geom_text(
#             # Use the calculated midpoint position (mid_prop) for y
#             aes(label = scales::percent(prop, accuracy = 1), y = mid_prop),
#             size = 3,
#             color = "black"
#         ) +
#         scale_y_continuous(
#             labels = scales::percent,
#             breaks = seq(0, 1, by = 0.25),
#             expand = expansion(mult = c(0, 0))
#         ) +
#         labs(
#             title = "Proportion of Up/Down",
#             x = "",
#             y = "Proportion",
#             fill = "Direction"
#         ) +
#         theme_minimal() +
#         coord_flip()
# 
#     plot_widths <- c(2, 1)
#     g4_combined <- (g4_overlap | g4b_overlap) +
#         plot_layout(
#             guides = "collect",
#             widths = plot_widths)
# 
#     f_name <- paste0("overaping_", suffix, "_FDR", FDR, ".pdf")
#     ggsave(here(plotDir, f_name),
#            g4_combined, width = 10, height = 7)
# 
#     message("All plots done!")
# 
# }
# 
# 
# #### Get the counts per cell type for each strategy ####
# 
# # Pairwise (DAR–Link edges)
# summarise_pairwise <- function(overlaps_df) {
#     overlaps_df |>
#         group_by(cell_type, direction) |>
#         summarise(n = n(), .groups = "drop")
# }
# 
# # Unique DAR peaks (per cell type & direction)
# summarise_unique_dars <- function(overlaps_df) {
#     overlaps_df |>
#         group_by(cell_type, direction) |>
#         summarise(n = n_distinct(peak_id.1), .groups = "drop")
# }
# 
# ## Define strategies
# strategies <- list(
#     "DAR-Link_pairwise" = summarise_pairwise,
#     "DAR-Link_unique"   = summarise_unique_dars,
#     "Link-DAR_unique_majority"  = function(df) summarise_unique_links(df, method = "majority"),
#     "Link-DAR_unique_exclusive" = function(df) summarise_unique_links(df, method = "exclusive")
#     # "Link-DAR_unique_any"     = function(df) summarise_unique_links(df, method = "any") # optional
# )
# 
# # Run all strategies
# walk2(strategies, names(strategies), function(fun, suffix) {
# 
#     # Apply summariser
#     summary_df <- fun(overlaps_df)
# 
#     # Prepare diverging data
#     summary_div <- prepare_data_to_plot(summary_df)
# 
#     # Make plots + save
#     make_div_prop_barplots(summary_df, summary_div, FDR, suffix)
# 
#     # Save raw summary
#     f_name <- here(processedDir, paste0("summary_", suffix, "_FDR", FDR, ".csv"))
#     write.csv(summary_df, f_name, row.names = FALSE)
# 
#     message("Finished: ", suffix)
# })




