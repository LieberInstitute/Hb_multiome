########################################################################
## Summaries/Prepare divergence plots for interpreting Peaks Overlapings
##
## Authors. CSC
## Date. Oct 02, 2025
## Recommended resources on interactive mode: srun --pty --mem=15GB --x11 bash
########################################################################

library("purrr")
library("patchwork")
library("tidyverse")
library("here")

#===============================================================================
# resolution_level = "Mid"      # 18 cell-types
#===============================================================================

## setup variable names

# resolution_level = "Mid" 
FDR = 0.1 # actual value to run the DAR-Links
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

## load full overlaps df, which is already FDR filtered
overlaps_df <- read_csv(
    here(inputCSV_Overlaps_Dir, "Overlaps_LinkPeak_DARs_FDR0.1_0.1.csv"),
    show_col_types = FALSE
)
# head(overlaps_df)


##==============================================================================
## Declare functions:

plt_shared_overlaps <- function(
        shared_overlaps,
        FDR,
        all_shared=TRUE
) {
    
    if (all_shared) { 
        title = "All 2 Cell Type Shared Overlaps"
    } else {
        title = "Only Hb Related"
    }
    
    color_vector <- c(
        "MHb/LHb Related" = "#ad1d8c", 
        "Other Cell Types" = "grey70"
    )
    
    g_shared_2ct <- shared_overlaps |>
        ggplot(aes(
            # Use .desc = FALSE (Ascending count) to achieve largest bar at the top of the flipped plot.
            x = fct_reorder(cell_types, cell_types, .fun = length, .desc = FALSE), 
            fill = bar_color_group
        )) + 
        geom_bar(color = "black") +
        geom_text(
            stat = "count", 
            aes(label = after_stat(count)), 
            hjust = -0.5, 
            size = 3
        ) +
        coord_flip() + # Flip coordinates for readable cell type labels
        scale_y_continuous(expand = expansion(mult = c(0, 0.25))) + # add 20% space on right
        scale_fill_manual(values = color_vector) +
        labs(
            title = title,
            subtitle = paste0("Overlaps: ", nrow(shared_overlaps), " | FDR thr = ", FDR),
            x = NULL,
            y = "Number of Overlaps"
        ) +
        theme_minimal() +
        theme(legend.position = "none")
    
    return(g_shared_2ct)
    
}

## create a single summary table for shared ct-levels
make_combined_summary_divergence_plots <- function( 
        all_overlaps_raw
){
    
    combined_summary <- all_overlaps_raw |>
        group_by(cell_type, peak_id, overlap_type) |>
        summarise(
            peak_accessibility = case_when(
                all(dar_logFC > 0) ~ "More",
                all(dar_logFC < 0) ~ "Less",
                TRUE ~ "Mixed" # Mark as Mixed if LinkPeak has both More and Less DARs
            ),
            .groups = "drop" 
        ) |>
        count(cell_type, overlap_type, peak_accessibility, name = "n") |>
        filter(peak_accessibility %in% c("More", "Less")) |> 
        
        # Add n_signed column for divergence (More = Negative, Less = Positive)
        mutate(
            n_signed = case_when(
                peak_accessibility == "More" ~ n, 
                peak_accessibility == "Less" ~ -n, 
                TRUE ~ 0
            ),
            # create a combined factor for filling/coloring/legend
            fill_group = factor(
                paste(overlap_type, peak_accessibility),
                levels = c("Shared Less", "Unique Less",
                           "Unique More", "Shared More")
            )
        )
    
    ## define cell order (based on total count magnitude)
    cell_order_combined <- combined_summary |>
        group_by(cell_type) |>
        summarise(total = sum(abs(n_signed)), .groups = "drop") |>
        arrange(total) |>
        pull(cell_type)
    
    combined_summary$cell_type <- factor(combined_summary$cell_type, levels = cell_order_combined)
    
    return(combined_summary)
    
}


## Single Stacked Diverging Plot 
plt_divergence_stacked <- function(
        combined_summary,
        n_unique_links,
        n_shared_links,
        FDR,
        ct_level,
        plotDir
){
    
    plot_title = paste("LinkPeaks-DARs: Unique vs. ", ct_level, "Shared Cell-Types")
    color_palette <- c(
        # Unique Group (Blue Tones) - to fit with Hb pilot colors 
        "Unique More" = "#1f78b4",  # Medium Blue (Original)
        "Unique Less" = "#AEC7E8",  # Light Blue (Lighter tint of #1f78b4)
        
        # Shared Group (Red/Teal Tones for better distinction)
        "Shared More" = "#ad1d8c",  # Rich Red/Crimson (Highly distinct from blue)
        "Shared Less" = "#FFB6C1"   # Light Pink/Rose (Lighter tint of the Rich Red)
    )
    
    g_divergence_stacked <- ggplot(combined_summary,
                                   aes(x = cell_type, y = n_signed, fill = fill_group)) +
        geom_col(color = "black", linewidth = 0.3) + 
        # geom_text(
        #     aes(label = n), # raw count to display
        #     position = position_stack(vjust = 0.5),
        #     size = 3
        # ) +
        coord_flip() +
        scale_fill_manual(
            values = color_palette,
            name = "Overlap - Accessibility"
        ) +
        scale_y_continuous(
            labels = abs, # Display only positive numbers on the axis
            expand = expansion(mult = c(0.15, 0.15))
        ) +
        labs(
            title = plot_title,
            subtitle = paste0("Unique = ", n_unique_links, " | Shared = ", n_shared_links, " | FDR thr = ", FDR),
            x = NULL,
            y = "Number of Overlaps"
        ) +
        theme_minimal() +
        theme(
            axis.text.y = element_text(size = 11), 
            axis.title.y = element_text(size = 12),
            axis.text.x = element_text(size = 11),
            axis.title.x = element_text(size = 12)
        )
    
    f_name <- paste0("overlaps_unique_and_", ct_level, "_shared_stackedbar_", FDR, ".pdf")
    ggsave(here(plotDir, f_name),
           g_divergence_stacked, width = 8, height = 8)
    
    #return(g_divergence_stacked)
    message("Divergence plots completed!")
    
}


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
f_name <- here(processedDir, paste0("overlaps_linkPeaks_unique_shared_detailed_FDR", FDR, ".csv"))
write.csv(unique_shared_overlaps, f_name, row.names = FALSE)


## Join classification back to cell-type level and count distinct LinkPeaks
plot_df_ct <- unique_shared_overlaps |>
    group_by(cell_type, overlap_type) |>
    summarise(
        n_links = n(),
        .groups = "drop"
    )
head(plot_df_ct)

## Calculate total counts for each overlap type
total_counts <- plot_df_ct |>
    group_by(overlap_type) |>
    summarise(total_n = sum(n_links))
# Create named vectors for colors and labels
overlap_colors <- c("Shared" = "#ad1d8c", "Unique" = "#1f78b4")
# Define the custom labels using the calculated totals
overlap_labels <- c(
    "Shared" = paste0("Shared (N=", total_counts$total_n[total_counts$overlap_type == "Shared"], ")"),
    "Unique" = paste0("Unique (N=", total_counts$total_n[total_counts$overlap_type == "Unique"], ")")
)

# Compute total_links first
plot_df_ct <- plot_df_ct |>
    group_by(cell_type) |>
    mutate(total_links = sum(n_links)) |>
    ungroup()

g_cell_type_overlap <- plot_df_ct |>
    ggplot(aes( # Order cell types by the total number of links (sum of Unique + Shared)
        x = fct_reorder(cell_type, n_links, .fun = sum), 
        y = n_links, 
        fill = overlap_type
    )) +
    # Use geom_col (or geom_bar(stat="identity")) for counts
    geom_col(position = position_stack(reverse = TRUE), color = "black") +
    geom_text( # Add total labels per cell type (outside bars)
        data = plot_df_ct |> distinct(cell_type, total_links),
        aes(
            x = cell_type, 
            y = total_links,
            label = total_links
        ),
        nudge_y = 150,   # Move text slightly to the right
        size = 4,
        inherit.aes = FALSE
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
        legend.background = element_rect(colour = "gray80", fill = "white"), 
        axis.text.y = element_text(size = 11), 
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.title.x = element_text(size = 12)
    )

f_name = paste0("overlaps_unique_shared_ct_FDR", FDR, ".pdf")
ggsave(here(plotDir, f_name),
       g_cell_type_overlap, width = 6, height = 8)


##==============================================================================
## Summarize unique LinkPeaks with n_DARs and n_cell_types

link_summary <- overlaps_df |>
    group_by(peak_id) |>
    summarise(
        n_DARs       = n(),                         # total DAR rows linked to this LinkPeak
        n_cell_types = n_distinct(cell_type),       # distinct od distinc cell types represented
        cell_types   = paste(unique(cell_type), collapse = "; "), # list cell types
        .groups = "drop"
    ) |>
    # LinkPeaks overlapping DARs from exactly one cell type (new defin ition, consistent with shared)
    mutate(overlap_type = ifelse(n_cell_types == 1, "Unique", "Shared"))

table(link_summary$overlap_type)
# Shared Unique 
# 3636   1967

f_name <- here(processedDir, paste0("overlaps_summary_linkPeak_DARs_unique_ct_FDR", FDR, ".csv"))
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
        values = c("Unique" = "#1f78b4", "Shared" = "#ad1d8c") 
    ) +
    labs(
        #title = "LinkPeak Overlap: Unique vs. Multi-Cell Type DARs",
        x = "Classification",
        y = "Number of LinkPeaks-DARs",
        fill = "Overlap Type"
    ) +
    theme_minimal() +
    theme(legend.position = "none")


## Distribution of Cell Type Multiplicity

g_multiplicity <- link_summary |>
    mutate(n_cell_types_factor = factor(n_cell_types)) |>
    # create the color classification column, I want first red
    mutate(
        bar_color = case_when(
            n_cell_types == 1 ~ "#1f78b4",  
            n_cell_types == 2 ~ "#ad1d8c",  
            TRUE ~ "grey"                   
        )
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
        x = "Cell-Type Overlapping Levels",
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
    )

f_name <- paste0("overlaps_combined_unique_distinct_shared_ct_FDR", FDR, ".pdf")
ggsave(here(plotDir, f_name),
       combined_plot, width = 7, height = 7)



##==============================================================================
## Make "Distribution of Cell Type-Specific LinkPeaks (Unique Overlaps)" plot

message("Creating barplot showing the frequency of each unique cell type .." )

# Filter the link_summary dataframe to include only "Unique" overlaps
uniques_df <- link_summary |>
    filter(overlap_type == "Unique") |>
    arrange(desc(n_DARs))
nrow(uniques_df) # 1967

g_uniques <- uniques_df |>
    ggplot(aes(
        # Use .desc = FALSE (Ascending count) to achieve largest bar at the top of the flipped plot.
        x = fct_reorder(cell_types, cell_types, .fun = length, .desc = FALSE), 
        fill = cell_types
    )) + 
    geom_bar(color = "black",  fill = "#1f78b4") +
    geom_text(
        stat = "count", 
        aes(label = after_stat(count)), 
        hjust = -0.5, 
        size = 4
    ) +
    coord_flip() + # Flip coordinates for readable cell type labels
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) + # add 10% space on right
    labs(
        title = "Cell Type-Specific LinkPeaks-DARs",
        subtitle = paste0("Total Overlaps: ", nrow(uniques_df), " | FDR thr = ", FDR),
        x = NULL,
        y = "Number of Unique Overlaps"
    ) +
    theme_minimal() +
    theme(
        #plot.title = element_text(hjust = 0.5),
        legend.position = "none"
    )

f_name <- paste0("overlaps_unique_ct_FDR", FDR, ".pdf")
ggsave(here(plotDir, f_name),
       g_uniques, width = 6, height = 8)



#===== Distribution of 2 Distinct Cell-Types Overlapping 

message("Creating barplot showing the frequency of 2 distinct cell-types .." )

table(link_summary$n_cell_types==2)
table(link_summary$n_cell_types)
#   1    2    3    4    5    6    7    8    9   10   11 
# 1967 1238  912  693  381  240  100   46   15    9    2 

# Filter the link_summary df to include shared overlaps with 2 shared cell-types
shared_2ct_df <- link_summary |>
    filter(overlap_type == "Shared" & n_cell_types == 2) |>
    arrange(desc(n_DARs))
nrow(shared_2ct_df) # 1238
head(shared_2ct_df)

# save significant DARs
f_name <- here(processedDir, paste0("overlaps_summary_linkPeak_DARs_shared_2ct_FDR", FDR, ".csv"))
write.csv(shared_2ct_df, f_name, row.names = FALSE)


# Prepare Data: add the classification column for all shared 2-ct overlaps
shared_2ct_df_colored <- shared_2ct_df |>
    mutate(
        is_hdb_related = grepl("MHb|LHb", cell_types),
        bar_color_group = case_when(
            is_hdb_related ~ "MHb/LHb Related",
            TRUE ~ "Other Cell Types"
        )
    )

## All shared `cell_types` 
g1_ove_all <- plt_shared_overlaps(shared_2ct_df_colored, FDR, all_shared=TRUE)
f_name = paste0("overlaps_shared_2ct_FDR", FDR, ".pdf")
ggsave(here(plotDir, f_name),
       g1_ove_all, width = 7, height = 8)

## only those shared `cell_types` that are MHb or LHb related 
shared_2_Hb_ct_df_colored <- shared_2ct_df_colored |>
    filter(is_hdb_related==TRUE)
#nrow(shared_2_Hb_ct_df_colored)
g2_ove_hb <- plt_shared_overlaps(shared_2_Hb_ct_df_colored, FDR, all_shared=FALSE)

combined_2ct_plot <- (g1_ove_all | g2_ove_hb) +
    plot_layout(
        guides = "collect",
        axis = "collect" )

f_name = paste0("overlaps_shared_2ct_FDR", FDR, "_combined.pdf")
ggsave(here(plotDir, f_name),
       combined_2ct_plot, width = 9, height = 8)




##==============================================================================
# Calculate and plot top 3 most frequent cell type combinations for each shared level. I have 10 levels 
# plot the results for the first 6 shared levels (e.g., n=1 through n=6).

message("Calculate the top 3 most frequent cell type combinations for each shared level")

# aggregate your link_summary data to find the top 3 most common cell type combinations 
# for each shared level (n_cell_types)

# Calculate up to 3 rows per shared level ( I have set 10 level )
head(link_summary)
top_combinations_key <- link_summary |>
    # Group by the shared level (the bar) and the specific cell type combination
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

f_name <- here(processedDir, paste0("overlaps_summary_linkPeak_ALL_shared_top3_freq_ct_", FDR, ".csv"))
write.csv(link_summary, f_name, row.names = FALSE)

message("Summary overlaps done!")

## plot shared peaks by shared cell-type level
g_top5_shared_combinations <- ggplot(
    #top_combinations_key,
    # Filter to include only the first 6 shared levels (n=1 to n=5)
    top_combinations_key |> filter(`Multiplicity (n)` <= 5),
    aes(x = reorder(display_text, n_peaks), y = n_peaks, fill = factor(`Multiplicity (n)`))
) +
    geom_col(show.legend = FALSE) +
    coord_flip() +
    facet_wrap(~`Multiplicity (n)`, scales = "free_y", ncol = 1) +
    labs(
        title = "Top 3 Cell-Type Combinations\nper Shared Level",
        subtitle = "Shared Cell Types",
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

f_name <- paste0("overlaps_shared_5ct_top3_FDR", FDR, ".pdf")
ggsave(here(plotDir, f_name),
       g_top5_shared_combinations, width = 8, height = 8)



#===========================================================================
## Process and Plot Accessibility UNIQUE LinkPeaks-DARs overlapping per cell type

message("Processing Open/Close LinkPeaks-DARs overlappings per cell type ... ")

## Define Accessibility  for all the unique overlaps
overlaps_df <- overlaps_df |>
    mutate(accessibility = case_when(
        dar_logFC > 0 ~ "More",
        dar_logFC < 0 ~ "Less",
        TRUE ~ "Neutral"
    )) |>
    filter(accessibility %in% c("More", "Less"))

nrow(overlaps_df) # 23691
nrow(uniques_df) # 1967

## verification
overlaps_df |>
    count(peak_id, sort = TRUE) |>
    filter(n > 1) |> head()

unique_peak_ids_to_filter <- uniques_df |>
    # We pull the peak_id column, which seems to correspond to the ATAC-seq peak identifier
    pull(peak_id) |>
    unique()
n_uniques_peak_ids <- length(unique_peak_ids_to_filter) # 1967

## Restrict overlaps to Unique LinkPeaks only
unique_overlaps <- overlaps_df |>
    filter(peak_id %in% unique_peak_ids_to_filter) 

nrow(unique_overlaps) # [1] 2904
head(unique_overlaps)

f_name <- here(processedDir, paste0("overlaps_unique_linkPeak_DARs_detail_FDR", FDR, ".csv"))
write.csv(unique_overlaps, f_name, row.names = FALSE)

message("Saved unique overlaps csv!")


## ==== Prepare data for divergence plots: Peak-level / Uniques vs shared-level 

unique_summary_collapsed <- unique_overlaps |>
    group_by(cell_type, peak_id) |>
    summarise(
        accessibility = case_when(
            all(dar_logFC > 0) ~ "More",
            all(dar_logFC < 0) ~ "Less",
            TRUE ~ "Mixed"   # at least one Up and one Down DAR
        ),
        .groups = "drop"
    ) |>
    count(cell_type, accessibility) |>
    group_by(cell_type) |>
    mutate(total = sum(n)) |>
    ungroup() |>
    mutate(
        # n_signed represent the count of unique DARs for each cell_type and accessibility group,
        # where the sign (+ or -) determines the direction in which the bar will be plotted
        n_signed = case_when(
            accessibility == "More"  ~ n, 
            accessibility == "Less"    ~ -n, 
            accessibility == "Mixed" ~ n   # show Mixed as positive
        )
    )
head(unique_summary_collapsed)
# cell_type  accessibility     n total n_signed
# <chr>      <chr>         <int> <int>    <int>
# 1 Astrocyte  Less             33    53       33
# 2 Astrocyte  More             20    53      -20

table(unique_summary_collapsed$accessibility)
# Less More 
# 10   15

cell_order1 <- unique_summary_collapsed |>
    group_by(cell_type) |>
    summarise(total = sum(abs(n_signed)), .groups = "drop") |>
    arrange(total) |>
    pull(cell_type)

unique_summary_collapsed$cell_type <- factor(unique_summary_collapsed$cell_type, levels = cell_order1)

unique_peak_ids <- uniques_df |> pull(peak_id)

# filter, label, and combine raw overlaps data
all_overlaps_raw <- overlaps_df |>
    filter(accessibility %in% c("More", "Less")) |> # Filter out Neutral
    mutate(
        overlap_type = ifelse(
            peak_id %in% unique_peak_ids, 
            "Unique", 
            "Shared"
        )
    )

## get total counts / 2-shared overlaps

n_unique_links <- all_overlaps_raw |> filter(overlap_type == "Unique") |> pull(peak_id) |> unique() |> length()
# 1967
n_shared_links <- all_overlaps_raw |> filter(overlap_type == "Shared") |> pull(peak_id) |> unique() |> length()
# 3636

## make summary table to divergence plot
combined_summary <- make_combined_summary_divergence_plots(
    all_overlaps_raw)

## verify
unique_summary_collapsed |> 
    group_by(accessibility) |> summarise(sum_signed = sum(n_signed))
# Expect: sum_signed(More) > 0, sum_signed(Less) < 0
combined_summary |> 
    group_by(peak_accessibility) |> summarise(sum_signed = sum(n_signed))
# Expect: sum_signed(More) > 0, sum_signed(Less) < 0

## make divergence plot
plt_divergence_stacked(
    combined_summary,
    n_unique_links,
    n_shared_links,
    FDR,
    "2",
    plotDir
)

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



