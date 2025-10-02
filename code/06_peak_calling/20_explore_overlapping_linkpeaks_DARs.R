########################################################################
## Summaries: three strategies (pairwise, unique DARs, unique Links)
##
## Authors. CSC
## Date. Sep 24, 2025
## Recommended resources on interactive mode: srun --pty --mem=30GB --x11 bash
## Note. Seurat objects were created with module load conda_R/4.3.x
########################################################################

library("GenomicRanges") 
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
resolution_level = "Mid" 
# test 2 thresholds
FDR = 0.2 # actual value to run the DAR-Links
FDR_thresh = c(0.1, 0.2)
# lfc_thresh = 0.2 # we do not use log FC for this exploratory analysis

## Set directory names
inputCSV_Links_Dir <- here(
    "processed-data",
    "06_peak_calling",
    "14_exploratory_pb_peak_scores_MACS2",
    "links_ct_merged"
)
inputCSV_DARs_Dir <- here(
    "processed-data",
    "06_peak_calling",
    "17_pseudobulk_DARs_MACS2_reduced_voomLmFit"
)
processedDir <- here(
    "processed-data",
    "06_peak_calling",
    "19_Linkage_DARs_analysis"
)
plotDir <- here(
    "plots",
    "06_peak_calling",
    "19_Linkage_DARs_analysis"
)

if (!dir.exists(processedDir)) { dir.create(processedDir) }
if (!dir.exists(plotDir)) { dir.create(plotDir) }

##==============================================================================

message("Loading  LinkPeaks results ...")

# List all files matching the specific clustering resolution level
lst_peak_files <- list.files(
    path = inputCSV_Links_Dir,
    pattern = paste0("^", resolution_level, ".*\\.csv$"),
)
#lst_peak_files = list.files(path = input_cvsDir)
message(length(lst_peak_files), " Link peak-genes files found ... ")
# lst_peak_files
## tier field is missing ?

# bind all files 
list_of_df <- lapply(here(inputCSV_Links_Dir, lst_peak_files), read.csv)
linkPeaks_results_all <- bind_rows(list_of_df)

message("Processing ", length(list_of_df), " LinkPeaks files for merged-peaks dataset")
message(nrow(linkPeaks_results_all), " total links") 

head(linkPeaks_results_all)


#===========================================================================

message("Ploting barplot of number of links by cellType .. ")

# get number of links (rows) for each cluster
cluster_counts <- linkPeaks_results_all |>
    count(cluster, sort = TRUE)

g1 <- ggplot(cluster_counts, aes(x = reorder(cluster, n), y = n)) +
    geom_col(fill = "steelblue") +
    geom_text(aes(label = n), hjust = -0.2, size = 3) +
    labs(
        title = "Links by Cell Type",
        subtitle = paste0("FDR thr = ", FDR),
        x = "Cell Type",
        y = "Number of Links"
    ) +
    theme_minimal() +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))  # add 10% space on right

f_name <- paste0("links_barplot_ct_FDR", FDR, ".pdf")
ggsave(here(plotDir, f_name),
       g1, width = 7, height = 7)



##==============================================================================

message("Loading DARs results ...")

all_DARs_files <- list.files(
    path = inputCSV_DARs_Dir,
    pattern = "voomlmFit_DAR_peaks.*\\.csv$",
    full.names = FALSE
)
# Subset the negation (files NOT containing "ALL_in_")
lst_DARs_files_cellType <- grep("voomlmFit_DAR_peaks_ALL_in_", all_DARs_files, value = TRUE, invert = TRUE)

message(length(lst_DARs_files_cellType), " DAR files found ... ")
lst_DARs_files_cellType


#===========================================================================

message("Parsing DARs for each cellType at ", FDR_thresh)

DARS_signif_df_lst = list()
# FDR_thresh = c(0.1, 0.2)

## extract/filters DARs from a given list of DARs at specific FDR thresh
filter_signific_DARs <- function(
        lst_DARs_files_cellType, # DARs csv list
        fdr_thresh,
        output_Dir
) {
    
    # testing
    # ct =  lst_DARs_files_cellType[2]
    # fdr_thresh = 0.2
    
    DARS_signif_df_lst <- list()
    
    for (ct in lst_DARs_files_cellType) {
        ct_name <- gsub("voomlmFit_DAR_peaks_|.csv", "", ct)
        message("Processing DARs for [", ct_name, "] with FDR=", fdr_thresh)
        
        DARs_csv <- here(inputCSV_DARs_Dir, ct)
        DARs_df <- read.csv(DARs_csv)
        
        message(nrow(DARs_df), " ", ct_name, " total DARs found ... ")
        
        DARs_significant <- DARs_df |>
            select(
                peak_id,
                logFC = paste0("logFC_", ct_name),
                fdr   = paste0("fdr_", ct_name)
            ) |>
            mutate(
                fdr = as.numeric(fdr),  # keep actual fdr values
                is_significant = if_else(fdr < fdr_thresh, "Significant", "Not Significant")
            ) |> 
            filter(fdr < fdr_thresh)
        
        DARS_signif_df_lst[[ct_name]] <- DARs_significant
        
        message(nrow(DARs_significant), " ", ct_name, " total significant DARs found ... ")
    
    }
    
    DARs_all <- bind_rows(DARS_signif_df_lst, .id = "cell_type")
    
    # save significant DARs
    f_name <- here(output_Dir, paste0("ALL_DARs_signif_FDR", gsub("\\.", "", as.character(fdr_thresh)), ".csv"))
    write.csv(DARs_all, f_name, row.names = FALSE)
    
    message(nrow(DARs_all), " total significant DARs found ... ")
    message("ALL DARs significant saved ... ")  
    
    return(DARs_all)

}

## for each FDR_thresh to test
DARs_results_all <- purrr::map_dfr(
    FDR_thresh, 
    ~ filter_signific_DARs(lst_DARs_files_cellType, .x, processedDir),
    .id = "FDR_threshold"
)

# Convert the .id column (1, 2) into meaningful labels
DARs_results_all <- DARs_results_all |>
    mutate(FDR_threshold = factor(FDR_threshold, 
                                  labels = paste0("FDR", FDR_thresh)))

unique(DARs_results_all$FDR_threshold)
DARs_results_all |> head()


#===========================================================================

message("Plotting barplot of DARs by cellType ...")

# Make FDR_threshold a factor for facet labels
DARs_results_all <- DARs_results_all |>
    mutate(FDR_threshold = as.factor(FDR_threshold))

# count DARs
DARs_counts <- DARs_results_all |>
    count(cell_type, FDR_threshold, sort = TRUE)

# Get unique thresholds
thresholds <- unique(DARs_counts$FDR_threshold)

# Loop through thresholds and plot one at a time
g2 <- lapply(thresholds, function(th) {
    
        df_sub <- filter(DARs_counts, FDR_threshold == th)
        
        ggplot(df_sub, aes(x = reorder(cell_type, n), y = n)) +
            geom_col(fill = "steelblue") +
            geom_text(aes(label = n), hjust = -0.2, size = 3) +
            labs(
                title = paste("DARs by Cell Type"),
                subtitle = paste("FDR thr =", str_split(th, "FDR")[[1]][2]),
                x = "Cell Type",
                y = "Number of DARs"
            ) +
            theme_minimal() +
            coord_flip() +
            scale_y_continuous(expand = expansion(mult = c(0, 0.2)))  # add 10% space on right
    })

walk2(
    .x = g2,           
    .y = thresholds,   # thresholds for filenames
    ~ ggsave(
        filename = file.path(plotDir, paste0("DARs_barplot_ct_FDR", .y, ".pdf")),
        plot = .x,
        width = 7,
        height = 7
    )
)

## overlay both thr for comparison purposes 

# Turn into a nicely formatted label
subtitle_label <- paste(
    names(table(DARs_results_all$FDR_threshold)),
    table(DARs_results_all$FDR_threshold),
    sep = ": ",
    collapse = " | "
)

g3 <- ggplot(DARs_counts, aes(x = reorder(cell_type, n), y = n, fill = FDR_threshold)) +
    geom_col(alpha = 0.5, position = "identity") +  # <- overlay instead of dodge
    labs(
        title = paste("DARs by Cell Type"),
        subtitle = subtitle_label,
        x = "Cell Type",
        y = "Number of DARs",
        fill = "FDR Threshold"
    ) +
    theme_minimal() +
    coord_flip()  +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2)))  # add 10% space on right

f_name <- paste0("DARs_barplot_overlay_FDRs.pdf")
ggsave(here(plotDir, f_name),
       g3, width = 7, height = 7)



#===========================================================================

message("Preparing GRanges for overlaping peaks ... ")

# inspect data
nrow(linkPeaks_results_all) # [1] 10948
colnames(linkPeaks_results_all)
head(linkPeaks_results_all$peak) # join id key
nrow(DARs_results_all) # [1] 582590
colnames(DARs_results_all)
head(DARs_results_all$peak_id) # join id key

# Split the peak column into coordinates. Note:
# seqnames / start / end: full chromatin accessibility window carried over in the aggregated assay.
# peak: a canonicalized peak ID string used for joins, overlaps, and downstream correlation reporting

link_split <- tidyr::separate(linkPeaks_results_all, peak,
                              into = c("seqnames","start","end"), sep = "-") |>
    mutate(start = as.integer(start),
           end   = as.integer(end))

# Convert LinkPeaks to GRanges to run findOverlaps()
gr_links <- 
    with(link_split, 
        GRanges(
             seqnames = seqnames,
             ranges   = IRanges(start, end),
             strand   = strand,
             # Metadata columns passed directly:
             peak_id_links = paste(seqnames, start, end, sep = "-"),
             CCscore  = score, # spearman CC
             gene_name  = gene,
             gene_id  = gene_id, 
             FDR_CC = FDR,
             cluster  = cluster,
             tss = tss,
             # the strand of the linked gene (+ or -)
             gene_strand = gene_strand, 
             distance = distance,
             distance_kb = distance_kb,
             # distance from peak center to TSS, 
             # keeping genomic sign (positive = downstream, negative = upstream)
             signed_distance = signed_distance, 
             # adjusts the sign of signed_distance to reflect the gene’s strand:
             # If the gene is on +, then upstream is negative, downstream positive
             # If the gene is on -, it flips accordingly
             signed_by_strand = signed_by_strand
        )
    )
length(gr_links) # [1] 10948
head(gr_links)


## ONLY use DARs at FDR02

nrow(DARs_results_all) # [1] 582590
DARs_FDRX <- DARs_results_all |>
    filter(FDR_threshold == paste0("FDR", FDR)) 
nrow(DARs_FDRX) # [1] 325580

# split the peak_id "chr-start-end" into seqnames, start, end
colnames(DARs_FDRX)
nrow(DARs_FDRX) # [1] 325580
head(DARs_FDRX)
# FDR_threshold cell_type            peak_id      logFC          fdr
# 1        FDR0.2 Astrocyte chr1-629811-630032  1.1115262 1.970054e-07
# 2        FDR0.2 Astrocyte chr1-633694-634122  1.1348873 8.347289e-09
DARs_split <- tidyr::separate(DARs_FDRX, peak_id, into = c("seqnames", "start", "end"), sep = "-") |>
    mutate(start = as.integer(start),
           end   = as.integer(end))
nrow(DARs_split) # [1] 325580
head(DARs_split)

# In GRanges(), anything other than seqnames, ranges, and strand is read as metadata
# Differential Accessibility Regions (DARs) are treated as unstranded 
gr_dars <- 
    with(DARs_split, 
        GRanges(
            seqnames = seqnames,
            ranges   = IRanges(start = start, end = end),
            strand   = "*",
            # Metadata columns passed directly:
            peak_id        = paste(seqnames, start, end, sep = "-"),
            cell_type      = cell_type,
            FDR_threshold  = FDR_threshold, # 0.1 / 0.2 
            # Add FC of chromatin accessibility:
            # - logFC > 0 → peak is more accessible (open) in the target cell type
            # - logFC < 0 → peak is less accessible (closed) in the target cell type
            logFC          = logFC,
            fdr_dars = fdr
            #is_significant = is_significant
    )
)

## inspect
length(gr_dars) # [1] 325580
head(gr_dars)



#===========================================================================

message("Finding overlaps ...") 

# ## testing:
# library("ggvenn")
# venn_list <- list(
#     DARs = unique(gr_dars$peak_id),
#     LinkPeaks = unique(gr_links$peak_id_links)
# )
# length(venn_list$LinkPeaks) # [1] 7098
# length(venn_list$DARs) # [1] 211558
# ggvenn(venn_list,
#        fill_color = c("skyblue", "orange"),
#        stroke_size = 0.5,
#        set_name_size = 4)

# findOverlaps() reports all pairs of ranges that overlap # “≥1 bp overlap”
# If one LinkPeaks region overlaps many DARs regions, you’ll get multiple rows for the same LinkPeak
#hits <- findOverlaps(gr_links, gr_dars)

# exact coordinate matches 
exact_hits <- findOverlaps(gr_links, gr_dars, type="equal") 
length(exact_hits) # 20650
## mini-test
# unique_equal <- unique(mcols(gr_links)$peak_id[queryHits(exact_hits)])
# length(unique_equal)

# Extract matched ranges
exact_links <- gr_links[queryHits(exact_hits)]
exact_dars  <- gr_dars[subjectHits(exact_hits)]
message("exact_links: ", length(exact_links))
message("exact_dars: ", length(exact_dars))

# # Compute percent overlap for both query (links) and subject (DARs)
# ov <- findOverlaps(gr_links, gr_dars)
# pi <- pintersect(gr_links_ct[queryHits(ov)], gr_dars_ct[subjectHits(ov)])
# prop_query  <- width(pi) / width(gr_links_ct[queryHits(ov)])
# prop_subject <- width(pi) / width(gr_dars_ct[subjectHits(ov)])
# # Keep only reciprocal overlaps ≥50% on both sides
# hits <- ov[prop_query >= 0.5 & prop_subject >= 0.5]

# Combine metadata from both sides
overlaps_df <- data.frame(
    as.data.frame(mcols(gr_links)[queryHits(exact_hits), ]),
    as.data.frame(mcols(gr_dars)[subjectHits(exact_hits), ])
)

nrow(overlaps_df) # [1] 20650
head(overlaps_df)

# Summarize LinkPeak overlaps
link_summary <- overlaps_df |>
    group_by(peak_id) |>
    summarise(
        n_DARs       = n(),                         # total matched DAR rows
        n_cell_types = n_distinct(cell_type),       # how many distinct DAR cell types
        cell_types   = paste(unique(cell_type), collapse = "; "), # list cell types
        .groups = "drop"
    ) |>
    mutate(overlap_type = ifelse(n_DARs == 1, "Unique", "Replicated"))

f_name <- here(processedDir, paste0("summary_LinkPeak_overlap_stats_FDR", FDR, ".csv"))
write.csv(link_summary, f_name, row.names = FALSE)

message("LinkPeak summary saved: ", f_name)

# ## Plot "Count of Unique vs. Replicated Overlaps"
# # Counts from the overlap_type column
# g_overlap_type <- link_summary |>
#     ggplot(aes(x = overlap_type, fill = overlap_type)) +
#     geom_bar(color = "black") +
#     geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 4) +
#     labs(
#         #title = "LinkPeak Overlap: Unique vs. Multi-Cell Type DARs",
#         x = "Overlap Classification",
#         y = "Number of LinkPeaks (Unique Regions)",
#         fill = "Overlap Type"
#     ) +
#     theme_minimal() +
#     theme(legend.position = "none")
# 
# ## Distribution of Cell Type Multiplicity
# g_multiplicity <- link_summary |>
#     mutate(n_cell_types_factor = factor(n_cell_types)) |>
#     ggplot(aes(x = n_cell_types_factor)) +
#     geom_bar(fill = "#56B4E9", color = "black") + # Use a distinct color
#     geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 3) +
#     labs(
#         #title = "Multiplicity of Cell Type Overlap (DARs per LinkPeak)",
#         x = "Number of Distinct Cell Types Overlapping the LinkPeak",
#         #y = "Number of LinkPeaks (Unique Regions)"
#     ) +
#     theme_minimal() +
#     theme(legend.position = "none")
# 
# # Adjust widths (1, 1.5) and keep a single y-axis title
# combined_plot <- g_overlap_type + g_multiplicity +
#     plot_layout(widths = c(1, 2)) &
#     theme(axis.title.y = element_text(size = 10))
# 
# # Add a shared y-axis label
# combined_plot <- combined_plot +
#     plot_annotation(
#         title = "LinkPeak-DARs Overlap",
#         subtitle = NULL,
#         caption = NULL,
#         theme = theme(
#             plot.margin = margin(5, 5, 5, 5),
#             axis.title.y = element_text(size = 12)
#         )
#     ) &
#     ylab("Number of LinkPeaks (Unique Regions)")
# 
# f_name <- paste0("Overlaps_Unique_vs_Replicated_FDR", FDR, ".pdf")
# ggsave(here(plotDir, f_name),
#        g_overlap_type, width = 7, height = 7)
# 
# 
# #===========================================================================
# ## Up vs Down DAR-Links count per cell type
# 
# # Define Up vs Down
# # peak_id → the LinkPeak ID (from gr_links)
# # peak_id.1 → the DAR peak ID (from gr_dars)
# overlaps_df <- overlaps_df |>
#     mutate(direction = case_when(
#         logFC > 0 ~ "Up",
#         logFC < 0 ~ "Down",
#         TRUE ~ "Neutral"
#     )) |>
#     filter(direction %in% c("Up", "Down"))
# 
# nrow(overlaps_df) # [1] 521984
# 
# #===========================================================================
# 
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
# # Define strategies
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



