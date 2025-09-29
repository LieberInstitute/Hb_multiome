########################################################################
##                   Compute LinkPeaks ∩ DARs by cellType 
## ##################### DAR → Links perspective ## ####################
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
## Load LinkPeaks results

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
## plot basic barplot of number of links by cellType

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

f_name <- paste0(resolution_level, "_level_pb_link_peaks_barplot_ct_FDR", FDR, ".pdf")
ggsave(here(plotDir, f_name),
       g1, width = 7, height = 7)



##==============================================================================
## Load DARs results

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
# parse DARs results for each cellType at specific FDR thr

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
                fdr = as.numeric(fdr),  # ✅ keep actual fdr values
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
## plot basic barplot to visualize number of DARs by cellType

message("Making intersection plots  DAR → Links perspective .... ")

# Make FDR_threshold is a factor for facet labels
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
        filename = file.path(plotDir, paste0(resolution_level, "_level_pb_DARs_", .y, ".pdf")),
        plot = .x,
        width = 7,
        height = 7
    )
)

## overlay both thr for comparison purposes 
g3 <- ggplot(DARs_counts, aes(x = reorder(cell_type, n), y = n, fill = FDR_threshold)) +
    geom_col(alpha = 0.5, position = "identity") +  # <- overlay instead of dodge
    labs(
        title = paste("DARs by Cell Type"),
        subtitle = "FDR thresholds",
        x = "Cell Type",
        y = "Number of DARs",
        fill = "FDR Threshold"
    ) +
    theme_minimal() +
    coord_flip()

f_name <- paste0(resolution_level, "_level_pb_DARs_overlay.pdf")
ggsave(here(plotDir, f_name),
       g3, width = 7, height = 7)



#===========================================================================
## Make intersection btw links and DARs

# inspect data
nrow(linkPeaks_results_all)
colnames(linkPeaks_results_all)
head(linkPeaks_results_all$peak) # join id key
nrow(DARs_results_all)
colnames(DARs_results_all)
head(DARs_results_all$peak_id) # join id key

## convert genomic ranges and do a proper overlap

# Convert LinkPeaks to GRanges
gr_links <- with(linkPeaks_results_all, 
                 GRanges(seqnames = seqnames,
                         ranges   = IRanges(start, end),
                         peak_id  = peak,
                         score    = score,
                         gene     = gene,
                         cluster  = cluster,
                         distance = distance,
                         tier     = tier))

# Convert DARs to GRanges

## only use DARs at FDR02
nrow(DARs_results_all) # [1] 582590
table(DARs_results_all$FDR_threshold)
DARs_FDRX <- DARs_results_all |>
    filter(FDR_threshold == paste0("FDR", FDR)) 
nrow(DARs_FDRX) # [1] 325580

# split the peak_id "chr-start-end" into seqnames, start, end
colnames(DARs_FDRX)
DARs_split <- tidyr::separate(DARs_FDRX, peak_id, into = c("seqnames", "start", "end"), sep = "-") |>
    mutate(start = as.integer(start),
           end   = as.integer(end))

gr_dars <- with(DARs_split,
                GRanges(seqnames = seqnames,
                        ranges   = IRanges(start, end),
                        peak_id  = paste(seqnames, start, end, sep = "-"),
                        cell_type = cell_type,
                        FDR_threshold = FDR_threshold,
                        logFC = logFC,
                        fdr = fdr,
                        is_significant = is_significant))

## inspect
head(gr_links)
head(gr_dars)

## find overlaps
# findOverlaps() reports all pairs of ranges that overlap
# If one LinkPeaks region overlaps many DARs regions, you’ll get multiple rows for the same LinkPeak
hits <- findOverlaps(gr_links, gr_dars)

# Combine metadata from both sides
overlaps_df <- data.frame(
    as.data.frame(mcols(gr_links)[queryHits(hits), ]),
    as.data.frame(mcols(gr_dars)[subjectHits(hits), ])
)

nrow(overlaps_df)
head(overlaps_df)


overlap_counts <- overlaps_df |>
    count(cell_type, FDR_thresh, sort = TRUE)

overlap_counts

#===========================================================================
## control overlaps: Unique overlapping LinkPeaks or DARs

# LinkPeaks overlap at least one DAR
unique_links <- unique(mcols(gr_links)$peak_id[queryHits(hits)])
length(unique_links)   # number of distinct LinkPeaks overlapping DARs

head(unique_links)

# only distinct DAR 
unique_dars <- unique(mcols(gr_dars)$peak_id[subjectHits(hits)])
length(unique_dars)    # number of distinct DARs overlapping LinkPeaks

## Count overlaps per region
link_counts <- as.data.frame(table(mcols(gr_links)$peak_id[queryHits(hits)]))
nrow(link_counts) # [1] 7089
head(link_counts)
# Var1 Freq
# 1 chr1-100213166-100213615   41
# 2 chr1-100231650-100231970   44
# 3 chr1-100265432-100266863    1


#===========================================================================
## Up vs Down DAR-Links count per cell type

# Define Up vs Down
# peak_id → the LinkPeak ID (from gr_links)
# peak_id.1 → the DAR peak ID (from gr_dars)
overlaps_df <- overlaps_df |>
    mutate(direction = case_when(
        logFC > 0 ~ "Up",
        logFC < 0 ~ "Down",
        TRUE ~ "Neutral"
    ))

########## Get the counts per cell type ##########

## This is a pairwise overlap count (DAR–Link edges)
# A DAR overlaps multiple LinkPeaks, it is counted once per overlap
DAR_link_summary <- overlaps_df |>
    group_by(cell_type, direction) |>
    summarise(n = n(), .groups = "drop")
DAR_link_summary

## counts unique DAR peaks per cell type & direction

# A DAR overlapping 5 Links is only counted once
summary_dar <- overlaps_df |>
    group_by(cell_type, direction) |>
    summarise(
        n_DARs = n_distinct(peak_id.1),   # distinct DAR peaks
        .groups = "drop"
    )
summary_dar
# A Link overlapping 5 DARs is only counted once
summary_links <- overlaps_df %>%
    group_by(cluster, direction) %>%
    summarise(
        n_Links = n_distinct(peak_id),   # distinct Link peaks
        .groups = "drop"
    )
summary_links

############

if (nrow(DAR_link_summary) > 0) { 
    
    f_name <- here(processedDir, paste0("summary_intersected_counts_per_ct_FDR", FDR, ".csv"))
    write.csv(DAR_link_summary, f_name, row.names = FALSE)
    message("Summary saved!")
    
} else (
    
    stop("None overlaps found!")
    
)


#===========================================================================

prepare_data_to_plot <- function(
        df_summary
) {
    # make a diverging plot (Up → right, Down → left), flip the sign of n for Down peaks
    df_summary_div <- df_summary |> 
        mutate(n_signed = ifelse(direction == "Down", -n, n))
    # order cell types by total number of overlaps
    cell_totals <- df_summary |>
        group_by(cell_type) |>
        summarise(total = sum(n), .groups = "drop")
    # Join back to signed summary
    df_summary_div <- df_summary |>
        mutate(n_signed = ifelse(direction == "Down", -n, n)) |> 
        left_join(cell_totals, by = "cell_type")
    
    return(df_summary_div)
    
}

DAR_link_summary_div <- function(DAR_link_summary)

# # make a diverging plot (Up → right, Down → left), flip the sign of n for Down peaks
# DAR_link_summary_div <- DAR_link_summary |> 
#     mutate(n_signed = ifelse(direction == "Down", -n, n))
# # order cell types by total number of overlaps
# cell_totals <- DAR_link_summary |>
#     group_by(cell_type) |>
#     summarise(total = sum(n), .groups = "drop")
# # Join back to signed summary
# DAR_link_summary_div <- DAR_link_summary |>
#     mutate(n_signed = ifelse(direction == "Down", -n, n)) |> 
#     left_join(cell_totals, by = "cell_type")

    
make_div_prop_barplots <- function(
        
) {
    
    g4_overlap <- ggplot(DAR_link_summary_div, 
                         aes(x = reorder(cell_type, total), y = n_signed, fill = direction)) +
        geom_col() +
        # add labels outside bars
        geom_text(aes(label = abs(n_signed),
                      hjust = ifelse(direction == "Down", 1.1, -0.1)), # left for Down, right for Up
                  position = position_identity(),
                  size = 3,
                  color = "#3D3936") +
        scale_y_continuous(labels = abs, expand = expansion(mult = c(0.15, 0.15))) + # add padding
        labs(
            title = "DAR-Links Overlaps by cell type",
            subtitle = paste0("FDR = ", FDR),
            x = "Cell Type",
            y = "Number of overlapping DAR-Links",
            fill = "Direction"
        ) +
        theme_minimal() +
        coord_flip()
    
    
    # define factor levels by total counts
    cell_totals <- DAR_link_summary |>
        group_by(cell_type) |>
        summarise(total = sum(n), .groups = "drop")
    
    cell_order <- cell_totals |>
        arrange(total) |>
        pull(cell_type)
    
    DAR_link_summary$cell_type <- factor(DAR_link_summary$cell_type, levels = cell_order)
    DAR_link_summary_div$cell_type <- factor(DAR_link_summary_div$cell_type, levels = cell_order)
    
    g4b_overlap <- DAR_link_summary|>
        group_by(cell_type) |>
        mutate(prop = n / sum(n)) |>
        ggplot(aes(x = cell_type, y = prop, fill = direction)) +
        geom_col() +
        scale_y_continuous(labels = scales::percent) +
        labs(
            title = "Proportion of Up/Down DAR-Links",
            x = "Cell Type",
            y = "Proportion",
            fill = "Direction"
        ) +
        theme_minimal() +
        coord_flip()
    
    g4_combined <- g4_overlap | g4b_overlap
    
    f_name <- paste0(resolution_level, "_level_pb_Links_DARs_overaping_FDR", FDR, ".pdf")
    ggsave(here(plotDir, f_name),
           g4_combined, width = 10, height = 7)
    
}

#===========================================================================
## Venn Diagram
# library("ggvenn")
# venn_list <- list(
#     LinkPeaks = unique(linkPeaks_results_all$peak),
#     DARs = unique(DARs_results_all$peak_id)
# )
# 
# ggvenn(venn_list,
#        fill_color = c("skyblue", "orange"),
#        stroke_size = 0.5,
#        set_name_size = 4)

