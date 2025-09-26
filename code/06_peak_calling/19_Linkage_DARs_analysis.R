########################################################################
## Compute LinkPeaks ∩ DARs by cellType 
## 
## Authors. CSC
## Date. Sep 24, 2025
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
## Note. Seurat objects were created with module load conda_R/4.3.x
########################################################################

# library("Seurat")
# library("Signac")
# library("BSgenome.Hsapiens.UCSC.hg38")
# library("GenomicRanges") 
library("purrr")
library("ggvenn")
library("tidyverse")
library("tidyr")
library("stringr")
library("here")

#===============================================================================
# resolution_level = "Mid"      # 18 cell-types
#===============================================================================

## setup variable names
p_met = "spearman"
w_size_num <- 5e5 
w_size = "5e5" 
resolution_level = "Mid" 

## FDR and signif scores
# FDR_thresh = 0.1 # In link peaks was used a FDR_thresh = 0.2 
lfc_thresh = 0.2

## Set directory names
inputCSV_Links_Dir <- here(
    "processed-data",
    "06_peak_calling",
    #"13_pseudobulk_LinkPeaks_MACS2_split_ct",
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
lst_peak_files


# bind all files 
list_of_df <- lapply(here(inputCSV_Links_Dir, lst_peak_files), read.csv)
# ## Testing - using all links:
# tmp_lst_peak_files = c(lst_peak_files[2], lst_peak_files[18])
# list_of_df <- lapply(here(inputCSV_Links_Dir, tmp_lst_peak_files), read.csv)

combined_data <- bind_rows(list_of_df)

message("Processing ", length(list_of_df), " LinkPeaks files for merged-peaks dataset")
message(nrow(combined_data), " total links") 

head(combined_data)


#===========================================================================
## plot basic barplot to visualize number of links by cellType

# get number of links (rows) for each cluster
cluster_counts <- combined_data |>
    count(cluster, sort = TRUE)

# 2. Create the bar plot using ggplot2
# - We use aes() to map the 'cluster' column to the x-axis and the 'n' (count) column to the y-axis.
# - geom_col() creates the bars.
# - labs() adds a title and improves the axis labels.
# - theme_minimal() provides a clean, minimalist plot theme.
# - coord_flip() flips the x and y axes for better readability, especially with long cluster names.
g1 <- ggplot(cluster_counts, aes(x = reorder(cluster, n), y = n)) +
    geom_col(fill = "steelblue") +
    geom_text(aes(label = n), hjust = -0.2, size = 3) +
    labs(
        title = "Number of Links by cell-type",
        x = "Cluster",
        y = "Number of Links"
    ) +
    theme_minimal() +
    coord_flip()

#print(g1)

f_name <- paste0(resolution_level, "_level_pb_link_peaks_barplot.pdf")
ggsave(here(plotDir, f_name),
       combined_plot, width = 7, height = 7)


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
# parse DARs results for each cellType

DARS_signif_df_lst = list()

# test 2 thresholds
FDR_thresh = c(0.1, 0.2)

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
DARs_results_all <- DARs_results_all %>%
    mutate(FDR_threshold = factor(FDR_threshold, 
                                  labels = paste0("FDR", FDR_thresh)))

unique(DARs_results_all$FDR_threshold)
DARs_results_all |> head()


#===========================================================================
## plot basic barplot to visualize number of DARs by cellType

# Make FDR_threshold is a factor for facet labels
DARs_results_all <- DARs_results_all |>
    mutate(FDR_threshold = as.factor(FDR_threshold))

# count DARs
DARs_counts <- DARs_results_all |>
    count(cell_type, FDR_threshold, sort = TRUE)

g2 <- ggplot(DARs_counts, aes(x = reorder(cell_type, n), y = n)) +
    geom_col(fill = "steelblue") +
    geom_text(aes(label = n), hjust = -0.2, size = 3) +
    facet_wrap(~ FDR_threshold, ncol = 1, scales = "free_y",
               labeller = labeller(FDR_threshold = label_both)) +
    labs(
        title = "Number of Significant DARs by Cell Type",
        x = "Cell Type",
        y = "Number of Significant DARs"
    ) +
    theme_minimal() +
    coord_flip()

g2

head(combined_data$peak)
head(DARs_results_all$peak_id)

#===========================================================================
# # Prepare as a named list
# venn_list <- list(
#     LinkPeaks = unique(combined_data$peak),
#     DARs = unique(DARs_results_all$peak_id)
# )
# 
# ggvenn(venn_list,
#        fill_color = c("skyblue", "orange"),
#        stroke_size = 0.5,
#        set_name_size = 4)


#===========================================================================
## Parse csv DAR files and plot Volcano and ViolinPlot for the 5 "Up/Down" DARs by cellType

# # empty list to store the plots
# barPlot_list <- list()
# 
# for (ct_DARs in lst_DARs_cvs) {
#     # ct_DARs = lst_DARs_cvs[2] # "Endo"
#     # ct_DARs = lst_DARs_cvs[18] 
#     
#     message("Processing:\n", ct_DARs)
#     
#     # load DARs
#     DAR_peaks_cvs <- here(input_cvsDir, ct_DARs)
#     if (file.exists(DAR_peaks_cvs)) {
#         DARs_df <- read.csv(DAR_peaks_cvs)
#         message("File loaded!")
#     } else {
#         stop(paste("File not found:", DAR_peaks_cvs))
#     }
#     #colnames(DARs_df)
#     
#     # get cluster name. ge. "Endo"
#     clust_name <- sub("^voomlmFit_DAR_peaks_ALL_in_(.*)\\.csv$", "\\1", ct_DARs)
#     
#     ## create plot for FDR_thr (s)     
#     purrr::map(FDR_thr, ~ {
#         create_volcano(
#             DARs_df, 
#             clust_name,
#             FDR_thr = .x,   # current FDR
#             lfc_thresh = lfc_thresh,
#             plot_Dir
#         )
#     })
#     
#     barPlot_list[[clust_name]] <- create_barPlot_significant_DARs(DARs_df, clust_name, FDR_thr, lfc_thresh)
#     
#     # Identify top up/down peaks based on FDR
#     fdr_cols <- grep("^fdr_", colnames(res_enrich), value = TRUE)
#     
#     # cluster-specific significant peaks
#     res_sig_ct <- res_enrich |>
#         # filter(.data[[paste0("fdr_", clus)]] < FDR_thr) |>
#         filter(.data[[paste0("fdr_", clus)]] < FDR_thr,
#                abs(.data[[paste0("logFC_", clus)]]) > lfc_thresh) |>
#         mutate(
#             logFC_ct = .data[[paste0("logFC_", clus)]],
#             direction = case_when(
#                 logFC_ct >  0 ~ "Up",    # opening
#                 logFC_ct <  0 ~ "Down",  # closing
#                 TRUE ~ "NS"              # should not occur if you filtered
#             )
#         )
#     
#     if (nrow(res_sig_ct) > 0) {
#         f_name <- here(output_Dir, paste0("voomlmFit_DAR_peaks_", clus, ".csv"))
#         write.csv(res_sig_ct, f_name, row.names = FALSE)
#         message("Enrichment statistics saved [", clus, "]")
#     }
#     
# }
# 
# 
# # Plot barPlots
# if (length(barPlot_list) > 0) {
#     
#     pdf(here(plot_Dir, "BarPlots_ALL_cellTypes_voomLmFit.pdf"), width = 10, height = 10)
#     
#     # Remove redundant y-axis labels
#     barPlot_list_clean <- lapply(barPlot_list, function(p) {
#         p + ylab(NULL) + theme(legend.position = "none")
#     })
#     
#     combined_plot <- wrap_plots(barPlot_list_clean, ncol = 5)
#     
#     title_name <- paste0("Number of Up/Down Enriched DARs by CellType (FDR < ", FDR_thr, ")") 
#     final_plot <- combined_plot +
#         plot_annotation(
#             title = title_name,
#             theme = theme(
#                 plot.title = element_text(size = 12, hjust = 0.5),
#                 axis.title.y = element_text(size = 9)
#             )
#         ) # +
#     #labs(tag = "Number of DARs") +
#     #theme(plot.tag = element_text(angle = 90), plot.tag.position = "left")
#     
#     print(final_plot)
#     
#     dev.off()
#     
# }
