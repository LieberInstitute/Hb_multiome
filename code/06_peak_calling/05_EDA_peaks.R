########################################################################
## Make visualizations for Peaks found on WNN clusters
##
## Authors. CSC
## Date. June 25, 2025
## Recommended resources on interactive mode: srun --pty --mem=40GB --x11 bash
## Resources:
## https://stuartlab.org/signac/articles/peak_calling
########################################################################

library("Signac")
library("ggplot2")
library("dplyr")
library("forcats")
library("here")

# Check/create directories
## clusters renamed for Spatial-Registration on Visium project

resolution_level = c("Fine", "Broad", "Mid")
# resolution_level = "Fine" # 42 clusters
# resolution_level = "Broad"  # 8 cell-types
# resolution_level = "Mid" # 8 cell-types

inputCSV_Dir <- here(
    "processed-data",
    "06_peak_calling",
    "01_call_peaks_MACS2"
)
plot_Dir <- here(
    "plots",
    "06_peak_calling",
    "05_EDA_peaks"
)
outputCSV_Dir  <- here(
    "processed-data",
    "06_peak_calling",
    "05_EDA_peaks"
)
if (!dir.exists(outputCSV_Dir)) { dir.create(outputCSV_Dir) }
if (!dir.exists(plot_Dir)) { dir.create(plot_Dir) }



for (clust_level in resolution_level) {

    message("Processing peaks for resolution:\n", clust_level)
    # suffix to save plots and csv
    f_sufix <- paste0(".resolution.", clust_level)
    f_sufix

    message("Loading all peaks found")
    
    f_name <- paste0("macs_peaks_", clust_level, "_resolution.csv")
    f_name <- here(inputCSV_Dir, f_name) 
    peaks_df <- read.csv(f_name, header = TRUE)
    
    message("Number of peaks found: ", nrow(peaks_df))
    
    message("Define target clusters: Habenula")
    
    # extract all unique cluster labels listed in the peak_called_in column of peaks_df
    all_clusters <- unlist(strsplit(peaks_df$peak_called_in, ","))
    unique_clusters <- sort(unique(trimws(all_clusters)))  # trimws removes any leading/trailing spaces
    unique_clusters
    
    message("Processing ", length(unique_clusters), " total clusters")
    
    target_clusters <- unique_clusters[grepl("MHb|LHb", unique_clusters)]
    
    message("Hb total: ", length(target_clusters))
    
    target_clusters
    # [1] "C.05.DD_LHb" "C.07.DD_MHb" "C.10.DD_MHb" "C.11.DD_MHb" "C.14.DD_MHb"
    # [6] "C.16.DD_MHb" "C.18.DD_LHb" "C.23.DD_LHb" "C.24.DD_LHb" "C.30.DD_LHb"
    # [11] "C.33.DD_LHb" "C.36.DD_MHb" "C.40.DD_LHb"
     
    # Match any target cluster in the peak_called_in field
    matches <- sapply(target_clusters, function(cl) {
        # use \\b (word boundary) to avoid partial matches, as the column contains ","
        grepl(paste0("\\b", cl, "\\b"), peaks_df$peak_called_in)
    })
    # counts how many clusters matched per peak to keep only peaks where at least one target cluster appears
    matched_rows <- rowSums(matches) > 0
    # Retains peak called in that cluster (TRUE)
    # Subset matched rows
    peaks_hb_subset <- peaks_df[matched_rows, ]
     
    message("Number of peaks in Hb clusters: ", nrow(peaks_hb_subset))
    # Number of peaks in Hb clusters: 100975
     
    head(peaks_hb_subset)[1:4]
    macs_peaks_Broad_resolution.csv
    f_name <- paste0("habenula_macs_peaks_", clust_level ,"_resolution.csv")
    write.csv(peaks_hb_subset, here(outputCSV_Dir, f_name), row.names = FALSE)
    
    message("Saved subsetted peaks for habenula clusters")
    
}

# # extract all unique cluster labels listed in the peak_called_in column of peaks_df
# all_clusters <- unlist(strsplit(peaks_df$peak_called_in, ","))
# unique_clusters <- sort(unique(trimws(all_clusters)))  # trimws removes any leading/trailing spaces
# unique_clusters
# length(unique_clusters)
# 
# target_clusters <- unique_clusters[grepl("MHb|LHb", unique_clusters)]
# target_clusters
# # [1] "C.05.DD_LHb" "C.07.DD_MHb" "C.10.DD_MHb" "C.11.DD_MHb" "C.14.DD_MHb"
# # [6] "C.16.DD_MHb" "C.18.DD_LHb" "C.23.DD_LHb" "C.24.DD_LHb" "C.30.DD_LHb"
# # [11] "C.33.DD_LHb" "C.36.DD_MHb" "C.40.DD_LHb"
# 
# # Match any target cluster in the peak_called_in field
# matches <- sapply(target_clusters, function(cl) {
#     # use \\b (word boundary) to avoid partial matches, as the column contains ","
#     grepl(paste0("\\b", cl, "\\b"), peaks_df$peak_called_in)
# })
# # counts how many clusters matched per peak to keep only peaks where at least one target cluster appears
# matched_rows <- rowSums(matches) > 0
# # Retains peak called in that cluster (TRUE)
# 
# # Subset matched rows
# peaks_hb_subset <- peaks_df[matched_rows, ]
# 
# message("Number of peaks in Hb clusters: ", nrow(peaks_hb_subset))
# # Number of peaks in Hb clusters: 100975
# 
# head(peaks_hb_subset)[1:4]

# write.csv(peaks_hb_subset, here(outputCSV_Dir, "hb_peaks.csv"), row.names = FALSE)
# message("Saved subsetted peaks for target clusters")

        
##==============================================================================

message("Make plots for EDA")

# Count Occurrences
cluster_peak_counts <- sort(table(all_clusters), decreasing = TRUE)
cluster_df <- as.data.frame(cluster_peak_counts)
colnames(cluster_df) <- c("cluster_ann", "n_peaks")

## (1) Plot all clusters by peak frequency 
p1 <- ggplot(cluster_df, aes(x = reorder(cluster_ann, -n_peaks), y = n_peaks)) +
    geom_bar(stat = "identity", fill = ifelse(cluster_df$cluster_ann %in% target_clusters, "#FF6F61", "grey")) +
    theme_minimal() +
    labs(title = "Number of peaks per cluster",
         #x = "Cluster",
         y = "Number of Peaks",
         caption = paste("Total peaks (raw):", nrow(peaks_hb_subset))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.6),
          axis.title.x = element_blank(),
          plot.caption = element_text(size = 10, hjust = 0),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA))

plot_name <- here(plot_Dir, "peaks_frequency_by_cluster.png")
ggsave(plot_name, plot = p1, width = 10, height = 6, dpi = 300)

message("Plot for number of peaks called per cluster done!")

##==============================================================================
# ## (2) Plot by regions of interest: Hb vs Others
# 
# # custom order, Hb clusters first, then the rest
# remaining_clusters <- setdiff(cluster_df$cluster_ann, target_clusters)
# custom_order <- c(target_clusters, sort(remaining_clusters))
# cluster_df$cluster_ann <- factor(cluster_df$cluster_ann, levels = custom_order)
# 
# p2 <- ggplot(cluster_df, aes(x = cluster_ann, y = n_peaks)) +
#     geom_bar(stat = "identity", fill = ifelse(cluster_df$cluster_ann %in% target_clusters, "#FF6F61", "grey")) +
#     theme_minimal() +
#     labs(title = "Number of peaks called per cluster and region",
#          y = "Number of Peaks",
#          caption = paste("Total peaks (raw):", nrow(peaks_hb_subset))) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 0.6),
#           axis.title.x = element_blank(),
#           plot.caption = element_text(size = 10, hjust = 0),
#           panel.background = element_rect(fill = "white", color = NA),
#           plot.background = element_rect(fill = "white", color = NA))
# 
# plot_name <- here(plot_Dir, "peaks_frequency_by_region.png")
# ggsave(plot_name, plot = p2, width = 10, height = 6, dpi = 300)
# 
# message("Plot for number of peaks called per cluster and region done!")


##==============================================================================
## (3) Plot by cluster sub-region: MHb, LHb, Others

# Assign region
cluster_df <- cluster_df |>
    mutate(region = case_when(
        grepl("MHb", cluster_ann) ~ "MHb",
        grepl("LHb", cluster_ann) ~ "LHb",
        TRUE ~ "Other"
    ))
# Assign panel width weights by number of bars
cluster_df <- cluster_df |>
    mutate(panel = case_when(
        region == "MHb" ~ "Hb",
        region == "LHb" ~ "Hb",
        TRUE ~ "Other"
    ))
unique(cluster_df$region)
# [1] "Other" "LHb"   "MHb"  
unique(cluster_df$panel)
# [1] "Other" "Hb" 

# Create separate order vectors for Hb and Other panels to avoid interleved bars 
cluster_df <- cluster_df |>
    mutate(
        cluster_ann = as.character(cluster_ann),
        order_rank = case_when(
            panel == "Hb" ~ match(cluster_ann, target_clusters),
            panel == "Other" ~ match(cluster_ann, sort(remaining_clusters))
        )
    ) |>
    arrange(panel, order_rank) |>
    mutate(cluster_ann = factor(cluster_ann, levels = unique(cluster_ann)))

p3 <- ggplot(cluster_df, aes(x = cluster_ann, y = n_peaks)) +
    geom_bar(stat = "identity", aes(fill = region)) +
    scale_fill_manual(values = c("LHb" = "darkblue", "MHb" = "#FF7F50", "Other" = "grey")) +
    facet_grid(. ~ panel, scales = "free_x", space = "free_x") +
    theme_minimal() +
    labs(title = "Number of peaks called per cluster and Hb sub-region",
         y = "Number of Peaks",
         caption = paste("Number of peaks found:", nrow(peaks_hb_subset))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.6),
          axis.title.x = element_blank(),
          legend.position = "none",
          plot.caption = element_text(size = 10, hjust = 0),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA)
    )
#p3
plot_name <- here(plot_Dir, "peaks_frequency_by_subregion.png")
ggsave(plot_name, plot = p3, width = 10, height = 6, dpi = 300)

message("Plot for number of peaks called per cluster and Hb region done!")


# library("slurmjobs")
# job_single(
#     "05_EDA_peaks",
#     create_shell = TRUE,
#     partition = "katun",
#     memory = "60G",
#     cores = 2,
#     logdir = "logs",
#     command = "05_EDA_peaks.R",
#     #create_logdir = TRUE
# )


##==============================================================================
# ## Preparing to compute pair peaks comparison in the Habenula clusters
# 
# all_cluster_IDs <- levels(SeuratOBJ)
# ## extract Hb clusters
# Hb_cluster_IDs <- cluster_IDs[grepl("MHb|LHb", all_cluster_IDs)]
# Hb_cluster_IDs
# # [1] "C.05.DD_LHb" "C.07.DD_MHb" "C.10.DD_MHb" "C.11.DD_MHb" "C.14.DD_MHb"
# # [6] "C.16.DD_MHb" "C.18.DD_LHb" "C.23.DD_LHb" "C.24.DD_LHb" "C.30.DD_LHb"
# # [11] "C.33.DD_LHb" "C.36.DD_MHb" "C.40.DD_LHb"
# pairwise_combinations <- combn(Hb_cluster_IDs, 2, simplify = FALSE)
# pairwise_combinations
# pairwise_df <- do.call(rbind, pairwise_combinations)
# colnames(pairwise_df) <- c("cluster_1", "cluster_2")
# pairwise_df <- as.data.frame(pairwise_df)
# pairwise_df
# 
# 
# for (clust_p in pairwise_df) {
#     da_peaks <- FindMarkers(
#         object = seurat_atac,
#         ident.1 = clust_p["cluster1"],
#         ident.2 = clust_p["cluster2"],
#         test.use = 'LR',
#         min.pct = 0.05
#     )
# }

