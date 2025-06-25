########################################################################
## Anlayze Hb Peaks
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

inputRDS_Dir <- here(
    "processed-data",
    "06_peak_calling"
)
plot_Dir <- here(
    "plots",
    "06_peak_calling",
    "05_analyze_hb_peaks"
)
if (!dir.exists(plot_Dir)) { dir.create(plot_Dir) }


message("Loading all peaks found")

f_name <- here(inputRDS_Dir, "all_peaks_by_cluster.csv") 
peaks_df <- read.csv(f_name, header = TRUE)

# extract all unique cluster labels listed in the peak_called_in column of peaks_df
all_clusters <- unlist(strsplit(peaks_df$peak_called_in, ","))
unique_clusters <- sort(unique(trimws(all_clusters)))  # trimws removes any leading/trailing spaces
unique_clusters
length(unique_clusters)

# Define target clusters: Habenula
target_clusters <- unique_clusters[grepl("MHb|LHb", unique_clusters)]
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

# Preview
message("Number of peaks in Hb clusters: ", nrow(peaks_hb_subset))
# Number of peaks in Hb clusters: 100975

head(peaks_hb_subset)
# seqnames  start    end width strand
#     1      chr1 191217 191620   404      *
#     3      chr1 629810 630397   588      *
#     5      chr1 633694 634120   427      *
#     6      chr1 778314 779306   993      *
#     9      chr1 819789 820144   356      *
#     10     chr1 827078 827686   609      *
#     peak_called_in
# 1                                                                                                                                                                                                C.26.DD_OPC,C.25.undetermined,C.19.DD_Inhib.Thal,C.02.DD_Oligo,C.16.DD_MHb,C.24.DD_LHb,C.13.no-match,C.09.undetermined,C.07.DD_MHb,C.30.DD_LHb,C.01.undetermined,C.05.DD_LHb,C.03.undetermined,C.38.DD_Inhib.Thal,C.22.undetermined,C.08.undetermined,C.20.DD_Astrocyte,C.06.DD_Excit.Thal,C.18.DD_LHb,C.28.DD_Inhib.Thal,C.12.undetermined,C.23.DD_LHb,C.04.undetermined,C.17.DD_Excit.Thal,C.10.DD_MHb,C.41.DD_Microglia,C.14.DD_MHb,C.21.DD_Astrocyte,C.34.DD_Oligo
# 3                                                                                              C.20.DD_Astrocyte,C.30.DD_LHb,C.16.DD_MHb,C.08.undetermined,C.15.DD_Excit.Thal,C.02.DD_Oligo,C.01.undetermined,C.23.DD_LHb,C.06.DD_Excit.Thal,C.19.DD_Inhib.Thal,C.32.undetermined,C.10.DD_MHb,C.18.DD_LHb,C.33.DD_LHb,C.21.DD_Astrocyte,C.14.DD_MHb,C.12.undetermined,C.13.no-match,C.04.undetermined,C.26.DD_OPC,C.05.DD_LHb,C.40.DD_LHb,C.17.DD_Excit.Thal,C.25.undetermined,C.39.DD_Inhib.Thal,C.07.DD_MHb,C.11.DD_MHb,C.36.DD_MHb,C.37.undetermined,C.28.DD_Inhib.Thal,C.22.undetermined,C.24.DD_LHb,C.29.DD_Endo,C.34.DD_Oligo,C.09.undetermined,C.27.DD_Microglia
# 5  C.21.DD_Astrocyte,C.29.DD_Endo,C.13.no-match,C.06.DD_Excit.Thal,C.41.DD_Microglia,C.20.DD_Astrocyte,C.11.DD_MHb,C.37.undetermined,C.02.DD_Oligo,C.30.DD_LHb,C.25.undetermined,C.28.DD_Inhib.Thal,C.01.undetermined,C.05.DD_LHb,C.32.undetermined,C.33.DD_LHb,C.36.DD_MHb,C.14.DD_MHb,C.26.DD_OPC,C.27.DD_Microglia,C.34.DD_Oligo,C.15.DD_Excit.Thal,C.23.DD_LHb,C.17.DD_Excit.Thal,C.19.DD_Inhib.Thal,C.08.undetermined,C.18.DD_LHb,C.09.undetermined,C.16.DD_MHb,C.24.DD_LHb,C.10.DD_MHb,C.12.undetermined,C.22.undetermined,C.07.DD_MHb,C.40.DD_LHb,C.31.DD_Excit.Thal,C.38.DD_Inhib.Thal,C.35.undetermined,C.39.DD_Inhib.Thal,C.04.undetermined,C.03.undetermined
# 6                                                                                                                               C.26.DD_OPC,C.24.DD_LHb,C.19.DD_Inhib.Thal,C.07.DD_MHb,C.15.DD_Excit.Thal,C.08.undetermined,C.05.DD_LHb,C.02.DD_Oligo,C.09.undetermined,C.20.DD_Astrocyte,C.04.undetermined,C.01.undetermined,C.03.undetermined,C.06.DD_Excit.Thal,C.21.DD_Astrocyte,C.28.DD_Inhib.Thal,C.17.DD_Excit.Thal,C.16.DD_MHb,C.12.undetermined,C.13.no-match,C.37.undetermined,C.25.undetermined,C.18.DD_LHb,C.22.undetermined,C.23.DD_LHb,C.30.DD_LHb,C.10.DD_MHb,C.14.DD_MHb,C.27.DD_Microglia,C.11.DD_MHb,C.31.DD_Excit.Thal,C.33.DD_LHb,C.41.DD_Microglia


##==============================================================================
## previews plots - EDA 

# Count Occurrences
cluster_peak_counts <- sort(table(all_clusters), decreasing = TRUE)
cluster_df <- as.data.frame(cluster_peak_counts)
colnames(cluster_df) <- c("cluster_ann", "n_peaks")

## (1) Plot all clusters by peak frequency 
p1 <- ggplot(cluster_df, aes(x = reorder(cluster_ann, -n_peaks), y = n_peaks)) +
    geom_bar(stat = "identity", fill = ifelse(cluster_df$cluster_ann %in% target_clusters, "#FF6F61", "grey")) +
    theme_minimal() +
    labs(title = "Number of peaks called per cluster",
         #x = "Cluster",
         y = "Number of Peaks",
         caption = paste("Number of peaks found:", nrow(peaks_hb_subset))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.6),
          axis.title.x = element_blank(),
          plot.caption = element_text(size = 10, hjust = 0))
#p1
plot_name <- here(plot_Dir, "peaks_frequency_by_cluster.png")
ggsave(plot_name, plot = p1, width = 10, height = 6, dpi = 300)




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

