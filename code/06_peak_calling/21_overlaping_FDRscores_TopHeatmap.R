########################################################################
## Explore FDR-Scores and Top overlapping (Heatmap)
## - (1) Heatmap summarize accessibility (logFC, FDR, directionality) 
##    → Subset by cell_type, because visual goal centers on cell-type–specific regulatory accessibility
##
## Authors. CSC
## Date. Oct 14, 2025
## Recommended resources on interactive mode: srun --pty --mem=20GB --x11 bash
########################################################################

library("pheatmap")
library("reshape2")
library("dplyr")
library("purrr")
library("ggplot2")
library("patchwork")
library("ggrepel")
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
    "06_peak_calling"
    #"20_explore_overlapping_linkpeaks_DARs"
    #"19_Linkage_DARs_analysis"
)
processedDir <- here(
    "processed-data",
    "06_peak_calling",
    "21_overlaping_FDRscores_TopHeatmap"
)
plotDir <- here(
    "plots",
    "06_peak_calling",
    "21_overlaping_FDRscores_TopHeatmap"
)

if (!dir.exists(processedDir)) { dir.create(processedDir) }
if (!dir.exists(plotDir)) { dir.create(plotDir) }

##==============================================================================

message("Loading Unique-Overlap with directionality ...")

## load raw overlaps and 2-shared ct overlaps 
overlap_df <- read.csv(here(inputCSV_Overlaps_Dir, "19_Linkage_DARs_analysis", 
                             paste0("Overlaps_LinkPeak_DARs_FDR", FDR, ".csv")))
shared_2ct_df <- read.csv(here(inputCSV_Overlaps_Dir, "20_explore_overlapping_linkpeaks_DARs", 
                             paste0("overlaps_summary_linkPeak_DARs_unique_ct_FDR", FDR, ".csv")))
                            

## filter uniques and 2-shared ct overlaps
all_unique_df <- overlap_df |>
    distinct(peak_id_links, .keep_all = TRUE)
nrow(all_unique_df) # 5603
head(all_unique_df)

## first, filter uniques
peaks_unique_cell_types <- shared_2ct_df |>
    filter(n_cell_types == 1) |>
    pull(peak_id)
length(peaks_unique_cell_types) # 1967

## Keep only those 1,967 unique LinkPeaks
unique_df <- all_unique_df |>
    filter(peak_id_links %in% peaks_unique_cell_types) |>
    mutate(
        overlap_type = "Unique",
        n_cell_types = 1,
    )

length(unique_df$peak_id_links) # 1967
table(unique_df$cell_type)
head(unique_df)
summary(unique_df)


## filter 2ct shared
peaks_shared_cell_types <- shared_2ct_df |>
    filter(n_cell_types == 2) |>
    pull(peak_id)
length(peaks_shared_cell_types) # 1238

# Keep only unique and 2-shared overlaps
unique_df2 <- all_unique_df |>
    filter(peak_id_links %in% peaks_shared_cell_types) |>
    mutate(
        overlap_type = "Shared",
        n_cell_types = 2,
    )

length(unique_df2$peak_id_links) 
table(unique_df2$cell_type)
head(unique_df2)
summary(unique_df2)

unique_df <- bind_rows(unique_df, unique_df2)
head(unique_df)
table(unique_df$overlap_type)
nrow(unique_df) # 3205

# save detailed unique and 2-shared overlaps with meta-data
f_name <- here(processedDir, paste0("overlaps_unique_2shared_ct_FDR", FDR, ".csv"))
write.csv(unique_df, f_name, row.names = FALSE)


## =============================================================================
## subset by cell_type and cluster, focus only on cis-links where both accessibility and expression specificity occur in the same cluster
# Conceptually: 
# cell_type → From the ATAC side: where the chromatin accessibility change occurs (DARs)
# cluster → From the RNA side: where the correlated gene was expressed or correlated in the LinkPeaks model


all_clusters <- unique(unique_df$cell_type)
# sort by MHh and LHb first 
clusters_sorted <- c(
    sort(grep("MHb|LHb", all_clusters, value = TRUE)), 
    sort(grep("MHb|LHb", all_clusters, value = TRUE, invert = TRUE))
)
clusters_sorted

# ## subset df by cell-Type / check point
# unique_df |> filter(cell_type=="LHb.2.7") |> nrow() #203
    
subset_cell_type <- function(unique_df, cluster_specific) {
    message("Subsetting [", cluster_specific, "] cell_type")
    subset_uniques <- unique_df |> filter(cell_type == cluster_specific)
    #message("Rows found: ", nrow(subset_uniques))
    message("Unique peaks: ", length(unique(subset_uniques$peak_id_links)))
    return(subset_uniques)
}


## =============================================================================
## Generate scatterplot

make_scattered_plot_dars_cc_real <- function(
        plot_data,
        top_genes,
        thr_CC,    # Thresh for Peak-Gene correlation
        thr_DAR,   # Thresh for Differential Accessibility Regions (DARs)
        thr_logFC, # Thresh for DARs Accessibility Direction
        fdr_cutoff = 0.2 # FDR used in both ds to define significant peaks 
) {

    clus_name <- unique(plot_data[["cell_type"]])
    
    message("Building plot for [", clus_name, "] (real FDR values)")
    
    top_hits <- plot_data |>
        filter(category %in% c("Active cCRE (+)", "Repressive cCRE (-)", "Neutral cCRE", "Linked OCR")) |>
        arrange(desc(abs(CCscore))) |> 
        head(top_genes)

    x_max <- min(1, max(plot_data$FDR_CC, na.rm = TRUE) * 1.05)
    y_max <- min(1, max(plot_data$fdr_dars, na.rm = TRUE) * 1.05)
    # define tick units
    axis_breaks <- seq(0, max(x_max, y_max), by = 0.05)
    axis_labels <- sprintf("%.2f", axis_breaks)
    
    g1 <- ggplot(plot_data, aes(x = CCscore, y = logFC, color = category)) +
        #geom_hline(yintercept = 0, linetype = "solid", color = "grey70") +
        #geom_hline(yintercept = c(-log2(1 + thr_DAR), log2(1 + thr_DAR)), linetype = "dashed", color = "darkgrey") +
        geom_hline(yintercept = c(-thr_logFC, thr_logFC), linetype = "dashed", color = "darkgrey") +
        geom_vline(xintercept = 0, linetype = "solid", color = "grey70") +
        geom_vline(xintercept = c(-thr_CC, thr_CC), linetype = "dashed", color = "darkgrey") +
        geom_point(alpha = 0.7, size = 1.6) +
        geom_text_repel(
            data = top_hits,
            aes(label = gene_name),
            size = 3,
            max.overlaps = top_genes
        ) +
        scale_color_manual(values = c(
            "Active cCRE (+)"  = "#E64B35FF",
            "Repressive cCRE (-)" = "#4DBBD5FF",
            "Neutral cCRE"     = "#CFB53B",
            "Linked OCR"      = "#00A087FF",
            "Unlinked OCR"    = "#3C5488FF",
            "Non-significant" = "lightgrey"
        )) +
        labs(
            title = paste(clus_name, " | Peak-Gene Correlation vs DARs"),
            subtitle = expression(
                paste("Spearman CC |", rho, "| > ", thr_CC,
                      " | FDR-DARs < ", thr_DAR, " & log2FC ", phantom(), pm, " ", thr_logFC)
            ),
            x = expression("Spearman CC |" * rho * "|"),
            y = expression(log[2] * "FC (DARs)")
        ) +
        theme_minimal(base_size = 12) +
        theme(
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold"),
            legend.position = "bottom",
            legend.title = element_blank()
        ) +
        plot_annotation(caption = paste("Uniques & 2-Shared CellTypes\nTop:", top_genes, "genes"))
    
    return(g1)
    
}

## subset the overlaps by cell_type
subsetted_df <- purrr::map_df(
    clusters_sorted,
    ~ subset_cell_type(unique_df, .x)
)

## verification
table(subsetted_df$cell_type)
summary(subsetted_df$cell_type)
head(subsetted_df)
subsetted_list_df <- split(subsetted_df, subsetted_df$cell_type)
names(subsetted_list_df)
length(subsetted_list_df) # 15 ct
#table(subsetted_list_df[[8]]["overlap_type"])


## make scattered plots using real values (non-normalized)

## Categorize peaks

## Set thresholds for LinkPeaks (CC) and DARs
# cluster_specific = "LHb.4"
top_genes_scattered_plt = 20
thr_CC = 0.3  # correlation strength
thr_DAR = 0.1 # fdr_dars / accessibility significance
thr_logFC = 0
thr_fdr = 0.2

plot_data_list <- purrr::map(subsetted_list_df, ~ .x |> 
    mutate(
        sig_CC  = abs(CCscore) > thr_CC,   # significant correlation
        sig_DAR = fdr_dars < thr_DAR,      # significant accessibility
        category = case_when(
            sig_CC & sig_DAR & logFC > thr_logFC  ~ "Active cCRE (+)",      # positively correlated & accessible
            sig_CC & sig_DAR & logFC < -thr_logFC  ~ "Repressive cCRE (-)", # negatively correlated & less accessible
            # Neural category was used to identify "Significantly linked, Significant DAR, but logFC is too small", when set log_FC=0 we do not need it any more
            # sig_CC & sig_DAR                ~ "Neutral cCRE",        # when other conditions not met  
            sig_CC & !sig_DAR               ~ "Linked OCR",          # correlated, not DAR
            !sig_CC & sig_DAR               ~ "Unlinked OCR",        # DAR, no correlation
            TRUE                            ~ "Non-significant"      # everything else
        )
    )
)

## Build a summary 

plot_data_full_df <- bind_rows(plot_data_list, .id = "source_df")

message("========= Summary of candidate RE by category ============\n")
# head(plot_data_full_df)
table(plot_data_full_df$cell_type)
as.data.frame(table(plot_data_full_df$category))
# Var1 Freq
# 1     Active cCRE (+)  706
# 2          Linked OCR  420
# 3     Non-significant  354
# 4 Repressive cCRE (-)  657
# 5        Unlinked OCR 1068

sum(table(plot_data_full_df$category)) #3205

hb_related_df <- plot_data_full_df |>
    filter(grepl("MHb|LHb", cell_type))

total_hb_related <- nrow(hb_related_df)
table(hb_related_df$category)

hb_related_signif_df <- hb_related_df |>
    filter(category %in% c("Active cCRE (+)", "Repressive cCRE (-)", "Linked OCR", "Unlinked OCR"))

total_hb_related_signif <- nrow(hb_related_signif_df) 

message("Total Hb related [thr_CC=", thr_CC, 
        " & thr_DAR=", thr_DAR,
        " & thr_logFC=", thr_logFC, 
        "]: ", total_hb_related)

message("Total Hb related significant: ", total_hb_related_signif)

# save overlapping with classification
f_name <- here(processedDir, paste0("overlaps_linkPeak_DARs_classified_thr_CC", thr_CC, "_thr_DAR", thr_DAR, ".csv"))
write.csv(plot_data_full_df, f_name, row.names = FALSE)

message("Saved linkPeak_DARs categories!")

scattered_plt_cell_type_real_values <- purrr::map(
    plot_data_list,
    ~ make_scattered_plot_dars_cc_real(
        .x, 
        top_genes_scattered_plt,
        thr_CC, thr_DAR, thr_logFC, 
        thr_fdr)
)

f_name = here(plotDir, paste0("ScatteredPlots_categories_2sharedCT_DARs_criteria_FDR", FDR, ".pdf"))
pdf(f_name, width = 8, height = 6)
walk(scattered_plt_cell_type_real_values, print)
dev.off()


## =============================================================================
## Heatmaps:

# (1) For summarizing accessibility (logFC, FDR, directionality) → Subset by cell_type
# - the goal centers on cell-type–specific DARs

# (2) For summarizing correlation strength (CCscore) →
# - if comparing co-regulated gene modules across clusters

# Filter plot_data_list to only consider those genes classified as: 
# categories:
# "Active CRE (+)"         # positively correlated & accessible
# “Repressive CRE (-)"     # negatively correlated & less accessible
# "Neutral CRE"            # New category for highly significant elements with marginal gene change
# "Shared CRE"             # correlated, not DAR
# "Unlinked OCR”.          # DAR, no correlation
# "Non-significant"        # everything else

# Active CRE (+)        Neutral CRE    Non-significant Repressive CRE (-) 
# 230                136                157                 13 
# Shared CRE       Unlinked OCR 
# 170                230


## =============================================================================
## Define functions to plot all the defined cRE classes 

make_heatmap_cRE <- function(
        heatmap_mat_ordered,
        subtitle_content
) {
    h1 <- pheatmap(
        heatmap_mat_ordered, 
        color = colorRampPalette(c("blue", "white", "red"))(100),
        cluster_rows = FALSE, # FALSE to respect row_order_final
        cluster_cols = FALSE,
        main = paste0( 
            "Top ", top_genes_heatmap, " Genes per Cell Type by Max |CCscore|",   
            "\n",
            subtitle_content
        ),
        fontsize_row = 7,
        fontsize_col = 10,
        show_rownames = TRUE,
        show_colnames = TRUE
    )
    return(h1)
    
}

## ====

## filter to consider those genes classified as Active CRE (+)"

top_genes_heatmap = 5
names(plot_data_list)
head(plot_data_list[[1]])


prepare_top_genes_heatmap <- function(data_list,
                                      category_peaks,
                                      top_genes_heatmap) 
    {
    
    # filtering within each cell type the specific category. ge. "Active CRE (+)"
    top_genes <- map(
        data_list,
        ~ .x |>
            filter(category == category_peaks) |> 
            arrange(desc(abs(CCscore))) |>
            slice_head(n = top_genes_heatmap) |>
            pull(gene_name) |>
            unique()
    )
    
    return(top_genes)

    }

##  ============================================================================


categories_to_plot <- c("Active cCRE (+)", "Repressive cCRE (-)", "Neutral cCRE", "Linked OCR")


for (cat_cRE in categories_to_plot) {
    # cat_cRE = "Repressive CRE (-)"
    
    top_genes_list <- prepare_top_genes_heatmap(
        plot_data_list,
        cat_cRE,
        top_genes_heatmap
        )
    
    names(top_genes_list) <- names(plot_data_list)
    top_genes <- unique(unlist(top_genes_list))
    message("Total unique top genes: ", length(top_genes))
    
    # Prepare the data for dcast (ensure no duplicates and correct type)
    wide_data_clean <- plot_data_full_df |>
        dplyr::filter(gene_name %in% top_genes) |>
        dplyr::select(gene_name, cell_type, CCscore) |>
        dplyr::group_by(gene_name, cell_type) |>
        # Keep the row with the largest absolute CCscore for each gene/cell pair
        dplyr::slice_max(order_by = abs(CCscore), n = 1, with_ties = FALSE) |>
        dplyr::ungroup()
    
    heatmap_mat <- reshape2::dcast(
        data = wide_data_clean, 
        formula = gene_name ~ cell_type, 
        value.var = "CCscore"
    )  |> dplyr::as_tibble()  # fun.aggregate is no longer needed
    
    # Set row names and remove the gene_name column
    heatmap_mat_clean <- heatmap_mat |> 
        # Move the 'gene_name' column to be the row names
        tibble::column_to_rownames(var = "gene_name") |>
        as.matrix()
    
    # Explicitly replace any lingering NAs with 0
    heatmap_mat_final <- replace(heatmap_mat_clean, is.na(heatmap_mat_clean), 0)
    
    # build diagonal
    max_col_index <- apply(abs(heatmap_mat_final), 1, which.max)
    max_col_name <- colnames(heatmap_mat_final)[max_col_index]
    # create a df for sorting
    sort_df <- data.frame(
        gene_name = rownames(heatmap_mat_final),
        max_cell_type = max_col_name,
        max_abs_score = apply(abs(heatmap_mat_final), 1, max) # include the max |CCscore| to break ties within the same cell type
    ) |>
        dplyr::arrange(max_cell_type, desc(max_abs_score))
    # Row order
    row_order_final <- sort_df$gene_name
    col_order_final <- unique(sort_df$max_cell_type) 
    heatmap_mat_ordered <- heatmap_mat_final[row_order_final, col_order_final]
    
    # make the plot
    f_name <- paste0("heatmap_", stringr::word(cat_cRE, 1), "_CCscore_top", top_genes_heatmap, "genes.pdf")
    f_name = here(plotDir, f_name)
    pdf(f_name, width = 5, height = 8)
    
    subtitle_content <- paste0(
        cat_cRE, "\n", 
        "Total Genes Plotted: ", nrow(heatmap_mat_ordered) 
    )
    
    # plot the heatmap with cRE
    tmp_h1 <- make_heatmap_cRE(
        heatmap_mat_ordered,
        subtitle_content)
    print(tmp_h1)
    
    dev.off()

}



## =============================================================================
## Classify links as within or cross: stacked bar (proportions and counts)

# Each row in links_summary is a peak–gene pair labeled according to whether
# the peak accessibility (cell_type) and correlated gene expression (cluster) belong to the same cluster.
links_summary <- uniques_df |>
    mutate(link_type = if_else(cell_type == cluster, "within_cluster", "cross_cluster"))
head(links_summary)

## Summarize counts per cell type
summary_counts <- links_summary |>
    count(cell_type, link_type) |>
    group_by(cell_type) |>
    mutate(
        total = sum(n),
        proportion = n / total
    )
summary_counts


## stacked bar (proportions) show, for each cell_type what fraction of its significant LinkPeaks 
## - connect to genes in the same cluster (green) vs other clusters

g1 <- ggplot(summary_counts, aes(x = cell_type, y = proportion, fill = link_type)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("within_cluster" = "#1b9e77", "cross_cluster" = "#d95f02")) +
    labs(
        title = "Within vs Cross-Cluster LinkPeaks-DARs overlaps by Cell Type",
        x = "Cell Type",
        y = "Proportion of Links",
        fill = "Link Type"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

f_name <- paste0("Cross-Cluster_proportions_FDR", FDR, ".pdf")
ggsave(here(plotDir, f_name),
       g1, width = 8, height = 8)


## stacked bar (absolute counts)
g1 <- ggplot(summary_counts, aes(x = cell_type, y = n, fill = link_type)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(y = "Number of Links") +
labs(
    title = "Within- vs Cross-Cluster LinkPeaks-DARs overlaps by Cell Type",
    x = "Cell Type",
    y = "Number of Links",
    fill = "Link Type"
) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

f_name <- paste0("Cross-Cluster_abs_counts_FDR", FDR, ".pdf")
ggsave(here(plotDir, f_name),
       g1, width = 8, height = 8)


## =============================================================================

## Highlighting reciprocal relationships
# matrix heatmap of regulatory interactions between cell-type pairs
# directionality and strength of inter-cluster links.

# heatmap h1:
# How many cross-cluster connections exist? 
# h1 - Excludes diagonal (cell_type != cluster), showing only cross-cluster links
# Shows how much cross-regulation occurs between clusters

pair_counts <- links_summary |>
    count(cell_type, cluster) |>
    filter(cell_type != cluster)
head(pair_counts)

f_name = here(plotDir, paste0("heatmap_CrossCluster_LinkCounts_h1_FDR", FDR,".pdf"))
pdf(f_name, width = 7, height = 6)
pheatmap::pheatmap(
    pivot_wider(pair_counts, names_from = cluster, values_from = n, values_fill = 0)[,-1],
    cluster_rows = TRUE, cluster_cols = TRUE,
    main = "Cross-Cluster Link Counts (peaks→genes)"
)
dev.off()

# ======

# heatmap h2:
# How strongly each cell type favors certain gene clusters?
# h2 Keeps all pairs, including within-cluster (diagonal)
# Highlights preferential link patterns (who regulates whom more strongly) across all clusters

pair_counts <- links_summary |>
    count(cell_type, cluster)
mat <- pair_counts |>
    pivot_wider(names_from = cluster, values_from = n, values_fill = 0)
mat_norm <- mat |>
    column_to_rownames("cell_type") |>
    as.matrix()
# normalize by row totals
mat_prop <- mat_norm / rowSums(mat_norm)
# # center/scale each row
mat_z <- t(scale(t(mat_norm)))  

# Y-axis (rows)	ATAC-defined cell_type	Where the open chromatin peak is located (e.g., LHb.4).
# X-axis (columns)	RNA-defined cluster	Where the correlated gene expression occurs (e.g., Excit.Thal).
# Color	Strength or proportion of links	Fraction (or Z-score) of LHb.4 peaks linked to genes in each expression cluster.
# plot
f_name = here(plotDir, paste0("heatmap_active_repressed_CCscore_top", top_genes_heatmap, "genes_FDR", thr_DAR,".pdf")) # FIX 4: Use thr_DAR
pdf(f_name, width = 8, height = 8)

subtitle_content <- paste0(
    "Genes Selected from 'Active CRE (+)' Category (CCscore > ", thr_CC, ", FDR-DAR < ", thr_DAR, ", logFC > ", thr_logFC, ")\n", # Recommended to be more specific
    "Total Genes Plotted: ", nrow(heatmap_mat_ordered) # Use ordered matrix size
)

pheatmap(
    heatmap_mat_ordered, # FIX 1: Plot the ordered matrix
    color = colorRampPalette(c("blue", "white", "red"))(100),
    cluster_rows = FALSE, # FIX 2: Set to FALSE to respect row_order_final
    cluster_cols = FALSE,
    main = paste0( # FIX 3: Clean up the main title
        "Top ", top_genes_heatmap, " Genes per Cell Type by Max |CCscore|",
        "\n",
        subtitle_content
    ),
    fontsize_row = 7,
    fontsize_col = 10,
    show_rownames = TRUE,
    show_colnames = TRUE
)

dev.off()

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

# ## compute log counts 
# make_scattered_plot_dars_cc <- function(subset_uniques) {
#     thr_CC <- 0.1
#     thr_DAR <- 0.1
#     thr_val <- thr_DAR
#     clus_name <- unique(subset_uniques[["cell_type"]])
#     
#     message("Building plot for [", clus_name, "]")
#     message("Using FDR threshold =", thr_val)
#     
#     plot_data <- subset_uniques |>
#         mutate(
#             neglog_FDR_CC = -log10(FDR_CC),
#             neglog_fdr_dars = -log10(fdr_dars),
#             sig_CC = FDR_CC < thr_CC,
#             sig_DAR = fdr_dars < thr_DAR,
#             category = case_when(
#                 sig_CC & sig_DAR & logFC > 0  ~ "Active CRE (+)",
#                 sig_CC & sig_DAR & logFC < 0  ~ "Repressive CRE (-)",
#                 sig_CC & !sig_DAR              ~ "Shared CRE",
#                 !sig_CC & sig_DAR              ~ "Unlinked OCR",
#                 TRUE                           ~ "Non-significant"
#             )
#         )
#     
#     top_hits <- plot_data |>
#         filter(category %in% c("Active CRE (+)", "Repressive CRE (-)")) |>
#         arrange(FDR_CC) |>
#         head(20)
#     
#     g1 <- ggplot(plot_data, aes(x = neglog_FDR_CC, y = neglog_fdr_dars, color = category)) +
#         geom_point(alpha = 0.7, size = 1.6) +
#         geom_text_repel(
#             data = top_hits,
#             aes(label = gene_name),
#             size = 3,
#             max.overlaps = 15
#         ) +
#         geom_vline(xintercept = -log10(thr_CC), linetype = "dashed", color = "darkgrey") +
#         geom_hline(yintercept = -log10(thr_DAR), linetype = "dashed", color = "darkgrey") +
#         scale_color_manual(values = c(
#             "Active CRE (+)" = "#E64B35FF",
#             "Repressive CRE (-)" = "#4DBBD5FF",
#             "Shared CRE" = "#00A087FF",
#             "Unlinked OCR" = "#3C5488FF",
#             "Non-significant" = "lightgrey"
#         )) +
#         +
#         scale_x_continuous(limits = c(0, 0.5)) +
#         scale_y_continuous(limits = c(0, 0.5)) +
#         labs(
#             title = "Peak–Gene Correlation vs. Differential Accessibility",
#             subtitle = paste0(clus_name, " | FDR = ", thr_val),
#             x = expression(-log[10](FDR[CC])),
#             y = expression(-log[10](FDR[DARs]))
#         ) +
#         theme_minimal(base_size = 12)
#     
#     return(g1)
# }

# # log10
# scattered_plt_cell_type <- purrr::map(
#     plot_data_list,
#     ~ make_scattered_plot_dars_cc(.x)
# )

