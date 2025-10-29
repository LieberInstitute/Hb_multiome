########################################################################
## Explore FDR-Scores and Top overlapping (Heatmap)
## - (1) Build a class of overlapping peaks
## - (2) Heatmap summarize accessibility (logFC, FDR, directionality) 
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

FDR = 0.2 # actual value to run the DAR-Links. We do not use log FC for this exploratory analysis

## Set directory names
inputCSV_Overlaps_Dir <- here(
    "processed-data",
    "06_peak_calling"
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
    subset_uniques <- unique_df |> filter(cell_type == cluster_specific)
    message("Subsetting [", cluster_specific, "]    ",
            length(unique(subset_uniques$peak_id_links))," cell_type")
    return(subset_uniques)
}


## =============================================================================
## Generate scatterplot

make_scattered_plot_dars_cc_real <- function(
        plot_data,
        categories_to_plot,
        top_genes,
        thr_CC,    # Thresh for Peak-Gene correlation
        thr_DAR,   # Thresh for Differential Accessibility Regions (DARs)
        thr_logFC, # Thresh for DARs Accessibility Direction
        fdr_cutoff = 0.2 # FDR used in both ds to define significant peaks 
) {

    message("cats_to_plot: ", paste(categories_to_plot, collapse = ", "))
    message("present in data: ", paste(unique(plot_data$category), collapse = ", "))
    
    clus_name <- unique(plot_data[["cell_type"]])
    message("Building plot for [", clus_name, "] (real FDR values)")
    
    #if (is.list(categories_to_plot)) { categories_to_plot <- unlist(categories_to_plot)}
    # Ensure categories_to_plot is a character vector
    if (is.null(categories_to_plot)) categories_to_plot <- character(0)
    if (is.list(categories_to_plot)) categories_to_plot <- unlist(categories_to_plot, use.names = FALSE)
    categories_to_plot <- as.character(categories_to_plot)
    # Use only categories present in the data
    cats_present <- intersect(categories_to_plot, unique(plot_data$category))
    
    # If nothing matches, keep all rows
    if (length(cats_present) == 0L) stop("Categories missed!")
    
    top_hits <- plot_data |>
        filter(category %in% cats_present) |>
        arrange(desc(abs(CCscore))) |> 
        head(top_genes)

    x_max <- min(1, max(plot_data$FDR_CC, na.rm = TRUE) * 1.05)
    y_max <- min(1, max(plot_data$fdr_dars, na.rm = TRUE) * 1.05)
    # define tick units
    axis_breaks <- seq(0, max(x_max, y_max), by = 0.05)
    axis_labels <- sprintf("%.2f", axis_breaks)
    
    subtitle_text <- paste0(
        "Spearman CC |p| > ", thr_CC,
        " | FDR-DARs < ", thr_DAR,
        " & log2FC ± ", thr_logFC
    )
    
    # Define the full palette (master color map) - we have 2 categories
    color_map <- c(
        "cell-specific cCRE (+)"   = "#E64B35FF",
        "cell-specific cCRE (-)"   = "#800080",
        "Linked_DAR (+) enriched"  = "#E64B35FF",  # reuse similar red for CSC+Nick version
        "Linked_DAR (-) depleted"  = "#800080",    # reuse purple for CSC+Nick version
        "(-) Linked (+) DAR enriched" ="#0424DB",  # rare concordance
        "Linked OCR"               = "#00A087FF",
        "Unlinked DAR"             = "#3C5488FF",
        "Non-significant"          = "lightgrey"
    )
    
    desired_order <- c(
        "cell-specific cCRE (+)",
        "cell-specific cCRE (-)",
        "Linked_DAR (+) enriched",
        "Linked_DAR (-) depleted",
        "(-) Linked (+) DAR enriched", 
        "Linked OCR",
        "Unlinked DAR",
        "Non-significant"
    )
    # Restrict to categories actually present in our subset
    present_categories <- intersect(desired_order, unique(plot_data$category))
    
    g1 <- ggplot(plot_data, aes(x = CCscore, y = logFC, color = category)) +
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
        scale_color_manual(
            values = color_map[present_categories],
            breaks = present_categories,
            drop = FALSE
        ) +
        labs(
            title = paste(clus_name, " | Peak-Gene Correlation vs DARs"),
            subtitle = subtitle_text,
            x = expression("Spearman CC " * "|" * rho * "|"),
            y = expression("log"[2] * "FC (DARs)")
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
summary(subsetted_df$cell_type)
head(subsetted_df)
subsetted_list_df <- split(subsetted_df, subsetted_df$cell_type)
message("Peaks found for ", length(subsetted_list_df), " cell types")
names(subsetted_list_df)
#table(subsetted_list_df[[8]]["overlap_type"])


## make scattered plots using real values (non-normalized)

## Set thresholds for LinkPeaks (CC) and DARs

top_genes_scattered_plt = 20
thr_CC = 0.3  # correlation strength
thr_DAR = 0.1 # fdr_dars / accessibility significance
thr_logFC = 0
thr_fdr = 0.2


## Categorize peaks overlaps ===================================================

# (1) Standard classification (Cynthia): classification-1 - removed

# (2) Leo's classification adapted to our current analysis: classification-2
peaks_classification2 = c("cell-specific cCRE (+)", "cell-specific cCRE (-)", "Linked OCR", "Unlinked DAR", "Non-significant")
        
# (3) Cynthia+Nick classification adapted to our current analysis: classification-3
peaks_classification3 = c("Linked_DAR (+) enriched", "Linked_DAR (-) depleted", "Linked OCR", "(-) Linked (+) DAR enriched", "Unlinked DAR", "Non-significant")

# ## Add peaks-classification name
# subsetted_list_df <- purrr::map(subsetted_list_df, ~ 
#                                     .x |> mutate(type_classification = "classification-2")  # you can change dynamically later
# )

## (A) Classification-2 logic  =================================================
plot_data_list_2 <- purrr::map(subsetted_list_df, ~ .x |> 
    mutate(
        sig_CC  = abs(CCscore) > thr_CC,   # significant correlation
        sig_DAR = fdr_dars < thr_DAR,      # significant accessibility
        category = case_when(
            sig_CC & sig_DAR & logFC > thr_logFC  ~ "cell-specific cCRE (+)",      # positively correlated & accessible
            sig_CC & sig_DAR & logFC < -thr_logFC  ~ "cell-specific cCRE (-)",     # negatively correlated & less accessible
            # Neural category was used to identify "Significantly linked, Significant DAR, but logFC is too small", when set log_FC=0 we do not need it any more
            # sig_CC & sig_DAR                ~ "Neutral cCRE",        # when other conditions not met  
            sig_CC & !sig_DAR               ~ "Linked OCR",          # correlated, not DAR
            !sig_CC & sig_DAR               ~ "Unlinked DAR",        # DAR, no correlation
            TRUE                            ~ "Non-significant"      # everything else
        ),
        type_classification = "classification-2"
    )
)

## (B) Classification-3 logic  =================================================
plot_data_list_3 <- purrr::map(subsetted_list_df, ~ .x |> 
   mutate(
       sig_CC  = abs(CCscore) > thr_CC,
       sig_DAR = fdr_dars < thr_DAR,
       category = case_when(
           # Primary interest: Concordant Linked DAR Categories plus secondary interest categories
           sig_CC & sig_DAR & logFC > thr_logFC   ~ "Linked_DAR (+) enriched",   # positive correlation + open chromatin
           sig_CC & sig_DAR & logFC < -thr_logFC  ~ "Linked_DAR (-) depleted",   # negative correlation + closed chromatin
           # This "New" category of Non-Concordant Categories (Negative/Positive logFC mismatch with CCscore sign)
           # - is rare and usually analyzed with co-accessibility or gene-silencing experiments
           # - suggest Distant Regulation / These DARs may be driving the differential expression of non-coding RNAs 
           (sig_CC < 0) & (sig_DAR & logFC > thr_logFC) ~ "(-) Linked (+) DAR enriched", 
           # Secondary interest categories
           sig_CC & !sig_DAR                      ~ "Linked OCR",
           !sig_CC & sig_DAR                      ~ "Unlinked DAR",
           # Nobody cares 
           TRUE                                   ~ "Non-significant"
       ),
       type_classification = "classification-3"
   )
)


## ============/

## Merge both classification versions for comparison
plot_data_full_df <- bind_rows(
    bind_rows(plot_data_list_2, .id = "source_df"),
    bind_rows(plot_data_list_3, .id = "source_df")
)
plot_data_full_df$type_classification <- factor(
    plot_data_full_df$type_classification,
    levels = c("classification-2", "classification-3")
)

## Build a summary 

message("========= Summary of candidate RE by category ============\n")

message("========= classification-2 ============\n")
plt_tmp <- plot_data_full_df |> filter(type_classification=="classification-2")
plt_tmp |> count(category, name = "n")
plt_tmp |> count("n") #1967
message("========= classification-3 ============\n")
plt_tmp2 <- plot_data_full_df |> filter(type_classification=="classification-3")
plt_tmp2 |> count(category, name = "n")

# category                    n
# 1              Linked OCR  420
# 2 Linked_DAR (+) enriched  706
# 3 Linked_DAR (-) depleted  657
# 4         Non-significant  354
# 5            Unlinked DAR 1068

#table(plot_data_full_df$cell_type) # includes both categories

## testing in one category only
hb_related_df <- plt_tmp2 |>
    filter(grepl("MHb|LHb", cell_type))
total_hb_related <- nrow(hb_related_df)
hb_related_signif_df <- hb_related_df |> 
    filter(category %in% peaks_classification3[peaks_classification3 != "Non-significant"])
total_hb_related_signif <- nrow(hb_related_signif_df)

message("Total Hb related [thr_CC=", thr_CC,
        " & thr_DAR=", thr_DAR,
        " & thr_logFC=", thr_logFC, "]: ", total_hb_related)
message("Total Hb related significant: ", total_hb_related_signif)
message("========================================================\n")


# save overlapping with classification
f_name <- here(processedDir, paste0("overlaps_linkPeak_DARs_classified_thr_CC", thr_CC, "_thr_DAR", thr_DAR, ".csv"))
write.csv(plot_data_full_df, f_name, row.names = FALSE)

message("Saved linkPeak_DARs overlapings with categories!")

## Plot both classifications
# plot_data_list_2 → list of 15 data frames 
# peaks_classification3 → character vectors (5 elements each)

scattered_plt_cell_type_real_values_2 <- purrr::map(
    plot_data_list_2,
    ~ make_scattered_plot_dars_cc_real(
        plot_data = .x,
        categories_to_plot = peaks_classification2,
        top_genes = top_genes_scattered_plt,
        thr_CC = thr_CC,
        thr_DAR = thr_DAR,
        thr_logFC = thr_logFC,
        fdr_cutoff = thr_fdr
    )
)
#scattered_plt_cell_type_real_values_2[2]
f_name = here(plotDir, paste0("ScatteredPlots_2sharedCT_class2_FDR", FDR, ".pdf"))
pdf(f_name, width = 8, height = 6)
walk(scattered_plt_cell_type_real_values_2, print)
dev.off()

scattered_plt_cell_type_real_values_3 <- purrr::map(
    plot_data_list_3,
    ~ make_scattered_plot_dars_cc_real(
        plot_data = .x,
        categories_to_plot = peaks_classification3,
        top_genes = top_genes_scattered_plt,
        thr_CC = thr_CC,
        thr_DAR = thr_DAR,
        thr_logFC = thr_logFC,
        fdr_cutoff = thr_fdr
    )
)
#scattered_plt_cell_type_real_values_3[2]
f_name = here(plotDir, paste0("ScatteredPlots_2sharedCT_class3_FDR", FDR, ".pdf"))
pdf(f_name, width = 8, height = 6)
walk(scattered_plt_cell_type_real_values_3, print)
dev.off()

## =============================================================================
## Heatmaps:

# (1) For summarizing accessibility (logFC, FDR, directionality) → Subset by cell_type
# - the goal centers on cell-type–specific DARs

# (2) For summarizing correlation strength (CCscore) →
# - if comparing co-regulated gene modules across clusters

## =============================================================================
## Define functions to plot heatmaps for specific categories 

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
            #"Top ", top_genes_heatmap, " Genes per Cell Type by Max |CCscore|",   
            #"\n",
            subtitle_content
        ),
        fontsize_row = 7,
        fontsize_col = 8,
        show_rownames = TRUE,
        show_colnames = TRUE
    )
    return(h1)
    
}


## function to plot heatmaps
prepare_top_genes_heatmap <- function(data_list,
                                      category_peaks,
                                      top_genes_heatmap) 
    {
    
    # filtering within each cell type the specific category. ge. "Active CRE (+)"
    top_genes <- map(
        data_list,
        ~ .x |>
            filter(category == !!category_peaks) |> 
            arrange(desc(abs(CCscore))) |>
            slice_head(n = top_genes_heatmap) |>
            select(gene_name, cell_type, CCscore) |>
            distinct() |>
            pivot_wider(
                names_from = cell_type,
                values_from = CCscore,
                values_fill = 0
            ) |>
            tibble::column_to_rownames("gene_name") |>
            as.matrix()
    )
    
    return(top_genes)

}

##  ============================================================================

## Plot classifiction-3 

top_genes_heatmap = 5
data_heatmap <- bind_rows(plot_data_list_3)
message("Rows in data_heatmap: ", nrow(data_heatmap))
print(table(data_heatmap$type_classification))
print(table(data_heatmap$category))

categories_to_plot <- c(
    "Linked_DAR (+) enriched",
    "Linked_DAR (-) depleted",
    "Linked OCR",
    "Unlinked DAR"
)

#categories_to_plot <- c("cell-specific cCRE (+)", "cell-specific cCRE (-)", "Linked OCR", "Unlinked DAR")
nrow(plot_data_full_df)
table(plot_data_full_df$type_classification)
for (cat_cRE in categories_to_plot) {
    # cat_cRE = "Linked_DAR (+) enriched"
    
    top_genes_list <- prepare_top_genes_heatmap(
        plot_data_list_3,
        cat_cRE,
        top_genes_heatmap
    )
    
    names(top_genes_list) <- names(plot_data_list)
    top_genes <- unique(unlist(top_genes_list))
    message("[", cat_cRE,"] Total unique top genes: ", length(top_genes))
    
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
        "Total Genes Plotted for: ", nrow(heatmap_mat_ordered) 
    )
    
    # plot the heatmap with cRE
    tmp_h1 <- make_heatmap_cRE(
        heatmap_mat_ordered,
        subtitle_content)
    print(tmp_h1)
    
    dev.off()

}


message("All plots done!!!")

# library("slurmjobs")
# job_single(
#   "21_overlaping_FDRscores_TopHeatmap",
#   create_shell = TRUE,
#   partition = "katun",
#   memory = "30G",
#   cores = 2,
#   logdir = "logs",
#   command = "Rscript 21_overlaping_FDRscores_TopHeatmap.R",
#   create_logdir = FALSE
# )

## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
