########################################################################
## Explore FDR-Scores and Top overlapping (Heatmap)
## - (1) Form a single CSV of DARs, linked peaks, and their overlaps.
##       Categorize all peaks
## - (2) Heatmap summarize accessibility (logFC, FDR, directionality) 
##    → Subset by cell_type, because visual goal centers on cell-type–specific regulatory accessibility
##
## Authors. CSC
## Date. Oct 14, 2025
## Recommended resources on interactive mode: srun --pty --mem=20GB --x11 bash
########################################################################

library("pheatmap")
library("reshape2")
library("patchwork")
library("ggrepel")
library("tidyverse")
library("here")

#===============================================================================
# resolution_level = "Mid"      # 18 cell-types
#===============================================================================

top_genes_scattered_plt = 20

overlap_path = here(
    "processed-data", "06_peak_calling", "19_Linkage_DARs_analysis",
    "Overlaps_LinkPeak_DARs_FDR0.1_0.1.csv"
)
links_dir = here(
    "processed-data", "06_peak_calling", "14_exploratory_pb_peak_scores_MACS2",
    "links_ct_merged"
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

color_map <- c(
    "Linked DAR enriched"   = "#E64B35FF",
    "Linked DAR depleted"   = "#800080",
    "Linked OCR enriched"   = "#00A087FF",
    "Linked OCR depleted"   = "#3B7000FF",
    "Discordant Linked DAR" = "#0424DB",  
    "Unlinked DAR"          = "#3C5488FF"
)

if (!dir.exists(processedDir)) { dir.create(processedDir) }
if (!dir.exists(plotDir)) { dir.create(plotDir) }

################################################################################
#   Form a single tibble of DARs, linked peaks, and their overlaps. Categorize
#   peaks
################################################################################

#   At this point we have a CSV of DARs, of linked peaks, and of their overlaps.
#   Merge into a single tibble and categorize them

#-------------------------------------------------------------------------------
#   Overlaps
#-------------------------------------------------------------------------------

overlap_df <- read_csv(overlap_path, show_col_types = FALSE)

#   Start with overlaps where the DAR is unique to one cell type

#   There are only 2 peaks where differential accessibility was found in 2
#   cell types. We might as well drop them, which allows us to claim all
#   DARs are unique to one cell type
overlap_df |>
    group_by(peak_id) |>
    filter(n() > 1) |>
    ungroup() |>
    select(peak_id, cell_type, link_gene_name) |>
    print()

overlap_df = overlap_df |>
    group_by(peak_id) |>
    filter(n() == 1) |>
    ungroup() |>
    mutate(
        category = case_when(
            (dar_logFC > 0) & (link_cc_score > 0) ~ "Linked DAR enriched",
            (dar_logFC < 0) & (link_cc_score < 0) ~ "Linked DAR depleted",
            TRUE ~ "Linked DAR discordant"
        )
    )

#-------------------------------------------------------------------------------
#   Linked peaks
#-------------------------------------------------------------------------------

link_files <- list.files(
    path = links_dir, pattern = ("^Mid.*\\.csv$"), full.names = TRUE
)
link_df = lapply(link_files, read_csv, show_col_types = FALSE) |>
    bind_rows() |>
    dplyr::rename(
        peak_id = peak,
        cell_type = cluster,
        link_cc_score = score,
        link_gene_id = gene_id,
        link_gene_name = gene,
        link_gene_strand = strand,
        link_fdr_cc = FDR
    ) |>
    mutate(
        dar_fdr = NA_real_,
        dar_logFC = NA_real_,
        category = ifelse(
            link_cc_score > 0, "Linked OCR enriched", "Linked OCR depleted"
        )
    ) |>
    select(all_of(colnames(overlap_df)))
    
stopifnot(setequal(colnames(overlap_df), colnames(link_df)))

## =============================================================================

all_clusters <- unique(overlap_df$cell_type)
# sort by MHh and LHb first 
clusters_sorted <- c(
    sort(grep("MHb|LHb", all_clusters, value = TRUE)), 
    sort(grep("MHb|LHb", all_clusters, value = TRUE, invert = TRUE))
)
    
subset_cell_type <- function(unique_df, cluster_specific) {
    subset_uniques <- unique_df |> filter(cell_type == cluster_specific)
    message("Subsetting [", cluster_specific, "] cell_type. Total unique OCR found: ",
            length(unique(subset_uniques$peak_id)))
    return(subset_uniques)
}


## =============================================================================
## Generate scatterplot

make_scattered_plot_dars_cc_real <- function(
        plot_data,
        categories_to_plot,
        top_genes,
        thr_CC = 0.3,    # Thresh for Peak-Gene correlation
        thr_DAR = 0.1,   # Thresh for Differential Accessibility Regions (DARs)
        fdr_cutoff = 0.1 # FDR used in both ds to define significant peaks 
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
        "Spearman CC |p| > ", thr_CC, " | FDR-DARs < ", thr_DAR
    )
    
    # Restrict to categories actually present in our subset
    present_categories <- intersect(
        names(color_map), unique(plot_data$category)
    )
    
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

## Categorize peaks overlaps ===================================================

## (B) Classification-3 logic  =================================================
plot_data_list_3 <- purrr::map(subsetted_list_df, ~ .x |> 
   mutate(
       sig_CC  = abs(CCscore) > thr_CC,
       sig_DAR = fdr_dars < thr_DAR,

       category = case_when(
           # ============================================================
           # 1. Concordant Linked DAR (+) enriched & depleted
           # ============================================================
           # Primary interest: Concordant Linked DAR Categories plus secondary interest categories
           sig_CC & sig_DAR &
               (CCscore > thr_CC) & (logFC > thr_logFC)  ~ "Linked_DAR (+) enriched", 
           sig_CC & sig_DAR &
               (CCscore < -thr_CC) & (logFC < thr_logFC) ~ "Linked_DAR (-) depleted", 
           
           # ============================================================
           # 2. Discordant Linked DAR (opposite signs!)
           # ============================================================
           # Non-Concordant Category (Negative/Positive logFC mismatch with CCscore sign)
           # - is rare and usually analyzed with co-accessibility or gene-silencing experiments
           # - suggest Distant Regulation or these DARs may be driving the differential expression of non-coding RNAs 
           sig_CC & sig_DAR &
               ((CCscore > thr_CC) & (logFC < thr_logFC) |
                    (CCscore < -thr_CC) & (logFC > thr_logFC)) ~ "Discordant Linked DAR",
           
           # ============================================================
           # 3. Linked OCR (+) enriched & depleted (gene strength)
           #    correlated but NOT a DAR
           # ============================================================
           # Asjusted categories correlated AND none DARs.  It will replace "Linked OCR"
           sig_CC & !sig_DAR &
               (CCscore > thr_CC) ~ "Linked OCR (+) enriched",
           sig_CC & !sig_DAR &
               (CCscore < -thr_CC) ~ "Linked OCR (-) depleted",
           
           # ============================================================
           # 4. Unlinked DAR (DAR but no correlation)
           # ============================================================
           # Secondary interest categories
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
    # bind_rows(plot_data_list_2, .id = "source_df"),
    bind_rows(plot_data_list_3, .id = "source_df")
)
plot_data_full_df$type_classification <- factor(
    plot_data_full_df$type_classification,
    levels = "classification-3"
    #levels = c("classification-2", "classification-3")
)

## Build a summary 

message("========= Summary of candidate RE by category ============\n")
message("========= classification-3 ============\n")
plt_tmp2 <- plot_data_full_df |> filter(type_classification=="classification-3")
plt_tmp2 |> count(category, name = "n")

## testing in one category only
hb_related_df <- plt_tmp2 |>
    filter(grepl("MHb|LHb", cell_type))
total_hb_related <- nrow(hb_related_df)

hb_related_df |> count(category)

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
#scattered_plt_cell_type_real_values_3[1]
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
            subtitle_content
        ),
        fontsize_row = 7,
        fontsize_col = 8,
        show_rownames = TRUE,
        show_colnames = TRUE,
        silent = TRUE 
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
            # Aggregate to find the single best (highest abs CCscore) link per gene per cell type
            group_by(gene_name, cell_type) |>
            reframe(
                CCscore = CCscore[which.max(abs(CCscore))]
            ) |>
            ungroup() |>
            arrange(desc(abs(CCscore))) |>
            slice_head(n = top_genes_heatmap) |>

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

## Plot Heatmap with classifiction-3 

top_genes_heatmap = 5
plot_data_full_df <- plot_data_full_df |> filter(type_classification=="classification-3")
table(plot_data_full_df$type_classification)
table(plot_data_full_df$category)

categories_to_plot <- c(
    "Linked_DAR (+) enriched",
    "Linked_DAR (-) depleted",
    "Linked OCR (+) enriched",
    "Linked OCR (-) depleted",
    "Discordant Linked DAR",
    "Unlinked DAR",
    "Non-significant"
)


for (cat_cRE in categories_to_plot) {
    # cat_cRE = "Linked_DAR (+) enriched"
    
    top_genes_list <- prepare_top_genes_heatmap(
        plot_data_list_3,
        cat_cRE,
        top_genes_heatmap
    )
    
    
    names(plot_data_list_3) <- names(plot_data_list_3)
    top_genes <- unique(unlist(purrr::map(top_genes_list, rownames)))
    
    if (length(top_genes) == 0) {
        message("Skipping [", cat_cRE, "] because no top genes were found across all cell types.")
        next
    } else {
        message("[", cat_cRE,"] Total unique top genes: ", length(top_genes))
    }
        
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
    cat_name <- case_when(
        cat_cRE == "Linked_DAR (+) enriched" ~ "Linked_DAR_enriched",
        cat_cRE == "Linked_DAR (-) depleted" ~ "Linked_DAR_depleted",
        cat_cRE == "Linked OCR (+) enriched" ~ "Linked_OCR_enriched",
        cat_cRE == "Linked OCR (-) depleted" ~ "Linked_OCR_depleted",
        cat_cRE == "Discordant Linked DAR" ~ "Discordant_Linked_DAR",
        cat_cRE == "Unlinked DAR" ~ "Unlinked_DAR",
        cat_cRE == "Non-significant" ~ "Non_significant"
    )
    f_name <- paste0("heatmap_", cat_name, "_CCscore_top", top_genes_heatmap, "genes.pdf")
    f_name = here(plotDir, f_name)
    pdf(f_name, width = 5, height = 8)
    
    subtitle_content <- paste0(
        cat_cRE, "/n", 
        "Plotted ", nrow(heatmap_mat_ordered), " OCRs" 
    )
    
    # plot the heatmap with cRE
    tmp_h1 <- make_heatmap_cRE(
        heatmap_mat_ordered,
        subtitle_content)
    # print(tmp_h1)
    # this is the correct print method for pheatmap
    grid::grid.newpage()
    grid::grid.draw(tmp_h1$gtable)
    
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
