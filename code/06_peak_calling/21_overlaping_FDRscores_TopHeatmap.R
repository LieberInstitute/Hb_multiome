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
top_genes_heatmap = 5

overlap_path = here(
    "processed-data", "06_peak_calling", "19_Linkage_DARs_analysis",
    "Overlaps_LinkPeak_DARs_FDR0.1_0.1.csv"
)
link_path = here(
    "processed-data", "06_peak_calling", "19_Linkage_DARs_analysis",
    "ALL_LinkPeaks_signif_FDR01.csv"
)
dars_path = here(
    "processed-data", "06_peak_calling", "19_Linkage_DARs_analysis",
    "ALL_DARs_signif_FDR01.csv"
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
    "Linked DAR discordant" = "#0424DB",
    "Unlinked DAR enriched" = "#3C5488FF",
    "Unlinked DAR depleted" = "#F39B7FFF"
)

if (!dir.exists(processedDir)) { dir.create(processedDir) }
if (!dir.exists(plotDir)) { dir.create(plotDir) }

################################################################################
#   Functions
################################################################################

make_scattered_plot_dars_cc_real <- function(
        plot_data,
        top_genes,
        #   The following are only used for plot text. Subsetting has already
        #   been done upstream
        thr_CC = 0.3,    # Thresh for Peak-Gene correlation
        thr_DAR = 0.1,   # Thresh for Differential Accessibility Regions (DARs)
        fdr_cutoff = 0.1 # FDR used in both ds to define significant peaks 
) {    
    clus_name <- unique(plot_data[["cell_type"]])
    
    top_hits <- plot_data |>
        arrange(desc(abs(link_cc_score))) |> 
        head(top_genes)

    x_max <- min(1, max(plot_data$link_fdr_cc, na.rm = TRUE) * 1.05)
    y_max <- min(1, max(plot_data$dar_fdr, na.rm = TRUE) * 1.05)
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
    
    g1 <- ggplot(plot_data, aes(x = link_cc_score, y = dar_logFC, color = category)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
        geom_vline(xintercept = 0, linetype = "solid", color = "grey70") +
        geom_vline(xintercept = c(-thr_CC, thr_CC), linetype = "dashed", color = "darkgrey") +
        geom_point(alpha = 0.7, size = 1.6) +
        geom_text_repel(
            data = top_hits,
            aes(label = link_gene_name),
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
        plot_annotation(caption = paste("Uniques CellTypes\nTop:", top_genes, "genes"))
    
    return(g1)
}

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
            # Aggregate to find the single best (highest abs link_cc_score) link per gene per cell type
            group_by(link_gene_name, cell_type) |>
            reframe(
                link_cc_score = link_cc_score[which.max(abs(link_cc_score))]
            ) |>
            ungroup() |>
            arrange(desc(abs(link_cc_score))) |>
            slice_head(n = top_genes_heatmap) |>

            pivot_wider(
                names_from = cell_type,
                values_from = link_cc_score,
                values_fill = 0
            ) |>
            tibble::column_to_rownames("link_gene_name") |>
            as.matrix()
    )
    
    return(top_genes)
}

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
    filter(n_distinct(cell_type) > 1) |>
    ungroup() |>
    select(peak_id, cell_type, link_gene_name) |>
    print()

overlap_df = overlap_df |>
    group_by(peak_id) |>
    filter(n_distinct(cell_type) == 1) |>
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

link_df = read_csv(link_path, show_col_types = FALSE) |>
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
    select(all_of(colnames(overlap_df))) |>
    filter(
        !(
            (peak_id %in% overlap_df$peak_id) &
            (cell_type %in% overlap_df$cell_type)
        )
    )
    
stopifnot(setequal(colnames(overlap_df), colnames(link_df)))

#-------------------------------------------------------------------------------
#   DARs
#-------------------------------------------------------------------------------

dars_df = read_csv(dars_path, show_col_types = FALSE) |>
    dplyr::rename(dar_fdr = fdr, dar_logFC = logFC) |>
    mutate(
        category = ifelse(
            dar_logFC > 0, "Unlinked DAR enriched", "Unlinked DAR depleted"
        )
    ) |>
    filter(
        !(
            (peak_id %in% overlap_df$peak_id) &
            (cell_type %in% overlap_df$cell_type)
        )
    )

for (other_colname in setdiff(colnames(overlap_df), colnames(dars_df))) {
    dars_df[[other_colname]] <- NA
}

stopifnot(setequal(colnames(overlap_df), colnames(dars_df)))

#   Merge the three
peak_df = bind_rows(overlap_df, link_df, dars_df)

message("Peak categories after merging:")
table(peak_df$category)

message("Category counts by cell type:")
table(peak_df$cell_type, peak_df$category)

write_csv(peak_df, here(processedDir, "all_peaks_categorized.csv.gz"))

################################################################################
#   Plots
################################################################################

split_peak_df <- split(peak_df, peak_df$cell_type)

scattered_plt_cell_type_real_values_3 <- purrr::map(
    split_peak_df,
    ~ make_scattered_plot_dars_cc_real(
        plot_data = .x, top_genes = top_genes_scattered_plt
    )
)
f_name = here(plotDir, "ScatteredPlots.pdf")
pdf(f_name, width = 8, height = 6)
walk(scattered_plt_cell_type_real_values_3, print)
dev.off()

## =============================================================================
## Heatmaps:

# (1) For summarizing accessibility (logFC, FDR, directionality) → Subset by cell_type
# - the goal centers on cell-type–specific DARs

# (2) For summarizing correlation strength (link_cc_score) →
# - if comparing co-regulated gene modules across clusters

for (cat_cRE in names(color_map)) {
    # cat_cRE = "Linked_DAR (+) enriched"
    
    top_genes_list <- prepare_top_genes_heatmap(
        split_peak_df,
        cat_cRE,
        top_genes_heatmap
    )
    
    top_genes <- unique(unlist(purrr::map(top_genes_list, rownames)))
    
    if (length(top_genes) == 0) {
        message("Skipping [", cat_cRE, "] because no top genes were found across all cell types.")
        next
    } else {
        message("[", cat_cRE,"] Total unique top genes: ", length(top_genes))
    }
        
    # Prepare the data for dcast (ensure no duplicates and correct type)
    wide_data_clean <- peak_df |>
        dplyr::filter(link_gene_name %in% top_genes) |>
        dplyr::select(link_gene_name, cell_type, link_cc_score) |>
        dplyr::group_by(link_gene_name, cell_type) |>
        # Keep the row with the largest absolute link_cc_score for each gene/cell pair
        dplyr::slice_max(order_by = abs(link_cc_score), n = 1, with_ties = FALSE) |>
        dplyr::ungroup()
    
    heatmap_mat <- reshape2::dcast(
        data = wide_data_clean, 
        formula = link_gene_name ~ cell_type, 
        value.var = "link_cc_score"
    )  |> dplyr::as_tibble()  # fun.aggregate is no longer needed
    
    # Set row names and remove the gene_name column
    heatmap_mat_clean <- heatmap_mat |> 
        # Move the 'link_gene_name' column to be the row names
        tibble::column_to_rownames(var = "link_gene_name") |>
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
        max_abs_score = apply(abs(heatmap_mat_final), 1, max) # include the max |link_cc_score| to break ties within the same cell type
    ) |>
        dplyr::arrange(max_cell_type, desc(max_abs_score))
    # Row order
    row_order_final <- sort_df$gene_name
    col_order_final <- unique(sort_df$max_cell_type) 
    heatmap_mat_ordered <- heatmap_mat_final[row_order_final, col_order_final]
    
    # make the plot
    cat_name <- gsub(" ", "_", cat_cRE)
    f_name <- paste0("heatmap_", cat_name, "_link_cc_score_top", top_genes_heatmap, "genes.pdf")
    f_name = here(plotDir, f_name)
    pdf(f_name, width = 5, height = 8)
    
    subtitle_content <- paste0(
        cat_cRE, "/n", 
        "Plotted ", nrow(heatmap_mat_ordered), " peaks" 
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
