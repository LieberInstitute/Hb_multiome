########################################################################
## Compare width peaks distributions between CellRanger Peaks and Signac::CallPeakss() outputs
##
## Authors. CSC
## Date. Aug 28, 2025
## Recommended resources on interactive mode: srun --pty --mem=30GB --x11 bash
########################################################################

library("EnsDb.Hsapiens.v86")           # Gene annotation (GTF-style)
library("ggplot2")
library("patchwork")
library("tidyverse")
library("dplyr")
library("scales")
library("here")


# for testing
p_met = "spearman"
# w_size = "1e5"
w_sizes = c("1e5", "5e4", "2.5e4")

here::here()

# Check/create directories
input_cvsDir <- here(
    "processed-data",
    "06_peak_calling",
    "00_link_peaks"
)
plotDir <- here(
    "plots",
    "06_peak_calling",
    "09_compare_peak_distributions"
)
# csvDir <- here(
#     "processed-data",
#     "06_peak_calling",
#     "06_exploratory_peak_scores"
# )

## Check directories
if (!dir.exists(plotDir)) {
    dir.create(plotDir)
}

## =============================================================================
## load CellRanger-ARC peaks counts and compute TSS stats 

## list of data frames with peak coordinates, central values and peak width
lst_peaks_df <- list()

message("Processing peak width histograms ...")


for (ws in w_sizes) {
    # ws = w_sizes[3]
    print(ws)    
    if (length(p_met)) {
        message(
            "Processing job for peak-method:\n",
            p_met,
            "\nWindow-size\n",
            ws
        )
        f_sufix <- paste0(".", p_met, ".", ws, ".cells_filtered_2perc")
        message("Processing: ", f_sufix)
    } else {
        message("Input arguments missed")
        stop()
    }
    
    # load link peak-gene csv
    gene_peaks_csv <- here(input_cvsDir, paste0("all_peak_gene_links", f_sufix, ".csv"))
    if (file.exists(gene_peaks_csv)) {
        link_df <- read.csv(gene_peaks_csv)
        message("File loaded!")
    } else {
        stop(paste("File not found:", gene_peaks_csv))
    }
    
    ##==========================================================================
    ## prepare df with peak ranges and width
    
    colnames(link_df)
    head(link_df, n=3)
    nrow(link_df)
    
    link_df <- link_df |>
        mutate(
            gene     = trimws(as.character(gene)),
            seqnames = as.character(seqnames),
            start    = as.numeric(start),
            end      = as.numeric(end)
        )
    head(link_df)
    
    message("Computing distance between peaks and TSS ...")
    
    ## Build TSS (strand-aware) GRanges table and compute the distance between each peak and its linked gene's TSS
    gene_coords <- genes(EnsDb.Hsapiens.v86)
    # tss_coords <- resize(gene_coords, width = 1, fix = "start")
    # tss_coords <- keepStandardChromosomes(tss_coords, pruning.mode = "coarse")
    # head(tss_coords)
    
    # Build TSS table (strand-aware)
    # extracts promoter regions from those gene_coords, upstream = 0 (don’t include any bases before the TSS, 
    # downstream = 1 (take exactly one base downstream from the TSS)
    tss_coords  <- promoters(gene_coords, upstream = 0, downstream = 1) %>%   # 1bp TSS, respects strand
        keepStandardChromosomes(pruning.mode = "coarse")
    # get a 1 bp range that pinpoints the TSS position for every gene, regardless of strand orientation
    
    # match UCSC-style peaks (chr1, chr2, etc.)
    seqlevelsStyle(tss_coords) <- "UCSC"
    
    tss_raw <- as.data.frame(tss_coords)
    has_biotype <- "gene_biotype" %in% colnames(tss_raw)
    
    # then take first per (gene_name, chr) - this avoid 1:many associations
    tss_df <- tss_raw %>%
        mutate(gene_biotype = if (has_biotype) gene_biotype else NA_character_) %>%
        #arrange(desc(gene_biotype == "protein_coding")) %>%  # prefer protein-coding where available
        group_by(gene_name, seqnames) %>%
        slice_head(n = 1) %>%
        ungroup() %>%
        transmute(
            gene_name,
            seqnames = as.character(seqnames),
            tss      = as.numeric(start),   # rename for clarity
            gene_strand = as.character(strand),
            gene_id
        )

    ## check duplicates
    dup_pairs <- tss_df %>%
        count(gene_name, seqnames, name = "n") %>%
        filter(n > 1)
    
    if (nrow(dup_pairs) > 0) {
        print(head(dup_pairs, 10))
        stop("Non-unique (gene_name, seqnames) in tss_df: ", nrow(dup_pairs), " duplicates.")
    }
    
    # Join by gene + chromosome to avoid many-to-many 
    colnames(link_df)
    colnames(tss_df)
    link_df2 <- link_df %>%
        left_join(tss_df, by = c("gene" = "gene_name", "seqnames" = "seqnames"))
    head(link_df2, n = 3)

    # drop rows with no TSS match
    n_before <- nrow(link_df2)
    link_df2 <- link_df2 %>% filter(!is.na(tss))
    message("Dropped ", n_before - nrow(link_df2), " rows with no TSS match.")
    
    nrow(link_df2)
    
    message("Distance between peaks and TSS added ...")
    
    message("Computing Peak center and distance to TSS ...")
    
    link_df2 <- link_df2 %>%
        mutate(
            peak_center       = (start + end) / 2,
            distance          = abs(peak_center - tss),
            signed_distance   = peak_center - tss,                       # genomic sign
            signed_by_strand  = ifelse(gene_strand == "-", -signed_distance, signed_distance),
            distance_kb       = distance / 1000
        )
    head(link_df2, n=2)

    message("Link gene-peak scores with TSS:")
    #table(link_df2$gene_strand, useNA = "ifany")

    message("Peak center and distance to TSS added ...")
    
    #===============================================================================
    # Parse true peak coordinates from the `peak` column
    # peak format assumed: "chrX-start-end"
    
    colnames(link_df2)
    head(link_df$peak)
    
    link_df2_parsed <- link_df2 %>%
        tidyr::separate(peak, into = c("p_chr","p_start","p_end"), sep = "-", remove = FALSE, convert = TRUE) %>%
        mutate(
            peak_width_bp = as.numeric(p_end) - as.numeric(p_start) + 1,
            peak_width_kb = peak_width_bp / 1000
        )
    
    head(link_df2_parsed)
    summary(link_df2_parsed)
    total_peaks <- nrow(link_df2_parsed)
    
    name_df = paste0("ARC_peaks_", p_met, "_ws_", ws)

    if (length(lst_peaks_df) == 0) {
        
        # initialize list with a named element
        lst_peaks_df <- list(link_df2_parsed)
        names(lst_peaks_df) <- name_df
        
    } else {
        
        # append with a name
        lst_peaks_df[[name_df]] <- link_df2_parsed
    }
    
    message("New dataframe added")
    
}

message("Peaks Dataframes: ")
names(lst_peaks_df)

# Combine all dataframes in lst_peaks_df, tagging each with its list name
combined_df <- bind_rows(
    lapply(names(lst_peaks_df), function(nm) {
        df <- lst_peaks_df[[nm]]
        df$group <- nm
        return(df)
    })
)

# covert grp to factor and check
combined_df$group <- factor(combined_df$group, levels = names(lst_peaks_df))
table(combined_df$group)

# Compute group-specific median and mean widths
group_stats <- combined_df |>
    group_by(group) |>
    summarise(
        median_width = median(peak_width_bp, na.rm = TRUE),
        mean_width = mean(peak_width_bp, na.rm = TRUE)
    )

p_hist <- ggplot(combined_df, aes(x = peak_width_bp, fill = group, color = group)) +
    geom_histogram(aes(y = after_stat(density)), bins = 100, alpha = 0.4, position = "identity") +
    geom_density(linewidth = 0.7,  alpha = 0.05) +
    # Median and mean lines per group
    geom_vline(data = group_stats, aes(xintercept = median_width, color = group),
               linetype = "dashed", linewidth = 0.7, show.legend = FALSE) +
    geom_vline(data = group_stats, aes(xintercept = mean_width, color = group),
               linetype = "dotted", linewidth = 0.7, show.legend = FALSE) +
    # Axes and theme
    scale_x_log10(labels = scales::label_number(scale_cut = scales::cut_si("b"))) +
    #coord_flip() +
    theme_minimal(base_size = 10) +
    theme(
        legend.position = "top",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        panel.background = element_rect(fill = "gray95", color = NA),
        plot.background = element_rect(fill = "gray98", color = NA)
    ) +
    labs(title = "Overlaid Peak Width Distributions - CellRanger ATAC Assay", 
         x = "Peak width (bp)", 
         y = "Density",
         fill = "Peak dataset:",
         color = "Peak dataset:")


f_name <- paste0("peak_width_histogram_cellrangerARC.pdf")
ggsave(here::here(plotDir, f_name), p_hist, width = 10, height = 6, dpi = 300)




