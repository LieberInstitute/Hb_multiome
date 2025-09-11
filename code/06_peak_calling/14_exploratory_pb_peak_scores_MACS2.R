########################################################################
## Explore/Evaluate Signac::LinkPeaks() from pseudobulk assays
## - Make distribution plots to evaluate Peak scores
##
## Authors. CSC
## Date. Sep 08, 2025
## Recommended resources mem=30GB
########################################################################

library("EnsDb.Hsapiens.v86")  # Gene annotation (GTF-style), gene names, positions, TSSs, chr locations, etc.
library("ggplot2")
library("patchwork")
library("tidyverse")
library("stringr")
library("dplyr")
library("scales")
library("here")


# Testing spearman at 5e4 on macs2 peaks 
resolution_level = "Mid"
p_met = "spearman"
w_size = "5e5"
FDR_thresh = 0.2
score_thresh = 0.2

if (length(resolution_level)) {
    message("Processing job for peak-method:\n",
            p_met,
            "\nWindow-size\n",
            w_size)
} else {
    message("Input arguments missed")
    stop()
}

# Check/create directories
input_cvsDir <- here(
    "processed-data",
    "06_peak_calling",
    "13_pseudobulk_LinkPeaks_MACS2_split_ct"
)
plotDir <- here(
    "plots",
    "06_peak_calling",
    "14_exploratory_pb_peak_scores_MACS2"
)
csvDir <- here(
    "processed-data",
    "06_peak_calling",
    "14_exploratory_pb_peak_scores_MACS2"
    
)


## Check directories
if (!dir.exists(plotDir)) {
    dir.create(plotDir)
}
if (!dir.exists(csvDir)) {
    dir.create(csvDir)
}

# List all files matching the specific clustering resolution level
pattern = paste0("^", resolution_level, ".*\\.csv$")
lst_peak_files <- list.files(
    path = input_cvsDir,
    pattern = pattern,
)
#lst_peak_files = list.files(path = input_cvsDir)
message("Link peak-genes files found:")
lst_peak_files

# for testing at Mid resolution:
# lst_peak_files <- lst_peak_files[5]

#===========================================================================

message("Build TSS (strand-aware) GRanges table...")

gene_coords <- genes(EnsDb.Hsapiens.v86)

# Build TSS table (strand-aware)
# extracts promoter regions from those gene_coords, upstream = 0 (don’t include any bases before the TSS, 
# downstream = 1 (take exactly one base downstream from the TSS)
tss_coords  <- promoters(gene_coords, upstream = 0, downstream = 1) %>%   # 1bp TSS, respects strand
    keepStandardChromosomes(pruning.mode = "coarse")
# get a 1 bp range that pinpoints the TSS position for every gene, regardless of strand orientation

# match UCSC-style peaks (chr1, chr2, etc.)
seqlevelsStyle(tss_coords) <- "UCSC"
tss_raw <- as.data.frame(tss_coords)
head(tss_raw)
has_biotype <- "gene_biotype" %in% colnames(tss_raw)
head(has_biotype)
message("Biotypes included:")
unique(tss_raw$gene_biotype)

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
head(tss_df)
nrow(tss_df) # [1] 56747

## check duplicates
dup_pairs <- tss_df %>%
    count(gene_name, seqnames, name = "n") %>%
    filter(n > 1)

if (nrow(dup_pairs) > 0) {
    print(head(dup_pairs, 10))
    stop("Non-unique (gene_name, seqnames) in tss_df: ", nrow(dup_pairs), " duplicates.")
}

message("Build TSS completed...")


#===========================================================================
# add function to make plots related with peak's width

make_width_plots <- function(
        link_df2,
        resolution_lev,
        ct_name,
        plotDir
) {
    
    print(paste0("Processing plots for ", resolution_lev, " for ", ct_name, " cell_type"))
    
    # Parse true peak coordinates from the `peak` column
    # peak format assumed: "chrX-start-end"
    colnames(link_df2)
    link_df2_parsed <- link_df2 %>%
        tidyr::separate(peak, into = c("p_chr","p_start","p_end"), sep = "-", remove = FALSE, convert = TRUE) %>%
        mutate(
            peak_width_bp = as.numeric(p_end) - as.numeric(p_start) + 1,
            peak_width_kb = peak_width_bp / 1000
        )
    
    head(link_df2_parsed)
    summary(link_df2_parsed)
    total_peaks <- nrow(link_df2_parsed)
    
    #=========================================
    
    # Overall peak width distribution (log scale)
    # Compute summary stats
    median_width <- median(link_df2_parsed$peak_width_bp, na.rm = TRUE)
    mean_width   <- mean(link_df2_parsed$peak_width_bp, na.rm = TRUE)
    
    suffix_subtitle <- paste(resolution_lev, "resolution. ", ct_name)
    
    p_hist <- ggplot(link_df2_parsed, aes(x = peak_width_bp)) +
        # histogram as horizontal bars
        geom_histogram(
            aes(y = after_stat(density)),  # normalize for density overlay
            bins = 100, fill = "steelblue", color = "white", alpha = 0.6
        ) +
        # density curve
        geom_density(color = "darkred", linewidth = 1) +
        geom_vline(xintercept = median_width, color = "black", linetype = "dashed", linewidth = 0.8) +
        geom_vline(xintercept = mean_width, color = "orange", linetype = "dotted", linewidth = 0.8) +
        # log scale for widths
        scale_x_log10(labels = label_number(scale_cut = cut_si("b"))) +
        labs(
            title = paste("Distribution of LinkPeak widths"),
            subtitle = suffix_subtitle,
            x = "Peak width (bp, log scale)",
            y = "Density",
            caption = paste("Dashed = median (", round(median_width), 
                            "bp), dotted = mean (", round(mean_width), "bp)\n", 
                            paste(total_peaks, "total peaks"))
        ) +
        theme_minimal(base_size = 12) +
        # theme(
        #     panel.background = element_rect(fill = "gray95", color = NA),
        #     plot.background = element_rect(fill = "gray98", color = NA)
        # ) +
        coord_flip()
    
    # Peak width vs. correlation score
    # log10 creates NaN/Inf, remove those rows to avoid warnings
    df <- link_df2_parsed %>%
        mutate(score = as.numeric(score),
               peak_width_bp = as.numeric(peak_width_bp)) %>%
        filter(is.finite(score), is.finite(peak_width_bp), peak_width_bp > 0)
    
    p_scatter_score <- ggplot(df, aes(x = score, y = peak_width_bp)) +
        geom_point(alpha = 0.25, size = 0.8, color = "grey30") +
        geom_smooth(method = "loess", se = FALSE, color = "darkred") +
        #scale_x_continuous(trans = pseudo_log_trans(base = 10, sigma = 0.01)) +
        scale_y_log10() +
        labs(title = "LinkPeak width vs correlation scores",
             subtitle = suffix_subtitle,
             x = "Correlation Score (pseudo-log scaled)",
             y = "Peak width (bp, log scale)",
             caption = paste(total_peaks, "total local peaks")) +
        theme_minimal()
    
    f_name <- paste0(resolution_lev, "_", ct_name, "_pb_linkpeaks_width_histogram.pdf")
    ggsave(here::here(plotDir, f_name), p_hist, width = 8, height = 8)
    f_name <- paste0(resolution_lev, "_", ct_name, "_pb_linkpeaks_width_vs_score.pdf")
    ggsave(here::here(plotDir, f_name), p_scatter_score, width = 8, height = 6)
    
    print(paste0("Width related plots for ", ct_name, " done!"))    
    
}


#===========================================================================

make_exploratory_plots <- function(
        link_df2,
        resolution_lev,
        ct_name,
        plotDir,
        suffix_name = ""
){
    
    # adding exploratory scores
    # define high-confidence
    # High: score ≥ 0.3
    # Moderate: 0.2 ≤ score < 0.3
    # Exploratory: 0.1 ≤ score < 0.2
    
    ## testing: (link_df2, resolution_level, ct_name, plotDir)
    # link_df2 = link_df2
    # resolution_lev = resolution_level
    # ct_name = ct_name
    # plotDir = plotDir
    
    print(paste0("Processing plots for ", resolution_lev, " for ", ct_name, " cell_type"))
    suffix_subtitle <- paste(resolution_lev, "resolution. ", ct_name)
    total_peaks <- nrow(link_df2)
    
    # set tiers in 4 thresholds - all links found
    link_df2 <- link_df2 %>%
        mutate(tier = case_when(
            score >= 0.3 ~ "High (>=0.3)",
            score >= 0.2 ~ "Moderate (0.2–0.3)",
            score >= 0.1 ~ "Exploratory (0.1–0.2)",
            TRUE ~ "Discarded"
        ))
    # use plain ASCII hyphens
    # link_df$tier <- gsub("\u2013", "-", link_df2$tier)
    table(link_df2$tier)
    df_plot <- link_df2  
    
    # Count total peaks and how many are below 0.1 to plot on discarted zone
    count_below_01 <- sum(df_plot$score < 0.1, na.rm = TRUE)
    count_below_02 <- sum((df_plot$score < 0.2 & df_plot$score > 0.1), na.rm = TRUE)
    count_below_03 <- sum((df_plot$score < 0.3 & df_plot$score > 0.2), na.rm = TRUE)
    count_above_03 <- sum((df_plot$score > 0.3), na.rm = TRUE)
    
    g1 <- ggplot(df_plot, 
                 aes(x = signed_by_strand/1000, y = score, color = tier),
                 method = "loess", se = FALSE, span = 0.8, color = "black", linewidth = 0.9) +
        geom_point(alpha = 0.5, size = 0.5) +
        # Threshold lines
        geom_hline(yintercept = 0.3, linetype = "dashed", color = "red") +
        geom_hline(yintercept = 0.2, linetype = "dashed", color = "darkorange") +
        geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey50") +
        # Labels for thresholds
        annotate("label", x = max(df_plot$distance/1000)*1.02, y = 0.4, 
                 label = paste(">0.3 (", count_above_03, " links)"), 
                 size=3, hjust = 0.8, vjust = +1, color = "black", fill = "white") +
        annotate("label", x = max(df_plot$distance/1000)*1.02, y = 0.3, 
                 label = paste(">0.2 and <0.3 (", count_below_03, " links)"), 
                 size=3, hjust = 0.8, vjust = +1, color = "black", fill = "white") +
        annotate("label", x = max(df_plot$distance/1000)*1.02, y = 0.2, 
                 label = paste(">0.1 and <0.2 (", count_below_02, " links)"), 
                 size=3, hjust = 0.8, vjust = +1, color = "black", fill = "white") +
        annotate("label", x = max(df_plot$distance/1000)*1.02, y = 0.1, 
                 label = paste("<0.1 (", count_below_01, " links)"), 
                 size=3, hjust = 0.8, vjust = +1, color = "black", fill = "white") +
        # Custom legend with count
        scale_color_manual(
            values = c(
                "High (>=0.3)" = "red",
                "Moderate (0.2–0.3)" = "deeppink4",
                "Exploratory (0.1–0.2)" = "darkorange",
                "Discarded" = "grey50"
            ) 
        )  +
        labs(
            title = "LinkPeak scores by tier",
            subtitle = paste(suffix_subtitle, " - spearman at 1e5", suffix_name),
            x = "Distance from TSS (kb)",
            y = "Correlation score",
        ) +
        theme_minimal() + 
        theme(legend.position = "none")
    # +
    # theme(legend.position = "bottom",
    #       legend.text = element_text(size = 8),
    #       legend.title = element_text(size = 9)) +
    # guides(color = guide_legend(override.aes = list(size = 3)))
    
    f_name <- paste0(resolution_lev, "_", ct_name, "_pb_linkpeaks_tier_scores", suffix_name, ".pdf")
    ggsave(here(plotDir, f_name), g1, width = 8, height = 5, device = cairo_pdf)
    
    message("Exploratory plots for ", ct_name, " done!")
    
}


#===============================================================================

## Check number of linked peaks per gene and vice-verse
# Number of linked peaks per gene
make_peak_gene_histograms <- funcion(
    link_df2 = link_df2,
    resolution_lev,
    ct_name,
    plotDir 
){
    
    
    ## Histogram TSS Scores
    pdf(file = here(plotDir, 
                    paste0(resolution_level, "_", ct_name, "_pb_histogram_TSS_dist.pdf")), 
        width = 7, height = 5)
    
    hist(link_df2$distance / 1000, breaks = 100,
         main = "LinkPeaks distance from Peaks to TSS",
         xlab = "Distance (kb)",
         col = "lightblue")
    dev.off()
    
    peaks_per_gene <- link_df2 %>%
        count(gene, name = "n_peaks") %>%
        arrange(desc(n_peaks))
    
    g1 <- ggplot(peaks_per_gene, aes(x = n_peaks)) +
        geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
        scale_x_continuous(breaks = scales::pretty_breaks()) +
        labs(
            title = "Local LinkPeaks by cell-type (pseudobulk)", 
            subtitle = "Peaks per gene",
            x = "Number of linked peaks per gene",
            y = "Number of genes"
        ) +
        theme_minimal()
    
    # Number of linked genes per peak
    genes_per_peak <- link_df2 %>%
        count(peak, name = "n_genes") %>%
        arrange(desc(n_genes))
    
    g2 <- ggplot(genes_per_peak, aes(x = n_genes)) +
        geom_histogram(binwidth = 1, fill = "firebrick", color = "white") +
        scale_x_continuous(breaks = scales::pretty_breaks()) +
        labs(
            subtitle = "Genes per peak",
            x = "Number of linked genes per peak",
            y = "Number of peaks"
        ) +
        theme_minimal()
    
    combined_plot <- g1 + g2
    combined_plot
    
    f_name <- paste0(resolution_lev, "_", ct_name, "_pb_link_peak_gene_histograms.pdf")
    ggsave(here(plotDir, f_name),
           combined_plot, width = 7, height = 7)
    
    message("Gene-Peak plots for ", ct_name, " done!")
    
}


#===========================================================================
# parse Linked peak-gene tables for each cell-type

message("Making plots for ", length(lst_peak_files), " cell-types")

for (ct in lst_peak_files) {
    # testing
    # ct =  lst_peak_files[6]
    
    gene_peaks_csv <- here(input_cvsDir, ct)
    
    message("Processing: ", basename(gene_peaks_csv))
    
    if (file.exists(gene_peaks_csv)) {
        link_df <- read.csv(gene_peaks_csv)
        message("File loaded!")
    } else {
        stop(paste("File not found:", gene_peaks_csv))
    }
    
    # Extract text between first and second "_"
    ct_name <- sub("^[^_]*_([^_]*)_.*", "\\1", ct) #ge. [1] "LHb.2.7"
    
    #=========================================
    # prepare df
    message("Processing ", nrow(link_df), " peaks found in ", ct_name, " ...")
    
    link_df <- link_df |>
        mutate(
            gene     = trimws(as.character(gene)),
            seqnames = as.character(seqnames),
            start    = as.numeric(start),
            end      = as.numeric(end)
            #peak_called_in = peak_called_in
        )
    head(link_df)
    
    #=========================================
    
    message("Computing distance between peaks and TSS ...")
    
    # Join by gene + chromosome to avoid many-to-many 
    colnames(link_df)
    colnames(tss_df)
    link_df2 <- link_df %>%
        left_join(tss_df, by = c("gene" = "gene_name", "seqnames" = "seqnames"))
    colnames(link_df2)
    head(link_df2, n = 3)
    
    # drop rows with no TSS match
    n_before <- nrow(link_df2)
    link_df2 <- link_df2 %>% filter(!is.na(tss))
    message("Dropped ", n_before - nrow(link_df2), " rows with no TSS match.")
    
    message("Computing Peak center and distance to TSS ...")
    
    link_df2 <- link_df2 %>%
        mutate(
            peak_center       = (start + end) / 2,
            distance          = abs(peak_center - tss),
            signed_distance   = peak_center - tss,                       # genomic sign
            signed_by_strand  = ifelse(gene_strand == "-", -signed_distance, signed_distance),
            distance_kb       = distance / 1000
        )
    head(link_df2, n=3)

    message("Link gene-peak scores with TSS:")
    summary(link_df2)
    table(link_df2$gene_strand, useNA = "ifany")
    
    message("Peak center and distance to TSS added ...")
    
    #===============================================================================
    
    f_name <- here(csvDir, paste0(resolution_level, "_", ct_name,  "_links_with_TSS_CC_FDR_spearman_5e5.csv"))
    write.csv(link_df2, file = f_name, row.names = FALSE)
    
    message("Link Gene-Peak table with TSS distances and CC scores saved!")
    
    #=========================================
    
    make_exploratory_plots(link_df2, resolution_level, ct_name, plotDir)
    make_peak_gene_histograms(link_df2, resolution_level, ct_name, plotDir)
    make_width_plots(link_df2, resolution_level, ct_name, plotDir)
    
    ## filter pRE links
    filtered_links <- link_df2 |>
        dplyr::filter(FDR <= FDR_thresh, score >= score_thresh)
        #dplyr::filter(FDR <= 0.1, abs(score) >= 0.3)
    nrow(filtered_links)    
    head(filtered_links)    
    
    suffix = paste0("FDR", FDR_thresh, "_score", score_thresh)
    f_name <- here(csvDir, paste0(resolution_level, "_", ct_name,
                                           "_links_FDR", suffix, ".csv"))
    write.csv(filtered_links, f_name, row.names = FALSE)
    message("Links passing threshold saved: ", nrow(filtered_links))
    
    make_exploratory_plots(filtered_links, resolution_level, ct_name, plotDir, suffix)
    
    
}


message("Plots done!!!")



## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
