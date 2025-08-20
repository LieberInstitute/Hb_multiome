########################################################################
## Explore/Evaluate Signac::LinkPeaks() output = FROM Signac::CallPeaks() --MACS2 Ooutput
## - Make several visualization to evaluate Peak scores
## - Make table with several confidence Peak scores
##
## Authors. CSC
## Date. Aug 20, 2025
## Recommended resources on interactive mode: srun --pty --mem=30GB --x11 bash
## Note. Seurat objects were created with module load conda_R/4.4.x
########################################################################

## I use EnsDb.Hsapiens.v86 for extracting gene names, positions, TSSs, chr locations, etc.
library("EnsDb.Hsapiens.v86")           # Gene annotation (GTF-style)
library("ggplot2")
library("patchwork")
library("tidyverse")
library("stringr")
library("dplyr")
library("scales")
library("here")


## read input arguments
args = commandArgs(trailingOnly = TRUE)
resolution_level <- args[2]
# resolution_level = "Fine" # 42 clusters
# resolution_level = "Broad"  # 8 cell-types
# resolution_level = "Mid" # 8 cell-types

# for testing ( Note only spearman at 5e4 was tested on macs2 peaks )
# resolution_level = "Broad"
# p_met = "spearman"
# w_size = "5e4"

if (length(resolution_level)) {
    
    if (length(p_met) && length(w_size)) {
        message(
            "Processing job for peak-method:\n",
            p_met,
            "\nWindow-size\n",
            w_size
        )
        f_sufix <- paste0(".resolution.", resolution_level, ".", p_met, ".", w_size, ".cells_filtered_2perc")
    } else {
        message("Methodology or Window-Size arguments missed")
        stop()
    }
    
} else {
    
    message("Resolution input arguments missed")
    stop()
    
}
f_sufix

# Check/create directories
input_cvsDir <- here(
    "processed-data",
    "06_peak_calling",
    "02_link_peaks_MACS2"
)
plotDir <- here(
    "plots",
    "06_peak_calling",
    "08_exploratory_peak_scores_MACS2"
)
csvDir <- here(
    "processed-data",
    "06_peak_calling",
    "08_exploratory_peak_scores_MACS2"
)


## Check directories
if (!dir.exists(plotDir)) {
  dir.create(plotDir)
}
if (!dir.exists(csvDir)) {
    dir.create(csvDir)
}

# List all files matching the specific clustering resolution level
# List only files that start with "LHb"
lst_peak_files <- list.files(
    path = input_cvsDir,
    pattern = resolution_level,     # ^ = beginning of string
)
#lst_peak_files = list.files(path = input_cvsDir)
message("Link peak-genes files found:")
lst_peak_files

# for testing at Mid resolution:
#lst_peak_files <- lst_peak_files[5]

#===========================================================================

message("Build TSS (strand-aware) GRanges table...")

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
#     gene_name seqnames       tss gene_strand gene_id        
#     <chr>     <chr>        <dbl> <chr>       <chr>          
# 1 5S_rRNA   chr1     143439605 +           ENSG00000252830
# 2 5S_rRNA   chr11    102057854 +           ENSG00000274097
# 3 5S_rRNA   chr17     37940790 -           ENSG00000277488
nrow(tss_df)
# [1] 56747

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
    
    suffix_subtitle <- paste("Spearman / 5e5d at ", resolution_lev, "resolution. ", ct_name)
    
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
            title = paste("Distribution of peak widths"),
            subtitle = suffix_subtitle,
            x = "Peak width (bp, log scale)",
            y = "Density",
            caption = paste("Dashed = median (", round(median_width), 
                            "bp), dotted = mean (", round(mean_width), "bp)\n", 
                            paste(total_peaks, "total peaks"))
        ) +
        theme_minimal(base_size = 12) +
        theme(
            panel.background = element_rect(fill = "gray95", color = NA),
            plot.background = element_rect(fill = "gray98", color = NA)
        ) +
        coord_flip()
    
    # Peak width vs. distance to TSS
    p_scatter_dist <- ggplot(link_df2_parsed, aes(x = distance_kb, y = peak_width_bp)) +
        geom_point(alpha = 0.25, size = 0.8, color = "grey30") +
        geom_smooth(method = "loess", se = FALSE, color = "darkred") +
        scale_x_log10(labels = label_number(scale_cut = cut_si("b"))) +
        labs(
            title = "Peak width vs distance to TSS",
            subtitle = suffix_subtitle,
            x = "Distance from TSS (kb)",
            y = "Peak width (bp, log scale)",
            caption = paste(total_peaks, "total local peaks")
        ) +
        theme_minimal()
    
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
        labs(title = "Peak width vs correlation score",
             subtitle = suffix_subtitle,
             x = "Correlation Score (pseudo-log scaled)",
             y = "Peak width (bp, log scale)",
             caption = paste(total_peaks, "total local peaks")) +
        theme_minimal()
    
    f_name <- paste0(resolution_lev, "_", ct_name, "_peak_width_histogram_spearman_5e5.pdf")
    ggsave(here::here(plotDir, f_name), p_hist, width = 8, height = 8)
    f_name <- paste0(resolution_lev, "_", ct_name, "_peak_width_vs_distance_spearman_5e5.pdf")
    ggsave(here::here(plotDir, f_name), p_scatter_dist, width = 8, height = 6)
    f_name <- paste0(resolution_lev, "_", ct_name, "_peak_width_vs_score_spearman_5e5.pdf")
    ggsave(here::here(plotDir, f_name), p_scatter_score, width = 8, height = 6)

    print(paste0("Width related plots for ", ct_name, " done!"))    

}


#===========================================================================
# parse Linked peak-gene tables for each cell-type

message("Making plots for ", length(lst_peak_files), " cell-types")

for (ct in lst_peak_files) {
    
    gene_peaks_csv <- here(input_cvsDir, ct)
    
    message("Processing: ", basename(gene_peaks_csv))
    
    if (file.exists(gene_peaks_csv)) {
        link_df <- read.csv(gene_peaks_csv)
        message("File loaded!")
    } else {
        stop(paste("File not found:", gene_peaks_csv))
    }
    
    # Extract text between first and second "_"
    ct_name <- sub("^[^_]*_([^_]*)_.*", "\\1", ct)
    
    #=========================================
    # preapare df
    message(nrow(link_df), " peaks found on ", ct_name, " ...")
    print(head(link_df))
    
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
    
    message("Distance between peaks and TSS added ...")
    
    #=========================================
    
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
    
    write.csv(link_df2,
              file = here(csvDir, 
                          paste0(resolution_level, "_", ct_name,  "_peak_gene_links_with_TSS_and_CC_spearman_5e5.csv")),
              row.names = FALSE)
    
    message("Link Gene-Peak table with TSS distances and CC scores saved!")
    
    #=========================================
    
    make_width_plots(link_df2, resolution_level, ct_name, plotDir)
    
}





#===============================================================================

message("Building plots ...")

## Histogram TSS Scores
pdf(file = here(plotDir, 
                paste0("peak_histogram_distance_TSS", f_sufix, ".pdf")), 
    width = 7, height = 5)

hist(link_df2$distance / 1000, breaks = 100,
     main = "Distance from Peaks to TSS",
     xlab = "Distance (kb)",
     col = "lightblue")
dev.off()


#===============================================================================
# adding exploratory scores
# define high-confidence
# High: score ≥ 0.30 & FDR < 0.05
# Moderate: 0.20 ≤ score < 0.30 & FDR < 0.10
# Exploratory: 0.10 ≤ score < 0.20 & FDR < 0.10 (treat as hypotheses)

colnames(link_df2)
## add adjusted p-value using the Benjamini–Hochberg correction
link_df2 <- link_df2 %>%
    mutate(FDR = p.adjust(pvalue, method = "BH"))

# set tiers due we have confidente peaks < 0.2 
link_df2 <- link_df2 %>%
    mutate(tier = case_when(
        score >= 0.30 & FDR < 0.05 ~ "High (>=0.30, FDR<0.05)",
        score >= 0.20 & FDR < 0.10 ~ "Moderate (0.20–0.30, FDR<0.10)",
        score >= 0.10 & FDR < 0.10 ~ "Exploratory (0.10–0.20, FDR<0.10)",
        TRUE ~ "Discarded"
    ))
# use plain ASCII hyphens
link_df$tier <- gsub("\u2013", "-", link_df2$tier)
head(link_df2)
table(link_df2$tier)

# quick view by distance (kb)
# Keep all data, no filtering of "Discarded" on the plot for visualization purposes
df_plot <- link_df2  

# Count total peaks and how many are below 0.1 to plot on discarted zone
count_below_01 <- sum(df_plot$score < 0.1, na.rm = TRUE)
count_below_02 <- sum((df_plot$score < 0.2 & df_plot$score > 0.1), na.rm = TRUE)
count_below_03 <- sum((df_plot$score < 0.3 & df_plot$score > 0.2), na.rm = TRUE)

g1 <- ggplot(df_plot, aes(x = distance/1000, y = score, color = tier)) +
    geom_point(alpha = 0.5, size = 0.8) +
    # trend over ALL tested links
    geom_smooth(
        data = df_plot,
        aes(x = distance_kb, y = score),
        method = "loess", se = FALSE, span = 0.8, color = "black", linewidth = 0.9) +
    # Threshold lines
    geom_hline(yintercept = 0.3, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 0.2, linetype = "dashed", color = "orange") +
    geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey50") +
    # Labels for thresholds
    annotate("text", x = max(df_plot$distance/1000)*1.02, y = 0.3, 
             label = paste("<0.3 (", count_below_03, " peaks)"), hjust = 0.8, vjust = -0.5, color = "red") +
    annotate("text", x = max(df_plot$distance/1000)*1.02, y = 0.2, 
             label = paste("<0.2 (", count_below_02, " peaks)"), hjust = 0.8, vjust = -0.5, color = "orange") +
    annotate("text", x = max(df_plot$distance/1000)*1.02, y = 0.1, 
             label = paste("<0.1 (", count_below_01, " peaks)"), hjust = 0.8, vjust = -0.5, color = "grey50") +
    # Custom legend with count
    scale_color_manual(
        values = c(
            "High (>=0.30, FDR<0.05)"          = "#b2182b",
            "Moderate (0.20–0.30, FDR<0.10)"   = "#ef8a62",
            "Exploratory (0.10–0.20, FDR<0.10)" = "#67a9cf",
            "Discarded"                        = "grey80"
        ),
        name = paste0("Tier (Count < 0.1: ", count_below_01, ")")
    )  +
    labs(
        x = "Distance from TSS (kb)",
        y = "Correlation score",
        title = "Peak–gene links by tier",
        subtitle = paste(p_met, 
                         format(w_size, scientific = TRUE), "filtered genes < 2%")
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

ggsave(here(plotDir, 
            paste0("peak_exploratory_scores_high_confidence", f_sufix, ".pdf")),
       g1, width = 8, height = 5,
       device = cairo_pdf)


#===============================================================================

## Check number of linked peaks per gene and viceverce
# Number of linked peaks per gene
peaks_per_gene <- link_df2 %>%
    count(gene, name = "n_peaks") %>%
    arrange(desc(n_peaks))

g1 <- ggplot(peaks_per_gene, aes(x = n_peaks)) +
    geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    labs(
        title = paste(p_met, format(w_size, scientific = TRUE)),
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

ggsave(here(plotDir, 
            paste0("link_peak_gene_histograms", f_sufix, ".pdf")),
       combined_plot, width = 8, height = 5)


message("Plots done!!!")


# library("slurmjobs")
# job_single(
#   "06_exploratory_peak_scores",
#   create_shell = TRUE,
#   partition = "katun",
#   memory = "30G",
#   cores = 2,
#   logdir = "logs",
#   command = "Rscript 06_exploratory_peak_scores.R",
#   create_logdir = TRUE
# )

## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
