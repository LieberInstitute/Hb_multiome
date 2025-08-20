########################################################################
## Explore/Evaluate Signac::LinkPeaks() output = FROM CELLRANGER OUTPUT
## - Make several visualization to evaluate Peak scores
## - Make table with several confidence Peak scores
##
## Authors. CSC
## Date. Aug 11, 2025
## Recommended resources on interactive mode: srun --pty --mem=30GB --x11 bash
## Note. Seurat objects were created with module load conda_R/4.4.x
########################################################################

## I use EnsDb.Hsapiens.v86 for extracting gene names, positions, TSSs, chr locations, etc.
library("EnsDb.Hsapiens.v86")           # Gene annotation (GTF-style)
library("ggplot2")
library("patchwork")
library("tidyverse")
library("dplyr")
library("scales")
library("here")


## read input arguments
args = commandArgs(trailingOnly = TRUE)
p_met <- args[2]
w_size <- args[4]
# 0: pearson, 1e5
# 1: pearson, 5e4
# 2: spearman, 1e5
# 3: spearman, 2.5e4

# for testing
# p_met = "spearman"
# w_size = "2.5e4"

## p_met:
# pearson -> peak-scores<0.2: likely due scATAC counts are ultra‑sparse; scRNA is zero‑inflated. Pearson r’s of 0.05–0.2 are common even for real links
# spearman -> as enhancer → gene relationships aren’t strictly linear; Pearson seems to underestimates. I will try spearman, more robust to nonlinearity/zeros

if (length(p_met) && length(w_size)) {
    message(
        "Processing job for peak-method:\n",
        p_met,
        "\nWindow-size\n",
        w_size
    )
    # ## Use numeric comparison first, then assign string labels
    # w_size_label <- case_when(
    #     isTRUE(all.equal(w_size, 25000))  ~ "2.5e4",
    #     isTRUE(all.equal(w_size, 50000))  ~ "5e4",
    #     isTRUE(all.equal(w_size, 100000)) ~ "1e5",
    #     TRUE                              ~ "00"
    # )
    #f_sufix <- paste0(".", p_met, ".", format(w_size_label, scientific = TRUE), ".cells_filtered_2perc")
    f_sufix <- paste0(".", p_met, ".", w_size, ".cells_filtered_2perc")
    message("Processing: ", f_sufix)
} else {
    message("Input arguments missed")
    stop()
}


# Check/create directories
input_cvsDir <- here(
    "processed-data",
    "06_peak_calling",
    "00_link_peaks"
)
plotDir <- here(
  "plots",
  "06_peak_calling",
  "06_exploratory_peak_scores"
)
csvDir <- here(
    "processed-data",
    "06_peak_calling",
    "06_exploratory_peak_scores"
)


## Check directories
if (!dir.exists(plotDir)) {
  dir.create(plotDir)
}
if (!dir.exists(csvDir)) {
    dir.create(csvDir)
}

# load link peak-gene csv
gene_peaks_csv <- here(input_cvsDir, paste0("all_peak_gene_links", f_sufix, ".csv"))
if (file.exists(gene_peaks_csv)) {
    link_df <- read.csv(gene_peaks_csv)
    message("File loaded!")
} else {
    stop(paste("File not found:", gene_peaks_csv))
}

## for testing: ================================================================
#link_df <- read.csv(file = here(input_cvsDir, "all_peak_gene_links.spearman.5e4_test.csv"))
## for testing: ================================================================

colnames(link_df)
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

## Build TSS (strand-aware) GRanges table and 
#  compute the distance between each peak and its linked gene's TSS
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
#unique(tss_raw$gene_biotype)

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

# Join by gene + chromosome to avoid many-to-many 
colnames(link_df)
colnames(tss_df)
link_df2 <- link_df %>%
    left_join(tss_df, by = c("gene" = "gene_name", "seqnames" = "seqnames"))
head(link_df2, n = 3)
#     seqnames  start     end width strand      score  gene               peak
# 1     chr1 921198 1001138 79941      * 0.06070581 ISG15 chr1-920766-921629
# 2     chr1 960688 1000172 39485      * 0.07694616  HES4 chr1-960318-961058
# 3     chr1 960688 1001138 40451      * 0.06605429 ISG15 chr1-960318-961058
#     zscore     pvalue     tss gene_strand         gene_id
# 1 2.052715 0.02005013 1001138           + ENSG00000187608
# 2 1.758599 0.03932286 1000172           - ENSG00000188290
# 3 2.016967 0.02184947 1001138           + ENSG00000187608

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
# seqnames  start     end width strand      score  gene               peak
# 1     chr1 921198 1001138 79941      * 0.06070581 ISG15 chr1-920766-921629
# 2     chr1 960688 1000172 39485      * 0.07694616  HES4 chr1-960318-961058
#   zscore     pvalue     tss gene_strand         gene_id peak_center distance
# 1 2.052715 0.02005013 1001138           + ENSG00000187608      961168    39970
# 2 1.758599 0.03932286 1000172           - ENSG00000188290      980430    19742
# signed_distance signed_by_strand distance_kb
# 1          -39970           -39970      39.970
# 2          -19742            19742      19.742
message("Link gene-peak scores with TSS:")
summary(link_df2)
#sum(link_df2$distance > 1e5)  # should be ~0 if you used LinkPeaks(..., distance=1e5)
table(link_df2$gene_strand, useNA = "ifany")
    # -    + 
    # 2714 2786 

message("Peak center and distance to TSS added ...")


#===============================================================================
## plot peak width

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

# compare to existing 'width' column
# - This will likely be FALSE for many rows; that's expected here
# - table(link_df2_parsed$width == link_df2_parsed$peak_width_bp, useNA = "ifany")

# Overall peak width distribution (log scale)
# Compute summary stats
median_width <- median(link_df2_parsed$peak_width_bp, na.rm = TRUE)
mean_width   <- mean(link_df2_parsed$peak_width_bp, na.rm = TRUE)

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
        title = paste("Distribution of peak widths", f_sufix),
        subtitle = paste(total_peaks, "total global peaks"),
        x = "Peak width (bp, log scale)",
        y = "Density",
        caption = paste("Dashed = median (", round(median_width), 
                        "bp), dotted = mean (", round(mean_width), "bp)")
    ) +
    theme_minimal(base_size = 12) +
    theme(
        panel.background = element_rect(fill = "gray95", color = NA),
        plot.background = element_rect(fill = "gray98", color = NA)
    ) +
    coord_flip()

f_name <- paste0("peak_width_histogram", f_sufix, ".png")
ggsave(here::here(plotDir, f_name), p_hist, width = 8, height = 8, dpi = 300)


# Peak width vs. distance to TSS
p_scatter_dist <- ggplot(link_df2_parsed, aes(x = distance_kb, y = peak_width_bp)) +
    geom_point(alpha = 0.25, size = 0.8, color = "grey30") +
    geom_smooth(method = "loess", se = FALSE, color = "darkred") +
    scale_x_log10(labels = label_number(scale_cut = cut_si("b"))) +
    labs(
        title = paste("Peak width vs distance to TSS", f_sufix),
        subtitle = paste(total_peaks, "total global peaks"),
        x = "Distance from TSS (kb)",
        y = "Peak width (bp, log scale)"
    ) +
    theme_minimal()

f_name <- paste0("peak_width_vs_distance", f_sufix, ".pdf")
ggsave(here::here(plotDir, f_name), p_scatter_dist, width = 8, height = 6)


#===============================================================================

message("Building plots ...")

## Histogram TSS Scores
pdf(file = here(plotDir, 
                paste0("histogram_scores", f_sufix, ".pdf")), 
    width = 7, height = 5)

hist(link_df2$distance / 1000, breaks = 100,
     main = "Distance from Peaks to TSS",
     xlab = "Distance (kb)",
     col = "lightblue")
dev.off()

## Correlation vs Distance with smoothing
g1 <- ggplot(link_df2, aes(distance_kb, score)) +
    geom_point(alpha = 0.3, color = "steelblue") +
    geom_hline(yintercept = 0.2, linetype = "dashed", color = "red") +
    labs(
        x = "Distance from TSS (kb)",
        y = paste("Correlation Score", p_met),
        title = "Peak-Gene Correlation vs.Distance"
    ) + geom_smooth(method = "loess", se = FALSE, color = "darkred") +
    theme_minimal()

ggsave(here(plotDir, 
            paste0("distribution_scores", f_sufix, ".pdf")),
            g1, width = 8, height = 5)


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
            paste0("exploratory_scores_high_confidence", f_sufix, ".pdf")),
       g1, width = 8, height = 5,
       device = cairo_pdf)

#===============================================================================


## filtered peaks
write.csv(link_df2,
          file = here(csvDir, paste0("peak_gene_links_with_TSS_and_CC", f_sufix, ".csv")),
          row.names = FALSE)

message("Link Gene-Peak table with TSS distances and CC scores saved!")


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
            paste0("histograms_link_peak_genes_peaks", f_sufix, ".pdf")),
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
