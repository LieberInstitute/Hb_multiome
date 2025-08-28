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
# w_size = "1e5"
w_sizes = c("1e5", "2.5e4")


for (ws in w_sizes) {
    print(ws)    
}

if (length(p_met) && length(w_size)) {
    message(
        "Processing job for peak-method:\n",
        p_met,
        "\nWindow-size\n",
        w_size
    )
    f_sufix <- paste0(".", p_met, ".", w_size, ".cells_filtered_2perc")
    message("Processing: ", f_sufix)
} else {
    message("Input arguments missed")
    stop()
}

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
        title = paste("Distribution of peak widths - ", f_sufix),
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




f_name <- paste0("peak_width_histogram", f_sufix, ".pdf")
ggsave(here::here(plotDir, f_name), p_hist, width = 8, height = 8)




