########################################################################
## Compare width peaks
## Input: LinkPeaks() results using CallPeaks (macs2) -- local
## Output a combined histogram
##
## Authors. CSC
## Date. Aug 28, 2025
## Recommended resources on interactive mode: srun --pty --mem=30GB --x11 bash
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


resolution_level = "Mid"

# Check/create directories
input_cvsDir <- here(
    "processed-data",
    "06_peak_calling",
    "02_link_peaks_MACS2", 
    "old_Mid" # <= temporal dir for testing / should be removed
)
plotDir <- here(
    "plots",
    "06_peak_calling",
    "10_compare_link_peak_gene_distributions_MACS2"
)
output_cvsDir <- here(
    "processed-data",
    "06_peak_calling",
    "10_compare_link_peak_gene_distributions_MACS2"
)

## Check directories
if (!dir.exists(plotDir)) {
    dir.create(plotDir)
}
if (!dir.exists(output_cvsDir)) {
    dir.create(output_cvsDir)
}


## load linked peaks from Signac::CallPeaks() and compute TSS stats 

# List all files matching the specific clustering resolution level
pattern = paste0("^", resolution_level, ".*\\.csv$")
lst_link_files <- list.files(
    path = input_cvsDir,
    pattern = pattern
)

message("Link peak-genes files found:")
lst_link_files

## list of file names with link peaks-to-genes information
lst_links_df <- list()

message("Processing peak width histograms ...")

for (ct_links in lst_link_files) {
    # ct_links = lst_link_files[1]

    message("Processing ", ct_links)
    
    # load link peak-gene csv
    gene_peaks_csv <- here(input_cvsDir, ct_links)
    if (file.exists(gene_peaks_csv)) {
        links_df <- read.csv(gene_peaks_csv)
        message("File loaded!")
    } else {
        stop(paste("File not found:", gene_peaks_csv))
    }
    
    ##==========================================================================
    ## prepare df with peak ranges and width
    
    colnames(links_df)
    head(links_df, n=3)
    nrow(links_df)
    
    links_df <- links_df |>
        mutate(
            gene     = trimws(as.character(gene)),
            seqnames = as.character(seqnames),
            start    = as.numeric(start),
            end      = as.numeric(end)
        )
    head(links_df)
    
    message("Computing distance between peaks and TSS ...")
    
    ## Build TSS (strand-aware) GRanges table and compute the distance between each peak and its linked gene's TSS
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
    colnames(links_df)
    colnames(tss_df)
    links_df2 <- links_df %>%
        left_join(tss_df, by = c("gene" = "gene_name", "seqnames" = "seqnames"))
    head(links_df2, n = 3)
    
    # drop rows with no TSS match
    n_before <- nrow(links_df2)
    links_df2 <- links_df2 %>% filter(!is.na(tss))
    message("Dropped ", n_before - nrow(links_df2), " rows with no TSS match.")
    
    nrow(links_df2)
    
    message("Distance between peaks and TSS added ...")
    
    message("Computing Peak center and distance to TSS ...")
    
    links_df2 <- links_df2 %>%
        mutate(
            peak_center       = (start + end) / 2,
            distance          = abs(peak_center - tss),
            signed_distance   = peak_center - tss,                       # genomic sign
            signed_by_strand  = ifelse(gene_strand == "-", -signed_distance, signed_distance),
            distance_kb       = distance / 1000
        )
    head(links_df2, n=2)
    
    message("Link gene-peak scores with TSS:")
    #table(links_df2$gene_strand, useNA = "ifany")
    
    message("Peak center and distance to TSS added ...")
    
    #===============================================================================
    # Parse true peak coordinates from the `peak` column
    # peak format assumed: "chrX-start-end"
    
    colnames(links_df2)
    head(links_df$peak)
    
    links_df2_parsed <- links_df2 %>%
        tidyr::separate(peak, into = c("p_chr","p_start","p_end"), sep = "-", remove = FALSE, convert = TRUE) %>%
        mutate(
            peak_width_bp = as.numeric(p_end) - as.numeric(p_start) + 1,
            peak_width_kb = peak_width_bp / 1000
        )
    
    head(links_df2_parsed, n=2)
    total_peaks <- nrow(links_df2_parsed)
    
    # extract name for named list    
    name_df = sub("^([^_]+_[^_]+)_.*", "\\1", ct_links)
    # build lst of df
    lst_links_df[[name_df]] <- links_df2_parsed
    
    message("New dataframe added")
    
}


message("Peaks Dataframes: ")
names(lst_links_df)

# Combine all dataframes in lst_links_df, tagging each with its list name
combined_df <- bind_rows(
    lapply(names(lst_links_df), function(nm) {
        df <- lst_links_df[[nm]]
        df$group <- nm
        return(df)
    })
)

# covert grp to factor and check
combined_df$group <- factor(combined_df$group, levels = names(lst_links_df))
table(combined_df$group)

# Compute group-specific median and mean widths
group_stats <- combined_df |>
    group_by(group) |>
    summarise(
        median_width = median(peak_width_bp, na.rm = TRUE),
        mean_width = mean(peak_width_bp, na.rm = TRUE)
    )

p_hist <- ggplot(combined_df, aes(x = peak_width_bp, fill = group, color = group)) +
    # geom_histogram(
    #     aes(y = after_stat(density)), 
    #     bins = 50, alpha = 0.4, position = "identity") +
    geom_density(
        aes(y = after_stat(density), group = group),
        linewidth = 0.7,  alpha = 0.1) +
    # Median and mean lines per group
    geom_vline(data = group_stats, aes(xintercept = median_width, color = group),
               linetype = "dashed", linewidth = 0.3, show.legend = FALSE) +
    # geom_vline(data = group_stats, aes(xintercept = mean_width, color = group),
    #           linetype = "dotted", linewidth = 0.7, show.legend = FALSE) +
    # Axes and theme
    scale_x_log10(
        labels = scales::label_number(scale_cut = scales::cut_si("b")),
        expand = expansion(mult = c(0.01, 0.01))
    ) +
    coord_cartesian(xlim = c(20, NA)) + # correspond with left tail cut off
    theme_minimal(base_size = 10) +
    theme(
        #legend.position = "top",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        panel.background = element_rect(fill = "gray95", color = NA),
        plot.background = element_rect(fill = "gray98", color = NA)
    ) +
    labs(title = "Overlaid Peak Width Distributions - Raw LinkPeaks (macs2)",
         subtitle = paste0(resolution_level, " cell-types"),
         x = "Peak width (bp, log scale)", 
         y = "Density",
         fill = "Correlation Window",
         color = "Correlation Window",
         caption = paste0(
             "Dashed lines = median peak width per group\n",
             "Total peaks: ", scales::comma(nrow(combined_df))
         ))


f_name <- paste0(resolution_level, "_link_peak_width_histogram_macs2.pdf")
ggsave(here(plotDir, f_name), p_hist, width = 6, height = 6, dpi = 300)


message("Plots done!")

summary_stats <- combined_df |>
    group_by(group) |>
    summarise(
        n = n(),
        median_width_bp = round(median(peak_width_bp, na.rm = TRUE), 1),
        Q1 = round(quantile(peak_width_bp, 0.25, na.rm = TRUE), 1),
        Q3 = round(quantile(peak_width_bp, 0.75, na.rm = TRUE), 1),
        IQR = Q3 - Q1
    ) |>
    arrange(median_width_bp)

summary_stats

f_name <- paste0(resolution_level, "_link_peak_width_summary_macs2.csv")
write.csv(summary_stats, here(output_cvsDir, f_name), row.names = FALSE)

message("Summary done!")

## Define Enhancer Class (Promoter vs Distal)

# use abs(signed_by_strand) = Distance from peak center to TSS, it is strand-aware
# This metric tells how far the linked peak is from the TSS of the gene
# note, peak_width_kb is the length of the peak itself
combined_df <- combined_df |>
    mutate(
        enhancer_class = case_when(
            abs(signed_by_strand) <= 1000 ~ "Promoter (<1kb)",
            abs(signed_by_strand) <= 10000 ~ "Proximal (1–10kb)",
            TRUE ~ "Distal (>10kb)"
        )
    )
## sort by distance to TSS
combined_df$enhancer_class <- factor(
    combined_df$enhancer_class,
    levels = c("Promoter (<1kb)", "Proximal (1–10kb)", "Distal (>10kb)")
)

p1 <- ggplot(combined_df, aes(x = enhancer_class, y = peak_width_bp, fill = enhancer_class)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.2) +
    geom_boxplot(width = 0.1, outlier.size = 0.3, alpha = 0.5) +  # optional overlay
    #scale_y_log10() +
    #facet_wrap(~ group) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "none"
    ) +
    labs(
        title = "Peak width by enhancer class",
        x = "Enhancer category (distance to TSS)",
        y = "Peak width (bp)"
    )

f_name <- paste0(resolution_level, "_link_peak_width_by_enhancer_class_macs2.pdf")
ggsave(here::here(plotDir, f_name), p1, width = 8, height = 6, dpi = 300)


# library("slurmjobs")
# job_single(
#   "06_exploratory_peak_scores",
#   create_shell = TRUE,
#   partition = "katun",
#   memory = "30G",
#   cores = 2,
#   logdir = "logs",
#   command = "Rscript 06_exploratory_peak_scores",
#   create_logdir = TRUE
# )

## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
