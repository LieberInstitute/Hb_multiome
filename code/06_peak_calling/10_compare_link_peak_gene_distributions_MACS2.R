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


resolution_level = "Broad"

# Check/create directories
input_cvsDir <- here(
    "processed-data",
    "06_peak_calling",
    "02_link_peaks_MACS2"
)
plotDir <- here(
    "plots",
    "06_peak_calling",
    "10_compare_peak_distributions_MACS2"
)

## Check directories
if (!dir.exists(plotDir)) {
  dir.create(plotDir)
}

## load linked peaks from Signac::CallPeaks() and compute TSS stats 

# List all files matching the specific clustering resolution level
lst_link_files <- list.files(
    path = input_cvsDir,
    pattern = paste0("^", resolution_level)
)

#lst_peak_files = list.files(path = input_cvsDir)
message("Link peak-genes files found:")
lst_link_files

## list of data frames with peak coordinates, central values and peak width
lst_links_df <- list()

message("Processing peak width histograms ...")

macs2_resolutions = c("Broad")  


for (clust_res in macs2_resolutions) {
    # clust_res = macs2_resolutions[1]

    print(clust_res)    
    # load link peak-gene csv
    gene_peaks_csv <- here(input_cvsDir, clust_res)
    if (file.exists(gene_peaks_csv)) {
        peaks_df <- read.csv(gene_peaks_csv)
        message("File loaded!")
    } else {
        stop(paste("File not found:", gene_peaks_csv))
    }
    
    ##==========================================================================
    ## prepare df with peak ranges and width
    
    colnames(peaks_df)
    head(peaks_df, n=3)
    nrow(peaks_df)
    
    peaks_df <- peaks_df |>
        mutate(
            gene     = trimws(as.character(gene)),
            seqnames = as.character(seqnames),
            start    = as.numeric(start),
            end      = as.numeric(end)
        )
    head(peaks_df)
    
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
    colnames(peaks_df)
    colnames(tss_df)
    peaks_df2 <- peaks_df %>%
        left_join(tss_df, by = c("gene" = "gene_name", "seqnames" = "seqnames"))
    head(peaks_df2, n = 3)
    
    # drop rows with no TSS match
    n_before <- nrow(peaks_df2)
    peaks_df2 <- peaks_df2 %>% filter(!is.na(tss))
    message("Dropped ", n_before - nrow(peaks_df2), " rows with no TSS match.")
    
    nrow(peaks_df2)
    
    message("Distance between peaks and TSS added ...")
    
    message("Computing Peak center and distance to TSS ...")
    
    peaks_df2 <- peaks_df2 %>%
        mutate(
            peak_center       = (start + end) / 2,
            distance          = abs(peak_center - tss),
            signed_distance   = peak_center - tss,                       # genomic sign
            signed_by_strand  = ifelse(gene_strand == "-", -signed_distance, signed_distance),
            distance_kb       = distance / 1000
        )
    head(peaks_df2, n=2)
    
    message("Link gene-peak scores with TSS:")
    #table(peaks_df2$gene_strand, useNA = "ifany")
    
    message("Peak center and distance to TSS added ...")
    
    #===============================================================================
    # Parse true peak coordinates from the `peak` column
    # peak format assumed: "chrX-start-end"
    
    colnames(peaks_df2)
    head(peaks_df$peak)
    
    peaks_df2_parsed <- peaks_df2 %>%
        tidyr::separate(peak, into = c("p_chr","p_start","p_end"), sep = "-", remove = FALSE, convert = TRUE) %>%
        mutate(
            peak_width_bp = as.numeric(p_end) - as.numeric(p_start) + 1,
            peak_width_kb = peak_width_bp / 1000
        )
    
    head(peaks_df2_parsed)
    summary(peaks_df2_parsed)
    total_peaks <- nrow(peaks_df2_parsed)
    
    name_df = paste0("ARC_peaks_", p_met, "_res_", clust_res)
    
    if (length(lst_links_df) == 0) {
        
        # initialize list with a named element
        lst_links_df <- list(peaks_df2_parsed)
        names(lst_links_df) <- name_df
        
    } else {
        
        # append with a name
        lst_links_df[[name_df]] <- peaks_df2_parsed
    }
    
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
