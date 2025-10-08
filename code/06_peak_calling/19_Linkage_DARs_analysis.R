########################################################################
## Find findOverlaps() between Unique LinkPeaks and DARs
## Plot:
## - barplot of number of links by cellType
## - barplot of DARs by cellType (FDR=0.1 and FDR=0.2)
## 
## Authors. CSC
## Date. Sep 24, 2025
## Recommended resources on interactive mode: srun --pty --mem=30GB --x11 bash
########################################################################

library("GenomicRanges") 
library("dplyr")
library("purrr")
library("ggplot2")
library("ggvenn")
library("patchwork")
library("tidyverse")
library("tidyr")
library("stringr")
library("here")

#===============================================================================
# resolution_level = "Mid"      # 18 cell-types
#===============================================================================

## setup variable names
resolution_level = "Mid" 
# test 2 thresholds
FDR = 0.2 # actual value to run the DAR-Links
FDR_thresh = c(0.1, 0.2)
# lfc_thresh = 0.2 # we do not use log FC for this exploratory analysis

## Set directory names
inputCSV_Links_Dir <- here(
    "processed-data",
    "06_peak_calling",
    "14_exploratory_pb_peak_scores_MACS2",
    "links_ct_merged"
)
inputCSV_DARs_Dir <- here(
    "processed-data",
    "06_peak_calling",
    "17_pseudobulk_DARs_MACS2_reduced_voomLmFit"
)
processedDir <- here(
    "processed-data",
    "06_peak_calling",
    "19_Linkage_DARs_analysis"
)
plotDir <- here(
    "plots",
    "06_peak_calling",
    "19_Linkage_DARs_analysis"
)

if (!dir.exists(processedDir)) { dir.create(processedDir) }
if (!dir.exists(plotDir)) { dir.create(plotDir) }

##==============================================================================

message("Loading  LinkPeaks results ...")

# List all files matching the specific clustering resolution level
lst_peak_files <- list.files(
    path = inputCSV_Links_Dir,
    pattern = paste0("^", resolution_level, ".*\\.csv$"),
)
#lst_peak_files = list.files(path = input_cvsDir)
message(length(lst_peak_files), " Link peak-genes files found ... ")

# bind all files 
list_of_df <- lapply(here(inputCSV_Links_Dir, lst_peak_files), read.csv)
linkPeaks_results_all <- bind_rows(list_of_df)

message("Processing ", length(list_of_df), " LinkPeaks files for merged-peaks dataset")
message(nrow(linkPeaks_results_all), " total links") 

head(linkPeaks_results_all)


#===========================================================================

message("Ploting barplot of number of links by cellType .. ")

# get number of links (rows) for each cluster
cluster_counts <- linkPeaks_results_all |>
    count(cluster, sort = TRUE)

g1 <- ggplot(cluster_counts, aes(x = reorder(cluster, n), y = n)) +
    geom_col(fill = "steelblue") +
    geom_text(aes(label = n), hjust = -0.2, size = 3) +
    labs(
        title = "Links by Cell Type",
        subtitle = paste0("FDR thr = ", FDR),
        x = "Cell Type",
        y = "Number of Links"
    ) +
    theme_minimal() +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))  # add 10% space on right

f_name <- paste0("links_barplot_ct_FDR", FDR, ".pdf")
ggsave(here(plotDir, f_name),
       g1, width = 7, height = 7)



##==============================================================================

message("Loading DARs results ...")

lst_DARs_files_cellType <- list.files(
    path = inputCSV_DARs_Dir,
    pattern = "voomlmFit_DAR_peaks_ALL_in_.*\\.csv$",
    full.names = FALSE
)

message(length(lst_DARs_files_cellType), " DAR files found ... ")
lst_DARs_files_cellType


#===========================================================================

message("Parsing DARs for each cellType at FDR ", paste("FDR=", FDR_thresh, " "))

#DARS_signif_df_lst = list()
# FDR_thresh = c(0.1, 0.2)

## extract/filters DARs from a given list of DARs at specific FDR thresh
filter_signific_DARs <- function(
        lst_DARs_files_cellType, # DARs csv list
        fdr_thresh,
        output_Dir
) {
    
    # testing
    # ct =  lst_DARs_files_cellType[2]
    # fdr_thresh = 0.2
    
    DARS_signif_df_lst <- list()
    
    for (ct in lst_DARs_files_cellType) {
        
        ct_name <- gsub("voomlmFit_DAR_peaks_ALL_in_|.csv", "", ct)
        message("Processing DARs for [", ct_name, "] with FDR=", fdr_thresh)
        
        DARs_csv <- here(inputCSV_DARs_Dir, ct)
        DARs_df <- read.csv(DARs_csv)
        
        message(nrow(DARs_df), " ", ct_name, " total DARs found ... ")
        
        DARs_significant <- DARs_df |>
            select(
                peak_id,
                logFC = paste0("logFC_", ct_name),
                fdr   = paste0("fdr_", ct_name)
            ) |>
            mutate(
                fdr = as.numeric(fdr),  # keep actual fdr values
                is_significant = if_else(fdr < fdr_thresh, "Significant", "Not Significant")
            ) |> 
            filter(fdr < fdr_thresh)
        
        DARS_signif_df_lst[[ct_name]] <- DARs_significant
        
        message(nrow(DARs_significant), " ", ct_name, " total significant DARs found ... ")
    
    }
    
    DARs_all <- bind_rows(DARS_signif_df_lst, .id = "cell_type")
    
    # save significant DARs
    f_name <- here(output_Dir, paste0("ALL_DARs_signif_FDR", gsub("\\.", "", as.character(fdr_thresh)), ".csv"))
    write.csv(DARs_all, f_name, row.names = FALSE)
    
    message(nrow(DARs_all), " total significant DARs found ... ")
    message("ALL DARs significant saved ... ")  
    
    return(DARs_all)

}

## for each FDR_thresh to test
DARs_results_all <- purrr::map_dfr(
    FDR_thresh, 
    ~ filter_signific_DARs(lst_DARs_files_cellType, .x, processedDir),
    .id = "FDR_threshold"
)

## Convert the .id column (1, 2) into meaningful labels
DARs_results_all <- DARs_results_all |>
    mutate(FDR_threshold = factor(FDR_threshold, 
                                  labels = paste0("FDR", FDR_thresh)))

table(DARs_results_all$FDR_threshold)
# FDR0.1 FDR0.2 
# 257010 325580 
DARs_results_all |> head()
nrow(DARs_results_all[DARs_results_all$cell_type == "LHb.4", ])

#===========================================================================

message("Plotting barplot of DARs by cellType ...")

# Make FDR_threshold a factor for facet labels
DARs_results_all <- DARs_results_all |>
    mutate(FDR_threshold = as.factor(FDR_threshold))

# count DARs
DARs_counts <- DARs_results_all |>
    count(cell_type, FDR_threshold, sort = TRUE)

# Get unique thresholds
thresholds <- unique(DARs_counts$FDR_threshold)

# Loop through thresholds and plot one at a time
g2 <- lapply(thresholds, function(th) {
    
        df_sub <- filter(DARs_counts, FDR_threshold == th)
        
        ggplot(df_sub, aes(x = reorder(cell_type, n), y = n)) +
            geom_col(fill = "steelblue") +
            geom_text(aes(label = n), hjust = -0.2, size = 3) +
            labs(
                title = paste("DARs by Cell Type"),
                subtitle = paste("FDR thr =", str_split(th, "FDR")[[1]][2]),
                x = "Cell Type",
                y = "Number of DARs"
            ) +
            theme_minimal() +
            coord_flip() +
            scale_y_continuous(expand = expansion(mult = c(0, 0.2)))  # add 10% space on right
    })

walk2(
    .x = g2,           
    .y = thresholds,   # thresholds for filenames
    ~ ggsave(
        filename = file.path(plotDir, paste0("DARs_barplot_ct_", .y, ".pdf")),
        plot = .x,
        width = 7,
        height = 7
    )
)


## overlay both thr for comparison purposes 

# Turn into a nicely formatted label
subtitle_label <- paste(
    names(table(DARs_results_all$FDR_threshold)),
    table(DARs_results_all$FDR_threshold),
    sep = ": ",
    collapse = " | "
)

g3 <- ggplot(DARs_counts, aes(x = reorder(cell_type, n), y = n, fill = FDR_threshold)) +
    geom_col(alpha = 0.5, position = "identity") +  # <- overlay instead of dodge
    labs(
        title = paste("DARs by Cell Type"),
        subtitle = subtitle_label,
        x = "Cell Type",
        y = "Number of DARs",
        fill = "FDR Threshold"
    ) +
    theme_minimal() +
    coord_flip()  +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2)))  # add 10% space on right

f_name <- paste0("DARs_barplot_overlay_FDRs.pdf")
ggsave(here(plotDir, f_name),
       g3, width = 7, height = 7)



#===========================================================================

message("Preparing GRanges for overlaping peaks ... ")

# inspect data
nrow(linkPeaks_results_all) # [1] 10948
colnames(linkPeaks_results_all)
head(linkPeaks_results_all$peak) # join id key
nrow(DARs_results_all) # [1] 582590
colnames(DARs_results_all)
head(DARs_results_all$peak_id) # join id key

# Split the peak column into coordinates. Note:
# seqnames / start / end: full chromatin accessibility window carried over in the aggregated assay.
# peak: a canonicalized peak ID string used for joins, overlaps, and downstream correlation reporting

link_split <- tidyr::separate(linkPeaks_results_all, peak,
                              into = c("seqnames","start","end"), sep = "-") |>
    mutate(start = as.integer(start),
           end   = as.integer(end))

# Convert LinkPeaks to GRanges to run findOverlaps()
gr_links <- 
    with(link_split, 
        GRanges(
             seqnames = seqnames,
             ranges   = IRanges(start, end),
             strand   = strand,
             # Metadata columns passed directly:
             peak_id_links = paste(seqnames, start, end, sep = "-"),
             CCscore  = score, # spearman CC
             gene_name  = gene,
             gene_id  = gene_id, 
             FDR_CC = FDR,
             cluster  = cluster,
             tss = tss,
             # the strand of the linked gene (+ or -)
             gene_strand = gene_strand, 
             distance = distance,
             distance_kb = distance_kb,
             # distance from peak center to TSS, 
             # keeping genomic sign (positive = downstream, negative = upstream)
             signed_distance = signed_distance, 
             # adjusts the sign of signed_distance to reflect the gene’s strand:
             # If the gene is on +, then upstream is negative, downstream positive
             # If the gene is on -, it flips accordingly
             signed_by_strand = signed_by_strand
        )
    )
length(gr_links) # [1] 10948
head(gr_links)


##==============================================================================
## Filter DARs with FDR=0.2

DARs_FDRX <- DARs_results_all |>
    filter(FDR_threshold == paste0("FDR", FDR)) 
nrow(DARs_FDRX) < nrow(DARs_results_all) # [1] 325580
nrow(DARs_FDRX) # [1] 325580
head(DARs_FDRX)
# FDR_threshold cell_type            peak_id      logFC          fdr
# 1        FDR0.2 Astrocyte chr1-629811-630032  1.1115262 1.970054e-07
# 2        FDR0.2 Astrocyte chr1-633694-634122  1.1348873 8.347289e-09


# split the peak_id "chr-start-end" into seqnames, start, end

DARs_split <- tidyr::separate(DARs_FDRX, peak_id, into = c("seqnames", "start", "end"), sep = "-") |>
    mutate(start = as.integer(start),
           end   = as.integer(end))

nrow(DARs_split) # [1] 325580


# In GRanges(), anything other than seqnames, ranges, and strand is read as metadata
# Differential Accessibility Regions (DARs) are treated as unstranded 
gr_dars <- 
    with(DARs_split, 
        GRanges(
            seqnames = seqnames,
            ranges   = IRanges(start = start, end = end),
            strand   = "*",
            # Metadata columns passed directly:
            peak_id        = paste(seqnames, start, end, sep = "-"),
            cell_type      = cell_type,
            FDR_threshold  = FDR_threshold, # 0.1 / 0.2 
            # Add FC of chromatin accessibility:
            # - logFC > 0 → peak is more accessible (open) in the target cell type
            # - logFC < 0 → peak is less accessible (closed) in the target cell type
            logFC          = logFC,
            fdr_dars = fdr
            #is_significant = is_significant
    )
)

## inspect
length(gr_dars) # [1] 325580
head(gr_dars)



#===========================================================================

message("Finding overlaps ...")

## Unique overlaps venn diagram
venn_list <- list(
    DARs = unique(gr_dars$peak_id),
    LinkPeaks = unique(gr_links$peak_id_links)
)
length(venn_list$LinkPeaks) # [1] 7098
length(venn_list$DARs) # [1] 211558
g_venn_overlaps <- ggvenn(venn_list,
       fill_color = c("skyblue", "orange"),
       stroke_size = 0.5,
       set_name_size = 4)

f_name <- paste0("Overlaps_venn_diagram_FDR", FDR, ".pdf")
ggsave(here(plotDir, f_name),
       g_venn_overlaps, width = 5, height = 5)


# findOverlaps() reports all pairs of ranges that overlap # “≥1 bp overlap”
# If one LinkPeaks region overlaps many DARs regions, you’ll get multiple rows for the same LinkPeak
# hits <- findOverlaps(gr_links, gr_dars)

# exact coordinate matches 
exact_hits <- findOverlaps(gr_links, gr_dars, type="equal") 
length(exact_hits) # 20650
## mini-test
# unique_equal <- unique(mcols(gr_links)$peak_id[queryHits(exact_hits)])
# length(unique_equal)

# Extract matched ranges
exact_links <- gr_links[queryHits(exact_hits)]
exact_dars  <- gr_dars[subjectHits(exact_hits)]
message("exact_links: ", length(exact_links))
message("exact_dars: ", length(exact_dars))

# # Compute percent overlap for both query (links) and subject (DARs)
# ov <- findOverlaps(gr_links, gr_dars)
# pi <- pintersect(gr_links_ct[queryHits(ov)], gr_dars_ct[subjectHits(ov)])
# prop_query  <- width(pi) / width(gr_links_ct[queryHits(ov)])
# prop_subject <- width(pi) / width(gr_dars_ct[subjectHits(ov)])
# # Keep only reciprocal overlaps ≥50% on both sides
# hits <- ov[prop_query >= 0.5 & prop_subject >= 0.5]

# Combine metadata from both sides
overlaps_df <- data.frame(
    as.data.frame(mcols(gr_links)[queryHits(exact_hits), ]),
    as.data.frame(mcols(gr_dars)[subjectHits(exact_hits), ])
)

nrow(overlaps_df) # [1] 20650
head(overlaps_df)

f_name <- here(processedDir, paste0("Overlaps_LinkPeak_DARs_FDR", FDR, ".csv"))
write.csv(overlaps_df, f_name, row.names = FALSE)


message("All done!!!")


# library("slurmjobs")
# job_single(
#   "19_Linkage_DARs_analysis",
#   create_shell = TRUE,
#   partition = "katun",
#   memory = "30G",
#   cores = 2,
#   logdir = "logs",
#   command = "Rscript 19_Linkage_DARs_analysis.R",
#   create_logdir = FALSE
# )

## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

