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

FDR_links = 0.1
FDR_dars = 0.1

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
linkPeaks_results_all <- bind_rows(list_of_df) |>
    as_tibble() |>
    filter(FDR < FDR_links)

message("Processing ", length(list_of_df), " LinkPeaks files for merged-peaks dataset")
message(nrow(linkPeaks_results_all), " total links") 

head(linkPeaks_results_all)
summary(linkPeaks_results_all$score)
nrow(linkPeaks_results_all) # 10948
# save LinkPeaks table 
f_name <- here(processedDir, paste0("ALL_LinkPeaks_signif_FDR", gsub("\\.", "", as.character(FDR_links)), ".csv"))
write.csv(linkPeaks_results_all, f_name, row.names = FALSE)

message("LinkPeaks saved for FDR", gsub("\\.", "", as.character(FDR_links)), ".csv")

#===========================================================================

message("Ploting barplot of number of links by cellType .. ")

# get number of links (rows) for each cluster
cluster_counts <- linkPeaks_results_all |>
    count(cluster, sort = TRUE)

g1 <- ggplot(cluster_counts, aes(x = reorder(cluster, n), y = n)) +
    geom_col(fill = "grey") +
    geom_text(aes(label = n), hjust = -0.2, size = 3) +
    labs(
        title = "LinkPeaks by Cell Type",
        subtitle = sprintf(
            "%d unique | FDR = %s", length(unique(linkPeaks_results_all$peak)),
            FDR_links
        ),
        x = NULL, # "Cell Type",
        y = "Number of Links"
    ) +
    theme_minimal() +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))  # add 10% space on right

f_name <- paste0("links_barplot_ct_FDR", FDR_links, ".pdf")
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

message("Parsing DARs for each cellType at ", paste("FDR=", FDR_dars, " "))

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
            mutate(fdr = as.numeric(fdr)) |> 
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
DARs_results_all <- filter_signific_DARs(
    lst_DARs_files_cellType, FDR_dars, processedDir
) |>
    as_tibble() |>
    #   keep peaks that are differentially accessible in at most 2 cell types
    group_by(peak_id) |>
    filter(n_distinct(cell_type) <= 2) |>
    ungroup()

#===========================================================================

message("Plotting barplot of DARs by cellType ...")

# count DARs
DARs_counts <- DARs_results_all |>
    count(cell_type, sort = TRUE)

# Extract cluster-order used in the LinkPeaks plot to set plots side to side in the same order
cluster_order <- cluster_counts |>
    pull(cluster)

# Clusters unique to each dataset
DAR_clust     <- setdiff(DARs_results_all$cell_type, cluster_order)
Linked_clust  <- setdiff(cluster_order, DARs_results_all$cell_type)

if (length(DAR_clust) > 0) {
    message("Clusters found in DARs but missing in LinkPeaks: ",
            paste(DAR_clust, collapse = ", "))
    cluster_order <- unique(c(cluster_order, DAR_clust))
} else if (length(Linked_clust) > 0) {
    message("Clusters found in LinkPeaks but missing in DARs: ",
            paste(Linked_clust, collapse = ", "))
    cluster_order <- unique(c(cluster_order, Linked_clust))
} else {
    message("Both LinkPeaks and DARs contain the same clusters.")
}

df_sub <- DARs_counts

# This prevents issues with missing levels and ensures we only order the visible data.
current_cluster_order <- cluster_order[cluster_order %in% unique(df_sub$cell_type)]
reversed_order <- rev(current_cluster_order)

g2 = ggplot(
        df_sub, aes(x = factor(cell_type, levels = reversed_order), y = n)
    ) +
    geom_col(fill = "grey") +
    geom_text(aes(label = n), hjust = -0.2, size = 3) +
    labs(
        title = paste("DARs by Cell Type"),
        subtitle = sprintf(
            "%d unique | FDR = %s", length(unique(DARs_results_all$peak_id)),
            FDR_dars
        ),
        x = NULL,
        y = "Number of DARs"
    ) +
    theme_minimal() +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2)))  # add 10% space on right

ggsave(
    filename = file.path(plotDir, paste0("DARs_barplot_ct_", FDR_dars, ".pdf")),
    plot = g2,
    width = 7,
    height = 7
)

combined_LinkPeaks_DARs <- g1 + g2
f_name <- sprintf("LinkPeaks_DARs_combined_barplot_FDR%s_FDR%s.pdf", FDR_links, FDR_dars)
ggsave(here(plotDir, f_name),
       combined_LinkPeaks_DARs, width = 7, height = 7)

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
             FDR_CC = FDR_links,
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

DARs_FDRX <- DARs_results_all
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
            FDR_dars  = FDR_dars,
            # Add FC of chromatin accessibility:
            # - logFC > 0 → peak is more accessible (open) in the target cell type
            # - logFC < 0 → peak is less accessible (closed) in the target cell type
            logFC          = logFC
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

f_name <- sprintf("Overlaps_venn_diagram_FDR%s_%s.pdf", FDR_links, FDR_dars)
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

f_name <- here(
    processedDir, sprintf("Overlaps_LinkPeak_DARs_FDR%s_%s.csv", FDR_links, FDR_dars)
)
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
