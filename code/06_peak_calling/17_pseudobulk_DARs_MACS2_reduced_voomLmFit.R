########################################################################
## Differential Accessibility (DA) using voomLmFit / pseudobulk multiome assays with merged peaks
##
## Authors. CSC
## Date. Sep 16, 2025
## Recommended resources mem=30GB
########################################################################

suppressPackageStartupMessages({
    library("edgeR")     # for DGEList, filterByExpr, voomLmFit
    library("limma")     # for contrasts.fit, eBayes, topTable, plot helpers
    library("Matrix")    # if counts are sparse; otherwise base matrix is fine
    library("data.table")
})
library("ggplot2")
library("ggrepel") 
library("future")
library("future.apply")
library("Seurat")
library("Signac")
library("dplyr")
library("here")


# Testing spearman at 5e4 on macs2 peaks 
resolution_level = "Mid"
p_met = "spearman"
w_size = "5e5"
sig_thresh <- 0.1    # FDR cutoff
lfc_thresh <- 0.25   # logFC cutoff


if (length(resolution_level)) {
    message("Processing job for peak-method:\n",
            p_met,
            "\nWindow-size\n",
            w_size)
} else {
    message("Input arguments missed")
    stop()
}

# Inputs:
# counts:   matrix of raw counts, rows = peaks, cols = pseudobulk samples
#           (e.g., aggregated by cluster x donor)
# coldata:  data.frame with sample metadata, nrow = ncol(counts)
#           must include at least: sample_id, cluster_id, group, donor
# peaks_df: data.frame/GRanges of peaks for annotation (optional)

# Check/create directories
inputRDS_Dir <- here(
    "processed-data",
    "06_peak_calling",
    "12_pseudobulk_MACS2"
)
Seurat_base_name <- "Mid_pseudobulk.spearman.5e5_merged_peaks.rds"

output_Dir <- here(
    "processed-data",
    "06_peak_calling",
    "17_pseudobulk_DARs_MACS2_reduced_voomLmFit"
)
plotDir <- here(
    "plots",
    "06_peak_calling",
    "17_pseudobulk_DARs_MACS2_reduced_voomLmFit"
)


## Check directories
if (!dir.exists(plotDir)) {
    dir.create(plotDir)
}
if (!dir.exists(output_Dir)) {
    dir.create(output_Dir)
}

# List all files matching the specific clustering resolution level
lst_peak_files <- list.files(
    path = inputRDS_Dir,
    pattern = "*.rds",
)
#lst_peak_files = list.files(path = input_cvsDir)
message("Link peak-genes files found:")
lst_peak_files
# [1] "mtx_merged_peaks_cell_level_Mid_resolution.rds"   
# [2] "Seurat_peaks_merged_cell_level_Mid_resolution.rds"


########################################################################
## Differential Accessibility Analysis using voomLmFit
## Input: Seurat object with pseudobulk RNA+ATAC assay with merged peaks
########################################################################

##==============================================================================
## Load Seurat / macs peaks / filtered genes. And make verification

seurat_name <- here(inputRDS_Dir, Seurat_base_name)
SeuratOBJ_pb <- readRDS(here(inputRDS_Dir, Seurat_base_name))

# check-ins
SeuratOBJ_pb
# An object of class Seurat 
# 366933 features across 169 samples within 2 assays 
# Active assay: ATAC_macs2_merged_pseudo (351037 features, 333643 variable features)
# 2 layers present: counts, data
# 1 other assay present: RNA
# 1 dimensional reduction calculated: lsi

# ## get rna counts
# DefaultAssay(SeuratOBJ_pb) <- "RNA"
# rna_counts <- GetAssayData(SeuratOBJ_pb, assay="RNA", layer="data")
# length(rownames(rna_counts)) # [1] 36601

## =============================================================================

## (1) get matrix of raw counts: rows = peaks, cols = pseudobulk samples
PSEUDO_ATAC_ASSAY <- "ATAC_macs2_merged_pseudo"
atac_counts <- GetAssayData(SeuratOBJ_pb, assay = PSEUDO_ATAC_ASSAY, laye ="counts")
length(rownames(atac_counts)) # [1] 351037
head(atac_counts)

## (2) get data.frame with sample metadata, nrow = ncol(counts)
#  must include at least: sample_id, cluster_id, group, donor
meta <- SeuratOBJ_pb@meta.data
colnames(meta)
# [1] "orig.ident"                        "nCount_RNA"                       
# [3] "nFeature_RNA"                      "nCount_ATAC_macs2_merged_pseudo"  
# [5] "nFeature_ATAC_macs2_merged_pseudo"
head(meta)

## add meta-data in the expected format
meta$sample_id <- row.names(meta)
meta$donor <- sapply(strsplit(row.names(meta), "_", fixed = TRUE), 
                "[[", 2)
tail(meta)

stopifnot(all(colnames(atac_counts) == meta$sample_id))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Tweakables
min_cpm        <- 1            # expression filter threshold (logical CPM rule)
min_samples    <- 3            # require in >= this many samples overall or per-group
use_block      <- TRUE         # TRUE = account for repeated measures (e.g., donor)
use_samp_wts   <- TRUE         # TRUE = estimate sample quality weights
trend_ebayes   <- TRUE         # TRUE = eBayes(trend=TRUE) often good for counts
robust_ebayes  <- TRUE         # TRUE = robust empirical Bayes
fdr_cutoff     <- 0.10

## ============================================================================/



## Ensure we are contrasting by cell-type
if (!"cell_type" %in% colnames(SeuratOBJ_pb@meta.data)) {
    SeuratOBJ_pb$cell_type <- sub("^[^_]+_", "", SeuratOBJ_pb$orig.ident)  # keep part after first underscore
}
length(SeuratOBJ_pb$cell_type)
# [1] 169

Idents(SeuratOBJ_pb) <- "cell_type"
cluster_ids <- levels(SeuratOBJ_pb)

message("Pseudobulk groups (cell-types):")
cluster_ids

## Define ATAC assay and metadata grouping
PSEUDO_ATAC_ASSAY <- "ATAC_macs2_merged_pseudo"  
grouping_var <- "orig.ident"    
# Set default assay
DefaultAssay(SeuratOBJ_pb) <- PSEUDO_ATAC_ASSAY

## double-check depth differences between pseudobulk samples
## Check total counts per pseudobulk sample
col_sums <- colSums(GetAssayData(SeuratOBJ_pb, assay = PSEUDO_ATAC_ASSAY, layer = "counts"))
col_sums
# Astrocyte_S03-Hb-r  Astrocyte_S04-Hb-r  Astrocyte_S05-Hb-r  Astrocyte_S06-Hb-r 
# 2254001              420381             1300877              507894 
# Astrocyte_S07-Hb-r  Astrocyte_S08-Hb-r  Astrocyte_S09-Hb-r  Astrocyte_S10-Hb-r 
# 2123583             1832793              830177             1654334 
message("Library size summary:\n"); print(summary(col_sums))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 4217    97388   623469  2605095  2827754 46470237 

# ## provide as input reduced peaks, so, this set should be consistent DA across cell-types
# ## function give per-cluster DA peaks vs all others: cluster marker accessibility
# 
# message("Starting DA across cell-types ...")
# 
# all_da <- setNames(vector("list", length(cluster_ids)), cluster_ids)
# 
# # FindMarkers() internally supports parallelization / MAC
# future::plan("multicore", workers = 4)
# 
# 
# for (ct in cluster_ids) {
#     message("[", ct, "]")
#     # Skip if this cell-type has very few samples
#     n_samples_ct <- sum(SeuratOBJ_pb$cell_type == ct)
#     if (n_samples_ct < 3) {
#         warning("Skipping ", ct, " (only ", n_samples_ct, " pseudobulk samples).")
#         next
#     }
#     res <- tryCatch(
#         {
#             tmp <- FindMarkers(
#                 object   = SeuratOBJ_pb,
#                 ident.1  = ct,
#                 ident.2  = NULL,     # group vs. rest
#                 only.pos = TRUE,
#                 test.use = "LR"      # used for binary accessibility
#                 #min.pct = 0.05, logfc.threshold = 0.25
#             )
#             tmp$cluster <- ct
#             tmp$peak    <- rownames(tmp)
#             tmp$FDR     <- p.adjust(tmp$p_val, method = "BH")
#             tmp
#         },
#         error = function(e) {
#             warning("FindMarkers failed for ", ct, ": ", conditionMessage(e))
#             NULL
#         }
#     )
#     all_da[[ct]] <- res
# }
# 
# ## Filter() keeps only elements of a list that return TRUE. Removes NULL entries
# da_results <- bind_rows(Filter(Negate(is.null), all_da))
# # Save combined
# write.csv(da_results, file.path(output_Dir, "DA_all_clusters.csv"), row.names = FALSE)
# # milestone
# # da_results <- read.csv(file.path(output_Dir, "DA_all_clusters.csv"))
# table(da_results$cluster)
# nrow(da_results)
# 
# # indiv. tests by cell-type
# for (cl in names(all_da)) {
#     if (is.null(all_da[[cl]])) next
#     output_path <- here(output_Dir, paste0("DA_", cl, ".csv"))
#     write.csv(all_da[[cl]], output_path, row.names = FALSE)
# }
# 
# message("Ends DA across cell-types!")
# 
# 
# ################################################################################
# ## summarize per-cluster stats and extract top enriched peaks
# 
# message("Starting summary ...")
# 
# da_sig <- da_results |>
#     filter(FDR < sig_thresh & abs(avg_log2FC) > lfc_thresh)
# nrow(da_sig)
# # Tue Sep 16 09:08:48 2025 ------------------------------
# tail(da_sig)
# table(da_sig$cluster)
# 
# # summary counts per cluster
# summary_table <- da_sig |>
#     group_by(cluster) |>
#     summarise(
#         n_sig_peaks = n(),
#         n_up   = sum(avg_log2FC > 0),
#         n_down = sum(avg_log2FC < 0),
#         top_peak = peak[which.max(abs(avg_log2FC))],
#         max_log2FC = max(abs(avg_log2FC))
#     ) |>
#     arrange(desc(n_sig_peaks))
# 
# # save summary
# f_name <- here(output_Dir, "DA_summary_per_cluster.csv")
# write.csv(summary_table, f_name, row.names = FALSE)
# 
# # extract top 10 enriched peaks per cluster
# top10_per_cluster <- da_sig |>
#     group_by(cluster) |>
#     arrange(desc(avg_log2FC)) |>
#     slice_head(n = 10)
# 
# table(top10_per_cluster$cluster)
# 
# f_name <- here(output_Dir, "DA_top10_per_cluster.csv")
# write.csv(top10_per_cluster, f_name, row.names = FALSE)
# 
# message("Summary and top10 tables saved in: ", output_Dir)
# 
# 
# ########################################################################
# 
# message("Starting summary Volcano plots per cluster")
# 
# f_name <- here(plotDir, "Volcano_all_clusters.pdf")
# pdf(f_name, width = 7, height = 6)  
# 
# for (ct in unique(da_results$cluster)) {
#     # ct = unique(da_results$cluster[1])
#     df <- da_results |> filter(cluster == ct)
#     
#     # mark significant points
#     df$signif <- with(df, ifelse(FDR < sig_thresh & abs(avg_log2FC) > lfc_thresh, "significant", "not"))
#     colnames(df)
#     
#     # select top 5 peaks by FDR
#     top5 <- df |>
#         arrange(FDR) |>
#         slice_head(n = 5)
#     
#     p <- ggplot(df, aes(x = avg_log2FC, y = -log10(FDR))) +
#         geom_point(aes(color = signif), alpha = 0.6, size = 1.2) +
#         scale_color_manual(values = c("significant" = "red", "not" = "grey70")) +
#         geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", color = "black") +
#         geom_hline(yintercept = -log10(sig_thresh), linetype = "dashed", color = "black") +
#         geom_text_repel(
#             data = top5,
#             aes(label = peak),
#             size = 3,
#             box.padding = 0.3,
#             point.padding = 0.2,
#             max.overlaps = 10
#         ) +
#         labs(
#             title = paste0("Volcano plot - ", ct),
#             x = "log2 Fold Change (1vsALL)",
#             y = "-log10(FDR)"
#         ) +
#         theme_bw() +
#         theme(legend.position = "bottom")
#     
#     print(p) 
#     
# }
# 
# dev.off()
# 
# message("Plots done!")

# library("slurmjobs")
# job_single(
#   "16_pseudobulk_DARs_MACS2_reduced",
#   create_shell = TRUE,
#   partition = "katun",
#   memory = "30G",
#   cores = 2,
#   logdir = "logs",
#   command = "Rscript 16_pseudobulk_DARs_MACS2_reduced.R",
#   create_logdir = FALSE
# )


## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
