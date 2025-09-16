########################################################################
## Differential Accessibility (DA) using Seurat/Signac from pseudobulk assays on meged peaks
## - Peaks were computed with MACS2
##
## Authors. CSC
## Date. Sep 15, 2025
## Recommended resources mem=30GB
########################################################################

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

# Check/create directories
inputRDS_Dir <- here(
    "processed-data",
    "06_peak_calling",
    "12_pseudobulk_MACS2"
)
output_Dir <- here(
    "processed-data",
    "06_peak_calling",
    "16_pseudobulk_DARs_MACS2_reduced"
)
plotDir <- here(
    "plots",
    "06_peak_calling",
    "16_pseudobulk_DARs_MACS2_reduced"
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
## Differential Accessibility Analysis in Seurat/Signac
## Input: Seurat object with pseudobulk RNA+ATAC assay and reduced peaks
########################################################################

##==============================================================================
## Load Seurat / macs peaks / filtered genes. And make verification

Seurat_base_name <- "Mid_pseudobulk.spearman.5e5_merged_peaks.rds"
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

meta <- SeuratOBJ_pb@meta.data
colnames(meta)
# [1] "orig.ident"                        "nCount_RNA"                       
# [3] "nFeature_RNA"                      "nCount_ATAC_macs2_merged_pseudo"  
# [5] "nFeature_ATAC_macs2_merged_pseudo"
head(meta)

## Ensure we are contrasting by cell-type
if (!"cell_type" %in% colnames(SeuratOBJ_pb@meta.data)) {
    SeuratOBJ_pb$cell_type <- sub("^[^_]+_", "", SeuratOBJ_pb$orig.ident)  # keep part after first underscore
}
length(SeuratOBJ_pb$cell_type)
# [1] 169

Idents(SeuratOBJ_pb) <- "cell_type"
cluster_ids <- levels(SeuratOBJ_pb)
unique(Idents(SeuratOBJ_pb))

message("Pseudobulk groups (cell-types):")
cluster_ids

## Define assay and metadata grouping
PSEUDO_ATAC_ASSAY <- "ATAC_macs2_merged_pseudo"  
grouping_var <- "orig.ident"    
# Set default assay
DefaultAssay(SeuratOBJ_pb) <- PSEUDO_ATAC_ASSAY

## double-check depth differences between pseudobulk samples
## Check total counts per pseudobulk sample
col_sums <- colSums(GetAssayData(SeuratOBJ_pb, assay = PSEUDO_ATAC_ASSAY, layer = "counts"))
col_sums
message("Library size summary:\n"); print(summary(col_sums))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 4217    97388   623469  2605095  2827754 46470237 


## Find differentially accessible peaks between cell-types
## Loop over multiple identities

# testing 
# da_res <- FindMarkers(
#     object = SeuratOBJ_pb,
#     ident.1 = "LHb.1",
#     ident.2 = NULL,              # NULL = rest of cells
#     only.pos = FALSE,            # get both up and down
#     min.pct = 0.1,               # require peak present in 10% of cells
#     test.use = "LR"              # likelihood ratio test (works well for binary peak presence)
#     # latent.vars = "nCount_ATAC"  # correct for sequencing depth if needed => this is pb data
# )


## provide as input reduced peaks, so, this set should be consistent DA across cell-types
## function give per-cluster DA peaks vs all others: cluster marker accessibility

message("Starting DA across cell-types ...")

all_da <- setNames(vector("list", length(cluster_ids)), cluster_ids)

# FindMarkers() internally supports parallelization / MAC
future::plan("multicore", workers = 4)


for (ct in cluster_ids) {
    message("[", ct, "]")
    # Skip if this cell-type has very few samples
    n_samples_ct <- sum(SeuratOBJ_pb$cell_type == ct)
    if (n_samples_ct < 3) {
        warning("Skipping ", ct, " (only ", n_samples_ct, " pseudobulk samples).")
        next
    }
    res <- tryCatch(
        {
            tmp <- FindMarkers(
                object   = SeuratOBJ_pb,
                ident.1  = ct,
                ident.2  = NULL,     # group vs. rest
                # only.pos = TRUE,
                test.use = "LR"      # used for binary accessibility
            )
            tmp$cluster <- ct
            tmp$peak    <- rownames(tmp)
            tmp$FDR     <- p.adjust(tmp$p_val, method = "BH")
            tmp
        },
        error = function(e) {
            warning("FindMarkers failed for ", ct, ": ", conditionMessage(e))
            NULL
        }
    )
    all_da[[ct]] <- res
}

## Filter() keeps only elements of a list that return TRUE. Removes NULL entries
da_results <- bind_rows(Filter(Negate(is.null), all_da))
# Save combined
write.csv(da_results, file.path(output_Dir, "DA_all_clusters.csv"), row.names = FALSE)
# milestone
# da_results <- read.csv(file.path(output_Dir, "DA_all_clusters.csv"))
table(da_results$cluster)
nrow(da_results)

# indiv. tests by cell-type
for (cl in names(all_da)) {
    if (is.null(all_da[[cl]])) next
    output_path <- here(output_Dir, paste0("DA_", cl, ".csv"))
    write.csv(all_da[[cl]], output_path, row.names = FALSE)
}

message("Ends DA across cell-types!")


################################################################################
## summarize per-cluster stats and extract top enriched peaks

message("Starting summary ...")

head(da_results)

## Keeps all statistically significant peaks, whether up (positive log2FC) or down (negative log2FC).
da_sig <- da_results |>
    filter(FDR < sig_thresh) #& abs(avg_log2FC) > lfc_thresh)
    #filter(FDR < sig_thresh & abs(avg_log2FC) > lfc_thresh)
nrow(da_sig)
# [1] 21849

tail(da_sig)
table(da_sig$cluster)

# summary counts per cluster / respect the fold-change cutoff when defining up vs down DARs
summary_table <- da_sig |>
    group_by(cluster) |>
    summarise(
        n_sig_peaks = n(),
        n_up   = sum(avg_log2FC >  lfc_thresh),   # significantly up 0.25
        n_down = sum(avg_log2FC < -lfc_thresh),   # significantly down 0.25
        top_up_peak   = peak[which.max(avg_log2FC)],        # most upregulated
        top_down_peak = peak[which.min(avg_log2FC)],        # most downregulated
        max_log2FC = max(avg_log2FC),             # give the extreme values per cluster.
        min_log2FC = min(avg_log2FC)
    ) |>
    arrange(desc(n_sig_peaks))


# save summary
f_name <- here(output_Dir, "DA_summary_per_cluster.csv")
write.csv(summary_table, f_name, row.names = FALSE)

# extract top 10 enriched peaks per cluster
top10_per_cluster <- da_sig |>
    group_by(cluster) |>
    arrange(desc(avg_log2FC)) |>
    slice_head(n = 10)

table(top10_per_cluster$cluster)

f_name <- here(output_Dir, "DA_top10_per_cluster.csv")
write.csv(top10_per_cluster, f_name, row.names = FALSE)

message("Summary and top10 tables saved in: ", output_Dir)


########################################################################

message("Starting summary Volcano plots per cluster")

f_name <- here(plotDir, "Volcano_all_clusters.pdf")
pdf(f_name, width = 7, height = 6)  


for (ct in unique(da_results$cluster)) {
    # ct = unique(da_results$cluster[1])
    df <- da_results |> filter(cluster == ct)
    
    # mark significant points
    df$signif <- with(df, ifelse(FDR < sig_thresh & abs(avg_log2FC) > lfc_thresh, "significant", "not"))
    
    # classify as Up / Down / Not
    df$signif <- case_when(
        df$FDR < sig_thresh & df$avg_log2FC >  lfc_thresh  ~ "Up",
        df$FDR < sig_thresh & df$avg_log2FC < -lfc_thresh  ~ "Down",
        TRUE                                               ~ "Not"
    )
    
    # select top 5 peaks by FDR
    top5 <- df |>
        arrange(FDR) |>
        slice_head(n = 5)
    
    p <- ggplot(df, aes(x = avg_log2FC, y = -log10(FDR))) +
        geom_point(aes(color = signif), alpha = 0.6, size = 1.2) +
        scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not" = "grey70")) +
        geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed") +
        geom_hline(yintercept = -log10(sig_thresh), linetype = "dashed") + 
        geom_text_repel(
            data = top5,
            aes(label = peak),
            size = 3,
            box.padding = 0.3,
            point.padding = 0.2,
            max.overlaps = 10
        ) +
        labs(
            title = paste0("Volcano plot - ", ct),
            x = "log2 Fold Change (1vsALL)",
            y = "-log10(FDR)"
        ) +
        theme_bw() +
        theme(legend.position = "bottom")
    
    print(p) 
    
}

dev.off()


## Barplot of Up vs Down DARs per cluster

barplot_name <- here(plotDir, "DA_barplot_UpDown_per_cluster.pdf")

# Prepare counts of Up and Down DARs
up_down_counts <- da_sig |>
    mutate(direction = case_when(
        avg_log2FC >  lfc_thresh ~ "Up",
        avg_log2FC < -lfc_thresh ~ "Down",
        TRUE ~ "Not"
    )) |>
    filter(direction != "Not") |>    # keep only true Up/Down DARs
    group_by(cluster, direction) |>
    summarise(n = n(), .groups = "drop")

# Order clusters by total number of DARs
cluster_order <- up_down_counts |>
    group_by(cluster) |>
    summarise(total = sum(n), .groups = "drop") |>
    arrange(desc(total)) |>
    pull(cluster)

up_down_counts$cluster <- factor(up_down_counts$cluster, levels = cluster_order)

# Plot
p1 <- ggplot(up_down_counts, aes(x = cluster, y = n, fill = direction)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Up" = "red", "Down" = "blue")) +
    labs(
        title = "Number of Up and Down DARs per cluster",
        x = "Cluster",
        y = "Number of DARs",
        fill = "Direction"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(barplot_name, p1, width = 8, height = 5, dpi = 300)

message("Barplot saved to: ", barplot_file)



message("Plots done!")

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
