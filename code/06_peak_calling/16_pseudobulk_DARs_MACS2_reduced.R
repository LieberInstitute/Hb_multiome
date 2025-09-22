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
lfc_thresh <- 0.1    # Default


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

## Ensure we are testing by cell-type
if (!"cellType" %in% colnames(SeuratOBJ_pb@meta.data)) {
    SeuratOBJ_pb$cellType <-  sub("_.*$", "", colnames(SeuratOBJ_pb)) # keep part after first underscore
}


Idents(SeuratOBJ_pb) <- "cellType"
cluster_ids <-levels(Idents(SeuratOBJ_pb))
#unique(Idents(SeuratOBJ_pb))
message("Pseudobulk groups (cell-types):")
cluster_ids

## Define assay and metadata grouping
PSEUDO_ATAC_ASSAY <- "ATAC_macs2_merged_pseudo"  
DefaultAssay(SeuratOBJ_pb) <- PSEUDO_ATAC_ASSAY

## add a library-size covariate for LR
covar_col <- paste0("nCount_", PSEUDO_ATAC_ASSAY)
stopifnot(covar_col %in% colnames(SeuratOBJ_pb@meta.data))

## Check total counts per pseudobulk sample
col_sums <- colSums(GetAssayData(SeuratOBJ_pb, assay = PSEUDO_ATAC_ASSAY, layer = "counts"))
col_sums
message("Library size summary:\n"); print(summary(col_sums))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 4217    97388   623469  2605095  2827754 46470237 


## Find differential peaks (markers) between cell-types

## provide as input peaks merged, so, this set should be consistent DA across cell-types
## function give per-cluster DA peaks vs all others: cluster marker accessibility

message("Starting DA across cell-types ...")

## Confirm peaks layer is non-empty
table(Idents(SeuratOBJ_pb)) |> print()
m <- GetAssayData(SeuratOBJ_pb, assay = PSEUDO_ATAC_ASSAY, layer = "counts") 
cat("nz features:", sum(Matrix::rowSums(m) > 0), " / ", nrow(m), "\n")

all_da <- setNames(vector("list", length(cluster_ids)), cluster_ids)

# FindMarkers() internally supports parallelization / for local use: MAC
# future::plan("multicore", workers = 16)


for (ct in cluster_ids) {
    # ct = "Endo"
    
    message("Processing [", ct, "]")
    n_samples_ct <- sum(SeuratOBJ_pb$cellType == ct)
    if (n_samples_ct < 3) {
        warning("Skipping ", ct, " (only ", n_samples_ct, " pseudobulk samples).")
        next
    }
    res <- tryCatch(
        {
            tmp <- FindMarkers(
                object   = SeuratOBJ_pb,
                ident.1  = ct,
                ident.2  = NULL,            # group vs. rest
                only.pos = FALSE,
                test.use = "LR",            # used for binary accessibility
                layer = "counts", 
                latent.vars = covar_col,    # control for library size
                min.pct = 0,                # turn off detection filtering while debugging
                logfc.threshold = 0         # no FC pre-filter; filter later by FDR
            )
            tmp$cluster <- ct
            tmp$peak    <- rownames(tmp)
            #tmp$FDR     <- p.adjust(tmp$p_val, method = "BH")
            dim(tmp); head(tmp)
        },
        error = function(e) {
            warning("FindMarkers failed for ", ct, ": ", conditionMessage(e))
            NULL
        }
    )
    all_da[[ct]] <- res
}

message("Ends DA!")


## Filter() keeps only elements of a list that return TRUE. Removes NULL entries
da_results <- bind_rows(Filter(Negate(is.null), all_da))

if (nrow(da_results) == 0) {
    message("No DA results produced. Check identities and FindMarkers settings.")
    q(save="no", status = 0)  # or stop("No DA results")
}

# Save combined
write.csv(da_results, file.path(output_Dir, "DA_all_clusters.csv"), row.names = FALSE)

# Save combinedper cellType
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
    filter(p_val_adj < sig_thresh & abs(avg_log2FC) > lfc_thresh)
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

# extract top upregulated per cluster
top10_up <- da_sig |>
    group_by(cluster) |>
    arrange(desc(avg_log2FC)) |>
    slice_head(n = 10)  |>
    mutate(direction = "Up")
# top 10 downregulated per cluster
top10_down <- da_sig |>
    group_by(cluster) |>
    arrange(avg_log2FC) |>
    slice_head(n = 10) |>
    mutate(direction = "Down")

top10_per_cluster <- bind_rows(top10_up, top10_down) |>
    arrange(cluster, desc(direction), desc(abs(avg_log2FC)))

f_name <- here(output_Dir, "DA_top10_up_down_per_cluster.csv")
write.csv(top10_per_cluster, f_name, row.names = FALSE)

message("Saved top 10 Up and Down DARs per cluster: ", output_Dir)


########################################################################

message("Starting Volcano plots per cluster ...")

f_name <- here(plotDir, "Volcano_all_clusters.pdf")
pdf(f_name, width = 7, height = 6)  


for (ct in unique(da_results$cluster)) {
    # ct = unique(da_results$cluster[1])
    
    df <- da_results |> filter(cluster == ct)
    if (nrow(df) == 0) next
    # classify as Up / Down / Not
    df$signif <- case_when(
        df$p_val_adj < sig_thresh & df$avg_log2FC >  lfc_thresh  ~ "Up",
        df$p_val_adj < sig_thresh & df$avg_log2FC < -lfc_thresh  ~ "Down",
        TRUE ~ "Not"
    )
    
    # select top 5 peaks by FDR
    top5 <- df |>
        arrange(p_val_adj) |>
        slice_head(n = 5)
    
    p <- ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
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

message("Volcano saved to: ", f_name)


## Barplot of Up vs Down DARs per cluster

message("Starting Barplot per cluster ...")

f_name <- here(plotDir, "DA_barplot_UpDown_per_cluster.pdf")

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

ggsave(f_name, p1, width = 8, height = 5, dpi = 300)

message("Barplot saved to: ", f_name)



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
