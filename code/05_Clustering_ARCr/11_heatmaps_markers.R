########################################################################
## Plot Heatmaps of WNN clustering
## - Top X DEG (FDR<5%) on aggregated cell-types normalized at z-scores
##
## Authors. CSC
## Date. May 08, 2024
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
## Note. Seurat objects were created with module load conda_R/4.3.x
########################################################################

library("Seurat")
library("Signac")
library("Matrix")
library("ComplexHeatmap")
library("circlize") # colorRamp2
library("tidyverse")
library("here")

# Check/create directories

## clusters renamed for Spatial-Registration on Visium project
inputRDS_Dir <- here(
    "processed-data",
    "05_Clustering_ARCr",
    "17_wnn_clustering_final_ct"
)
plotDir <- here(
    "plots",
    "05_Clustering_ARCr",
    "11_heatmaps_markers"
)
inputCVS_Dir <- here(
    "processed-data",
    "05_Clustering_ARCr",
    "02_Hb_celltypes_from_seurat_reanalyze_v3",
    "cvs_files_markers"
)
input_meanRatio_Dir <- here(
    "processed-data",
    "05_Clustering_ARCr",
    "13_wnn_geneExp_plt_mean_ratio_annotated"
)


## Check directories
if (!dir.exists(plotDir)) {
    dir.create(plotDir)
}

## Load Seurat with WNN

message("Loading Seurat ....")

# clusters renamed for sharing with Visium project(s)
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_HD.rds"
seurat_name <- here(inputRDS_Dir, Seurat_base_name)

# Load Seurat
SeuratOBJ <- readRDS(here(inputRDS_Dir, Seurat_base_name))
colnames(SeuratOBJ@meta.data)
DefaultAssay(SeuratOBJ) <- "RNA"
#levels(SeuratOBJ)
# [1] "C.04.LHb.4"      "C.05.LHb.2.7"    "C.06.LHb.4"      "C.07.MHb.2"     
# [5] "C.08.LHb.4"      "C.09.LHb.4"      "C.10.MHb.1"      "C.11.MHb.1.2"   
# [9] "C.13.LHb.4"      "C.14.MHb.1"      "C.16.MHb.1.2"    "C.18.LHb.1.3.4"  ...

# Aggregate average expression per cluster
avg_expr <- AggregateExpression(
    SeuratOBJ, 
    group.by = "cluster_ann", 
    return.seurat = FALSE
)$RNA
head(avg_expr)

base_name <- str_extract(seurat_name, regex("C\\.\\w+"))


## =============================================================================
## Read DEG file and prepare top genes with current cluster annotation

message("Reading and preparing DEG (Wilcox test - Seurat) ....")

DEG_file_name <- "WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_cellTypes_integrated_top50.csv"
DEG_file_name <- here(inputCVS_Dir, DEG_file_name)
df_markers_findALLSeurat <- read.csv(DEG_file_name)
df_markers_findALLSeurat <- df_markers_findALLSeurat |> drop_na(cell_type)
head(df_markers_findALLSeurat)
# p_val avg_log2FC pct.1 pct.2 p_val_adj cluster     gene            cell_type
# 1     0   4.336764 0.938 0.100         0       1 OTX2-AS1        DD_Inhib.Thal
# 2     0   3.849508 0.900 0.083         0       1      KIT        DD_Inhib.Thal
# 3     0   3.682782 0.927 0.122         0       1    MEIS2      LB_Thalamus/MDm

table(unique(SeuratOBJ@meta.data$C.leiden_wnn)==unique(SeuratOBJ@meta.data$seurat_clusters))
# TRUE 
# 40 
head(unique(SeuratOBJ@meta.data$cluster_ann))

# Replace numeric IDs clusters with annotated cluster names and remove clusters removed (NAs)
head(df_markers_findALLSeurat)
table(SeuratOBJ$seurat_clusters)
table(SeuratOBJ$cluster_ann)

# Build mapping from Seurat - Get numeric annotated name mapping
cluster_map <- SeuratOBJ@meta.data |>
    select(seurat_clusters, cluster_ann) |>
    distinct() |>
    mutate(seurat_clusters = as.integer(as.character(seurat_clusters)))  # ensures consistent ordering
#cluster_map
#                       seurat_clusters     cluster_ann
# S04_AAACATGCAGTAATAG-1               1 C.01.Inhib.Thal
# S04_AACCTTGCATTATGAC-1               2      C.02.Oligo
# S04_ATCTTTGGTGATGGCT-1               3 C.03.Excit.Thal
# S04_AAACAGCCAGCAAGGC-1               4      C.04.LHb.4

# clean DEG table (clusters removed) and join
df_markers_findALLSeurat_clean <- df_markers_findALLSeurat |>
    mutate(cluster = as.integer(as.character(cluster))) |>
    left_join(cluster_map, by = c("cluster" = "seurat_clusters")) |>
    rename(cluster_name = cluster_ann) |>
    filter(!is.na(cluster_name))

unique(df_markers_findALLSeurat_clean$cluster_name)
head(df_markers_findALLSeurat_clean[c("cluster", "cluster_name")])
# cluster    cluster_name
# 1       1 C.01.Inhib.Thal
# 2       1 C.01.Inhib.Thal
# 3       1 C.01.Inhib.Thal ...

# check cleaned marker table is compatible with your expression matrix
stopifnot(all(df_markers_findALLSeurat_clean$cluster_name %in% colnames(avg_expr)))


## =============================================================================
# function to aggregate and scaled mxt to plot 'Diagonal heatmap' by cluster and broad cell-type

message("Loading functions ....")

topGenes_mtx <- function(dge_annotated_clusters, avg_expr, top_genes = 3) {

    # Get top N genes per cluster (highest avg_log2FC)    
    top_markers <- dge_annotated_clusters |>
        group_by(cluster_name) |>
        top_n(n = top_genes, wt = avg_log2FC)

    # before sub-setting I force diagonal layout
    # order by cluster → avg_log2FC → one gene per row
    ordered_gene_cluster <- top_markers |>
        arrange(factor(cluster_name, levels = colnames(avg_expr)), desc(avg_log2FC)) |>
        distinct(gene, cluster_name)
    
    # Keep only genes that exist in avg_expr mtx
    ordered_gene_cluster <- ordered_gene_cluster |>
        filter(gene %in% rownames(avg_expr))
    #table(ordered_gene_cluster$cluster_name %in% colnames(avg_expr))

    # Subset matrix and define column order mtx: rows = genes, columns = clusters
    # unique cluster order (to avoid repeating columns in heatmap)
    unique_clusters <- unique(ordered_gene_cluster$cluster_name)
    mat_ordered <- avg_expr[ordered_gene_cluster$gene, unique_clusters]

    #  scale (row-wise z-score)
    mat_scaled <- t(scale(t(as.matrix(mat_ordered))))
    
    return(list(
        mat_scaled = mat_scaled,
        ordered_genes = ordered_gene_cluster
    ))
    
}

make_heatmap_row_group <- function(mat_scaled, row_cluster_group, top_anno, row_anno, 
                                   row_names_font_size, f_name, grid_ann, height_htm) {
    
    hm1 <- Heatmap(
        mat_scaled,
        name = "Z-score",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        top_annotation = top_anno,
        left_annotation = row_anno,
        row_split = row_cluster_group,
        border = TRUE,  # adds a horizontal line between row groups
        row_title_gp = gpar(fontsize = 10, fontface = "bold"),  # customize strip label
        row_title_rot = 0,
        gap = unit(1, "mm"),  # spacing between row blocks
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = row_names_font_size)
    )
    pdf(f_name, width = 12, height = height_htm)
    draw(hm1)
    grid::grid.text(
        grid_ann,
        x = unit(.9, "npc"),    # right-aligned
        y = unit(0.02, "npc"),  # distance from bottom
        gp = gpar(fontsize = 9, fontface = "italic")
    )
    dev.off()
    
    return(hm1)    
}

make_heatmap_column_group <- function(mat_scaled, row_cluster_group, column_cluster_group, top_anno, row_anno, 
                                      row_names_font_size, f_name, grid_ann, height_htm) {
    hm1 <- Heatmap(
        mat_scaled,
        name = "Z-score",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        top_annotation = top_anno,
        left_annotation = row_anno,
        row_split = row_cluster_group,
        column_split = column_cluster_group,  # adds visual column blocks
        border = TRUE,  # adds a horizontal line between row groups
        row_title_gp = gpar(fontsize = 9, fontface = "bold"),  # customize strip label
        column_title_gp = gpar(fontsize = 9, fontface = "bold"),  # customize strip label
        row_names_gp = gpar(fontsize = row_names_font_size),
        column_names_gp = gpar(fontsize = 8),
        row_title_rot = 0,
        column_title_rot = 45,
        gap = unit(1, "mm")  # spacing between row blocks
    )
    pdf(f_name, width = 12, height = height_htm)
    draw(hm1)
    grid::grid.text(
        grid_ann,
        x = unit(.9, "npc"),    # right-aligned
        y = unit(0.02, "npc"),  # distance from bottom
        gp = gpar(fontsize = 9, fontface = "italic")
    )
    dev.off()
    
    return(hm1)  
    
}

## =============================================================================
# process data and make heatmap of aggregated and scaled mtx by WNN clusterID

message("Start processing ...")

set.seed(7312025)

top_genes_number = c(3, 5)

## define colors for the group (columns)
group_colors <- c(
    LHb = "#1f78b4",
    MHb = "#ad1d8c",
    Oligo = "#384a08",
    Astrocyte = "#532222", 
    OPC = "#829454",
    Microglia = "#141b02",
    Endo = "#d95f02",
    Thal = "#4d55b7",
    Excit.Thal = "#6855A3",
    Inhib.Thal = "#805EE6"
    #Other = "black"
)


for (top in top_genes_number) {
    # top = 3
    message("Processing top ", top, " genes by cluster ...")
    
    topGenes_mtx_lst <- topGenes_mtx(df_markers_findALLSeurat_clean, avg_expr, top) 
    names(topGenes_mtx_lst)
    mat_scaled <- topGenes_mtx_lst[["mat_scaled"]]
    
    ## verification
    #str(mat_scaled)
    dimnames(mat_scaled)
    length(dimnames(mat_scaled)[[1]]) # ge. For top3 = 114 genes
    length(dimnames(mat_scaled)[[2]]) # ge. For top3 = 39 clusters
    #head(as.data.frame(as.matrix(mat_scaled[ , "C.37.Thal", drop = FALSE])))
    
    ## =============================================================================
    # Define unique cluster names and assign colors based on broad cell types in cluster names
    
    clusters <- colnames(mat_scaled)
    
    # define group membership at broad level
    merged_cluster <- sapply(clusters, function(cl) {
        case_when(
            grepl("LHb", cl) ~ "LHb",
            grepl("MHb", cl) ~ "MHb",
            grepl("Excit\\.Thal", cl) ~ "Excit.Thal",
            grepl("Inhib\\.Thal", cl) ~ "Inhib.Thal",
            grepl("Thal", cl) ~ "Thal",
            grepl("Astro", cl) ~ "Astrocyte",
            grepl("Oligo", cl) ~ "Oligo",
            grepl("OPC", cl) ~ "OPC",
            grepl("Microglia", cl) ~ "Microglia",
            grepl("Endo", cl) ~ "Endo",
            TRUE ~ "Other"
        )
    })
    
    # make it a named factor 
    merged_cluster <- factor(merged_cluster, 
                             levels = names(group_colors))
    unique(merged_cluster)                         
    names(merged_cluster) <- clusters
    
    top_anno <- HeatmapAnnotation(
        Region = merged_cluster,  # name shown in legend
        col = list(Region = group_colors),
        annotation_name_side = "left"
    )
    
    ## =============================================================================
    # Grouped row strips by cluster
    ordered_gene_cluster <- topGenes_mtx_lst[["ordered_genes"]]
    row_cluster <- ordered_gene_cluster$cluster_name
    names(row_cluster) <- ordered_gene_cluster$gene
    
    # Group each row (gene) by high-level cluster group (LHb, MHb, Thal, Other)
    row_cluster_group <- sapply(row_cluster, function(cl) {
        case_when(
            grepl("LHb", cl) ~ "LHb",
            grepl("MHb", cl) ~ "MHb",
            grepl("Excit\\.Thal", cl) ~ "Excit.Thal",
            grepl("Inhib\\.Thal", cl) ~ "Inhib.Thal",
            grepl("Thal", cl) ~ "Thal",
            grepl("Astro", cl) ~ "Astrocyte",
            grepl("Oligo", cl) ~ "Oligo",
            grepl("OPC", cl) ~ "OPC",
            grepl("Microglia", cl) ~ "Microglia",
            grepl("Endo", cl) ~ "Endo",
            TRUE ~ "Other"
        )
    })
    
    row_cluster_group <- factor(row_cluster_group, 
                                levels = names(group_colors))
    unique(row_cluster_group)                            
    names(row_cluster_group) <- names(row_cluster)  # Ensure names = genes
    
    row_anno <- rowAnnotation(
        Region = row_cluster_group,
        col = list(Region = group_colors),
        show_annotation_name = FALSE,
        show_legend = FALSE
    )
    
    ## =============================================================================
    
    message("Preparing row cluster-group heatmap ...")
    
    if (top == 3) { row_names_font_size <- 7.5; 6.8 }
    if (top == 3) { height_htm <- 12; 22 }
    
    # make ands save heatmap
    f_name <- paste0("heatmap_top", top,"_genes-row_grouped-Broad_res.pdf")
    grid_ann <- paste0("Top ", top," DEG (FDR < 5%)\nWilcoxon test (Seurat)")
    
    hm1 <- make_heatmap_row_group(mat_scaled, row_cluster_group, top_anno, row_anno, 
                                  row_names_font_size, f_name, grid_ann, height_htm) # file name and bottom ann
    
    message("Heatmap with row cluster-group done and saved for top ", top, " genes")
    
    ## =============================================================================
    # Alternative version: blocks column-grouped by brain region (e.g., LHb, MHb, Thal, etc.)
    
    message("Processing top ", top, " genes by cluster ...")
    
    # Assign each column (cluster) to a region
    column_cluster_group <- sapply(colnames(mat_scaled), function(cl) {
        case_when(
            grepl("LHb", cl) ~ "LHb",
            grepl("MHb", cl) ~ "MHb",
            grepl("Excit\\.Thal", cl) ~ "Excit.Thal",
            grepl("Inhib\\.Thal", cl) ~ "Inhib.Thal",
            grepl("Thal", cl) ~ "Thal",
            grepl("Astro", cl) ~ "Astrocyte",
            grepl("Oligo", cl) ~ "Oligo",
            grepl("OPC", cl) ~ "OPC",
            grepl("Microglia", cl) ~ "Microglia",
            grepl("Endo", cl) ~ "Endo"
        )
    })
    
    message("Preparing column cluster-group heatmap ...")
    
    column_cluster_group <- factor(column_cluster_group, 
                                   levels = names(group_colors))
    
    names(column_cluster_group) <- colnames(mat_scaled)
    
    f_name <- paste0("heatmap_top", top,"_genes-column_grouped-Broad_res.pdf")
    
    hm1 <- make_heatmap_column_group(mat_scaled, row_cluster_group, column_cluster_group, top_anno, row_anno, 
                                     row_names_font_size, f_name, grid_ann, height_htm) # file name and bottom ann
    
    message("Plot done and saved for top ", top, " genes")
    

}



## =============================================================================
## Sub-setting data for Medial and Lateral Hb 

# define genes to plot by cluster
top_subset <- 10

## define colors for the group
group_colors_subset <- c(
    LHb = "#1f78b4",
    MHb = "#ad1d8c"
)

# build mtx
topGenes_mtx_lst <- topGenes_mtx(df_markers_findALLSeurat_clean, avg_expr, top_subset) 
names(topGenes_mtx_lst)
mat_scaled <- topGenes_mtx_lst[["mat_scaled"]]

# Subset clusters from scaled matrix
sel_clusters <- grep("MHb|LHb", colnames(mat_scaled), value = TRUE)
mat_LHb_MHb <- mat_scaled[, sel_clusters]
genes_LHb_MHb <- rownames(mat_LHb_MHb)

# Subset ordered gene-cluster mapping
ordered_gene_cluster <- topGenes_mtx_lst[["ordered_genes"]] |>
    filter(as.character(cluster_name) %in% sel_clusters)
# Filter to relevant clusters and genes
ordered_gene_cluster_sub <- topGenes_mtx_lst[["ordered_genes"]] |>
    filter(grepl("LHb|MHb", cluster_name), gene %in% rownames(mat_LHb_MHb)) |>
    distinct(gene, .keep_all = TRUE)  # remove duplicated genes

# cluster order from ordered_gene_cluster_sub
ordered_clusters <- unique(ordered_gene_cluster_sub$cluster_name)
ordered_clusters <- ordered_clusters[ordered_clusters %in% colnames(mat_LHb_MHb)]

# Reorder columns
ordered_clusters <- as.character(ordered_clusters)
mat_LHb_MHb <- mat_LHb_MHb[, ordered_clusters]

# Build gene-to-cluster mapping
row_cluster <- ordered_gene_cluster_sub$cluster_name
names(row_cluster) <- ordered_gene_cluster_sub$gene

# Group assignment (LHb/MHb)
row_cluster_group <- sapply(row_cluster, function(cl) {
    case_when(
        grepl("LHb", cl) ~ "LHb",
        grepl("MHb", cl) ~ "MHb",
        TRUE ~ "Other"  # just in case
    )
})

row_cluster_group <- factor(row_cluster_group, levels = c("LHb", "MHb"))
names(row_cluster_group) <- names(row_cluster)

# Ensure consistent gene list between matrix and cluster group
common_genes <- intersect(rownames(mat_LHb_MHb), names(row_cluster_group))

# Subset both to keep common genes only
mat_LHb_MHb <- mat_LHb_MHb[common_genes, , drop = FALSE]
row_split_sub <- row_cluster_group[common_genes]
# Re-check
stopifnot(all(rownames(mat_LHb_MHb) == names(row_split_sub)))

# Order genes by their cluster (as in the original topGenes_mtx())
gene_order <- ordered_gene_cluster_sub |>
    arrange(match(cluster_name, ordered_clusters)) |>
    pull(gene)

gene_order <- gene_order[gene_order %in% rownames(mat_LHb_MHb)]  # safety
# Reorder rows
mat_LHb_MHb <- mat_LHb_MHb[gene_order, , drop = FALSE]
row_split_sub <- row_split_sub[gene_order]

# Build row annotation
row_anno_sub <- rowAnnotation(
    Region = row_split_sub,
    col = list(Region = group_colors_subset),
    show_annotation_name = FALSE,
    show_legend = FALSE
)

# Broad group assignment for column clusters
merged_Hb_clusters <- sapply(colnames(mat_LHb_MHb), function(cl) {
    case_when(
        grepl("LHb", cl) ~ "LHb",
        grepl("MHb", cl) ~ "MHb",
        TRUE ~ "Other"
    )
}) |> factor(levels = c("LHb", "MHb"))

# Set names to match columns of matrix
names(merged_Hb_clusters) <- colnames(mat_LHb_MHb)
# build the annotation
top_anno_sub <- HeatmapAnnotation(
    Region = merged_Hb_clusters,
    col = list(Region = group_colors_subset),
    annotation_name_side = "left"
)
# Ensure column names match top_annotation
stopifnot(all(colnames(mat_LHb_MHb) %in% names(merged_Hb_clusters)))

f_name <- paste0("heatmap_subset_MHb_LHb_top", top_subset, "_genes-column_grouped-Broad_res.pdf")
grid_ann <- paste0("Top ", top_subset," DEG (FDR < 5%)\nWilcoxon test (Seurat)")

hm_LHb_MHb <- make_heatmap_column_group(mat_LHb_MHb, row_split_sub, merged_Hb_clusters, top_anno_sub, row_anno_sub, 
                                        row_names_font_size, f_name, grid_ann, 18) # file name and bottom ann

message("All done!!!")


# library("slurmjobs")
# job_single(
#     "11_heatmaps_markers",
#     cores = 2,
#     partition = "katun",
#     memory = "80G",
#     create_shell = TRUE
#     )


## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
