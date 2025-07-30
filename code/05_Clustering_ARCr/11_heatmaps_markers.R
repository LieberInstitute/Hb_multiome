########################################################################
## Plot Heatmaps of GEX on WNN clustering
## - Top3 and 5 DEG (FDR<5%) on aggregated cell-types normalized at z-scores
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
library("circlize")
library("bluster")
library("tidyverse")
library("stringr")
library("here")

## input directories

here()

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

## Check directories
if (!dir.exists(plotDir)) {
    dir.create(plotDir)
}

## Load Seurat with WNN

# clusters renamed for sharing with Visium project(s)
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_HD.rds"
# old: "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium.rds"
seurat_name <- here(inputRDS_Dir, Seurat_base_name)

# Load Seurat
SeuratOBJ <- readRDS(here(inputRDS_Dir, Seurat_base_name))
colnames(SeuratOBJ@meta.data)
DefaultAssay(SeuratOBJ) <- "RNA"
cluster_levels <- levels(SeuratOBJ)
# [1] "C.04.LHb.4"      "C.05.LHb.2.7"    "C.06.LHb.4"      "C.07.MHb.2"     
# [5] "C.08.LHb.4"      "C.09.LHb.4"      "C.10.MHb.1"      "C.11.MHb.1.2"   
# [9] "C.13.LHb.4"      "C.14.MHb.1"      "C.16.MHb.1.2"    "C.18.LHb.1.3.4" 
# [13] "C.23.LHb.1"      "C.24.LHb.4"      "C.30.LHb.7"      "C.31.LHb.4"     
# [17] "C.33.LHb.1.3"    "C.36.MHb.3"      "C.40.LHb.4"      "C.01.Inhib.Thal"
# [21] "C.02.Oligo"      "C.03.Excit.Thal" "C.12.Excit.Thal" "C.15.Excit.Thal"
# [25] "C.17.Excit.Thal" "C.19.Inhib.Thal" "C.20.Astrocyte"  "C.21.Astrocyte" 
# [29] "C.22.Oligo"      "C.25.Excit.Thal" "C.26.OPC"        "C.27.Microglia" 
# [33] "C.28.Inhib.Thal" "C.29.Endo"       "C.32.Excit.Thal" "C.35.Excit.Thal"
# [37] "C.37.Thal"       "C.38.Inhib.Thal" "C.39.Inhib.Thal" "C.41.Microglia" 

# Aggregate average expression per cluster
avg_expr <- AggregateExpression(
    SeuratOBJ, 
    group.by = "cluster_ann", 
    return.seurat = FALSE
)$RNA
head(avg_expr)

base_name <- str_extract(seurat_name, regex("C\\.\\w+"))

message("Processing Heatmap for ", base_name)


## =============================================================================
## Read DEG file and prepare top genes with current cluster annotation

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
cluster_map
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
# Prepare aggregated and scaled mxt to plot 'Diagonal heatmap' by cluster and broad cell-type

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
    table(ordered_gene_cluster$cluster_name %in% colnames(avg_expr))
    # TRUE 
    # 114 
    
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

## =============================================================================

# get and plot heatmap of aggregated and scaled mtx by cluster and broad cell-type

top_genes_number = c(3, 5)

## define colors for the group
group_colors <- c(
    LHb = "#1f78b4",
    MHb = "#ad1d8c",
    Oligo = "#384a08",
    Astrocyte = "#532222", 
    OPC = "#829454",
    Microglia = "#141b02",
    Endo = "#d95f02",
    Thal = "#4d55b7"
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
    length(dimnames(mat_scaled)[[1]]) # 188 genes
    length(dimnames(mat_scaled)[[2]]) # 40 clusters
    #head(as.data.frame(as.matrix(mat_scaled[ , "C.37.Thal", drop = FALSE])))
    
    ## =============================================================================
    # Define unique cluster names and assign colors based on broad cell types in cluster names
    
    clusters <- colnames(mat_scaled)
    
    # define group membership at broad level
    merged_cluster <- sapply(clusters, function(cl) {
        case_when(
            grepl("LHb", cl) ~ "LHb",
            grepl("MHb", cl) ~ "MHb",
            grepl("Thal", cl) ~ "Thal",
            grepl("Astro", cl) ~ "Astrocyte",
            grepl("Oligo", cl) ~ "Oligo",
            grepl("OPC", cl) ~ "OPC",
            grepl("Microglia", cl) ~ "Microglia",
            grepl("Endo", cl) ~ "Endo"
        )
    })
    
    # make it a named factor 
    merged_cluster <- factor(merged_cluster, 
                             levels = c("LHb", "MHb", "Thal", "Astrocyte", "Oligo","OPC", "Microglia", "Endo"))
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
            grepl("Thal", cl) ~ "Thal",
            grepl("Astro", cl) ~ "Astrocyte",
            grepl("Oligo", cl) ~ "Oligo",
            grepl("OPC", cl) ~ "OPC",
            grepl("Microglia", cl) ~ "Microglia",
            grepl("Endo", cl) ~ "Endo"
        )
    })
    
    row_cluster_group <- factor(row_cluster_group, 
                                levels = c("LHb", "MHb", "Thal", "Astrocyte", "Oligo","OPC", "Microglia", "Endo"))
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
    if (top == 3) { height_htm <- 12; 20 }
    
    # Add bottom annotation with text to describe basic stats to process the heatmap
    bottom_anno <- HeatmapAnnotation(
        annotation_label = anno_text(
            "Top 3 DEG (FDR < 5%) - Wilcoxon test (Seurat)",
            gp = gpar(fontsize = 9, fontface = "italic"),
            just = "center"
        ),
        annotation_name_side = "bottom",
        annotation_height = unit(1.2, "cm"),
        show_annotation_name = FALSE
    )

    make_heatmap <- function(mat_scaled, row_cluster_group, top_anno, row_anno) {
        
        hm <- Heatmap(
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
            row_names_gp = gpar(fontsize = row_names_font_size),
            bottom_annotation = bottom_anno
        )
    
        return(hm)    
    }
    
    hm1 <- make_heatmap(mat_scaled, row_cluster_group, top_anno, row_anno)
    hm1
    f_name <- paste0("heatmap_top", top,"_genes-row_grouped-Broad_res.pdf")
    
    pdf(here(plotDir, f_name), width = 12, height = height_htm)
    draw(hm1)
    dev.off()
    
    message("Heatmap with row cluster-group done and saved for top ", top, " genes")
    
    
    ## =============================================================================
    # Applied to the same mtx from topGenes_mtx_lst(), this is an alternative version
    #    to visualize multiple diagonal blocks grouped by brain region (e.g., LHb, MHb, Thal, etc.)
    
    message("Processing top ", top, " genes by cluster ...")
    
    # Assign each column (cluster) to a region
    column_cluster_group <- sapply(colnames(mat_scaled), function(cl) {
        case_when(
            grepl("LHb", cl) ~ "LHb",
            grepl("MHb", cl) ~ "MHb",
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
                                   levels = c("LHb", "MHb", "Thal", "Astrocyte", "Oligo","OPC", "Microglia", "Endo"))
    
    names(column_cluster_group) <- colnames(mat_scaled)
    
    make_heatmap <- function(mat_scaled, column_cluster_group, top_anno, row_anno) {
        Heatmap(
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
            gap = unit(1, "mm"),  # spacing between row blocks
            bottom_annotation = bottom_anno
        )
        
    }
    
    hm1 <- make_heatmap(mat_scaled, column_cluster_group, top_anno, row_anno)
    
    f_name <- paste0("heatmap_top", top,"_genes-column_grouped-Broad_res.pdf")
    
    pdf(here(plotDir, f_name), width = 12, height = height_htm)
    draw(hm1)
    dev.off()
    
    message("Plot done and saved for top ", top, " genes")
    

}



## =============================================================================
## Subsettig data for Medial and Lateral Hb 

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

# # row split 
# row_split_sub <- row_cluster_group[rownames(mat_LHb_MHb)]
# # check
# stopifnot(all(rownames(mat_LHb_MHb) == names(row_split_sub)))
# # Error: all(rownames(mat_LHb_MHb) == names(row_split_sub)) is not TRUE
# Ensure consistent gene list between matrix and cluster group
common_genes <- intersect(rownames(mat_LHb_MHb), names(row_cluster_group))
# Subset both to keep common genes only
mat_LHb_MHb <- mat_LHb_MHb[common_genes, , drop = FALSE]
row_split_sub <- row_cluster_group[common_genes]
# Re-check the sanity
stopifnot(all(rownames(mat_LHb_MHb) == names(row_split_sub)))

# Build row annotation
row_anno_sub <- rowAnnotation(
    Region = row_split_sub,
    col = list(Region = group_colors_subset),
    show_annotation_name = FALSE,
    show_legend = FALSE
)


# Ensure column names match top_annotation
stopifnot(all(colnames(mat_LHb_MHb) %in% names(merged_Hb_clusters)))





# # ## get thalamus clusters
# # cells_to_keep <- WhichCells(SeuratOBJ, idents = Thal)
# # SeuratOBJ_subset <- subset(SeuratOBJ, cells = cells_to_keep)
# # levels(SeuratOBJ_subset)
# # unique(Idents(SeuratOBJ_subset))
# 
# top_markers <- df_markers_findALLSeurat |>
#     group_by(cluster) |>
#     top_n(n = 5, wt = avg_log2FC)
# head(top_markers)
# 
# ## Filter genes for Thal clusters
# colnames(df_markers_findALLSeurat)
# thal_markers <- df_markers_findALLSeurat |>
#     filter(grepl("Thal", cell_type))
# 
# ## Subset the matrix 
# # Keep only Thal columns
# #colnames(mat)
# thal_cols <- grep("Thal", colnames(mat), value = TRUE)
# mat_thal <- mat[, thal_cols]
# 
# # Keep only Thal genes
# mat_thal <- mat_thal[intersect(thal_markers$gene, rownames(mat_thal)), ]
# # Scale
# mat_thal_scaled <- t(scale(t(as.matrix(mat_thal))))
# 
# # create merged_cluster annotation (Thal group only)
# merged_cluster_thal <- sapply(thal_cols, function(cl) {
#     if (grepl("LHb", cl)) {
#         "LHb"
#     } else if (grepl("MHb", cl)) {
#         "MHb"
#     } else if (grepl("Thal", cl)) {
#         "Thal"
#     } else {
#         "Other"
#     }
# })
# 
# top_anno_thal <- HeatmapAnnotation(
#     Region = factor(merged_cluster_thal, levels = names(group_colors)),
#     col = list(Region = group_colors),
#     annotation_name_side = "left",
#     show_legend = FALSE  # this disables only the annotation legend
# )
# 
# hm_thal <- make_heatmap(mat_thal_scaled, top_anno_thal)
# 
# pdf("thal_top3_marker_heatmap.pdf", width = 10, height = 8)
# draw(hm_thal)
# dev.off()
# 
# ## =============================================================================
# 
# 
# 
# 
# 
# 
# 
# ########## Gene markers lists. New function to join LB and DD gene markers lists
# 
# ## Source gene markers lists 
# 
# source(here("code", "04_DiffExpr_Clustering_seurat", "remote_DGE_marker_gene_lists.R"))   
# 
# markers.custom <- get_multiple_markers_genes_lst()
# tmp <- names(markers.custom)
# tmp <- paste(tmp, collapse=', ')
# message("Processing ", length(markers.custom), " categories of gene-markers list \n *****(", tmp, ")*****")
# 
# 
# ## merge DD markers for heatmap with 'Data-driven' markers and merge LB markers for heatmap with 'Literature-Based'markers
# 
# # prepare DD marker's list
# DD_markers_lst <- append(markers.custom$DD_MHb, markers.custom$DD_LHb) 
# # prepare LB marker's list
# LB_markers_lst <- markers.custom[grep("^LB", names(markers.custom))]
# LB_markers_lst <- as.vector(unlist(LB_markers_lst, recursive = FALSE))
# 
# ## scale expression once for fair comparison across all clusters
# 
# # Keep only genes found in the dataset
# dd_markers <- DD_markers_lst[DD_markers_lst %in% rownames(SeuratOBJ)]
# lb_markers <- LB_markers_lst[LB_markers_lst %in% rownames(SeuratOBJ)]
# 
# if (length(dd_markers)==0 || length(lb_markers)==0) {
#     stop()
# } else {
#     combined_named_list <- list(
#         Hb_DD_markers = dd_markers,
#         Hb_LB_markers = lb_markers
#     )
#         
# }
# names(combined_named_list)
# # [1] "Hb_DD_markers" "Hb_LB_markers"
# 
# 
# ## prepare Seurat object
# 
# # Set cell type as the identity class for grouping
# colnames(SeuratOBJ@meta.data)
# head(SeuratOBJ@meta.data$cluster_ann)
# unique(SeuratOBJ@meta.data$cluster_ann)
# #Idents(SeuratOBJ) <- "seurat_clusters"
# 
# # Multi-cell-type heatmap, scale all relevant genes once, across all cells
# genes_to_scale <- as.vector(unlist(append(combined_named_list[1], combined_named_list[2])))
# SeuratOBJ <- ScaleData(SeuratOBJ, features = genes_to_scale)
# 
# 
# 
# ######## Plot 1: plot the 2 gene-markers list (DD, LB) across all clusters in separate heatmaps
# 
# heatmap_list <- list()
# 
# for (ct in seq_along(combined_named_list)) {
# 
#     # ct=2
#     # grab cell type names and label
#     ct_name <- names(combined_named_list[ct])
#     print(paste0("Cell-type: ", ct_name))
#     markers_to_plt <- unlist(combined_named_list[ct])
#     # markers_to_plt=c("PRKCB"      "ADCY6"      "DCHS2"      "GLIS1"      "GREB1L")
#     title_label <- paste0("Cell type class: ", ct_name)
#     message("Preparing heatmap with class: ", ct_name)    
#     
#     # ####### this heatmap plot all Hb medial and lateral markers across all clusters - cluster size aware
#     # heatmap_plot <- DoHeatmap(SeuratOBJ,
#     #                           group.by = "cluster_ann",
#     #                           features = markers_to_plt, size = 2,
#     #                           disp.min = -2.5, disp.max = 2.5,
#     #                           group.bar = TRUE, # Omits color bar by identity class (from group.by)
#     #                           slot = "scale.data") +
#     #     scale_fill_gradientn(colors = c("blue", "white", "red")) +
#     #     ggtitle(title_label) +
#     #     theme(
#     #         plot.title = element_text(hjust = 0.5),
#     #         axis.text.y = element_text(size = 6),
#     #         legend.position = "none"
#     #     )
#     # # disp.min / disp.max = 2.5: gene-expr after scaling can have extreme values. Clipping keeps the heatmap visually interpretable
#     # print(heatmap_plot)
#     
#     
#     ####### this heatmap plot all Hb medial and lateral markers across all clusters - same cluster size
#     
#     # Downsample to equal cell numbers per group
#     # Number of cells per group you want (e.g., 50)
#     n_cells <- 50
#     
#     # Randomly sample equal number of cells from each group to control column size
#     # Add cell IDs as a column first (from rownames)
#     meta_df <- SeuratOBJ@meta.data
#     meta_df$cell_id <- rownames(meta_df)
#     # Sample cells evenly across clusters
#     cells_to_plot <- meta_df |>
#         group_by(cluster_ann) |>
#         sample_n(size = min(n_cells, n()), replace = FALSE) |>
#         arrange(cluster_ann) |>   # This sets a fixed order to remove dendogram manually
#         pull(cell_id)
#     length(cells_to_plot)
#     # Reorder markers if needed
#     # ordered_markers <- markers_to_plt[markers_to_plt %in% rownames(SeuratOBJ)]
#     
#     # Plot with fixed number of cells per group (uniform column width)
#     
#     temp_plot <- DoHeatmap(SeuratOBJ,
#                               features = markers_to_plt,
#                               group.by = "cluster_ann",
#                               cells = cells_to_plot,
#                               group.bar = TRUE,
#                               label = TRUE,
#                               size = 2,
#                               disp.min = -2.5, disp.max = 2.5, 
#                               slot = "scale.data") +
#         scale_fill_gradientn(colors = c("blue", "white", "red")) +
#         #ggtitle(title_label) +
#         theme(
#             plot.title = element_blank(),               # remove internal title
#             #plot.title = element_text(hjust = 0.5, size = 8),
#             plot.margin = margin(t = 20, r = 5, b = 30, l = 5),  # extra bottom space
#             axis.text.y = element_text(size = 5),
#             axis.text.x = element_blank(),        # hide x-axis labels (just in case)
#             axis.ticks.x = element_blank(),       # remove x-axis ticks
#             legend.position = "none"
#         )
#     print(temp_plot)
#     # Optionally add horizontal lines for this specific marker group
#     
#     if (ct_name == "Hb_LB_markers") {
#         # Note: y-axis is reversed — top gene = lowest y value
#         # make a vector with size of gene blocks to add them to heatmpas - seperate cell types by marker class
#         gene_blks <- markers.custom[grep("^LB", names(markers.custom))]
#         gene_blks_lng <- map_int(gene_blks, ~ length(.x) )    
#         gene_block_sizes <- as.vector(gene_blks_lng)
#         cumulative_positions <- cumsum(gene_block_sizes)
#         
#         # Add lines at block boundaries (skip last)
#         for (y_pos in cumulative_positions[-length(cumulative_positions)]) {
#             temp_plot <- temp_plot +
#                 geom_hline(yintercept = y_pos + 0.5, color = "black", linetype = "solid", linewidth = 0.3)
#         }
#     }
#     
#     # Store the final plot in the list
#     heatmap_list[[ct]] <- temp_plot
#     
# }
# 
# # Combine all heatmaps side by side
# 
# length(heatmap_list)
# 
# combined_plot <- wrap_plots(heatmap_list, ncol = length(heatmap_list)) +
#     plot_annotation(
#         title = "Habenula Data-Driven and Literature-Based gene markers side to side",
#         theme = theme(
#             plot.title = element_text(size = 10, hjust = 0.5, face = "bold")
#         )
#     )
# 
# f_name <- paste0(
#     base_name,
#     "_heatmap_all_reference_markers_width5.pdf"
# )
# pdf(file = here(plotDir, f_name), width = 4 * length(heatmap_list), height = 6)
# print(combined_plot)
# dev.off()
# 
# length(heatmap_list)
# #heatmap_list[1]
# f_name <- paste0(
#     base_name,
#     "_heatmap_all_reference_markers_width10.pdf"
# )
# combined_plot <- wrap_plots(heatmap_list, ncol = length(heatmap_list))
# pdf(file = here(plotDir, f_name), width = 10 * length(heatmap_list), height = 6)
# print(combined_plot)
# dev.off()
# 
# 
# 
# ######## Plot 2: plot the top-x genes (markers) used to annotate each cluster
# 
# # df_markers_findALLSeurat$cluster
# # df_markers_findALLSeurat$gene
# # df_markers_findALLSeurat$cell_type
# 
# # Downsample to equal cell numbers per group
# # Number of cells per group you want (e.g., 50)
# n_cells <- 50
# 
# # Randomly sample equal number of cells from each group to control column size
# # Add cell IDs as a column first (from rownames)
# meta_df <- SeuratOBJ@meta.data
# meta_df$cell_id <- rownames(meta_df)
# # Sample cells evenly across clusters
# cells_to_plot <- meta_df |>
#     group_by(cluster_ann) |>
#     sample_n(size = min(n_cells, n()), replace = FALSE) |>
#     arrange(cluster_ann) |>   # This sets a fixed order to remove dendrogram manually
#     pull(cell_id)
# 
# # extract top 10 genes per cluster
# top_markers <- df_markers_findALLSeurat %>%
#     group_by(cluster) %>%
#     top_n(n = 10, wt = avg_log2FC)
# 
# # Unique gene list
# marker_genes <- unique(top_markers$gene)
# marker_genes <- marker_genes[marker_genes %in% rownames(SeuratOBJ)]
# #SeuratOBJ <- ScaleData(SeuratOBJ, features = marker_genes, verbose = FALSE)
# 
# #Idents(SeuratOBJ) <- "cluster_ann"
# 
# plt <- DoHeatmap(SeuratOBJ,
#           features = marker_genes,
#           group.by = "cluster_ann",
#           cells = cells_to_plot,
#           group.bar = TRUE,
#           size = 3) +
#     scale_fill_gradientn(colors = c("blue", "white", "red")) +
#     #ggtitle("Top Marker Genes per Cluster") +
#     theme(#plot.title = element_text(hjust = 0.5, size = 8),
#           plot.margin = margin(t = 20, r = 5, b = 30, l = 5),  # extra bottom space
#           axis.text.y = element_text(size = 5),
#           axis.text.x = element_blank(),        # hide x-axis labels (just in case)
#           axis.ticks.x = element_blank(),       # remove x-axis ticks
#           legend.position = "none")
# 
# f_name <- paste0(
#     base_name,
#     "_heatmap_top10genes.pdf"
# )
# pdf(file = here(plotDir, f_name), width = 5 * length(heatmap_list), height = 6)
# print(plt)
# dev.off()
# 
# message("All plots done!")


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
