########################################################################
## Plot Heatmaps of GEX on WNN clustering
##
## Authors. CSC
## Date. May 08, 2024
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
## Note. Seurat objects were created with module load conda_R/4.3.x

## Note. this pipeline requires module load conda_R/4.3.x to keep the integrity of the seurat object

########################################################################

library("Seurat")
library("Signac")
library("ggplot2")
library("bluster")
# library("viridisLite")
library("patchwork")
# library("ggplotify")
# library("gridExtra")
library("purrr")
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
    "08_wnn_gene_expression_plts_renamed_idents"
)
plotDir <- here(
    "plots",
    "05_Clustering_ARCr",
    "08_wnn_gene_expression_plts_renamed_idents"
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

## Load input with RDS wnn to compare

# WNN clustering results of interest. To plot annotated or not annotatted clusters
# For inputRDS_Dir_not_annotated
# inputRDS_Dir_not_annotated <- here("processed-data", "05_Clustering_ARCr", "01_clustering_std_method")
# Seurat_base_name <- "seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.louvain_lsi_r1"
# For inputRDS_Dir_annotated
# inputRDS_Dir_annotated <- here("processed-data", "05_Clustering_ARCr", "05_rename_idents")
# Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2.rds"

# For inputRDS_Dir, clusters renamed for Spatial-Registration on Visium project
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium.rds"
seurat_name <- here(inputRDS_Dir, Seurat_base_name)

# Load Seurat
SeuratOBJ <- readRDS(here(inputRDS_Dir, Seurat_base_name))
levels(SeuratOBJ)
## Levels should be
# [1] "C.05 DD_LHb" "C.07 DD_MHb" "C.10 DD_MHb" "C.11 DD_MHb" "C.14 DD_MHb"
# [6] "C.16 DD_MHb" "C.18 DD_LHb" "C.23 DD_LHb" "C.24 DD_LHb" "C.30 DD_LHb"
# [11] "C.33 DD_LHb" "C.36 DD_MHb" "C.40 DD_LHb" "C.01"        "C.02"
# [16] "C.03"        "C.04"        "C.06"        "C.08"        "C.09"
# [21] "C.12"        "C.13"        "C.15"        "C.17"        "C.19"
# [26] "C.20"        "C.21"        "C.22"        "C.25"        "C.26"
# [31] "C.27"        "C.28"        "C.29"        "C.31"        "C.32"
# [36] "C.34"        "C.35"        "C.37"        "C.38"        "C.39"
# [41] "C.41"        "C.42"

DefaultAssay(SeuratOBJ) <- "RNA"

Seurat_base_name <- str_extract(seurat_name, regex("C\\.\\w+"))
Seurat_base_name
# C.leiden_lsi_r2_renamed_visium

## Read DEG file

DEG_file_name <- "WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_cellTypes_integrated_top50.csv"
DEG_file_name <- here(inputCVS_Dir, DEG_file_name)
df_cluster_names <- read.csv(DEG_file_name)
df_cluster_names <- df_cluster_names |> drop_na(cell_type)
head(df_cluster_names)
# p_val avg_log2FC pct.1 pct.2 p_val_adj cluster     gene            cell_type
# 1     0   4.336764 0.938 0.100         0       1 OTX2-AS1        DD_Inhib.Thal
# 2     0   3.849508 0.900 0.083         0       1      KIT        DD_Inhib.Thal
# 3     0   3.682782 0.927 0.122         0       1    MEIS2      LB_Thalamus/MDm



########## Gene markers lists. New function to join LB and DD gene markers lists

## Source gene markers lists 

source(here("code", "04_DiffExpr_Clustering_seurat", "remote_DGE_marker_gene_lists.R"))   

markers.custom <- get_multiple_markers_genes_lst()
tmp <- names(markers.custom)
tmp <- paste(tmp, collapse=', ')
message("Processing ", length(markers.custom), " categories of gene-markers list \n *****(", tmp, ")*****")

# ## Check duplicated marker genes
# names(markers.custom)
# x <- markers.custom
# length(unlist(x)) # 481
# # table(unname(unlist(x)))
# v_dup <- duplicated(unname(unlist(x)))
# dup_genes <- unname(unlist(x))[v_dup]


## merge DD markers for heatmap with 'Data-driven' markers and merge LB markers for heatmap with 'Literature-Based'markers

# prepare DD marker's list
DD_markers_lst <- append(markers.custom$DD_MHb, markers.custom$DD_LHb) 
# prepare LB marker's list
LB_markers_lst <- markers.custom[grep("^LB", names(markers.custom))]
LB_markers_lst <- as.vector(unlist(LB_markers_lst, recursive = FALSE))

## scale expression once for fair comparison across all clusters

# Keep only genes found in the dataset
dd_markers <- DD_markers_lst[DD_markers_lst %in% rownames(SeuratOBJ)]
lb_markers <- LB_markers_lst[LB_markers_lst %in% rownames(SeuratOBJ)]

if (length(dd_markers)==0 || length(lb_markers)==0) {
    stop()
} else {
    combined_named_list <- list(
        Hb_DD_markers = dd_markers,
        Hb_LB_markers = lb_markers
    )
        
}
names(combined_named_list)
# [1] "Hb_DD_markers" "Hb_LB_markers"


## prepare Seurat object

# Set cell type as the identity class for grouping
#colnames(SeuratOBJ@meta.data)
head(SeuratOBJ@meta.data$seurat_clusters)
Idents(SeuratOBJ) <- "seurat_clusters"
## add / tranfer cell-types if you want to plot by cell-type

# Multi-cell-type heatmap, scale all relevant genes once, across all cells
genes_to_scale <- as.vector(unlist(append(combined_named_list[1], combined_named_list[2])))
# SeuratOBJ <- NormalizeData(SeuratOBJ)
SeuratOBJ <- ScaleData(SeuratOBJ, features = genes_to_scale)

## plot all the list in separate heatmaps

heatmap_list <- list()

for (ct in seq_along(combined_named_list)) {

    # ct=1
    # grab cell type names and label
    ct_name <- names(combined_named_list[ct])
    print(paste0("Cell-type: ", ct_name))
    markers_to_plt <- unlist(combined_named_list[ct])
    # markers_to_plt=c("PRKCB"      "ADCY6"      "DCHS2"      "GLIS1"      "GREB1L")
    title_label <- paste0("Cell type class: ", ct_name)
    message("Preparing heatmap with class: ", ct_name)    
    
    # ####### this heatmap plot all Hb medial and lateral markers across all clusters - cluster size aware
    # heatmap_plot <- DoHeatmap(SeuratOBJ,
    #                           group.by = "seurat_clusters", 
    #                           features = markers_to_plt, size = 2,
    #                           disp.min = -2.5, disp.max = 2.5, 
    #                           group.bar = TRUE, # Omits color bar by identity class (from group.by)
    #                           slot = "scale.data") +
    #     scale_fill_gradientn(colors = c("blue", "white", "red")) +
    #     ggtitle(title_label) +
    #     theme(
    #         plot.title = element_text(hjust = 0.5),
    #         axis.text.y = element_text(size = 6),
    #         legend.position = "none"
    #     )
    # # disp.min / disp.max = 2.5: gene-expr after scaling can have extreme values. Clipping keeps the heatmap visually interpretable
    # print(heatmap_plot)
    
    
    ####### this heatmap plot all Hb medial and lateral markers across all clusters - same cluster size
    
    # Downsample to equal cell numbers per group
    # Number of cells per group you want (e.g., 50)
    n_cells <- 50
    
    # Randomly sample equal number of cells from each group to control column size
    # Add cell IDs as a column first (from rownames)
    meta_df <- SeuratOBJ@meta.data
    meta_df$cell_id <- rownames(meta_df)
    # Sample cells evenly across clusters
    cells_to_plot <- meta_df |>
        group_by(seurat_clusters) |>
        sample_n(size = min(n_cells, n()), replace = FALSE) |>
        arrange(seurat_clusters) |>   # This sets a fixed order to remove dendogram manually
        pull(cell_id)
    
    # Reorder markers if needed
    # ordered_markers <- markers_to_plt[markers_to_plt %in% rownames(SeuratOBJ)]
    
    # Plot with fixed number of cells per group (uniform column width)
    
    temp_plot <- DoHeatmap(SeuratOBJ,
                              features = markers_to_plt,
                              group.by = "seurat_clusters",
                              cells = cells_to_plot,
                              group.bar = TRUE,
                              label = TRUE,
                              size = 2,
                              disp.min = -2.5, disp.max = 2.5, 
                              slot = "scale.data") +
        scale_fill_gradientn(colors = c("blue", "white", "red")) +
        #ggtitle(title_label) +
        theme(
            plot.title = element_blank(),               # remove internal title
            #plot.title = element_text(hjust = 0.5, size = 8),
            plot.margin = margin(t = 20, r = 5, b = 30, l = 5),  # extra bottom space
            axis.text.y = element_text(size = 5),
            axis.text.x = element_blank(),        # hide x-axis labels (just in case)
            axis.ticks.x = element_blank(),       # remove x-axis ticks
            legend.position = "none"
        )

    # Optionally add horizontal lines for this specific marker group
    
    if (ct_name == "Hb_LB_markers") {
        # Note: y-axis is reversed — top gene = lowest y value
        # make a vector with size of gene blocks to add them to heatmpas - seperate cell types by marker class
        gene_blks <- markers.custom[grep("^LB", names(markers.custom))]
        gene_blks_lng <- map_int(gene_blks, ~ length(.x) )    
        gene_block_sizes <- as.vector(gene_blks_lng)
        cumulative_positions <- cumsum(gene_block_sizes)
        
        # Add lines at block boundaries (skip last)
        for (y_pos in cumulative_positions[-length(cumulative_positions)]) {
            temp_plot <- temp_plot +
                geom_hline(yintercept = y_pos + 0.5, color = "black", linetype = "solid", linewidth = 0.3)
        }
    }
    
    # Store the final plot in the list
    heatmap_list[[ct]] <- temp_plot
    
}

# Combine all heatmaps side by side

length(heatmap_list)

combined_plot <- wrap_plots(heatmap_list, ncol = length(heatmap_list)) +
    plot_annotation(
        title = "Habenula Data-Driven and Literature-Based gene markers side to side",
        theme = theme(
            plot.title = element_text(size = 10, hjust = 0.5, face = "bold")
        )
    )

f_name <- paste0(
    Seurat_base_name,
    "_heatmap_all_reference_markers_width5.pdf"
)
pdf(file = here(plotDir, f_name), width = 4 * length(heatmap_list), height = 6)
print(combined_plot)
dev.off()

length(heatmap_list)
#heatmap_list[1]
f_name <- paste0(
    Seurat_base_name,
    "_heatmap_all_reference_markers_width10.pdf"
)
combined_plot <- wrap_plots(heatmap_list, ncol = length(heatmap_list))
pdf(file = here(plotDir, f_name), width = 10 * length(heatmap_list), height = 6)
print(combined_plot)
dev.off()

## Subset Hb clusters. Use length of cluster ID as criteria
## extract clusters IDs and cluster label
## Subset Hb clusters. Use length of cluster ID as criteria
## extract clusters IDs and cluster label

# SeuOBJ_clusters <- Idents(SeuratOBJ)

# all_clusters <- levels(SeuOBJ_clusters)
# no_hb_clust <- all_clusters[nchar(all_clusters) <= 4]
# hb_clusters <- all_clusters[!all_clusters %in% c(no_hb_clust)]
# hb_clusters_ann <- hb_clusters
# hb_clusters_ann
# # [1] "C.05 DD_LHb" "C.07 DD_MHb" "C.10 DD_MHb" "C.11 DD_MHb" "C.14 DD_MHb"
# # [6] "C.16 DD_MHb" "C.18 DD_LHb" "C.23 DD_LHb" "C.24 DD_LHb" "C.30 DD_LHb"
# # [11] "C.33 DD_LHb" "C.36 DD_MHb" "C.40 DD_LHb"
# # length(hb_clusters)
# hb_clusters <- as.integer(substr(hb_clusters, 3, 4))
# hb_clusters
# 
# hb_df <- data.frame(
#     cluster = hb_clusters,
#     cluster_ann = hb_clusters_ann,
#     stringsAsFactors = FALSE 
# )
# hb_df

# ##########Heatmap 1: plot the top genes by cluster
# 
# f_name <- paste0(
#     Seurat_base_name,
#     "_heatmap_hb_top20_fdr5.pdf"
# )
# pdf(file = here(plotDir, f_name))
# 
# # Loop through all clusters in your list
# for (i in seq_along(hb_df$cluster_ann)) {
#     #i=2
#     # Get cluster ID and annotation
#     cluster_id <- hb_df$cluster[i]
#     cluster_label <- hb_df$cluster_ann[i]
#     
#     # Subset Seurat object to current cluster
#     seurat_subset <- subset(SeuratOBJ, idents = cluster_label)
#     
#     # Filter top 50 genes for this cluster
#     TopGenes <- df_cluster_names %>%
#         filter(cluster == cluster_id) %>%
#         top_n(n = 50, wt = avg_log2FC)
#     
#     # Filter genes that exist in Seurat object
#     TopGenes <- TopGenes %>% filter(gene %in% rownames(seurat_subset))
#     
#     # Plot heatmap
#     heatmap_plot <- DoHeatmap(seurat_subset, features = TopGenes$gene, size = 3) +
#         scale_fill_gradientn(colors = c("blue", "white", "red")) +
#         ggtitle(paste("Cluster", cluster_label))
#     print(heatmap_plot)
#     
# }
# dev.off()



## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
