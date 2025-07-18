########################################################################
## Plot Heatmaps of GEX on WNN clustering
## - Downsample to equal cell numbers per group (n=50)
## - Plot top50 DEG match cell-types by LB or DD marker genes 
## Authors. CSC
## Date. May 08, 2024
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
## Note. Seurat objects were created with module load conda_R/4.3.x

## Note. this pipeline requires module load conda_R/4.3.x to keep the integrity of the seurat object

########################################################################

library("Seurat")
library("Signac")
library("Matrix")
library("ComplexHeatmap")
library("circlize")
# library("ggplot2")
# library("bluster")
# library("patchwork")
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
    #"08_wnn_gene_expression_plts_renamed_idents"
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
DefaultAssay(SeuratOBJ) <- "RNA"
levels(SeuratOBJ)
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

base_name <- str_extract(seurat_name, regex("C\\.\\w+"))

message("Processing Heatmap for ", base_name)


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

# Find top gene per cluster (highest avg_log2FC or pct diff)
top_markers <- df_cluster_names |>
    group_by(cluster) |>
    top_n(n = 3, wt = avg_log2FC)
head(top_markers)

avg_expr <- AggregateExpression(
    SeuratOBJ, 
    group.by = "cluster_ann", 
    return.seurat = FALSE)$RNA
head(avg_expr)

# Subset only for top genes
mat <- avg_expr[unique(top_markers$gene), ]

# order genes per cluster order
ordered_genes <- top_markers |>
    arrange(factor(cluster, levels = colnames(mat)), desc(avg_log2FC))

mat_ordered <- mat[ordered_genes, ]

# Scale across rows (genes)
mat_scaled <- t(scale(t(as.matrix(mat_ordered))))

# Heatmap
Heatmap(
    mat_scaled,
    name = "Z-score",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
    row_names_gp = gpar(fontsize = 8)
)



########## Gene markers lists. New function to join LB and DD gene markers lists

## Source gene markers lists 

source(here("code", "04_DiffExpr_Clustering_seurat", "remote_DGE_marker_gene_lists.R"))   

markers.custom <- get_multiple_markers_genes_lst()
tmp <- names(markers.custom)
tmp <- paste(tmp, collapse=', ')
message("Processing ", length(markers.custom), " categories of gene-markers list \n *****(", tmp, ")*****")


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
colnames(SeuratOBJ@meta.data)
head(SeuratOBJ@meta.data$cluster_ann)
unique(SeuratOBJ@meta.data$cluster_ann)
#Idents(SeuratOBJ) <- "seurat_clusters"

# Multi-cell-type heatmap, scale all relevant genes once, across all cells
genes_to_scale <- as.vector(unlist(append(combined_named_list[1], combined_named_list[2])))
SeuratOBJ <- ScaleData(SeuratOBJ, features = genes_to_scale)



######## Plot 1: plot the 2 gene-markers list (DD, LB) across all clusters in separate heatmaps

heatmap_list <- list()

for (ct in seq_along(combined_named_list)) {

    # ct=2
    # grab cell type names and label
    ct_name <- names(combined_named_list[ct])
    print(paste0("Cell-type: ", ct_name))
    markers_to_plt <- unlist(combined_named_list[ct])
    # markers_to_plt=c("PRKCB"      "ADCY6"      "DCHS2"      "GLIS1"      "GREB1L")
    title_label <- paste0("Cell type class: ", ct_name)
    message("Preparing heatmap with class: ", ct_name)    
    
    # ####### this heatmap plot all Hb medial and lateral markers across all clusters - cluster size aware
    # heatmap_plot <- DoHeatmap(SeuratOBJ,
    #                           group.by = "cluster_ann",
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
        group_by(cluster_ann) |>
        sample_n(size = min(n_cells, n()), replace = FALSE) |>
        arrange(cluster_ann) |>   # This sets a fixed order to remove dendogram manually
        pull(cell_id)
    length(cells_to_plot)
    # Reorder markers if needed
    # ordered_markers <- markers_to_plt[markers_to_plt %in% rownames(SeuratOBJ)]
    
    # Plot with fixed number of cells per group (uniform column width)
    
    temp_plot <- DoHeatmap(SeuratOBJ,
                              features = markers_to_plt,
                              group.by = "cluster_ann",
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
    print(temp_plot)
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
    base_name,
    "_heatmap_all_reference_markers_width5.pdf"
)
pdf(file = here(plotDir, f_name), width = 4 * length(heatmap_list), height = 6)
print(combined_plot)
dev.off()

length(heatmap_list)
#heatmap_list[1]
f_name <- paste0(
    base_name,
    "_heatmap_all_reference_markers_width10.pdf"
)
combined_plot <- wrap_plots(heatmap_list, ncol = length(heatmap_list))
pdf(file = here(plotDir, f_name), width = 10 * length(heatmap_list), height = 6)
print(combined_plot)
dev.off()



######## Plot 2: plot the top-x genes (markers) used to annotate each cluster

# df_cluster_names$cluster
# df_cluster_names$gene
# df_cluster_names$cell_type

# Downsample to equal cell numbers per group
# Number of cells per group you want (e.g., 50)
n_cells <- 50

# Randomly sample equal number of cells from each group to control column size
# Add cell IDs as a column first (from rownames)
meta_df <- SeuratOBJ@meta.data
meta_df$cell_id <- rownames(meta_df)
# Sample cells evenly across clusters
cells_to_plot <- meta_df |>
    group_by(cluster_ann) |>
    sample_n(size = min(n_cells, n()), replace = FALSE) |>
    arrange(cluster_ann) |>   # This sets a fixed order to remove dendrogram manually
    pull(cell_id)

# extract top 10 genes per cluster
top_markers <- df_cluster_names %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)

# Unique gene list
marker_genes <- unique(top_markers$gene)
marker_genes <- marker_genes[marker_genes %in% rownames(SeuratOBJ)]
#SeuratOBJ <- ScaleData(SeuratOBJ, features = marker_genes, verbose = FALSE)

#Idents(SeuratOBJ) <- "cluster_ann"

plt <- DoHeatmap(SeuratOBJ,
          features = marker_genes,
          group.by = "cluster_ann",
          cells = cells_to_plot,
          group.bar = TRUE,
          size = 3) +
    scale_fill_gradientn(colors = c("blue", "white", "red")) +
    #ggtitle("Top Marker Genes per Cluster") +
    theme(#plot.title = element_text(hjust = 0.5, size = 8),
          plot.margin = margin(t = 20, r = 5, b = 30, l = 5),  # extra bottom space
          axis.text.y = element_text(size = 5),
          axis.text.x = element_blank(),        # hide x-axis labels (just in case)
          axis.ticks.x = element_blank(),       # remove x-axis ticks
          legend.position = "none")

f_name <- paste0(
    base_name,
    "_heatmap_top10genes.pdf"
)
pdf(file = here(plotDir, f_name), width = 5 * length(heatmap_list), height = 6)
print(plt)
dev.off()

message("All plots done!")


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
