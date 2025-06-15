########################################################################
## Plot Heatmap, DimPlots, Jaccard, VPlots
##
## Authors. CSC
## Date. Jan 24, 2024
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
## Note. Seurat objects were created with module load conda_R/4.3.x
########################################################################

library("Seurat")
library("Signac")
library("ggplot2")
library("pheatmap")
library("bluster")
library("viridisLite")
library("patchwork")
library("ggplotify")
library("gridExtra")
library("tidyverse")
library("stringr")
library("here")

## directories

## clusters renamed for Spatial-Registration on Visium project

inputRDS_Dir <- here(
  "processed-data",
  "05_Clustering_ARCr",
  "05_rename_idents"
)
inputCVS_Dir <- here(
  "processed-data",
  "05_Clustering_ARCr",
  "02_Hb_celltypes_from_seurat_reanalyze_v3",
  "cvs_files_markers"
)
plotDir <- here(
    "plots",
    "05_Clustering_ARCr",
    "08_wnn_geneExp_plts_annotated"
)

## Check directories
if (!dir.exists(plotDir)) {
  dir.create(plotDir)
}

## Load input with RDS wnn to compare
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium.rds"

seurat_name <- here(inputRDS_Dir, Seurat_base_name)
title_name <- str_extract(seurat_name, regex("C\\.\\w*\\_r2"))

# Load Seurat
SeuratOBJ <- readRDS(here(inputRDS_Dir, Seurat_base_name))
levels(SeuratOBJ)

DefaultAssay(SeuratOBJ) <- "RNA"
Seurat_base_name <- str_extract(seurat_name, regex("C\\.\\w+"))
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


## To plot subset Seurat object where identity (cluster name) contains "MHb" or "LHb"

Seurat_subset <- subset(SeuratOBJ, idents = grep("MHb|LHb", Idents(SeuratOBJ), value = TRUE))
levels(Seurat_subset)
# [1] "C.05.DD_LHb" "C.07.DD_MHb" "C.10.DD_MHb" "C.11.DD_MHb" "C.14.DD_MHb"
# [6] "C.16.DD_MHb" "C.18.DD_LHb" "C.23.DD_LHb" "C.24.DD_LHb" "C.30.DD_LHb"
# [11] "C.33.DD_LHb" "C.36.DD_MHb" "C.40.DD_LHb"
table(Idents(Seurat_subset))

hb_clusters <- levels(Seurat_subset)
hb_numeric_cluster <- as.numeric(sub("^C\\.(\\d+)\\..*$", "\\1", hb_clusters))
hb_df <- data.frame(
    clusterID = hb_numeric_cluster,
    cluster_ann = hb_clusters,
    stringsAsFactors = FALSE
)
head(hb_df)
#   clusterID cluster_ann
# 1         5 C.05.DD_LHb
# 2         7 C.07.DD_MHb
# 3        10 C.10.DD_MHb
# 4        11 C.11.DD_MHb
# 5        14 C.14.DD_MHb
# 6        16 C.16.DD_MHb


## =============================================================================
## Prepare Violin Plot on top 10 gene markers

# Get top 10 genes per habenula cluster
# df_cluster_names$cluster
top10 <- df_cluster_names |>
    filter(cluster %in% hb_numeric_cluster) |>
    group_by(cluster) |>
    top_n(n = 10, wt = avg_log2FC)
head(top10)

# Start PDF output
f_name <- paste0(Seurat_base_name, "_VPlot_hb_top10_fdr5_by_cluster.pdf")
pdf(file = here(plotDir, f_name), width = 8.5, height = 11)  # standard letter size

# Loop through each habenula cluster
for (clus in unique(top10$cluster)) {

    message("Processing habenula cluster: ", clus)

    top10_cluster <- top10 |>
        filter(cluster == clus)

    # Create list of plots
    vln_plots <- lapply(top10_cluster$gene, function(gen) {
        VlnPlot(
            object = SeuratOBJ,
            layer = "data",
            features = gen,
            pt.size = 0
        ) +
            labs(title = gen) +
            theme(
                text = element_text(size = 7),
                axis.text.x = element_text(size = 5),
                axis.text.y = element_text(size = 5),
                plot.title = element_text(hjust = 0.5, size = 8)
            ) +
            NoLegend()
    })

    # Combine plots using patchwork
    combined_plot <- wrap_plots(vln_plots, ncol = 2) +
        plot_annotation(
            title = paste0("Top 10 genes in cluster: ", clus),
            theme = theme(plot.title = element_text(hjust = 0.5, size = 12))
        )

    print(combined_plot)
}

dev.off()


# ## =============================================================================
# ## Prepare Violin Plot on top 10 gene markers
# 
# BiocManager::install("DeconvoBuddies")
# BiocManager::install("LieberInstitute/DeconvoBuddies")
# # ERROR: this R is version 4.3.2, package 'DeconvoBuddies' requires R >=  4.4.0
# library("SingleCellExperiment")
# 
# sce_your_data <- as.SingleCellExperiment(SeuratOBJ)
# # top 10 genes per habenula cluster
# top10
# 
# library(DeconvoBuddies)
# 
# plot_marker_express(
#     sce = sce_your_data,
#     stat = top10,
#     cellType_col = "cell_type",    # adjust to match your column name
#     cell_type = "DD_LHb",       # the cell type you want to visualize
#     gene_col = "gene"              # column with gene symbols
# )
# 

## =============================================================================
## Heatmap: plot top.x genes by cluster

f_name <- paste0(
    Seurat_base_name,
    "_heatmap_hb_top50_fdr5.pdf"
)
pdf(file = here(plotDir, f_name))
levels(Seurat_subset)
# Loop through all clusters in your list
for (i in seq_along(hb_df$cluster_ann)) {
    #i=2
    # Get cluster ID and annotation
    cluster_id <- hb_df$clusterID[i]
    cluster_label <- hb_df$cluster_ann[i]
    message("Processing cluster:", cluster_label)
    
    seurat_cluster <- subset(Seurat_subset, idents = cluster_label)
    TopGenes <- df_cluster_names %>%
        filter(cluster == cluster_id) %>%
        top_n(n = 50, wt = avg_log2FC)
    
    # Filter genes that exist in Seurat object
    TopGenes <- TopGenes %>% filter(gene %in% rownames(seurat_cluster))
    
    heatmap_plot <- DoHeatmap(seurat_cluster, features = TopGenes$gene, size = 3) +
        scale_fill_gradientn(colors = c("blue", "white", "red")) +
        ggtitle(paste("Cluster", cluster_label))
    print(heatmap_plot)

}
dev.off()


## =============================================================================
##  compuused Jaccard for RNA and ATAC

plot_wnn_rna_heatmap <- function(
        jaccard_matrix,
        title
) {
    pheatmap::pheatmap(
        mat = jaccard_matrix,
        main = title,
        color_palette = viridisLite::plasma(101),
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        show_numbers = TRUE,
        number_format = "%.1f",
        fontsize = 10,
        fontsize_number = 7,
        na_color = "black",
        angle_col = 90,
        legend = TRUE
    )
}

##  compute Jaccard for RNA
clust.wnn1 <- as.vector(SeuratOBJ$seurat_clusters)
clust.wnn2 <- as.vector(SeuratOBJ$C.leiden)
jacc.mat <- linkClustersMatrix(clust.wnn1, clust.wnn2)
# rownames(jacc.mat)
colnames(jacc.mat) <- paste0("RNA.C.", colnames(jacc.mat))
plt_wnn_rna <- plot_wnn_rna_heatmap(jaccard_matrix = jacc.mat, "Overlap between RNA and WNN cluster identities")

##  compute Jaccard for ATAC
clust.wnn2 <- as.vector(SeuratOBJ$C.leiden_atac)
jacc.mat <- linkClustersMatrix(clust.wnn1, clust.wnn2)
colnames(jacc.mat) <- paste0("ATAC.C.", colnames(jacc.mat))
plt_wnn_atac <- plot_wnn_rna_heatmap(jaccard_matrix = jacc.mat, "Overlap between ATAC and WNN cluster identities")

# arrange and save plots
plot_list <- list()
plot_list[['rna']] <- as.ggplot(plt_wnn_rna)
plot_list[['atac']] <- as.ggplot(plt_wnn_atac)
g <- grid.arrange(grobs = plot_list, ncol = 2)

tmp_png <- paste0(Seurat_base_name, "_Jaccard_WNN_RNA_ATAC.png")
ggsave(g, filename = here(plotDir, tmp_png), height = 6, width = 17)


## =============================================================================
## DimPlot ATAC, RNA and WNN annotated

seu_c <- as.vector(head(SeuratOBJ$seurat_clusters))
wnn_c <- as.vector(head(SeuratOBJ$`C.leiden_wnn`))
rna_c <- as.vector(head(SeuratOBJ$`C.leiden`))
atac_c <- as.vector(head(SeuratOBJ$`C.leiden_atac`))
tbl <- data.frame(rna_c, atac_c, wnn_c, seu_c)
#     rna_c atac_c wnn_c seu_c
# 1    24     11    25  C.25
# 2     4     11     4  C.04
# 3    15      4     9  C.09
# 4     1      9     4  C.04
# 5     2     11     1  C.01
# 6     2      8     1  C.01
SeuratOBJ@meta.data
clust_name = "seurat_clusters"
plt1 <- DimPlot(
  SeuratOBJ,
  reduction = "umap.integrated",
  group.by = clust_name,
  label = TRUE,
  label.size = 2.5,
  repel = TRUE
) +
  ggtitle(
    "RNA",
    subtitle = paste(
      " Leiden at res=2 knn=30; SNN Clusters=",
      length(table(SeuratOBJ[["C.leiden"]])),
      "\nAnnotated by WNN clusters"
    )
  ) &
  NoLegend()

plt2 <- DimPlot(
  SeuratOBJ,
  reduction = "umap.lsi.integrated",
  #group.by = clust_name_atac,
  group.by = clust_name,
  label = TRUE,
  label.size = 2.5,
  repel = TRUE
) +
  ggtitle(
    "ATAC",
    subtitle = paste(
      " Leiden at res=2 knn=30; SNN Clusters=",
      length(table(SeuratOBJ[["C.leiden_atac"]])),
      "\nAnnotated by WNN clusters"
    )
  ) &
  NoLegend()

plt3 <- DimPlot(
  SeuratOBJ,
  reduction = "wnn.umap",
  group.by = clust_name,
  label = TRUE,
  label.size = 2.5
) +
  ggtitle(
    "WNN",
    subtitle = paste(
      "Leiden at res=2 knn=30; WNN Clusters=",
      length(table(SeuratOBJ[[clust_name]]))
    )
  )

pltALL <- plt1 + plt2 + plt3 & theme(plot.title = element_text(hjust = 0.5))
tmp_png <- paste0(Seurat_base_name, "_RNA_ATAC_WNN_DimPlots.png")
ggsave(
  pltALL,
  filename = here(
    plotDir,
    tmp_png
  ),
  height = 7,
  width = 20
)


## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
