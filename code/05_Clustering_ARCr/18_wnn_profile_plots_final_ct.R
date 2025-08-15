########################################################################
## Script on process .... Plot Violin for 1vsALL and mean-ratio
##
## Authors. CSC
## Date. Jul 20, 2025
##
## Recommended resources on interactive mode: srun --pty --mem=80GB --x11 bash
########################################################################

library("Seurat")
library("purrr")
library("ggplot2")
library("colorspace") # make color gradients 
library("patchwork")
library("dplyr")
library("stringr")
library("here")

# directories

Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_HD.rds"
inputSeuratRDS <- here(
    "processed-data", 
    "05_Clustering_ARCr", 
    "17_wnn_clustering_final_ct", 
    Seurat_base_name)

processedDir <- here(
    "processed-data",
    "05_Clustering_ARCr",
    "18_wnn_profile_plots_final_ct"
)

plotDir <- here(
    "plots",
    "05_Clustering_ARCr",
    "18_wnn_profile_plots_final_ct"
)

## Check directories
if (!dir.exists(plotDir)) {
    dir.create(plotDir)
}
if (!dir.exists(processedDir)) {
    dir.create(processedDir)
}

#===============================================================================

message(Sys.time(), " - Load Seurat with WNN clusters")

SeuratOBJ <- readRDS(inputSeuratRDS)
SeuratOBJ
DefaultAssay(SeuratOBJ) <- "RNA"

## some verification
class(SeuratOBJ[["ATAC"]])
levels(SeuratOBJ)
colnames(SeuratOBJ@meta.data)

# check cell-types
table(SeuratOBJ$seurat_clusters)

table(sort(SeuratOBJ$cluster_ann))

table(SeuratOBJ$merged_cluster)
# Astrocyte       Endo Excit_Thal Inhib_Thal        LHb        MHb  Microglia 
# 2684        343       9738       6024      19673      10944        663 
# Oligo        OPC       Thal 
# 4882        638        111 

levels(SeuratOBJ)
# [1] "C.04.LHb.4"      "C.05.LHb.2.7"    "C.06.LHb.4"      "C.07.MHb.2"     
# [5] "C.08.LHb.4"      "C.09.LHb.4"      "C.10.MHb.1"      "C.11.MHb.1.2"   
# [9] "C.13.LHb.4"      "C.14.MHb.1"      "C.16.MHb.1.2"    "C.18.LHb.1.3.4" 
# [13] "C.23.LHb.1"      "C.24.LHb.4"      "C.30.LHb.7"      "C.31.LHb.4"     
# [17] "C.33.LHb.1.3"    "C.36.MHb.3"      "C.40.LHb.4"      "C.01.Inhib.Thal" ....

## =============================================================================

## merged clusters, I picked up Hex color codes similar to those used on human pilot
my_colors <- c(
    LHb = "#1f78b4",
    MHb = "#ad1d8c",
    Oligo = "#384a08",
    Astrocyte = "#532222", 
    OPC = "#829454",
    Microglia = "#141b02",
    Endo = "#d95f02",
    Inhib_Thal = "#9a9fe7",
    Excit_Thal = "#42467b",
    Thal = "#4d55b7"
)

## assign color gradients to fine resolution clusters based on Broad cell-types
# extract LHb and MHb clusters
cluster_levels <- levels(SeuratOBJ)
LHb_clusters <- grep("LHb", cluster_levels, value = TRUE)
MHb_clusters <- grep("MHb", cluster_levels, value = TRUE)
# Create tonal gradients for LHb and MHb
LHb_colors <- sequential_hcl(length(LHb_clusters), h = 210, c = 80, l = c(30, 80))
MHb_colors <- sequential_hcl(length(MHb_clusters), h = 320, c = 80, l = c(30, 80))
# Build full cluster color map
my_colors_fine <- setNames(rep("#bdbdbd", length(cluster_levels)), cluster_levels)
my_colors_fine[LHb_clusters] <- LHb_colors
my_colors_fine[MHb_clusters] <- MHb_colors
# Assign base color for other types from your existing palette
for (category in c("Oligo", "Astrocyte", "OPC", "Microglia", "Endo", "Inhib.Thal", "Excit.Thal", "Thal")) {
    matched <- grep(category, cluster_levels, value = TRUE)
    my_colors_fine[matched] <- my_colors[[gsub("\\.", "_", category)]]
}
#scales::show_col(my_colors_fine)

## =============================================================================


message("Processing Violin plots ...")

## Plot Hb canonical genes for merged_clusters

plot_violin_clusters <- function(seurat_obj, genes, group_col = "merged_cluster", colors = NULL) {
    
    plots <- purrr::map(genes, ~ {
        VlnPlot(
            object = SeuratOBJ,
            layer = "data",
            group.by = group_col, 
            features = .x,
            pt.size = 0.2,
            alpha = 0.1,
            cols = colors
        ) +
            labs(title = .x) +
            theme(
                text = element_text(size = 10),
                axis.text.x = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 1),
                axis.text.y = element_text(size = 10),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                plot.title = element_text(hjust = 0.5, size = 12)
            ) +
            coord_flip() +
            NoLegend()
    }) |> purrr::set_names(genes)
    
    return(plots)
    
}

## make violin plots for merged clusters - with solid color vector 

genes_to_plot <- c("GPR151", "POU4F1", "TAC3")

my_plots <- plot_violin_clusters(SeuratOBJ, genes_to_plot, colors = my_colors)
plt1 <- my_plots[["GPR151"]] + my_plots[["POU4F1"]] + my_plots[["TAC3"]] 
plt1[[3]]
ggsave(here(plotDir, "WNN_Vplots_Hb_canonical_broad_clusters.pdf"), plt1, width = 6, height = 7)

## make violin plots for all clusters detail - with gradient tonalities for MHb and LHb, other cell-types solid color
my_plots <- plot_violin_clusters(SeuratOBJ, group_col = "cluster_ann", genes_to_plot, colors = my_colors_fine)
plt1 <- my_plots[["GPR151"]] + my_plots[["POU4F1"]] + my_plots[["TAC3"]] 
plt1[[1]]

ggsave(here(plotDir, "WNN_Vplots_Hb_canonical_fine_clusters.pdf"), plt1, width = 6, height = 7)

message("Violin plots done!")

## =============================================================================


# library("slurmjobs")
# job_single(
#     "17_wnn_hierarchical_clustering_final_ct", 
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



