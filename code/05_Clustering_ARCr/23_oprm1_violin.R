#   Violin plots of OPRM1 across mid and fine multiome clusters

library(here)
library(Seurat)
library(Signac)
library(DeconvoBuddies)
library(tidyverse)

seur_path = here(
    'processed-data', '05_Clustering_ARCr', '22_add_mid_level_clustering',
    'seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_HD.rds'
)
plot_dir = here('plots', '05_Clustering_ARCr', '23_oprm1_violin')
cluster_cols = c('cluster_ann', 'mid_cluster')
cluster_res = c('fine', 'mid')

dir.create(plot_dir, showWarnings = FALSE)

seur = readRDS(seur_path)

for (i in seq_along(cluster_cols)) {
    #   Order Mhb then LHb then others
    temp = unique(seur@meta.data[[cluster_cols[i]]]) |> as.character()
    hb_temp = sub('^C\\.[0-9]+\\.', '', temp)
    names(hb_temp) = temp
    seur@meta.data[[cluster_cols[i]]] = factor(
        seur@meta.data[[cluster_cols[i]]],
        levels = c(
            sort(hb_temp[grepl('MHb', hb_temp)]),
            sort(hb_temp[grepl('LHb', hb_temp)]),
            sort(hb_temp[!grepl('[ML]Hb', hb_temp)])
        ) |> names()
    )

    p = VlnPlot(
            object = seur,
            layer = "data",
            group.by = cluster_cols[i], 
            features = "OPRM1",
            pt.size = 0.2,
            alpha = 0.1
        ) +
        labs(x = "Cluster") +
        NoLegend()
    png(
        file.path(plot_dir, sprintf('%s.png', cluster_res[i])),
        width = 30 * length(temp), height = 500
    )
    print(p)
    dev.off()
}
