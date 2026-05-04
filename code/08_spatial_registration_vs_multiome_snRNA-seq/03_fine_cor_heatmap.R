#   Correlation heatmap: registration of up-to-date fine multiome annotations
#   against Yalcinbas (fine)

library(here)
library(tidyverse)
library(spatialLIBD)
library(sessioninfo)

plot_dir = here('plots', '08_spatial_registration_vs_multiome_snRNA-seq')
model_path = here(
    'processed-data', '11_link_prep', '04_registration_wrapper',
    'model_results_fine.rds'
)
ref_path = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/05_snRNA-seq_model_stats/enrichment_final_Annotations.rds"
out_path = here(
    'processed-data', '08_spatial_registration_vs_multiome_snRNA-seq',
    'fine_vs_yalcinbas.rds'
)

dir.create(plot_dir, showWarnings = FALSE)
dir.create(dirname(out_path), showWarnings = FALSE)

#   Read in enrichment stats for each dataset
t_stats = readRDS(model_path)$enrichment
results_enrichment = list(enrichment = readRDS(ref_path))

this_cor = layer_stat_cor(
    t_stats, modeling_results = results_enrichment, model_type = "enrichment",
    top_n = 100
)

#   Annotate clusters
annotated_clusters = annotate_registered_clusters(
    this_cor, confidence_threshold = 0.25, cutoff_merge_ratio = 0.25
)

#   Make heatmaps
pdf(file.path(plot_dir, "fine_vs_yalcinbas.pdf"))
layer_stat_cor_plot(
    this_cor, annotation = annotated_clusters,
    heatmap_legend_param = list(title = "Cor", at = c(-1, 0, 1))
)
dev.off()

saveRDS(this_cor, file = out_path)

session_info()
