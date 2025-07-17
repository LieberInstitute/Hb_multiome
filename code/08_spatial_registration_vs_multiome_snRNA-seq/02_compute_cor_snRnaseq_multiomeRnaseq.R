########################################################################
## Compute Spatial-Registration between snRNAseq (human pilot) and multiome scRNAseq + scATACseq
##
## Notes:
## For 60 to 80k spots: $srun --pty --mem=60GB --x11 bash
##
## Authors. CSC
##
#################### Compute correlation of Fine clusters snRNAseq vs Multiome snRNAseq ##########################

library("here")
library("purrr")
library("spatialLIBD")
library("ComplexHeatmap")
library("grid") # plot title clipped by internal function of layer_stat_cor_plot()
library("sessioninfo")
library("tidyverse")


## Input / Output dirs
rds_input <- here(
  "processed-data",
  "08_spatial_registration_vs_multiome_snRNA-seq",
  "enrichment_snRNA-multiome_v5.rds"
  )
## Create output directories
dir_rdata <- here(
  "processed-data",
  "08_spatial_registration_vs_multiome_snRNA-seq"
  )
dir_plot <- here(
    "plots", 
    "08_spatial_registration_vs_multiome_snRNA-seq"
    )
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_plot, showWarnings = FALSE, recursive = TRUE)

## read input arguments (plot vertical or horizontal format)
args = commandArgs(trailingOnly = TRUE)

designF = args[2]
# design_format=("vertical" "horizontal")
# designF = "vertical"
## x-axis = snRNAseq cell-types
## y-axis = snRNAseq multiome cell-types

message("Compute correlations and Spatial-Registration in ", designF, " design")

## load snRNAseq t-stats enrichment data (fine resolution)

results_enrichment <- readRDS(here(dir_rdata, "enrichment_final_Annotations.rds"))
colnames(results_enrichment)

## filter only enrichment t-stats

modeling_res_enrichment <- results_enrichment[, grep(
  "^t_stat_",
  colnames(results_enrichment)
)]
colnames(modeling_res_enrichment) <- gsub(
  "^t_stat_",
  "",
  colnames(modeling_res_enrichment)
)
modeling_res_enrichment[1:3,]
#                 Astrocyte        Endo Excit.Thal Inhib.Thal      LHb.1
# ENSG00000238009  0.6873346  1.23174976  1.2820563   1.526893 -0.6502180
# ENSG00000241860  0.7626024 -2.42952913  1.4923036   2.026956  1.2102714
# ENSG00000237491 -0.2347861 -0.06842523  0.8264148   1.981078 -0.1328745


################################################################################
##  Compute correlation for FINE cluster annotations
################################################################################

if (designF == "vertical") {
  # designF="vertical"

  ## x-axis = snRNAseq cell-types
  ## y-axis = snRNAseq multiome cell-types

  results_enrichment_multiome <- readRDS(rds_input)$enrichment |>
    filter(!duplicated(ensembl))
  head(rownames(results_enrichment_multiome))
  class(results_enrichment_multiome) # [1] "data.frame"
  head(results_enrichment_multiome[5:10])
  #                 t_stat_C.05.DD_LHb t_stat_C.06 t_stat_C.07.DD_MHb t_stat_C.08
  # ENSG00000238009          1.2234449    1.759594         -0.5585619  2.26415183
  # ENSG00000241860          0.2151696    1.234082         -0.2209549  0.07914740
  # ENSG00000237491          0.4194972    1.935467         -0.7272817  1.36389851

  cor_fine <- layer_stat_cor(
    stats = results_enrichment_multiome, # data.frame / y-axis
    modeling_results = list(enrichment = results_enrichment), # list / x-axis
    model_type = "enrichment",
    top_n = 100
  )
  #head(cor_fine)
  annotated_clusters_fine <- annotate_registered_clusters(
    cor_fine,
    confidence_threshold = 0.25,
    cutoff_merge_ratio = 0.25
  )

  # ## Use annotation labels on the correlation matrices
  # rownames(cor_fine) <- paste0(rownames(cor_fine), " ~ ",
  #                              annotated_clusters_fine$layer_label[match(rownames(cor_fine), annotated_clusters_fine$cluster)])
  #
  # ## With default confidence and cutoff_merge_ratio
  # annotated_clusters_fine <- annotate_registered_clusters(cor_fine, confidence_threshold = 0.25, cutoff_merge_ratio = 0.25)

  head(cor_fine)
  #           C.01        C.02       C.03        C.04  C.05.DD_LHb        C.06
  # LHb.2 -0.1959626 -0.06687099 -0.1797187 -0.04842418  0.636548725  0.01036223
  # LHb.7 -0.1068896 -0.01307742 -0.1562915 -0.02322458  0.735296169 -0.01500998
  # LHb.6 -0.1717414  0.01985625 -0.1635756 -0.22551095  0.147718838 -0.22087233
} else {
  # designF=="horizontal"

  ## x-axis = snRNAseq multiome cell-types
  ## y-axis = snRNAseq cell-types

  ## load multiome snRNAseq t-stats enrichment data

  sn_multiome_data <- readRDS(rds_input)
  head(sn_multiome_data$enrichment[5:10])
  #                 t_stat_C.05.DD_LHb t_stat_C.06 t_stat_C.07.DD_MHb t_stat_C.08
  # ENSG00000238009          1.2234449    1.759594         -0.5585619  2.26415183
  # ENSG00000241860          0.2151696    1.234082         -0.2209549  0.07914740
  # ENSG00000237491          0.4194972    1.935467         -0.7272817  1.36389851

  cor_fine <- layer_stat_cor(
    stats = modeling_res_enrichment, # data.frame
    modeling_results = sn_multiome_data,
    model_type = "enrichment",
    top_n = 100
  )

  annotated_clusters_fine <- annotate_registered_clusters(
    cor_fine,
    confidence_threshold = 0.25,
    cutoff_merge_ratio = 0.25
  )
  head(modeling_res_enrichment)
}

## save t-stats rData and Plots

f_name <- paste0("cor_multiome_vs_snRNA-seq_top100_", designF, ".Rdata")
save(cor_fine, file = file.path(dir_rdata, f_name))

##   Make heatmaps fine clusters snRNAseq vs Multiome snRNAseq

plt_name <- paste0("cor_top100_registration_snMultiome_snRNAseq_v2_", designF, ".pdf")
# cor_top100_registration_snMultiome_snRNAseq_v2_vertical.pdf
pdf(here(dir_plot, plt_name))

layer_stat_cor_plot(
  cor_fine,
  annotation = annotated_clusters_fine,
  heatmap_legend_param = list(title = "Cor", at = c(-1, 0, 1)),
  column_names_gp = gpar(fontsize = 10),
  row_names_gp = gpar(fontsize = 10)
)

dev.off()


message("Spatial Registration for FINE cluster annotations done!")


################################################################################
##  Compute correlation for FINE cluster annotations
##  Subset HB cell-types of interest
################################################################################

# Prepare matrix 
# subset columns that contain "LHb", "MHb", or "Thal"? in the matrix 
colnames(cor_fine)
rownames(cor_fine)
cor_fine_subset <- cor_fine[, grep("LHb|MHb", colnames(cor_fine))]
colnames(cor_fine_subset)
rownames(cor_fine_subset)
# subset rows in the matrix that contain "LHb", "MHb", or "Thal"
rownames(cor_fine_subset)
cor_fine_subset <- cor_fine_subset[grep("LHb|MHb", rownames(cor_fine_subset)), ]
colnames(cor_fine_subset)
rownames(cor_fine_subset)

# arrange rownnames for visualization purposes, first "MHb" and then by number of cluster
rows <- rownames(cor_fine_subset)
is_mhb <- grepl("MHb", rows)
get_num <- function(x) as.numeric(sub("C\\.(\\d+)\\..*", "\\1", x))
cor_fine_subset <- cor_fine_subset[sorted_rows, ]
# extract each Hb group
mhb_rows <- rows[is_mhb]
lhb_rows <- rows[!is_mhb]
# apply numeric sort
mhb_sorted <- mhb_rows[order(get_num(mhb_rows))]
lhb_sorted <- lhb_rows[order(get_num(lhb_rows))]
# combine the final order
sorted_rows <- c(mhb_sorted, lhb_sorted)
# apply to the matrix
cor_fine_subset <- cor_fine_subset[sorted_rows, ]
colnames(cor_fine_subset)
rownames(cor_fine_subset)
head(cor_fine_subset)

# extract Hb annotations of interest from 'annotated_clusters_fine' 
annotated_clusters_fine_subset <- annotated_clusters_fine[grepl("LHb|MHb", annotated_clusters_fine$cluster), ]
head(annotated_clusters_fine_subset)
# cluster layer_confidence       layer_label
# 1  C.05.DD_LHb             good       LHb.7/LHb.2
# 2  C.18.DD_LHb             good LHb.1/LHb.3/LHb.4
# 3  C.23.DD_LHb             good             LHb.1
# 4  C.33.DD_LHb             good       LHb.3/LHb.1
# 33 C.16.DD_MHb             good             LHb.6
# 34 C.10.DD_MHb             good             MHb.1


# arrange rownnames for visualization purposes, first "MHb" and then by number of cluster
rows <- annotated_clusters_fine_subset$cluster
is_mhb <- grepl("MHb", rows)
get_num <- function(x) as.numeric(sub("C\\.(\\d+)\\..*", "\\1", x))
cor_fine_subset <- cor_fine_subset[sorted_rows, ]
# extract each Hb group
mhb_rows <- rows[is_mhb]
lhb_rows <- rows[!is_mhb]
# apply numeric sort
mhb_sorted <- mhb_rows[order(get_num(mhb_rows))]
lhb_sorted <- lhb_rows[order(get_num(lhb_rows))]
# combine the final order
sorted_rows <- c(mhb_sorted, lhb_sorted)
sorted_rows
# reorder the annotated_clusters_fine_subset
annotated_clusters_fine_subset <- annotated_clusters_fine_subset[
    match(sorted_rows, annotated_clusters_fine_subset$cluster),
]
annotated_clusters_fine_subset$cluster

## plot heatmap
plt_name <- "cor_top100_registration_snMultiome_snRNAseq_Habenula_clusters.pdf"
pdf(here(dir_plot, plt_name), width = 10, height = 10)

layer_stat_cor_plot(
    cor_fine_subset,
    annotation = annotated_clusters_fine_subset,
    heatmap_legend_param = list(title = "Cor", at = c(-1, 0, 1)),
    column_names_gp = gpar(fontsize = 14),
    row_names_gp = gpar(fontsize = 14),
    #cluster_rows = FALSE  # <-- turn off row clustering
) 
# hm + draw( # # Draw the heatmap with title
#     hm,
#     column_title = "Spatial-Registration: LHb, MHb, and Thal",
#     column_title_gp = gpar(fontsize = 16, fontface = "bold"))


dev.off()






# library("slurmjobs")
#
# ## A regular job with 10 cores on the 'imaginary' partition
# job_single("01_compute_cor_snRnaseq_multiomeRnaseq", cores = 2, partition = "katun", create_shell = TRUE)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# > session_info()
# [90mCRAN (R 4.4.0)
#     beachmat               2.22.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
#     beeswarm               0.4.0     2021-06-01 [2] CRAN (R 4.4.0)
#     benchmarkme            1.0.8     2022-06-12 [2] CRAN (R 4.4.0)
#     benchmarkmeData        1.0.4     2020-04-23 [2] CRAN (R 4.4.0)
#     Biobase              * 2.66.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
#     BiocFileCache          2.14.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
#     BiocGenerics         * 0.52.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
#     BiocIO                 1.16.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
#     BiocManager            1.30.25   2024-08-28 [2] CRAN (R 4.4.1)
#     BiocNeighbors          2.0.1     2024-11-28 [2] Bioconductor 3.20 (R 4.4.2)
#     BiocParallel           1.40.2    2025-04-10 [2] Bioconductor
#     BiocSingular           1.22.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
#     BiocVersion            3.20.0    2024-05-01 [2] Bioconductor 3.20 (R 4.4.0)
#     Biostrings             2.74.1    2024-12-16 [2] Bioconductor 3.20 (R 4.4.2)
#     bit                    4.6.0     2025-03-06 [2] CRAN (R 4.4.3)
#     bit64                  4.6.0-1   2025-01-16 [2] CRAN (R 4.4.2)
#     bitops                 1.0-9     2024-10-03 [2] CRAN (R 4.4.1)
#     blob                   1.2.4     2023-03-17 [2] CRAN (R 4.4.0)
#     bslib                  0.9.0     2025-01-30 [2] CRAN (R 4.4.2)
#     cachem                 1.1.0     2024-05-16 [2] CRAN (R 4.4.0)
#     Cairo                  1.6-2     2023-11-28 [2] CRAN (R 4.4.0)
#     circlize               0.4.16    2024-02-20 [2] CRAN (R 4.4.0)
#     cli                    3.6.5     2025-04-23 [2] CRAN (R 4.4.3)
#     clue                   0.3-66    2024-11-13 [2] CRAN (R 4.4.2)
#     cluster                2.1.8     2024-12-11 [3] CRAN (R 4.4.3)
#     codetools              0.2-20    2024-03-31 [3] CRAN (R 4.4.3)
#     colorspace             2.1-1     2024-07-26 [2] CRAN (R 4.4.1)
#     ComplexHeatmap       * 2.22.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
#     config                 0.3.2     2023-08-30 [2] CRAN (R 4.4.0)
#     cowplot                1.1.3     2024-01-22 [2] CRAN (R 4.4.0)
#     crayon                 1.5.3     2024-06-20 [2] CRAN (R 4.4.1)
#     curl                   6.2.2     2025-03-24 [2] CRAN (R 4.4.3)
#     data.table             1.17.2    2025-05-12 [2] CRAN (R 4.4.3)
#     DBI                    1.2.3     2024-06-02 [2] CRAN (R 4.4.0)
#     dbplyr                 2.5.0     2024-03-19 [2] CRAN (R 4.4.0)
#     DelayedArray           0.32.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
#     dichromat              2.0-0.1   2022-05-02 [2] CRAN (R 4.4.0)
#     digest                 0.6.37    2024-08-19 [2] CRAN (R 4.4.1)
#     doParallel             1.0.17    2022-02-07 [2] CRAN (R 4.4.0)
#     dplyr                * 1.1.4     2023-11-17 [2] CRAN (R 4.4.0)
#     DT                     0.33      2024-04-04 [2] CRAN (R 4.4.0)
#     edgeR                  4.4.2     2025-01-27 [2] Bioconductor 3.20 (R 4.4.2)
#     ExperimentHub          2.14.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
#     farver                 2.1.2     2024-05-13 [2] CRAN (R 4.4.0)
#     fastmap                1.2.0     2024-05-15 [2] CRAN (R 4.4.0)
#     filelock               1.0.3     2023-12-11 [2] CRAN (R 4.4.0)
#     forcats              * 1.0.0     2023-01-29 [2] CRAN (R 4.4.0)
#     foreach                1.5.2     2022-02-02 [2] CRAN (R 4.4.0)
#     generics               0.1.4     2025-05-09 [2] CRAN (R 4.4.3)
#     GenomeInfoDb         * 1.42.3    2025-01-27 [2] Bioconductor 3.20 (R 4.4.2)
#     GenomeInfoDbData       1.2.13    2024-10-01 [2] Bioconductor
#     GenomicAlignments      1.42.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
#     GenomicRanges        * 1.58.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
#     GetoptLong             1.0.5     2020-12-15 [2] CRAN (R 4.4.0)
#     ggbeeswarm             0.7.2     2023-04-29 [2] CRAN (R 4.4.0)
#     ggplot2              * 3.5.2     2025-04-09 [2] CRAN (R 4.4.3)
#     ggrepel                0.9.6     2024-09-07 [2] CRAN (R 4.4.1)
#     GlobalOptions          0.1.2     2020-06-10 [2] CRAN (R 4.4.0)
#     glue                   1.8.0     2024-09-30 [2] CRAN (R 4.4.1)
#     golem                  0.5.1     2024-08-27 [2] CRAN (R 4.4.1)
#     gridExtra              2.3       2017-09-09 [2] CRAN (R 4.4.0)
#     gtable                 0.3.6     2024-10-25 [2] CRAN (R 4.4.2)
#     here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.4.0)
#     hms                    1.1.3     2023-03-21 [2] CRAN (R 4.4.0)
#     htmltools              0.5.8.1   2024-04-04 [2] CRAN (R 4.4.0)
#     htmlwidgets            1.6.4     2023-12-06 [2] CRAN (R 4.4.0)
#     httpuv                 1.6.16    2025-04-16 [2] CRAN (R 4.4.3)
#     httr                   1.4.7     2023-08-15 [2] CRAN (R 4.4.0)
#     IRanges              * 2.40.1    2024-12-05 [2] Bioconductor 3.20 (R 4.4.2)
#     irlba                  2.3.5.1   2022-10-03 [2] CRAN (R 4.4.0)
#     iterators              1.0.14    2022-02-05 [2] CRAN (R 4.4.0)
#     jquerylib              0.1.4     2021-04-26 [2] CRAN (R 4.4.0)
#     jsonlite               2.0.0     2025-03-27 [2] CRAN (R 4.4.3)
#     KEGGREST               1.46.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
#     later                  1.4.2     2025-04-08 [2] CRAN (R 4.4.3)
#     lattice                0.22-6    2024-03-20 [3] CRAN (R 4.4.3)
#     lazyeval               0.2.2     2019-03-15 [2] CRAN (R 4.4.0)
#     lifecycle              1.0.4     2023-11-07 [2] CRAN (R 4.4.0)
#     limma                  3.62.2    2025-01-09 [2] Bioconductor 3.20 (R 4.4.2)
#     locfit                 1.5-9.12  2025-03-05 [2] CRAN (R 4.4.3)
#     lubridate            * 1.9.4     2024-12-08 [2] CRAN (R 4.4.2)
#     magick                 2.8.6     2025-03-23 [2] CRAN (R 4.4.3)
#     magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.4.0)
#     Matrix                 1.7-2     2025-01-23 [3] CRAN (R 4.4.3)
#     MatrixGenerics       * 1.18.1    2025-01-09 [2] Bioconductor 3.20 (R 4.4.2)
#     matrixStats          * 1.5.0     2025-01-07 [2] CRAN (R 4.4.2)
#     memoise                2.0.1     2021-11-26 [2] CRAN (R 4.4.0)
#     mime                   0.13      2025-03-17 [2] CRAN (R 4.4.3)
#     paletteer              1.6.0     2024-01-21 [2] CRAN (R 4.4.0)
#     pillar                 1.10.2    2025-04-05 [2] CRAN (R 4.4.3)
#     pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.4.0)
#     plotly                 4.10.4    2024-01-13 [2] CRAN (R 4.4.0)
#     png                    0.1-8     2022-11-29 [2] CRAN (R 4.4.0)
#     promises               1.3.2     2024-11-28 [2] CRAN (R 4.4.2)
#     purrr                * 1.0.4     2025-02-05 [2] CRAN (R 4.4.2)
#     R6                     2.6.1     2025-02-15 [2] CRAN (R 4.4.2)
#     rappdirs               0.3.3     2021-01-31 [2] CRAN (R 4.4.0)
#     RColorBrewer           1.1-3     2022-04-03 [2] CRAN (R 4.4.0)
#     Rcpp                   1.0.14    2025-01-12 [2] CRAN (R 4.4.2)
#     RCurl                  1.98-1.17 2025-03-22 [2] CRAN (R 4.4.3)
#     readr                * 2.1.5     2024-01-10 [2] CRAN (R 4.4.0)
#     rematch2               2.1.2     2020-05-01 [2] CRAN (R 4.4.0)
#     restfulr               0.0.15    2022-06-16 [2] CRAN (R 4.4.0)
#     rjson                  0.2.23    2024-09-16 [2] CRAN (R 4.4.1)
#     rlang                  1.1.6     2025-04-11 [2] CRAN (R 4.4.3)
#     rprojroot              2.0.4     2023-11-05 [2] CRAN (R 4.4.0)
#     Rsamtools              2.22.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
#     RSQLite                2.3.11    2025-05-04 [2] CRAN (R 4.4.3)
#     rsvd                   1.0.5     2021-04-16 [2] CRAN (R 4.4.0)
#     rtracklayer            1.66.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
#     S4Arrays               1.6.0     2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
#     S4Vectors            * 0.44.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
#     sass                   0.4.10    2025-04-11 [2] CRAN (R 4.4.3)
#     ScaledMatrix           1.14.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
#     scales                 1.4.0     2025-04-24 [2] CRAN (R 4.4.3)
#     scater                 1.34.1    2025-03-03 [2] Bioconductor 3.20 (R 4.4.3)
#     scuttle                1.16.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
#     sessioninfo          * 1.2.3     2025-02-05 [2] CRAN (R 4.4.2)
#     shape                  1.4.6.1   2024-02-23 [2] CRAN (R 4.4.0)
#     shiny                  1.10.0    2024-12-14 [2] CRAN (R 4.4.2)
#     shinyWidgets           0.9.0     2025-02-21 [2] CRAN (R 4.4.3)
#     SingleCellExperiment * 1.28.1    2024-11-10 [2] Bioconductor 3.20 (R 4.4.2)
#     SparseArray            1.6.2     2025-02-20 [2] Bioconductor 3.20 (R 4.4.3)
#     SpatialExperiment    * 1.16.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
#     spatialLIBD          * 1.21.5    2025-05-16 [1] Github (LieberInstitute/spatialLIBD@aff00db)
#     statmod                1.5.0     2023-01-06 [2] CRAN (R 4.4.0)
#     stringi                1.8.7     2025-03-27 [2] CRAN (R 4.4.3)
#     stringr              * 1.5.1     2023-11-14 [2] CRAN (R 4.4.0)
#     SummarizedExperiment * 1.36.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
#     tibble               * 3.2.1     2023-03-20 [2] CRAN (R 4.4.0)
#     tidyr                * 1.3.1     2024-01-24 [2] CRAN (R 4.4.0)
#     tidyselect             1.2.1     2024-03-11 [2] CRAN (R 4.4.0)
#     tidyverse            * 2.0.0     2023-02-22 [2] CRAN (R 4.4.0)
#     timechange             0.3.0     2024-01-18 [2] CRAN (R 4.4.0)
#     tzdb                   0.5.0     2025-03-15 [2] CRAN (R 4.4.3)
#     UCSC.utils             1.2.0     2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
#     vctrs                  0.6.5     2023-12-01 [2] CRAN (R 4.4.0)
#     vipor                  0.4.7     2023-12-18 [2] CRAN (R 4.4.0)
#     viridis                0.6.5     2024-01-29 [2] CRAN (R 4.4.0)
