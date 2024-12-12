########################################################################
## Compares pre-selected Seurats including WNN clusters with different methods, resolutions and knn-connectivity settings
## INPUT:
##      (1) First Seurat with WNN to compare
##      (2) Second Seurat with WNN to compare
##
## OUPUT:
##      1) Jaccard Index Heatmap
##
## Authors. CSC 
## Date. Dec 11, 2024
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
########################################################################

library("Seurat")
library("Signac")
library("scran")
library("bluster")
library("SingleCellExperiment")
#library("spatialLIBD")
library("dplyr")
library("bluster")
library("ComplexHeatmap")
library("ggplot2")
library("here")
library("viridisLite")

## input directories

# Check/create directories
inputRDS_Dir <- here("processed-data", "05_Clustering_ARCr", "01_clustering_std_method")
outputCVS_Dir <- here("processed-data", "05_Clustering_ARCr")
plotDir <- here("plots", "05_Clustering_ARCr", "03_jaccard")

## Check directories
if (!dir.exists(plotDir)) {dir.create(plotDir)}
if (!dir.exists(outputCVS_Dir)) {dir.create(outputCVS_Dir)}

## Load the Seurats with WNN clusters

message("Reading Seurat(s) to compute Jaccard Index on WNN clusters")

## Main arguments to specify which WNN clustering results to compare

Seurat_base_name = "seurat.norm_counts_Harmony_ARCr_QCed"
# knn = 20
knn = 30
# knn = 40
# methodWNN = "C.leiden_lsi_r"
# methodWNN = "C.SLM_lsi_r"
methodWNN = "C.louvain_lsi_r"
# methodWNN = "C.louvainM_lsi_r"
# res = 0.8 
resolution = 1
# res = 1.5
# res = 2


Seurat_base_name_1 = paste0(Seurat_base_name, "_WNN_k", knn, "_", methodWNN, resolution)
seurat_RDSname_1 = paste0(Seurat_base_name_1, ".rds")
if ( !length(list.files(inputRDS_Dir, pattern = seurat_RDSname_1)==1) ) { message("Seurat object missed!");  stop() }


## Load FIRST Seurat with WNN clustering 

seurat_RDSname_1 <- here(inputRDS_Dir, seurat_RDSname_1)
SeuratOBJ_1 <- readRDS(seurat_RDSname_1)
# length(Cells(x = SeuratOBJ))
n_clust <- nrow(unique(SeuratOBJ_1[["seurat_clusters"]]))
message("\nFirst Seurat with WNN loaded: `", Seurat_base_name_1, "`")
message(n_clust ," clusters")
table(SeuratOBJ_1[["seurat_clusters"]])
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 5519 3930 3647 3112 2653 2277 2272 2265 2130 2100 2049 2044 1739 1697 1640 1614 
# 16   17   18   19   20   21   22   23   24   25   26   27   28   29   30   31 
# 1518 1500 1497 1392 1268 1165  990  833  833  702  616  559  466  292  284  238 
# 32   33   34   35   36   37 
# 236  202  194  104   77   48


## Prepare arguments for SECOND WNN data clustering 

# Seurat_base_name = "seurat.norm_counts_Harmony_ARCr_QCed"
methodWNN = "C.leiden_lsi_r"

Seurat_base_name_2 = paste0(Seurat_base_name, "_WNN_k", knn, "_", methodWNN, resolution)
seurat_RDSname_2 = paste0(Seurat_base_name_2, ".rds")
if ( !length(list.files(inputRDS_Dir, pattern = seurat_RDSname_2)==1) ) { message("Seurat 2 object missed!");  stop() }

## Load SECOND Seurat with WNN clustering 

seurat_RDSname_2 <- here(inputRDS_Dir, seurat_RDSname_2)
SeuratOBJ_2 <- readRDS(seurat_RDSname_2)
# length(Cells(x = SeuratOBJ))
n_clust <- nrow(unique(SeuratOBJ_2[["seurat_clusters"]]))
message("\nSecond Seurat with WNN loaded: `", Seurat_base_name_2, "`")
message(n_clust ," clusters")
table(SeuratOBJ_2[["seurat_clusters"]])
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 5516 5240 3041 3023 2378 2330 2295 2261 2123 2077 2041 1971 1809 1758 1743 1700 
# 17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32 
# 1694 1641 1485 1396 1256 1193  968  892  831  706  615  556  450  238  202  195 
# 33 
# 78 

## Some fast checking

Reductions(SeuratOBJ_1) #umap.lovain
Reductions(SeuratOBJ_2) #umap.leiden
colnames(SeuratOBJ_1@meta.data)
colnames(SeuratOBJ_2@meta.data)
#tail(SeuratOBJ_1[["seurat_clusters"]], n=3)
tail(SeuratOBJ_1[["wsnn_res.1"]], n=3)
#tail(SeuratOBJ_2[["seurat_clusters"]], n=3)
tail(SeuratOBJ_2[["wsnn_res.1"]], n=3)

message("Starting Jaccard Index Processing ...")


## Plot approximate silhouette for evaluating cluster separation
## Identified and save closest neighboring cluster for each cell in each cluster 

plot_approxSilhouette <- function(sce, name_reduction, name_method, re, k){
  
  sil.approx <- approxSilhouette(reducedDim(sce, name_reduction), clusters=colData(sce)$seurat_clusters)
  #sil.approx <- approxSilhouette(reducedDim(sce.pbmc, "PCA"), clusters=colLabels(sce.pbmc))
  sil.approx
  sil.data <- as.data.frame(sil.approx)
  sil.data$closest <- factor(ifelse(sil.data$width > 0, colData(sce)$seurat_clusters, sil.data$other))
  sil.data$cluster <- colData(sce)$seurat_clusters
  
  ## identified the closest neighboring cluster for each cell in each cluster
  tbl_aprox_sil <- table(Cluster=colData(sce)$seurat_clusters, sil.data$closest)
  cvs_file <- paste0("Silhouette_closest_neighboring_cluster_tbl_", name_method, "_r", re, "_knn", k, ".cvs")
  cvs_file <- here(outputCVS_Dir, cvs_file)
  write.csv(tbl_aprox_sil, cvs_file)
  
  plt1 <- ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) +
    ggbeeswarm::geom_quasirandom(method="smiley") + labs(title = paste0(name_method, " at resolution = ", re, " with k.nn=", k))
    #+ labs(subtitle = "CellRangerARC-reanalyze Human Hb")
  ggsave(plt1, filename = here(plotDir, paste0("Silhouette_", name_method, "_r", re, "_knn", k, ".png")), height = 6, width = 10)
  
  return(plt1)
}

sce.1 <- as.SingleCellExperiment(SeuratOBJ_1, assay = "RNA")
reducedDimNames(sce.1)
# [1] "PCA"                "UMAP.UNINTEGRATED"  "INTEGRATED.CCA"    
# [4] "UMAP"               "INTEGRATED.HARMONY" "UMAP.LOVAIN"       
# [7] "LSI"                "UMAP.ATAC"          "WNN.UMAP"
spe1.red_name <-  "UMAP.LOVAIN"
rm("SeuratOBJ_1")

sce.2 <- as.SingleCellExperiment(SeuratOBJ_2, assay = "RNA")
reducedDimNames(sce.2)
spe2.red_name <-  "UMAP.LEIDEN"
rm("SeuratOBJ_2")

## Prepare Silhoutte plot
plt_lov <- plot_approxSilhouette(sce.1, spe1.red_name, substring(spe1.red_name, 6, nchar(spe1.red_name)), resolution, knn)
plt_leid <- plot_approxSilhouette(sce.2, spe2.red_name, substring(spe2.red_name, 6, nchar(spe2.red_name)), resolution, knn)
plt1 <- plt_lov / plt_leid
tmp_name <- paste0(substring(spe1.red_name, 6, nchar(spe1.red_name)), "_", substring(spe2.red_name, 6, nchar(spe2.red_name)))
tmp_name <- paste0("Silhouette_", tmp_name, "_r", resolution, "_knn", knn, ".png")
ggsave(plt1, filename = here(plotDir, tmp_name), height = 12, width = 10)


## Comparing different clusterings

# Sys.time()
# load(
#     here(
#         "processed-data",
#         "rdata",
#         "spe",
#         "01_build_spe",
#         "spe_filtered_final_with_clusters.Rdata"
#     ),
#     verbose = TRUE
# )
# Sys.time()

# ## Read the layers
# bayes_layers <-
#     get(load(
#         here(
#             "processed-data",
#             "rdata",
#             "spe",
#             "08_spatial_registration",
#             "bayesSpace_layer_annotations.Rdata"
#         )
#     )) |>
#     select(Annotation = bayesSpace, layer_long = cluster, layer_combo) |>
#     filter(Annotation %in% c("k09", "k16", "k28"))


## Compute the jaccard matrices, just like at
## https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/4b94e5bf1986df546bdb8624769e2ab746c23e70/code/05_explore_sce/06_explore_azimuth_annotations.R#L109

SeuratOBJ_1@meta.data

jacc.mat <-
    with(
        colData(spe),
        linkClustersMatrix(bayesSpace_harmony_9, bayesSpace_harmony_16)
    )

## Mark 0s as NAs
jacc.mat[jacc.mat == 0] <- NA

Sp09_order <- bayes_layers |>
    filter(Annotation == "k09") |>
    select(layer_combo) |>
    (`[[`)("layer_combo") |>
    as.character()

Sp16_order <- bayes_layers |>
    filter(Annotation == "k16") |>
    select(layer_combo) |>
    (`[[`)("layer_combo") |>
    as.character()

rownames(jacc.mat) <- sort(Sp09_order)
colnames(jacc.mat) <- sort(Sp16_order)


domain_colors_k09 <-
    setNames(
        Polychrome::palette36.colors(9),
        rownames(jacc.mat)
    )
row_ha <- rowAnnotation(
    df = data.frame(Sp09 = Sp09_order),
    col = list(Sp09 = domain_colors_k09[Sp09_order]),
    show_legend = c(FALSE)
)
domain_colors_k16 <-
    setNames(
        Polychrome::palette36.colors(16),
        colnames(jacc.mat)
    )
col_ha <-
    HeatmapAnnotation(
        df = data.frame(Sp16 = Sp16_order),
        col = list(Sp16 = domain_colors_k16[Sp16_order]),
        annotation_name_side = "left",
        show_legend = c(FALSE)
    )

pdf(here(plot_dir, "Sp09_vs_Sp16_complex.pdf"),
    height = 8,
    width = 12
)
Heatmap(
    jacc.mat[Sp09_order, ][, Sp16_order],
    name = "Correspondence",
    right_annotation = row_ha,
    bottom_annotation = col_ha,
    col = viridisLite::plasma(101),
    na_col = "black",
    cluster_rows = FALSE,
    cluster_columns = FALSE
)
dev.off()






## Repeat but for k09 vs k28
jacc.mat <-
    with(
        colData(spe),
        linkClustersMatrix(bayesSpace_harmony_9, bayesSpace_harmony_28)
    )

## Mark 0s as NAs
jacc.mat[jacc.mat == 0] <- NA

## Note that not all k28 domains are present: 18 and 21 are missing
present_col <- as.integer(colnames(jacc.mat))
Sp28_order <- bayes_layers |>
    filter(Annotation == "k28") |>
    select(layer_combo) |>
    (`[[`)("layer_combo") |>
    as.character()


rownames(jacc.mat) <- sort(Sp09_order)
colnames(jacc.mat) <- sort(Sp28_order)

domain_colors_k28 <-
    setNames(
        Polychrome::palette36.colors(28)[present_col],
        colnames(jacc.mat)
    )
col_ha <-
    HeatmapAnnotation(
        df = data.frame(Sp28 = Sp28_order),
        col = list(Sp28 = domain_colors_k28[Sp28_order]),
        annotation_name_side = "left",
        show_legend = c(FALSE)
    )

pdf(here(plot_dir, "Sp09_vs_Sp28_complex.pdf"),
    height = 8,
    width = 12
)
Heatmap(
    jacc.mat[Sp09_order, ][, Sp28_order],
    name = "Correspondence",
    right_annotation = row_ha,
    bottom_annotation = col_ha,
    col = viridisLite::plasma(101),
    na_col = "black",
    cluster_rows = FALSE,
    cluster_columns =
    )
dev.off()

## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.2.2 Patched (2022-12-14 r83473)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2022-12-19
#  pandoc   2.19.2 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/bin/pandoc
#
# ─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package                * version   date (UTC) lib source
#  AnnotationDbi            1.60.0    2022-11-01 [2] Bioconductor
#  AnnotationHub            3.6.0     2022-11-01 [2] Bioconductor
#  assertthat               0.2.1     2019-03-21 [2] CRAN (R 4.2.1)
#  attempt                  0.3.1     2020-05-03 [2] CRAN (R 4.2.1)
#  beachmat                 2.14.0    2022-11-01 [2] Bioconductor
#  beeswarm                 0.4.0     2021-06-01 [2] CRAN (R 4.2.1)
#  benchmarkme              1.0.8     2022-06-12 [2] CRAN (R 4.2.1)
#  benchmarkmeData          1.0.4     2020-04-23 [2] CRAN (R 4.2.1)
#  Biobase                * 2.58.0    2022-11-01 [1] Bioconductor
#  BiocFileCache            2.6.0     2022-11-01 [1] Bioconductor
#  BiocGenerics           * 0.44.0    2022-11-01 [2] Bioconductor
#  BiocIO                   1.8.0     2022-11-01 [2] Bioconductor
#  BiocManager              1.30.19   2022-10-25 [2] CRAN (R 4.2.2)
#  BiocNeighbors            1.16.0    2022-11-01 [2] Bioconductor
#  BiocParallel             1.32.4    2022-12-01 [1] Bioconductor
#  BiocSingular             1.14.0    2022-11-01 [2] Bioconductor
#  BiocVersion              3.16.0    2022-04-26 [2] Bioconductor
#  Biostrings               2.66.0    2022-11-01 [2] Bioconductor
#  bit                      4.0.5     2022-11-15 [2] CRAN (R 4.2.2)
#  bit64                    4.0.5     2020-08-30 [2] CRAN (R 4.2.1)
#  bitops                   1.0-7     2021-04-24 [2] CRAN (R 4.2.1)
#  blob                     1.2.3     2022-04-10 [2] CRAN (R 4.2.1)
#  bluster                * 1.8.0     2022-11-01 [2] Bioconductor
#  bslib                    0.4.2     2022-12-16 [2] CRAN (R 4.2.2)
#  cachem                   1.0.6     2021-08-19 [2] CRAN (R 4.2.1)
#  Cairo                    1.6-0     2022-07-05 [2] CRAN (R 4.2.1)
#  circlize                 0.4.15    2022-05-10 [2] CRAN (R 4.2.1)
#  cli                      3.4.1     2022-09-23 [2] CRAN (R 4.2.1)
#  clue                     0.3-63    2022-11-19 [2] CRAN (R 4.2.2)
#  cluster                  2.1.4     2022-08-22 [3] CRAN (R 4.2.2)
#  codetools                0.2-18    2020-11-04 [3] CRAN (R 4.2.2)
#  colorout                 1.2-2     2022-11-02 [1] Github (jalvesaq/colorout@79931fd)
#  colorspace               2.0-3     2022-02-21 [2] CRAN (R 4.2.1)
#  ComplexHeatmap         * 2.14.0    2022-11-01 [2] Bioconductor
#  config                   0.3.1     2020-12-17 [2] CRAN (R 4.2.1)
#  cowplot                  1.1.1     2020-12-30 [2] CRAN (R 4.2.1)
#  crayon                   1.5.2     2022-09-29 [2] CRAN (R 4.2.1)
#  curl                     4.3.3     2022-10-06 [2] CRAN (R 4.2.1)
#  data.table               1.14.6    2022-11-16 [2] CRAN (R 4.2.2)
#  DBI                      1.1.3     2022-06-18 [2] CRAN (R 4.2.1)
#  dbplyr                   2.2.1     2022-06-27 [2] CRAN (R 4.2.1)
#  DelayedArray             0.24.0    2022-11-01 [2] Bioconductor
#  DelayedMatrixStats       1.20.0    2022-11-01 [1] Bioconductor
#  desc                     1.4.2     2022-09-08 [2] CRAN (R 4.2.1)
#  digest                   0.6.31    2022-12-11 [2] CRAN (R 4.2.2)
#  doParallel               1.0.17    2022-02-07 [2] CRAN (R 4.2.1)
#  dotCall64                1.0-2     2022-10-03 [2] CRAN (R 4.2.1)
#  dplyr                  * 1.0.10    2022-09-01 [2] CRAN (R 4.2.1)
#  dqrng                    0.3.0     2021-05-01 [2] CRAN (R 4.2.1)
#  DropletUtils             1.18.1    2022-11-22 [2] Bioconductor
#  DT                       0.26      2022-10-19 [2] CRAN (R 4.2.1)
#  edgeR                    3.40.1    2022-12-14 [2] Bioconductor
#  ellipsis                 0.3.2     2021-04-29 [2] CRAN (R 4.2.1)
#  ExperimentHub            2.6.0     2022-11-01 [2] Bioconductor
#  fansi                    1.0.3     2022-03-24 [2] CRAN (R 4.2.1)
#  fastmap                  1.1.0     2021-01-25 [2] CRAN (R 4.2.1)
#  fields                   14.1      2022-08-12 [2] CRAN (R 4.2.1)
#  filelock                 1.0.2     2018-10-05 [2] CRAN (R 4.2.1)
#  foreach                  1.5.2     2022-02-02 [2] CRAN (R 4.2.1)
#  fs                       1.5.2     2021-12-08 [2] CRAN (R 4.2.1)
#  generics                 0.1.3     2022-07-05 [2] CRAN (R 4.2.1)
#  GenomeInfoDb           * 1.34.4    2022-12-01 [1] Bioconductor
#  GenomeInfoDbData         1.2.9     2022-09-29 [2] Bioconductor
#  GenomicAlignments        1.34.0    2022-11-01 [2] Bioconductor
#  GenomicRanges          * 1.50.2    2022-12-16 [2] Bioconductor
#  GetoptLong               1.0.5     2020-12-15 [2] CRAN (R 4.2.1)
#  ggbeeswarm               0.7.1     2022-12-16 [2] CRAN (R 4.2.2)
#  ggplot2                  3.4.0     2022-11-04 [2] CRAN (R 4.2.2)
#  ggrepel                  0.9.2     2022-11-06 [2] CRAN (R 4.2.2)
#  GlobalOptions            0.1.2     2020-06-10 [2] CRAN (R 4.2.1)
#  glue                     1.6.2     2022-02-24 [2] CRAN (R 4.2.1)
#  golem                    0.3.5     2022-10-18 [2] CRAN (R 4.2.1)
#  gridExtra                2.3       2017-09-09 [2] CRAN (R 4.2.1)
#  gtable                   0.3.1     2022-09-01 [2] CRAN (R 4.2.1)
#  HDF5Array                1.26.0    2022-11-01 [2] Bioconductor
#  here                   * 1.0.1     2020-12-13 [2] CRAN (R 4.2.1)
#  htmltools                0.5.4     2022-12-07 [2] CRAN (R 4.2.2)
#  htmlwidgets              1.6.0     2022-12-15 [2] CRAN (R 4.2.2)
#  httpuv                   1.6.7     2022-12-14 [2] CRAN (R 4.2.2)
#  httr                     1.4.4     2022-08-17 [2] CRAN (R 4.2.1)
#  igraph                   1.3.5     2022-09-22 [2] CRAN (R 4.2.1)
#  interactiveDisplayBase   1.36.0    2022-11-01 [2] Bioconductor
#  IRanges                * 2.32.0    2022-11-01 [2] Bioconductor
#  irlba                    2.3.5.1   2022-10-03 [2] CRAN (R 4.2.1)
#  iterators                1.0.14    2022-02-05 [2] CRAN (R 4.2.1)
#  jquerylib                0.1.4     2021-04-26 [2] CRAN (R 4.2.1)
#  jsonlite                 1.8.4     2022-12-06 [2] CRAN (R 4.2.2)
#  KEGGREST                 1.38.0    2022-11-01 [2] Bioconductor
#  knitr                    1.41      2022-11-18 [2] CRAN (R 4.2.2)
#  later                    1.3.0     2021-08-18 [2] CRAN (R 4.2.1)
#  lattice                  0.20-45   2021-09-22 [3] CRAN (R 4.2.2)
#  lazyeval                 0.2.2     2019-03-15 [2] CRAN (R 4.2.1)
#  lifecycle                1.0.3     2022-10-07 [2] CRAN (R 4.2.1)
#  limma                    3.54.0    2022-11-01 [1] Bioconductor
#  locfit                   1.5-9.6   2022-07-11 [2] CRAN (R 4.2.1)
#  magick                   2.7.3     2021-08-18 [2] CRAN (R 4.2.1)
#  magrittr                 2.0.3     2022-03-30 [2] CRAN (R 4.2.1)
#  maps                     3.4.1     2022-10-30 [2] CRAN (R 4.2.2)
#  Matrix                   1.5-3     2022-11-11 [2] CRAN (R 4.2.2)
#  MatrixGenerics         * 1.10.0    2022-11-01 [1] Bioconductor
#  matrixStats            * 0.63.0    2022-11-18 [2] CRAN (R 4.2.2)
#  memoise                  2.0.1     2021-11-26 [2] CRAN (R 4.2.1)
#  mime                     0.12      2021-09-28 [2] CRAN (R 4.2.1)
#  munsell                  0.5.0     2018-06-12 [2] CRAN (R 4.2.1)
#  paletteer                1.5.0     2022-10-19 [2] CRAN (R 4.2.1)
#  pillar                   1.8.1     2022-08-19 [2] CRAN (R 4.2.1)
#  pkgconfig                2.0.3     2019-09-22 [2] CRAN (R 4.2.1)
#  pkgload                  1.3.2     2022-11-16 [2] CRAN (R 4.2.2)
#  plotly                   4.10.1    2022-11-07 [2] CRAN (R 4.2.2)
#  png                      0.1-8     2022-11-29 [2] CRAN (R 4.2.2)
#  Polychrome               1.5.1     2022-05-03 [1] CRAN (R 4.2.2)
#  promises                 1.2.0.1   2021-02-11 [2] CRAN (R 4.2.1)
#  purrr                    0.3.5     2022-10-06 [2] CRAN (R 4.2.1)
#  R.methodsS3              1.8.2     2022-06-13 [2] CRAN (R 4.2.1)
#  R.oo                     1.25.0    2022-06-12 [2] CRAN (R 4.2.1)
#  R.utils                  2.12.2    2022-11-11 [2] CRAN (R 4.2.2)
#  R6                       2.5.1     2021-08-19 [2] CRAN (R 4.2.1)
#  rappdirs                 0.3.3     2021-01-31 [2] CRAN (R 4.2.1)
#  RColorBrewer             1.1-3     2022-04-03 [2] CRAN (R 4.2.1)
#  Rcpp                     1.0.9     2022-07-08 [2] CRAN (R 4.2.1)
#  RCurl                    1.98-1.9  2022-10-03 [2] CRAN (R 4.2.1)
#  rematch2                 2.1.2     2020-05-01 [2] CRAN (R 4.2.1)
#  restfulr                 0.0.15    2022-06-16 [2] CRAN (R 4.2.1)
#  rhdf5                    2.42.0    2022-11-01 [2] Bioconductor
#  rhdf5filters             1.10.0    2022-11-01 [2] Bioconductor
#  Rhdf5lib                 1.20.0    2022-11-01 [2] Bioconductor
#  rjson                    0.2.21    2022-01-09 [2] CRAN (R 4.2.1)
#  rlang                    1.0.6     2022-09-24 [2] CRAN (R 4.2.1)
#  rmote                    0.3.4     2022-11-02 [1] Github (cloudyr/rmote@fbce611)
#  roxygen2                 7.2.3     2022-12-08 [2] CRAN (R 4.2.2)
#  rprojroot                2.0.3     2022-04-02 [2] CRAN (R 4.2.1)
#  Rsamtools                2.14.0    2022-11-01 [2] Bioconductor
#  RSQLite                  2.2.19    2022-11-24 [2] CRAN (R 4.2.2)
#  rstudioapi               0.14      2022-08-22 [2] CRAN (R 4.2.1)
#  rsvd                     1.0.5     2021-04-16 [2] CRAN (R 4.2.1)
#  rtracklayer              1.58.0    2022-11-01 [2] Bioconductor
#  S4Vectors              * 0.36.1    2022-12-05 [1] Bioconductor
#  sass                     0.4.4     2022-11-24 [2] CRAN (R 4.2.2)
#  ScaledMatrix             1.6.0     2022-11-01 [1] Bioconductor
#  scales                   1.2.1     2022-08-20 [2] CRAN (R 4.2.1)
#  scater                   1.26.1    2022-11-13 [2] Bioconductor
#  scatterplot3d            0.3-42    2022-09-08 [1] CRAN (R 4.2.2)
#  scuttle                  1.8.3     2022-12-14 [2] Bioconductor
#  servr                    0.25      2022-11-04 [1] CRAN (R 4.2.2)
#  sessioninfo            * 1.2.2     2021-12-06 [2] CRAN (R 4.2.1)
#  shape                    1.4.6     2021-05-19 [2] CRAN (R 4.2.1)
#  shiny                    1.7.4     2022-12-15 [2] CRAN (R 4.2.2)
#  shinyWidgets             0.7.5     2022-11-17 [2] CRAN (R 4.2.2)
#  SingleCellExperiment   * 1.20.0    2022-11-01 [1] Bioconductor
#  spam                     2.9-1     2022-08-07 [2] CRAN (R 4.2.1)
#  sparseMatrixStats        1.10.0    2022-11-01 [2] Bioconductor
#  SpatialExperiment      * 1.8.0     2022-11-01 [2] Bioconductor
#  spatialLIBD            * 1.11.4    2022-12-17 [1] Github (LieberInstitute/spatialLIBD@1aecde8)
#  statmod                  1.4.37    2022-08-12 [2] CRAN (R 4.2.1)
#  stringi                  1.7.8     2022-07-11 [2] CRAN (R 4.2.1)
#  stringr                  1.5.0     2022-12-02 [2] CRAN (R 4.2.2)
#  SummarizedExperiment   * 1.28.0    2022-11-01 [2] Bioconductor
#  tibble                   3.1.8     2022-07-22 [2] CRAN (R 4.2.1)
#  tidyr                    1.2.1     2022-09-08 [2] CRAN (R 4.2.1)
#  tidyselect               1.2.0     2022-10-10 [2] CRAN (R 4.2.1)
#  usethis                  2.1.6     2022-05-25 [2] CRAN (R 4.2.1)
#  utf8                     1.2.2     2021-07-24 [2] CRAN (R 4.2.1)
#  vctrs                    0.5.1     2022-11-16 [2] CRAN (R 4.2.2)
#  vipor                    0.4.5     2017-03-22 [2] CRAN (R 4.2.1)
#  viridis                  0.6.2     2021-10-13 [2] CRAN (R 4.2.1)
#  viridisLite              0.4.1     2022-08-22 [2] CRAN (R 4.2.1)
#  withr                    2.5.0     2022-03-03 [2] CRAN (R 4.2.1)
#  xfun                     0.35      2022-11-16 [2] CRAN (R 4.2.2)
#  XML                      3.99-0.13 2022-12-04 [2] CRAN (R 4.2.2)
#  xml2                     1.3.3     2021-11-30 [2] CRAN (R 4.2.1)
#  xtable                   1.8-4     2019-04-21 [2] CRAN (R 4.2.1)
#  XVector                  0.38.0    2022-11-01 [2] Bioconductor
#  yaml                     2.3.6     2022-10-18 [2] CRAN (R 4.2.1)
#  zlibbioc                 1.44.0    2022-11-01 [1] Bioconductor
#
#  [1] /users/lcollado/R/4.2.x
#  [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/site-library
#  [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/library
#
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
