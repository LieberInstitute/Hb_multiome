
########################################################################
## Build JackStrawPlot for comparing the distribution of p-values for each PCA, CCA and Harmony integration
## INPUT: Seurat object with PCA, CCA and Harmony slots
## 
## OUPUT: JackStrawPlot plots
## 
## Date: April 9th, 2024 
########################################################################

# load libraries
library('Seurat')        
library('purrr')
library("here")
library("ggplot2")

here::here()



############        Initials      ############

count_mtx_type <-'data_counts'

## read directory with Seurat objects
if (count_mtx_type=='data_counts') { 
  Seurat_base_name <- 'seurat.combined.data_counts'
} else { 
  Seurat_base_name <- 'seurat.combined.norm_counts' 
}

## Select seurat objects with PCA, CCA and Harmony data
rds_namePCA <- here('processed-data/02_merge_seurats', paste0(Seurat_base_name, '_PCA.rds'))
rds_nameCCA <- here('processed-data/02_merge_seurats', paste0(Seurat_base_name, '_PCA_CCA.rds'))
rds_nameHarm <- here('processed-data/02_merge_seurats', paste0(Seurat_base_name, '_PCA_Harmony.rds'))

lst_seurats = list()
lst_sub_titles = list()

lst_seurats <- append(lst_seurats, rds_namePCA)
lst_sub_titles <- append(lst_sub_titles, 'PCA')

if (!is_empty(rds_nameCCA)) { 
  lst_seurats <- append(lst_seurats, rds_nameCCA)
  lst_sub_titles <- append(lst_sub_titles, 'CCA')
}
if (!is_empty(rds_nameHarm)) { 
  lst_seurats <- append(lst_seurats, rds_nameHarm) 
  lst_sub_titles <- append(lst_sub_titles, 'Harmony')
}

lst_seurats
lst_sub_titles

## Check if plots directory exists, if not create it

if (!dir.exists(here("plots/02_merge_seurats/"))) {
  dir.create(here("plots/02_merge_seurats/"))
}
rds_path_out <- here('plots/02_merge_seurats')




## Function to build JackStraw and Elbow plow before and after integration

plot_JackStraw <- function(g, t) {

  message('Loading seurat object ', g)
  
  #g <- lst_seurats[1]
  #t <- lst_sub_titles[1]
  
  SeuratOBJ <-readRDS(g) 
  
  if ( is.null(SeuratOBJ@reductions$pca) )  {
    
    print('Processing JackStraw scores ...')
    
    SeuratOBJ <- JackStraw(SeuratOBJ, dims = 20, prop.freq = 0.01, num.replicate = 100)
    SeuratOBJ <- ScoreJackStraw(SeuratOBJ, dims = 1:20)
    
    print('JackStraw scores calculated ... ')

    p1 <- JackStrawPlot(SeuratOBJ, dims = 1:20)  + ggtitle(label = Seurat_base_name, subtitle = t)
    p2 <- ElbowPlot(SeuratOBJ, ndims = 30)
    p3 <- p1 / p2
    
    print('JackStraw and elbow plots done!')
    
  }
  
}


## map the lists to build Violin plots for Hb, LHb and MHb

plot_list = list()

## test
#plot_list <- map2(lst_seurats[1], lst_sub_titles[1], 
#                  ~ plot_JackStraw(g = .x, t = .y))

if (length(lst_seurats) == length(lst_sub_titles)) {
  plot_list <- map2(lst_seurats, lst_sub_titles, 
                    ~ plot_JackStraw(g = .x, t = .y)) 
  }


print('JackStraw plots done')


dir_plt <- here('plots/02_merge_seurats/')
pdf_name <- file.path(dir_plt, paste0('JackStrawPlot_Elbow_',Seurat_base_name, '.pdf'))
pdf(file = pdf_name)
print(plot_list)
dev.off()


print('JackStraw plots saved!')


## slurm script reproducibility

# slurmjobs::job_loop(
#   loops = list(type_mtx = c("data_counts", "norm_counts")),
#   name = "04_JackStrawPlot",
#   cores = 2,
#   create_shell = TRUE
# )










############ Reproducibility information ####################

library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
proc.time()
options(width = 120)
session_info()



# > library("sessioninfo")
# > print('Reproducibility information:')
# [1] "Reproducibility information:"
# > # Last modification
#   > Sys.time()
# [1] "2024-04-10 13:21:31 EDT"
# > proc.time()
# user  system elapsed 
# 175.860  13.563 577.445 
# > options(width = 120)
# > session_info()
# 9m
# digest             0.6.34     2024-01-11 [2] CRAN (R 4.3.2)
# dotCall64          1.1-1      2023-11-28 [2] CRAN (R 4.3.2)
# dplyr              1.1.4      2023-11-17 [2] CRAN (R 4.3.2)
# ellipsis           0.3.2      2021-04-29 [2] CRAN (R 4.3.2)
# fansi              1.0.6      2023-12-08 [2] CRAN (R 4.3.2)
# fastDummies        1.7.3      2023-07-06 [2] CRAN (R 4.3.2)
# fastmap            1.1.1      2023-02-24 [2] CRAN (R 4.3.2)
# fitdistrplus       1.1-11     2023-04-25 [2] CRAN (R 4.3.2)
# future             1.33.1     2023-12-22 [2] CRAN (R 4.3.2)
# future.apply       1.11.1     2023-12-21 [2] CRAN (R 4.3.2)
# generics           0.1.3      2022-07-05 [2] CRAN (R 4.3.2)
# ggplot2          * 3.4.4      2023-10-12 [2] CRAN (R 4.3.2)
# ggrepel            0.9.5      2024-01-10 [2] CRAN (R 4.3.2)
# ggridges           0.5.6      2024-01-23 [2] CRAN (R 4.3.2)
# globals            0.16.2     2022-11-21 [2] CRAN (R 4.3.2)
# glue               1.7.0      2024-01-09 [2] CRAN (R 4.3.2)
# goftest            1.2-3      2021-10-07 [2] CRAN (R 4.3.2)
# gridExtra          2.3        2017-09-09 [2] CRAN (R 4.3.2)
# gtable             0.3.4      2023-08-21 [2] CRAN (R 4.3.2)
# here             * 1.0.1      2020-12-13 [2] CRAN (R 4.3.2)
# htmltools          0.5.7      2023-11-03 [2] CRAN (R 4.3.2)
# htmlwidgets        1.6.4      2023-12-06 [2] CRAN (R 4.3.2)
# httpuv             1.6.14     2024-01-26 [2] CRAN (R 4.3.2)
# httr               1.4.7      2023-08-15 [2] CRAN (R 4.3.2)
# ica                1.0-3      2022-07-08 [2] CRAN (R 4.3.2)
# igraph             2.0.1.9008 2024-02-09 [2] Github (igraph/rigraph@39158c6)
# irlba              2.3.5.1    2022-10-03 [2] CRAN (R 4.3.2)
# jsonlite           1.8.8      2023-12-04 [2] CRAN (R 4.3.2)
# KernSmooth         2.23-22    2023-07-10 [3] CRAN (R 4.3.2)
# later              1.3.2      2023-12-06 [2] CRAN (R 4.3.2)
# lattice            0.22-5     2023-10-24 [3] CRAN (R 4.3.2)
# lazyeval           0.2.2      2019-03-15 [2] CRAN (R 4.3.2)
# leiden             0.4.3.1    2023-11-17 [2] CRAN (R 4.3.2)
# lifecycle          1.0.4      2023-11-07 [2] CRAN (R 4.3.2)
# listenv            0.9.1      2024-01-29 [2] CRAN (R 4.3.2)
# lmtest             0.9-40     2022-03-21 [2] CRAN (R 4.3.2)
# lubridate          1.9.3      2023-09-27 [2] CRAN (R 4.3.2)
# magrittr           2.0.3      2022-03-30 [2] CRAN (R 4.3.2)
# MASS               7.3-60.0.1 2024-01-13 [3] CRAN (R 4.3.2)
# Matrix             1.6-5      2024-01-11 [3] CRAN (R 4.3.2)
# matrixStats        1.2.0      2023-12-11 [2] CRAN (R 4.3.2)
# mime               0.12       2021-09-28 [2] CRAN (R 4.3.2)
# miniUI             0.1.1.1    2018-05-18 [2] CRAN (R 4.3.2)
# munsell            0.5.0      2018-06-12 [2] CRAN (R 4.3.2)
# nlme               3.1-164    2023-11-27 [3] CRAN (R 4.3.2)
# parallelly         1.36.0     2023-05-26 [2] CRAN (R 4.3.2)
# patchwork          1.2.0      2024-01-08 [2] CRAN (R 4.3.2)
# pbapply            1.7-2      2023-06-27 [2] CRAN (R 4.3.2)
# pillar             1.9.0      2023-03-22 [2] CRAN (R 4.3.2)
# pkgconfig          2.0.3      2019-09-22 [2] CRAN (R 4.3.2)
# plotly             4.10.4     2024-01-13 [2] CRAN (R 4.3.2)
# plyr               1.8.9      2023-10-02 [2] CRAN (R 4.3.2)
# png                0.1-8      2022-11-29 [2] CRAN (R 4.3.2)
# polyclip           1.10-6     2023-09-27 [2] CRAN (R 4.3.2)
# progressr          0.14.0     2023-08-10 [2] CRAN (R 4.3.2)
# promises           1.2.1      2023-08-10 [2] CRAN (R 4.3.2)
# purrr            * 1.0.2      2023-08-10 [2] CRAN (R 4.3.2)
# R6                 2.5.1      2021-08-19 [2] CRAN (R 4.3.2)
# RANN               2.6.1      2019-01-08 [2] CRAN (R 4.3.2)
# RColorBrewer       1.1-3      2022-04-03 [2] CRAN (R 4.3.2)
# Rcpp               1.0.12     2024-01-09 [2] CRAN (R 4.3.2)
# RcppAnnoy          0.0.22     2024-01-23 [2] CRAN (R 4.3.2)
# RcppHNSW           0.6.0      2024-02-04 [2] CRAN (R 4.3.2)
# reshape2           1.4.4      2020-04-09 [2] CRAN (R 4.3.2)
# reticulate         1.35.0     2024-01-31 [2] CRAN (R 4.3.2)
# rlang              1.1.3      2024-01-10 [2] CRAN (R 4.3.2)
# ROCR               1.0-11     2020-05-02 [2] CRAN (R 4.3.2)
# rprojroot          2.0.4      2023-11-05 [2] CRAN (R 4.3.2)
# RSpectra           0.16-1     2022-04-24 [2] CRAN (R 4.3.2)
# Rtsne              0.17       2023-12-07 [2] CRAN (R 4.3.2)
# scales             1.3.0      2023-11-28 [2] CRAN (R 4.3.2)
# scattermore        1.2        2023-06-12 [2] CRAN (R 4.3.2)
# sctransform        0.4.1      2023-10-19 [2] CRAN (R 4.3.2)
# sessioninfo      * 1.2.2      2021-12-06 [2] CRAN (R 4.3.2)
# Seurat           * 5.0.1      2023-11-17 [2] CRAN (R 4.3.2)
# SeuratObject     * 5.0.1      2023-11-17 [2] CRAN (R 4.3.2)
# shiny              1.8.0      2023-11-17 [2] CRAN (R 4.3.2)
# slurmjobs          1.2.1      2024-03-20 [1] Github (LieberInstitute/slurmjobs@3fbddb2)
# sp               * 2.1-3      2024-01-30 [2] CRAN (R 4.3.2)
# spam               2.10-0     2023-10-23 [2] CRAN (R 4.3.2)
# spatstat.data      3.0-4      2024-01-15 [2] CRAN (R 4.3.2)
# spatstat.explore   3.2-6      2024-02-01 [2] CRAN (R 4.3.2)
# spatstat.geom      3.2-8      2024-01-26 [2] CRAN (R 4.3.2)
# spatstat.random    3.2-2      2023-11-29 [2] CRAN (R 4.3.2)
# spatstat.sparse    3.0-3      2023-10-24 [2] CRAN (R 4.3.2)
# spatstat.utils     3.0-4      2023-10-24 [2] CRAN (R 4.3.2)
# stringi            1.8.3      2023-12-11 [2] CRAN (R 4.3.2)
# stringr            1.5.1      2023-11-14 [2] CRAN (R 4.3.2)
# survival           3.5-7      2023-08-14 [3] CRAN (R 4.3.2)
# tensor             1.5        2012-05-05 [2] CRAN (R 4.3.2)
# tibble             3.2.1      2023-03-20 [2] CRAN (R 4.3.2)
# tidyr              1.3.1      2024-01-24 [2] CRAN (R 4.3.2)
# tidyselect         1.2.0      2022-10-10 [2] CRAN (R 4.3.2)
# timechange         0.3.0      2024-01-18 [2] CRAN (R 4.3.2)
# utf8               1.2.4      2023-10-22 [2] CRAN (R 4.3.2)
# uwot               0.1.16     2023-06-29 [2] CRAN (R 4.3.2)
# vctrs              0.6.5      2023-12-01 [2] CRAN (R 4.3.2)
# viridisLite        0.4.2      2023-05-02 [2] CRAN (R 4.3.2)
# withr              3.0.0      2024-01-16 [2] CRAN (R 4.3.2)
# xtable             1.8-4      2019-04-21 [2] CRAN (R 4.3.2)
# zoo                1.8-12     2023-04-13 [2] CRAN (R 4.3.2)
# 
# [1] /users/csoto/R/4.3.x
# [2] /jhpce/shared/community/core/conda_R/4.3.x/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.3.x/R/lib64/R/library