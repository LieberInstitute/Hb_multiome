########################################################################
## Find conserved genes across samples/conditions 
##         
## Authors. CSC/HT
## Implementation from https://hbctraining.github.io/scRNA-seq_online/lessons/09_merged_SC_marker_identification.html
## Date. Dic 5th, 2023 / last md. 
########################################################################

# load libraries
library(Seurat)                                 # 4.9.9.9045 2023-05-17 [1] Github (satijalab/seurat@7d1094c)
library(Signac)                                 # 1.9.0.9000 2023-05-08 [1] Github (stuart-lab/signac@cf31022)
library(dplyr)
library(purrr)

library(patchwork)
library(ggplot2)
library(ggpubr)

library(stringr)
library(here)


# To point to compatible versions
stopifnot(packageVersion("Seurat") >='4.9.9.9060')
stopifnot(packageVersion("Signac") >='1.11.9000')

# Paths in the JHPCE cluster
here::here()

# read directory
file_dir <- here("plots/06_Diff_expr_genes")


############  Initials ############

# Set default pcs as initial, but get custom PCs used to run UMAPs. Also this variable is used to run FindNeighbors() function (inciso 3)
umap_dims <- 1:30


###### 1. Read pre-existing Seurat Object ######

# Search the seurat pre-existing
sample_tmp <- commandArgs(trailingOnly = TRUE)
# Samples for testing
#sample_tmp <- '1_HPC_KDM,human'  # testing HUMAN tissue
#sample_tmp <- '3_HPC_KDM,mouse'  # testing MOUSE tissue 
#sample_tmp <- 'hippo42_1,human'    # testing human tissue -- smaller 615 cells

s_sample <- sample_tmp %>% str_split(",") %>%flatten_chr() %>% .[[1]]
s_tissue <- sample_tmp %>% str_split(",") %>%flatten_chr() %>% .[[2]]
message("Processing sample: ", s_sample, ' / ', s_tissue)

# Function to load Seurat object with dim reduction
load_seurat_object <- function(s_sample) {
  # Construct file pattern
  file_pattern <- paste0(s_sample, "_DE_gene_markers.rds")
  file_dir <- here("plots/06_Diff_expr_genes")
  # Read and return Seurat object
  message(paste0("Reading Seurat object for sample: ", s_sample))
  readRDS(file.path(file_path, file_pattern))
}

# Read command line arguments and load Seurat object
SeuratOBJ <- load_seurat_object(s_sample)

print(SeuratOBJ)


########  1. Integrate assays from different conditions

## seurat_integrated object

# View cluster levels
head(Idents(seurat_integrated))
# Get number of clusters to process 
cluster_idx <- length(levels(seurat_integrated@active.ident))

########  2. Find conserved markers from RNA Assay

DefaultAssay(seurat_integrated) <- "RNA"


FindConservedMarkers(seurat_integrated,
                     ident.1 = cluster,
                     grouping.var = "sample",
                     only.pos = TRUE,
                     min.diff.pct = 0.25,
                     min.pct = 0.25,
                     logfc.threshold = 0.25)

## test it out on one cluster to see how it works:
cluster0_conserved_markers <- FindConservedMarkers(seurat_integrated,
                                                   ident.1 = 0,
                                                   grouping.var = "sample",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)


## Adding Gene Annotations

annotations <- read.csv("data/annotation.csv")

# Combine markers with gene descriptions 
cluster0_ann_markers <- cluster0_conserved_markers %>% 
    rownames_to_column(var="gene") %>% 
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name"))

View(cluster0_ann_markers)


########### Save the Seurat with cell-markers

file_path <- "/processed-data/06_Diff_expr_genes"
file_name <-  paste0(s_sample,'_DE_gene_markers.rds')
filex <- here(file_path, file_name)
saveRDS(SeuratOBJ, file = filex)


message("Tasks completed successfully!")



################################################################################


library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
#"2023-04-04 12:42:26 EDT"
proc.time()
options(width = 120)
session_info()

