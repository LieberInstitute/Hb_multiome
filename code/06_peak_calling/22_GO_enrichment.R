########################################################################
## GO enrichment in cCRE and classified Linked OCR and DARs unlinked
##
## Authors. CSC
## Date. Oct 24, 2025
## Recommended resources on interactive mode: srun --pty --mem=20GB --x11 bash
########################################################################

#library("pheatmap")
#library("reshape2")
library("dplyr")
library("purrr")
library("ggplot2")
library("patchwork")
#library("ggrepel")
library("tidyverse")
library("tidyr")
library("stringr")
library("here")

library("getopt")
library("org.Hs.eg.db")
library("clusterProfiler")
library("rrvgo")
library("ComplexHeatmap")

#===============================================================================
# resolution_level = "Mid"      # 18 cell-types
#===============================================================================


## Set directory names
inputCSV_cCRE_ORC_Dir <- here(
    "processed-data",
    "06_peak_calling",
    "21_overlaping_FDRscores_TopHeatmap"
)
processedDir <- here(
    "processed-data",
    "06_peak_calling",
    "22_GO_enrichment"
)
plotDir <- here(
    "plots",
    "06_peak_calling",
    "22_GO_enrichment"
)

if (!dir.exists(processedDir)) { dir.create(processedDir) }
if (!dir.exists(plotDir)) { dir.create(plotDir) }


##==============================================================================

message("Loading cCRE and ORC file ...")

f_name = "overlaps_linkPeak_DARs_classified_thr_CC0.3_thr_DAR0.1.csv"

## load raw overlaps and 2-shared ct overlaps 
cCRE_df <- read.csv(here(inputCSV_cCRE_ORC_Dir, f_name))
nrow(cCRE_df) # 5603
head(cCRE_df)

## validation
unique(cCRE_df$type_classification)
as.data.frame(table(cCRE_df$category))

only_cCRE_df <- cCRE_df |> 
    filter(category %in% c("cell-specific cCRE (-)", "cell-specific cCRE (+)"))

table(only_cCRE_df$category)
# cell-specific cCRE (-) cell-specific cCRE (+) 
# 657                    706 
table(only_cCRE_df$overlap_type)
# Shared Unique 
# 508    855 

cluster_levels <- only_cCRE_df$cell_type |> unique()

only_cCRE_df |> count(cluster)

head(only_cCRE_df)

summary(only_cCRE_df$fdr_dars)

#DE_data |> filter(vlmf_adj.P.Val < 0.05) |> count(cluster)
only_cCRE_df |> filter(fdr_dars < 0.1) |> count(cell_type)


## =============================================================================

## ENTREZID look up
names(only_cCRE_df)
## gene symbol and gene entrez name
only_cCRE_df[c("gene_name", "gene_id")]

entrez_search <- bitr(only_cCRE_df$gene_id, 
                      fromType = "ENSEMBL", 
                      toType = "ENTREZID", 
                      OrgDb = "org.Hs.eg.db")
# Warning message:
#     In bitr(only_cCRE_df$gene_id, fromType = "ENSEMBL", toType = "ENTREZID",  :
#                 1.43% of input gene IDs are fail to map.

entrez_search |> count(ENSEMBL) |> count(n)

# DE_entrez <- DE_data |> 
#     left_join(entrez_search, by = c("gene_id" = "ENSEMBL"), relationship = "many-to-many") |>
#     filter(!is.na(ENTREZID)) |>
#     mutate(DE_class = case_when(vlmf_logFC > 0 & vlmf_adj.P.Val < 0.05 ~ "up",
#                                 vlmf_logFC < 0 & vlmf_adj.P.Val < 0.05 ~ "down",
#                                 TRUE ~ "None"),
#            DE_class_cluster = paste0(gsub("\\.", "-", cluster), "_",DE_class)) ## doesn't like .  in cluster names

DE_entrez <- only_cCRE_df |> 
    left_join(entrez_search, by = c("gene_id" = "ENSEMBL"), relationship = "many-to-many") |>
    filter(!is.na(ENTREZID)) |>
    mutate(DE_class = case_when(logFC > 0 & fdr_dars < 0.05 ~ "up",
                                    logFC < 0 & fdr_dars < 0.05 ~ "down",
                                TRUE ~ "None"),
           DE_class_cluster = paste0(gsub("\\.", "-", cell_type), "_",only_cCRE_df)) ## doesn't like .  in cluster names

DE_entrez |> count(cluster)
DE_entrez |> filter(DE_class != "None") |> count(DE_class_cluster)




