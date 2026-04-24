library(SingleCellExperiment)
library(tidyverse)
library(DeconvoBuddies)
library(here)
library(sessioninfo)
library(qs2)
library(rtracklayer)

sce_path = here("processed-data", "10_MAGMA", "RNA", "sce.rds")
out_dir = here("processed-data", "10_MAGMA", "RNA", "gene_sets")
mean_ratio_threshold = 1.05
max_num_genes = 200

dir.create(out_dir, showWarnings = FALSE)

export_set = function(sce, cell_type_col, file_tag) {    
    marker_stats = get_mean_ratio(
            sce = sce, assay_name = "logcounts",
            cellType_col = cell_type_col, gene_ensembl = "gene_name",
            gene_name = "gene_id"
        ) |>
        filter(MeanRatio > mean_ratio_threshold) |>
        dplyr::rename(
            set_id = cellType.target,
            gene_id = gene_name,
            gene_name = gene,
            mean_ratio = MeanRatio
        ) |>
        group_by(set_id) |>
        arrange(desc(mean_ratio), .by_group = TRUE) |>
        slice_head(n = max_num_genes) |>
        ungroup() |>
        select(set_id, gene_id, gene_name, mean_ratio) |>
        arrange(set_id, desc(mean_ratio))

    message(sprintf("Marker counts for %s resolution:", file_tag))
    print(table(marker_stats$set_id))

    write_tsv(marker_stats, file.path(out_dir, sprintf("%s.tsv", file_tag)))
}

# sce = readRDS(sce_path)
sce = qs_read('/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/05_03_annotation_adjustments/06_refined_annotations/refined_annotation_multiomeHab_SCE.qs2')

# add gene_id and gene_name to rowData
rowData(sce)$gene_name = rownames(rowData(sce))

reference_gtf = '/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz'
gtf = import(reference_gtf)
gtf = gtf[gtf$type == "gene"]

rowData(sce)$gene_id <- gtf$gene_id[match(rowData(sce)$gene_name, gtf$gene_name)]

table(colData(sce)$refined_mid_cluster)
table(colData(sce)$refined_cluster_ann)

keep <- !is.na(rowData(sce)$gene_id) & rowData(sce)$gene_id != ""
sce <- sce[keep, ]

# export gene sets for both broad and fine resolutions
export_set(sce, "refined_mid_cluster", "mid")
export_set(sce, "refined_cluster_ann", "fine")

session_info()


#     Astrocyte          Endo     Ependymal    Excit.Thal Inhib_LHb_4.1 
#            62           121           200           200           100 
# Inhib_LHb_4.2    Inhib.Thal     LHb.1.3.4       LHb.2.7         LHb.4 
#           200           200           200            94            11 
#         MHb.1       MHb.1.2         MHb.2         MHb.3     Microglia 
#            50            90            62            91            45 
#         Oligo           OPC 
#           200            57 


# C.01.Inhib.Thal      C.02.Oligo C.03.Excit.Thal      C.04.LHb.4    C.05.LHb.2.7 
#              10              96             200              34              26 
#      C.07.MHb.2      C.08.LHb.4      C.09.LHb.4      C.10.MHb.1    C.11.MHb.1.2 
#              63              83               1              56              13 
# C.12.Excit.Thal      C.13.LHb.4      C.14.MHb.1 C.15.Excit.Thal    C.16.MHb.1.2 
#             200             200               8              73              46 
# C.17.Excit.Thal  C.18.LHb.1.3.4 C.19.Inhib.Thal  C.20.Astrocyte  C.21.Ependymal 
#              34              30              48              46             200 
#      C.22.Oligo      C.23.LHb.1 C.25.Excit.Thal        C.26.OPC  C.27.Microglia 
#             115              61               1              38              41 
# C.28.Inhib.Thal       C.29.Endo      C.31.LHb.4 C.32.Excit.Thal    C.33.LHb.1.3 
#              31             100             101             200             145 
# C.35.Excit.Thal      C.36.MHb.3 C.38.Inhib.Thal  C.41.Microglia   Inhib_LHb_4.1 
#              75              36             121              17               3 
#   Inhib_LHb_4.2 
#               4 