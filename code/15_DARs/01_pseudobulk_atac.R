library(sessioninfo)
library(Seurat)
library(Signac)
library(tidyverse)
library(here)
library(qs2)
library(GenomicRanges)

resolution = c("broad", "mid", "fine")[
    as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
]

seur_in_path = here(
    'processed-data', '11_link_prep', '02_rebuild_atac_assay',
    'cell_level_seur.qs2'
)
seur_out_path = here(
    'processed-data', '15_DARs', '01_pseudobulk_atac',
    sprintf('seur_pb_%s.qs2', resolution)
)
task_map_out_path = here(
    'processed-data', '15_DARs', '01_pseudobulk_atac',
    'task_map.csv'
)

dir.create(dirname(seur_out_path), showWarnings = FALSE)

seur = qs_read(seur_in_path)

#   Assign the 3 cell-type resolutions. Shockingly, the 'seur@meta.data$' syntax
#   is required; assigning like 'seur$' causes unexpected NAs and in some cases,
#   failure to even assign values
seur@meta.data$broad_class = ifelse(
    grepl('[ML]Hb|Thal', seur$refined_mid_cluster),
    'Neuron',
    seur$refined_mid_cluster
)
seur@meta.data$mid_class = case_when(
    grepl('^MHb', seur$refined_mid_cluster) ~ 'MHb',
    grepl('LHb', seur$refined_mid_cluster) ~ 'LHb',
    grepl('Thal', seur$refined_mid_cluster) ~ 'Thalamus',
    TRUE ~ seur$refined_mid_cluster
)
seur@meta.data$fine_class = seur$refined_mid_cluster

#   Helper to make the DAR array job work cleanly. Note slurmjobs::job_loop()
#   doesn't support loops where different elements have different lengths in
#   the other variable (number of cell types depends on resolution)
if (resolution == 'fine') {
    rbind(
            tibble(resolution = "broad", cell_type = unique(seur$broad_class)),
            tibble(resolution = "mid", cell_type = unique(seur$mid_class)),
            tibble(resolution = "fine", cell_type = unique(seur$fine_class))
        ) |>
        arrange(resolution, cell_type) |>
        #   Since we're using the same underlying data across resolutions,
        #   there's no need to run the identical non-neuronal cell types 3
        #   times. We'll just do it for fine
        filter(
            str_detect(cell_type, 'Hb|^Neuron|Thal') | (resolution == 'fine')
        ) |>
        mutate(task_id = row_number()) |>
        write_csv(task_map_out_path)
}

#   Pseudobulk and rebuild
atac_mat = AggregateExpression(
    seur, assays = 'ATAC',
    group.by = c(paste0(resolution, '_class'), 'orig.ident'),
    return.seurat = FALSE, verbose = TRUE
)[['ATAC']]

seur_pb = CreateSeuratObject(
    counts = CreateChromatinAssay(
        counts = atac_mat, ranges = granges(seur[['ATAC']]),
        annotation = Annotation(seur[['ATAC']])
    ),
    assay = "ATAC"
)

#   Not sure why this was dropped
seur_pb@meta.data$donor = rownames(seur_pb@meta.data) |>
    str_extract('S[0-9]{2}-Hb-r$') |>
    str_replace_all('-', '_')

#   Drop 'data' layer (we only need 'counts' for DAR calculation)
seur_pb[["ATAC"]] = subset(
    seur_pb[["ATAC"]], cells = colnames(seur_pb), layer = "counts"
)

qs_save(seur_pb, seur_out_path)

session_info()
