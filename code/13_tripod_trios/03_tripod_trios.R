#   Using previously preprocessed objects, fit TRIPOD models and write
#   significant trios as a CSV for each type of test TRIPOD allows. Largely
#   following the vignette:
#   https://htmlpreview.github.io/?https://github.com/yuchaojiang/TRIPOD/blob/main/vignettes/TRIPOD_pbmc.html

library(TRIPOD)
library(Seurat)
library(Signac)
library(GenomicRanges)
library(tidyverse)
library(here)
library(sessioninfo)
library(qs2)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
library(BiocParallel)
library(dendextend)
library(duckplyr)

cell_types = c(
    "Astrocyte", "Endo", "Excit.Thal", "Inhib_LHb_4.1", "Inhib_LHb_4.2",
    "Inhib.Thal", "LHb.1.3.4", "LHb.2.7", "LHb.4", "MHb.1", "MHb.1.2",
    "MHb.2", "MHb.3", "Microglia", "Oligo", "OPC", "all"
)
this_cell_type = cell_types[as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))]

in_path = here(
    "processed-data", "13_tripod_trios", "02_tripod_preprocess",
    sprintf("preprocessed_objects_%s.qs2", this_cell_type)
)
out_dir = here("processed-data", "13_tripod_trios", "03_tripod_trios")
fdr_cutoff = 0.1

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(out_dir, showWarnings = FALSE)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (num_cores == 1) {
    BPPARAM = SerialParam()
} else {
    BPPARAM = MulticoreParam(num_cores)
}

#   Read preprocessed objects
pre_list = qs_read(in_path)
tripod_seur = pre_list$tripod_seur
seur = pre_list$seur
metacell_seur = pre_list$metacell_seur

#   Fit models
gene_vec = tripod_seur$transcripts.gr$gene_name
xy_mat_list = bplapply(
    gene_vec,
    getXYMatrices,
    ext.upstream = 5e5,
    ext.downstream = 5e5,
    transcripts.gr = tripod_seur$transcripts.gr,
    peaks.gr = tripod_seur$peaks.gr,
    metacell.rna = metacell_seur$rna,
    metacell.peak = metacell_seur$peak,
    peakxmotif = tripod_seur$peakxmotif,
    motifxTF = tripod_seur$motifxTF,
    BPPARAM = BPPARAM
)

#   See https://github.com/yuchaojiang/TRIPOD/issues/8. Basically TRIPOD can
#   generate output objects that are not valid in further steps (representing
#   genes without nearby peaks or TF binding sites). Drop those genes
valid_candidates = sapply(xy_mat_list, function(x) ncol(x[['Xt']]) > 0)
gene_vec = gene_vec[valid_candidates]
xy_mat_list = xy_mat_list[valid_candidates]
names(xy_mat_list) = gene_vec

#   match.by = "Xt": Conditioned on peak accessibility, is TF expression
#   correlated with each gene's expression?
#   match.by = "Yj": Conditioned on TF expression, is peak accessibility
#   correlated with gene expression?"
#
#   Stringency level 1 fixes the matching variable level and asks if the
#   other variable is correlated with gene expression. Stringency level 2 asks
#   if raising expression of the matching variable results in a stronger
#   relationship between the other variable and gene expression
result_list = list()
for (condition_type in c("Xt", "Yj")) {
    this_result <- bplapply(
        xy_mat_list,
        fitModel,
        model.name = "TRIPOD",
        match.by = condition_type,
        BPPARAM = BPPARAM
    )
    names(this_result) = gene_vec
    
    for (stringency_level in c(1, 2)) {
        result_list[[length(result_list) + 1]] = getTrios(
                xymats.list = this_result,
                fdr.thresh = fdr_cutoff,
                sign = "positive",
                model.name = "TRIPOD",
                level = stringency_level
            ) |>
            as_tibble() |>
            mutate(
                condition_on = condition_type,
                stringency_level = stringency_level
            )
    }
}

#   Write the main results
result_list |>
    bind_rows() |>
    compute_parquet(
        file.path(out_dir, sprintf("trios_%s.parquet", this_cell_type))
    )

#   It also appears from the vignette that the fit model lists are also required
#   for certain visualizations, so we'll save those too
qs_save(
    xy_mat_list,
    file.path(out_dir, sprintf("fit_models_%s.qs2", this_cell_type))
)

session_info()
