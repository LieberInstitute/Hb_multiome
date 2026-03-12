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

in_path = here(
    "processed-data", "12_new_peaks", "09_tripod_preprocess",
    "preprocessed_objects.qs2"
)
out_dir = here("processed-data", "12_new_peaks", "10_tripod_trios")
fdr_cutoff = 0.05

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

seur@meta.data$mid_cluster = factor(seur@meta.data$mid_cluster)
color_seur = getColors(object = seur, celltype.col.name = "mid_cluster")

#   Fit models
gene_vec = tripod_seur$transcripts.gr$gene_name[1:100]
xy_mat_list = bplapply(
    gene_vec,
    getXYMatrices,
    ext.upstream = 1e5,
    transcripts.gr = tripod_seur$transcripts.gr,
    peaks.gr = tripod_seur$peaks.gr,
    metacell.rna = metacell_seur$rna,
    metacell.peak = metacell_seur$peak,
    peakxmotif = tripod_seur$peakxmotif,
    motifxTF = tripod_seur$pbmc.motifxTF,
    metacell.celltype = color_seur$metacell$celltype,
    metacell.celltype.col = color_seur$metacell$color,
    BPPARAM = BPPARAM
)
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
    result_list[[condition_type]] = this_result
    
    for (stringency_level in c(1, 2)) {
        result_list[[length(result_list) + 1]] = getTrios(
                xymats.list = this_result,
                fdr.thresh = fdr_cutoff,
                sign = "positive",
                model.name = "TRIPOD",
                level = stringency_level
            ) |>
            as_tibble() |>
            mutate(condition_on = condition_type, stringency_level = stringency_level)
    }
}

#   Write the main results
result_list |>
    bind_rows() |>
    write_csv(file.path(out_dir, "trios.csv.gz"))

#   It also appears from the vignette that the fit model lists are also required
#   for certain visualizations, so we'll save those too
qs_save(result_list, file.path(out_dir, "fit_models.qs2"))

session_info()
