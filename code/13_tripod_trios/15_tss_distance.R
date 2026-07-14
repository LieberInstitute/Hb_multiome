library(sessioninfo)
library(tidyverse)
library(here)
library(duckplyr)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(GenomicRanges)

trio_path = here(
    "processed-data", "13_tripod_trios", "09_gather_trios",
    "filtered_trios_fine.parquet"
)
plot_dir = here("plots", "13_tripod_trios", "15_tss_distance")

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(plot_dir, showWarnings = FALSE)

trio_df = read_parquet_duckdb(trio_path, prudence = "stingy") |>
    dplyr::filter(is_intersect, stringency_level == 1) |>
    dplyr::select(peak, gene, TF, cell_type, coef, adj) |>
    collect()

tss_gr = genes(EnsDb.Hsapiens.v86)
tss_gr = promoters(tss_gr, upstream = 0, downstream = 1)
tss_gr = keepStandardChromosomes(tss_gr, pruning.mode = "coarse")
seqlevelsStyle(tss_gr) = "UCSC"

tss_df = as.data.frame(tss_gr) |>
    as_tibble() |>
    transmute(
        gene = gene_name,
        peak_chr = as.character(seqnames),
        tss = start,
        gene_strand = as.character(strand)
    ) |>
    distinct(gene, peak_chr, .keep_all = TRUE)

trio_distance_df = trio_df |>
    separate_wider_delim(
        peak, delim = "-", names = c("peak_chr", "peak_start", "peak_end")
    ) |>
    mutate(
        peak_start = as.numeric(peak_start),
        peak_end = as.numeric(peak_end),
        peak_center = (peak_start + peak_end) / 2
    ) |>
    left_join(tss_df, by = c("gene", "peak_chr")) |>
    mutate(
        distance_bp = abs(peak_center - tss),
        distance_kb = distance_bp / 1000
    )

p = trio_distance_df |>
    ggplot(aes(x = distance_kb)) +
    geom_density() +
    labs(
        x = "Distance: Peak Center to Gene TSS (kb)",
        y = "Density"
    ) +
    theme_bw(base_size = 20)

ggsave(
    filename = here(plot_dir, "trio_peak_center_to_tss_density.pdf"),
    plot = p
)

session_info()
