library(tidyverse)
library(here)
library(ComplexHeatmap)
library(RColorBrewer)
library(sessioninfo)

cell_type_groups = c('broad', 'semi_broad', 'mid')
gwas_groups = c('MDD', 'MDD2019', 'panic', 'SCZ', 'SUD2020')
gene_set_paths = here(
    'processed-data', '10_MAGMA', 'gene_sets',
    sprintf('%s.tsv', cell_type_groups)
)
gene_stat_paths = here(
    'processed-data', '10_MAGMA', gwas_groups,
    sprintf('%s.genes.out', gwas_groups)
)
out_path = here('processed-data', '10_MAGMA', 'gene_level_results.csv')
reference_gtf = '/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz'
names(gene_set_paths) = cell_type_groups
names(gene_stat_paths) = gwas_groups
sig_cutoff = 0.05

gene_df_list = list()
for (gwas in names(gene_stat_paths)) {
    gene_stat_df = read_table(
            gene_stat_paths[[gwas]], show_col_types = FALSE
        ) |>
        mutate(gwas = gwas) |>
        select(GENE, P, gwas)
    for (cell_type_group in cell_type_groups) {
        gene_df = read_table(
                gene_set_paths[[cell_type_group]], show_col_types = FALSE
            ) |>
            left_join(gene_stat_df, by = c('gene_id' = 'GENE')) |>
            mutate(
                cell_type = str_extract(set_id, '^[^_]+'),
                peak_category = str_extract(set_id, '(?<=_).+')
            ) |>
            rename(p = P) |>
            select(gene_id, cell_type, peak_category, gwas, p)
        
        #   Cell-type resolutions only differ in how they treat habenula types,
        #   which means there would be plenty of duplicated information if we
        #   didn't drop the non-Hb cell types in all but the broad resolution
        if (cell_type_group == 'broad') {
            gene_df_list[[paste(gwas, cell_type_group, sep = '_')]] = gene_df
        } else {
            gene_df_list[[paste(gwas, cell_type_group, sep = '_')]] = gene_df |>
                filter(grepl('^[ML]Hb', cell_type))
        }
    }
}

gene_df = bind_rows(gene_df_list)

#   Warn about percentage of genes missing MAGMA stats
for (gwas in names(gene_stat_paths)) {
    message(
        sprintf(
            'Dropping %d%% of genes for GWAS %s missing MAGMA stats',
            round(100 * mean(is.na(gene_df$p[gene_df$gwas == gwas]))),
            gwas
        )
    )
}

gene_df |>
    filter(!is.na(p)) |>
    group_by(gwas, cell_type, peak_category) |>
    arrange(p) |>
    slice_head(n = 5) |>
    ungroup() |>
    write_csv(out_path)

session_info()
