library(tidyverse)
library(here)
library(rtracklayer)
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
set_stat_paths = here('processed-data', '10_MAGMA', '%s', '%s.gsa.out')
out_path = here('processed-data', '10_MAGMA', 'top_genes.csv')
reference_gtf = '/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz'
names(gene_set_paths) = cell_type_groups
names(gene_stat_paths) = gwas_groups
sig_cutoff = 0.05
gwas_renaming = c(
    'MDD2019' = 'MDD',
    'panic' = 'Panic Disorder',
    'SCZ' = 'SCZ',
    'SUD2020' = 'OUD'
)

#   MAGMA set-level outputs have a variable amount of header lines. Auto-detect
#   the header length and read in dynamically
read_table_auto_skip = function(path, check_lines = 100) {
    n_skip = sum(grepl('^#', readLines(path, n = check_lines)))
    clean_df = read_table(path, skip = n_skip, show_col_types = FALSE)
    return(clean_df)
}

gene_df_list = list()
for (gwas in names(gene_stat_paths)) {
    gene_stat_df = read_table(
            gene_stat_paths[[gwas]], show_col_types = FALSE
        ) |>
        mutate(gwas = gwas) |>
        select(GENE, P, gwas)
    for (cell_type_group in cell_type_groups) {
        #   Read in set-level stats to identify significant sets
        set_df = sprintf(set_stat_paths, gwas, cell_type_group) |>
            read_table_auto_skip() |>
            mutate(
                cell_type = str_extract(FULL_NAME, '^[^_]+'),
                peak_category = str_extract(FULL_NAME, '(?<=_).+'),
                set_is_sig = P < sig_cutoff
            ) |>
            select(cell_type, peak_category, set_is_sig)

        #   Read in gene sets themselves and merge with gene- and set-level
        #   stats
        gene_df = read_table(
                gene_set_paths[[cell_type_group]], show_col_types = FALSE
            ) |>
            left_join(gene_stat_df, by = c('gene_id' = 'GENE')) |>
            mutate(
                cell_type = str_extract(set_id, '^[^_]+'),
                peak_category = str_extract(set_id, '(?<=_).+')
            ) |>
            left_join(set_df, by = c('cell_type', 'peak_category')) |>
            dplyr::rename(p = P) |>
            select(gene_id, cell_type, peak_category, gwas, p, set_is_sig)
        
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
            'Dropping %d%% of genes for %s GWAS missing MAGMA stats',
            round(100 * mean(is.na(gene_df$p[gene_df$gwas == gwas]))),
            gwas
        )
    )
}

#   Read in the GTF to get gene symbols
gtf = import(reference_gtf)
gtf = gtf[gtf$type == 'gene'] |>
    as.data.frame() |>
    as_tibble() |>
    select(gene_id = gene_id, gene_name = gene_name)

#   Export final gene sets, only including genes where the set
#   as a whole was significant
gene_df |>
    filter(!is.na(p), p < sig_cutoff, set_is_sig, gwas != 'MDD') |>
    mutate(gwas = gwas_renaming[gwas]) |>
    arrange(gwas, cell_type, peak_category, p) |>
    left_join(gtf, by = 'gene_id') |>
    select(gwas, cell_type, peak_category, gene_id, gene_name, p) |>
    write_csv(out_path)

session_info()
