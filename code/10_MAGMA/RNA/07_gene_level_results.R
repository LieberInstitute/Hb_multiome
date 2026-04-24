library(tidyverse)
library(here)
library(rtracklayer)
library(sessioninfo)

cell_type_groups = c('broad', 'mid', 'fine')
gwas_groups = list(
    substance = c(
        'SUD2020', 'AUD', 'CUD', 'ext_cannabis', 'lifetime_cannabis', 'OUD',
        'SUD2', 'SUD3'
    ),
    non_substance = c(
        'MDD', 'panic', 'SCZ', 'compulsive', 'internalizing', 'neurodev',
        'p_factor', 'SCZ_BPD'
    )
)
gene_set_paths = here(
    "processed-data", "10_MAGMA", "RNA",
    'gene_sets', sprintf('%s.tsv', cell_type_groups)
)
gene_stat_paths = here(
    "processed-data", "10_MAGMA", "RNA",
    unlist(gwas_groups), sprintf('mean_ratio_%s.genes.out', unlist(gwas_groups))
)
set_stat_paths = here(
    'processed-data', '10_MAGMA', 'RNA',
    '%s', 'mean_ratio', '%s.gsa.out'
)
out_path = here(
    'processed-data', '10_MAGMA', 'RNA',
    'top_genes.csv'
)
reference_gtf = '/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz'
names(gene_set_paths) = cell_type_groups
names(gene_stat_paths) = unlist(gwas_groups)
sig_cutoff = 0.05
gwas_renaming = c(
    'MDD' = 'MDD',
    'panic' = 'Panic Disorder',
    'compulsive' = 'Compuls. Dis.',
    'SCZ' = 'SCZ',
    'SCZ_BPD' = 'SCZ/BPD',
    'AUD' = 'AUD',
    'CUD' = 'CUD',
    'ext_cannabis' = 'Ext. Cannabis',
    'lifetime_cannabis' = 'Life. Cannabis',
    'SUD2020' = 'OUD 1',
    'OUD' = 'OUD 2',
    'SUD2' = 'SUD 1',
    'SUD3' = 'SUD 2',
    'internalizing' = 'Intern. Disorders',
    'neurodev' = 'Neurodev.',
    'p_factor' = 'P Factor'
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
            dplyr::rename(cell_type = VARIABLE) |>
            mutate(
                set_is_sig = P < sig_cutoff,
                cell_type_res = cell_type_group
            ) |>
            select(cell_type, set_is_sig, cell_type_res)

        #   Read in gene sets themselves and merge with gene- and set-level
        #   stats
        gene_df_list[[paste(gwas, cell_type_group, sep = '_')]] = read_table(
                gene_set_paths[[cell_type_group]], show_col_types = FALSE
            ) |>
            left_join(gene_stat_df, by = c('gene_id' = 'GENE')) |>
            dplyr::rename(cell_type = set_id) |>
            left_join(set_df, by = 'cell_type') |>
            dplyr::rename(p = P) |>
            select(gene_id, cell_type, gwas, p, set_is_sig, cell_type_res)
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

gene_df = gene_df |>
    filter(!is.na(p), p < sig_cutoff, set_is_sig) |>
    mutate(gwas = gwas_renaming[gwas])

#   Do we have enough genes for meaningful testing?
gene_df |>
    group_by(cell_type, gwas, cell_type_res) |>
    summarize(n = n()) |>
    arrange(cell_type_res, cell_type) |>
    print(n = Inf)

#   Read in the GTF to get gene symbols
gtf = import(reference_gtf)
gtf = gtf[gtf$type == 'gene'] |>
    as.data.frame() |>
    as_tibble() |>
    select(gene_id = gene_id, gene_name = gene_name)

#   Export final gene sets, only including genes where the set
#   as a whole was significant
gene_df |>
    arrange(cell_type_res, gwas, cell_type, p) |>
    left_join(gtf, by = 'gene_id') |>
    select(cell_type_res, gwas, cell_type, gene_id, gene_name, p) |>
    write_csv(out_path)

session_info()