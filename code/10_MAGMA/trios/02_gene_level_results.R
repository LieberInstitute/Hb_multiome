library(tidyverse)
library(here)
library(rtracklayer)
library(sessioninfo)

gwas_name_path = here('processed-data', '10_MAGMA', 'RNA', 'gwas_info.csv')
gwas_map = read_csv(gwas_name_path, show_col_types = FALSE) |>
    select(nickname, manuscript_name)
these_gwas_names = ifelse(
    gwas_map$nickname == 'MDD2019', 'MDD', gwas_map$nickname
)

cell_type_groups = c('broad', 'fine')
gene_set_paths = here(
    "processed-data", "13_tripod_trios", "11_GO",
    sprintf("gene_sets_%s.tsv", cell_type_groups)
)
gene_stat_paths = here(
    "processed-data", "10_MAGMA", "RNA",
    these_gwas_names, sprintf('mean_ratio_%s.genes.out', these_gwas_names)
)
set_stat_paths = here('processed-data', '10_MAGMA', 'trios', '%s', '%s.gsa.out')
out_path = here('processed-data', '10_MAGMA', 'trios', 'top_genes.csv')
out_low_genes_path = here(
    'processed-data', '10_MAGMA', 'trios', 'low_gene_sets.csv'
)
reference_gtf = '/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz'
names(gene_set_paths) = cell_type_groups
names(gene_stat_paths) = these_gwas_names
sig_cutoff = 0.05
min_genes_per_set = 10

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
    filter(!is.na(p)) |>
    mutate(gwas = ifelse(gwas == 'MDD', 'MDD2019', gwas)) |>
    left_join(gwas_map, by = c('gwas' = 'nickname'))

#   Do we have enough genes for meaningful testing? Export sets with too few
#   genes, so we can note them in the heatmaps
gene_df |>
    group_by(cell_type, gwas, cell_type_res) |>
    filter(n() < min_genes_per_set) |>
    ungroup() |>
    distinct(cell_type_res, cell_type, manuscript_name) |>
    write_csv(out_low_genes_path)

#   Read in the GTF to get gene symbols
gtf = import(reference_gtf)
gtf = gtf[gtf$type == 'gene'] |>
    as.data.frame() |>
    as_tibble() |>
    select(gene_id = gene_id, gene_name = gene_name)

#   Export final gene sets, only including genes where the set
#   as a whole was significant
gene_df |>
    filter(p < sig_cutoff, set_is_sig) |>
    #   Require sets to have a minimum number of genes (to accurately determine
    #   set-level significance)
    group_by(cell_type, gwas, cell_type_res) |>
    filter(n() >= min_genes_per_set) |>
    ungroup() |>
    arrange(cell_type_res, gwas, cell_type, p) |>
    left_join(gtf, by = 'gene_id') |>
    select(cell_type_res, manuscript_name, cell_type, gene_id, gene_name, p) |>
    dplyr::rename(gwas = manuscript_name) |>
    write_csv(out_path)

session_info()
