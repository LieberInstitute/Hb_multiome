library(tidyverse)
library(here)
library(rtracklayer)
library(ComplexHeatmap)
library(RColorBrewer)
library(sessioninfo)

cell_type_groups = c('broad', 'semi_broad', 'mid')
gwas_groups = list(
    substance = c(
        'SUD2020', 'AUD', 'CUD', 'ext_cannabis', 'lifetime_cannabis', 'OUD',
        'SUD2', 'SUD3'
    ),
    non_substance = c(
        'MDD2019', 'panic', 'SCZ', 'compulsive', 'internalizing', 'neurodev',
        'p_factor', 'SCZ_BPD'
    )
)
gene_set_paths = here(
    'processed-data', '10_MAGMA', 'gene_sets',
    sprintf('%s.tsv', cell_type_groups)
)
gene_stat_paths = here(
    'processed-data', '10_MAGMA', unlist(gwas_groups),
    sprintf('%s.genes.out', unlist(gwas_groups))
)
set_stat_paths = here('processed-data', '10_MAGMA', '%s', '%s.gsa.out')
out_path = here('processed-data', '10_MAGMA', 'top_genes.csv')
plot_dir = here('plots', '10_MAGMA')
reference_gtf = '/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz'
names(gene_set_paths) = cell_type_groups
names(gene_stat_paths) = unlist(gwas_groups)
sig_cutoff = 0.05
gwas_renaming = c(
    'MDD2019' = 'MDD',
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

dir.create(plot_dir, showWarnings = FALSE)

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
                set_is_sig = P < sig_cutoff,
                cell_type_res = cell_type_group
            ) |>
            select(cell_type, peak_category, set_is_sig, cell_type_res)

        #   Read in gene sets themselves and merge with gene- and set-level
        #   stats
        gene_df = read_table(
                gene_set_paths[[cell_type_group]], show_col_types = FALSE
            ) |>
            left_join(gene_stat_df, by = c('link_gene_id' = 'GENE')) |>
            mutate(
                cell_type = str_extract(set_id, '^[^_]+'),
                peak_category = str_extract(set_id, '(?<=_).+')
            ) |>
            left_join(set_df, by = c('cell_type', 'peak_category')) |>
            dplyr::rename(p = P, gene_id = link_gene_id) |>
            select(
                gene_id, cell_type, peak_category, gwas, p, set_is_sig,
                cell_type_res
            )
        
        #   Cell-type resolutions only differ in how they treat habenula types,
        #   which means there would be plenty of duplicated information if we
        #   didn't drop the non-Hb cell types in all but the mid resolution
        if (cell_type_group == 'mid') {
            gene_df_list[[paste(gwas, cell_type_group, sep = '_')]] = gene_df
        } else {
            gene_df_list[[paste(gwas, cell_type_group, sep = '_')]] = gene_df |>
                filter(grepl('^([ML])*Hb', cell_type))
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
gene_df = gene_df |>
    filter(!is.na(p), p < sig_cutoff, set_is_sig) |>
    mutate(gwas = gwas_renaming[gwas])

gene_df |>   
    arrange(gwas, cell_type, peak_category, p) |>
    left_join(gtf, by = 'gene_id') |>
    select(gwas, cell_type, peak_category, gene_id, gene_name, p) |>
    write_csv(out_path)

#   Grab unique genes per substance-use-related GWAS and cell type
gene_df = gene_df |>
    filter(gwas %in% gwas_renaming[gwas_groups[['substance']]]) |>
    group_by(gwas, cell_type, cell_type_res) |>
    filter(!duplicated(gene_id)) |>
    ungroup() |>
    select(gene_id, cell_type, gwas, cell_type_res) |>
    mutate(cell_type_res = factor(cell_type_res, levels = cell_type_groups))

#   Within a given cell-type resolution and cell type, how many significant
#   genes are unique to each GWAS (taking the union of gene sets across peak
#   categories first)?
unique_df_list = list()
for (gwas in unique(gene_df$gwas)) {
    for (cell_type_group in cell_type_groups)
        for (cell_type in unique(gene_df[gene_df$cell_type_res == cell_type_group, ]$cell_type)) {
            these_genes = gene_df |>
                filter(
                    gwas == !!gwas, cell_type == !!cell_type,
                    cell_type_res == !!cell_type_group
                ) |>
                pull(gene_id)
            other_genes = gene_df |>
                filter(
                    gwas != !!gwas, cell_type == !!cell_type,
                    cell_type_res == !!cell_type_group
                ) |>
                pull(gene_id)

            unique_df_list[[length(unique_df_list) + 1]] = tibble(
                num_unique = length(setdiff(these_genes, other_genes)),
                cell_type = cell_type,
                gwas = gwas,
                cell_type_res = cell_type_group
            )
        }
}
unique_df = bind_rows(unique_df_list)

p = unique_df |>
    ggplot(aes(x = gwas, y = cell_type, fill = num_unique)) +
        geom_tile() +
        scale_fill_viridis_c() +
        facet_grid(
            rows = vars(cell_type_res), scales = "free_y", space  = "free_y"
        ) +
        theme_bw(base_size = 20) +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            strip.text.y.right = element_text(angle = 0)
        ) +
        labs(x = "GWAS", y = "Cell Type", fill = "Num Unique\nGenes")
pdf(
    file = file.path(plot_dir, "shared_heatmap_substance_num_unique.pdf"),
    width = 9, height = 6
)
print(p)
dev.off()

#   Get the original size of the union across peak categories per cell type
#   which we'll use to scale the number of unique genes from above (for
#   easier interpretation)
gene_set_df_list = list()
for (cell_type_group in cell_type_groups) {
    gene_set_df_list[[cell_type_group]] = read_table(
            gene_set_paths[[cell_type_group]], show_col_types = FALSE
        ) |>
        mutate(
            cell_type = str_extract(set_id, '^[^_]+'),
            cell_type_res = cell_type_group
        )
}
gene_set_df = bind_rows(gene_set_df_list) |>
    group_by(cell_type, cell_type_res) |>
    summarize(num_genes = length(unique(link_gene_id))) |>
    ungroup()

p = unique_df |>
    left_join(gene_set_df, by = c('cell_type', 'cell_type_res')) |>
    mutate(
        perc_unique = 100 * num_unique / num_genes
    ) |>
    ggplot(aes(x = gwas, y = cell_type, fill = perc_unique)) +
        geom_tile() +
        scale_fill_viridis_c() +
        facet_grid(
            rows = vars(cell_type_res), scales = "free_y", space  = "free_y"
        ) +
        theme_bw(base_size = 20) +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            strip.text.y.right = element_text(angle = 0)
        ) +
        labs(x = "GWAS", y = "Cell Type", fill = "% Unique\nGenes")
pdf(
    file = file.path(plot_dir, "shared_heatmap_substance_perc_unique.pdf"),
    width = 9, height = 6
)
print(p)
dev.off()

session_info()
      