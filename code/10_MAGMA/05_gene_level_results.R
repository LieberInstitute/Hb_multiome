library(tidyverse)
library(here)
library(ComplexHeatmap)
library(RColorBrewer)
library(sessioninfo)

cell_type_groups = c('broad', 'semi_broad', 'mid')
gene_set_paths = here(
    'processed-data', '10_MAGMA', 'gene_sets',
    sprintf('%s.tsv', cell_type_groups)
)
names(gene_set_paths) = cell_type_groups
sig_cutoff = 0.05

#   Paths to gene-level statistics
hb_pilot_dir = '/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/13_MAGMA/GWAS'
gene_stat_paths = c(
    MDD = file.path(hb_pilot_dir, 'MDD/MDD.phs001672.pha005122.genes.out'),
    MDD2019 = file.path(
        hb_pilot_dir, 'mdd2019edinburgh/PGC_UKB_depression_genome-wide.genes.out'
    ),
    panic = file.path(hb_pilot_dir, 'panic2019/pgc-panic2019.genes.out'),
    SCZ = file.path(
        hb_pilot_dir,
        'scz2022/PGC3_SCZ_wave3.european.autosome.public.v3.ensembl.genes.out'
    ),
    SUD2020 = file.path(hb_pilot_dir, 'sud2020op/opi.DEPvEXP_EUR.noAF.genes.out')
)

#   Why are there NAs in 'p' and 'gwas' (probably after the left join)?





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
            select(gene_id, cell_type, peak_category, p, gwas) |>
            group_by(cell_type, peak_category) |>
            arrange(p) |>
            slice_head(n = 5) |>
            ungroup()
        
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

gene_df_list |>
    bind_rows()