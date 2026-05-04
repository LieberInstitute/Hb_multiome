library(tidyverse)
library(here)
library(viridis)
library(sessioninfo)

multiome_path = here('processed-data', '10_MAGMA', 'RNA', 'heatmap_data.csv')
extra_hd_path = '/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/09_HD_cell_level/no_secondary/MAGMA/extracellular/heatmap_data.csv'
cell_hd_path = '/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/09_HD_cell_level/no_secondary/MAGMA/heatmap_data.csv'
ficture_cell_hd_path = '/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/new_samples2/ficture_harmony/MAGMA/heatmap_data.csv'
ficture_extra_hd_path = '/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/cell_environment/MAGMA/heatmap_data.csv'
plot_dir = here('plots', '10_MAGMA', 'RNA')
cell_type_order = list(
    fine = c(
        'MHb.1', 'MHb.1.2', 'MHb.2', 'MHb.3', 'Excit_LHb', 'LHb.1.3.4',
        'LHb.2.7', 'LHb.4', 'Inhib_LHb_4.1', 'Inhib_LHb_4.2',
        'Excit.Thal/Inhib_LHb_4.2', 'Excit.Thal', 'Inhib.Thal', 'Astrocyte',
        'Endo', 'Endo/microglia', 'Microglia', 'Oligo', 'OPC', 'Ependymal',
        'Subependymal', paste0('X', 0:9)
    ),
    mid = c(
        'MHb', 'LHb', 'Excit.Thal/Inhib_LHb_4.2', 'Excit.Thal', 'Inhib.Thal',
        'Astrocyte', 'Endo', 'Endo/microglia', 'Microglia', 'Oligo', 'OPC',
        'Ependymal', 'Subependymal'
    ),
    broad = c(
        'Hb', 'Excit.Thal/Inhib_LHb_4.2', 'Excit.Thal', 'Inhib.Thal',
        'Astrocyte', 'Endo', 'Endo/microglia', 'Microglia', 'Oligo', 'OPC',
        'Ependymal', 'Subependymal'
    )
)
gwas_order = c(
    'MDD_Howard', 'internalizing_F4_Grotzinger', 'compulsive_F1_Grotzinger',
    'externalizing_Linnér', 'neurodev_F3_Grotzinger', 'p_factor_Grotzinger',
    'SCZ_Trubetskoy', 'SCZ/BPD_F2_Grotzinger', 'panic_Forster', 'AUD_Zhou',
    'CUD_Johnson', 'CUD_Pasman', 'SUD_Hotoum', 'SUD_Polimanti',
    'SUD_F5_Grotzinger', 'OUD_Deak'
)

################################################################################
#   Functions
################################################################################

p_val_heatmap = function(results_df, f_name) {
    p = ggplot(
            results_df,
            aes(
                x = cell_type, y = gwas_group, fill = neg_log_p, label = p_label
            )
        ) +
        geom_tile() +
        geom_text(size = 6) +
        scale_fill_viridis_c() +
        facet_grid(gwas_category ~ dataset, scales = "free", space = "free") +
        theme_bw(base_size = 20) +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            strip.text.y.right = element_text(angle = 0)
        ) +
        labs(x = "Cell Type", y = "GWAS Trait", fill = "-log10(p)")
    pdf(
        file.path(plot_dir, f_name), width = 20, height = 8
    )
    print(p)
    dev.off()
}

################################################################################
#   Main
################################################################################

multiome_df = read_csv(multiome_path, show_col_types = FALSE) |>
    mutate(dataset = "Multiome (RNA)")
extra_hd_df = read_csv(extra_hd_path, show_col_types = FALSE) |>
    mutate(dataset = "Extra HD")
cell_hd_df = read_csv(cell_hd_path, show_col_types = FALSE) |>
    mutate(dataset = "Cellular HD")
ficture_extra_hd_df = read_csv(ficture_extra_hd_path, show_col_types = FALSE) |>
    mutate(dataset = "Ficture Extra HD", cell_type_group = "fine")
ficture_cell_hd_df = read_csv(ficture_cell_hd_path, show_col_types = FALSE) |>
    mutate(dataset = "Ficture Cellular HD", cell_type_group = "fine")
results_df = bind_rows(
        multiome_df, extra_hd_df, cell_hd_df, ficture_extra_hd_df,
        ficture_cell_hd_df
    ) |>
    filter(cell_type != 'LHb.4/Inhib_LHb_4.2') |>
    mutate(
        gwas_category = case_when(
            grepl('^[ACOS]UD', gwas_group) ~ 'Substance Use',
            gwas_group == 'p_factor_Grotzinger' ~ 'P-Factor',
            TRUE ~ 'Psychiatric'
        )
    )

stopifnot(setequal(results_df$gwas_group, gwas_order))
for (this_res in names(cell_type_order)) {
    this_results_df = results_df |>
        filter(cell_type_group == this_res) |>
        mutate(
            cell_type = factor(cell_type, levels = cell_type_order[[this_res]]),
            gwas_group = factor(gwas_group, levels = gwas_order)
        )
  
    stopifnot(setequal(this_results_df$cell_type, cell_type_order[[this_res]]))
  
    p_val_heatmap(this_results_df, sprintf("heatmap_%s.pdf", this_res))
}

session_info()
