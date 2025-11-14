#   Prepare MAGMA inputs for the OUD GWAS, which appeared incompletely processed
#   in the habenula pilot repo

library(tidyverse)
library(here)
library(sessioninfo)
library(data.table)

gwas_path = '/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/13_MAGMA/GWAS/OUD/OUD.phs001672.pha004954.txt'
gwas_colnames = c('SNP ID', 'P-value', 'Chr ID', 'Chr Position')
snp_out_path = here('processed-data', '10_MAGMA', 'OUD', 'SNPs.tsv')
p_out_path = here('processed-data', '10_MAGMA', 'OUD', 'p_values.tsv')

gwas_df = fread(gwas_path, skip = 17, select = gwas_colnames) |>
    as_tibble() |>
    dplyr::rename(
        SNP = `SNP ID`, P = `P-value`, CHR = `Chr ID`, BP = `Chr Position`
    )

gwas_df |>
    select(SNP, CHR, BP) |>
    write_tsv(snp_out_path)

gwas_df |>
    select(SNP, P) |>
    write_tsv(p_out_path)

session_info()
