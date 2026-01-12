#   Prepare SNP and p-value input files for MAGMA for each GWAS. At times, we
#   also need to lift from hg38 to hg19; all analysis is done with LD estimated
#   in European ancestry and using hg19 coordinates.

library(here)
library(rtracklayer)
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(sessioninfo)

chain_path = here('processed-data', '10_MAGMA', 'hg38ToHg19.over.chain')
mdd_path = '/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/13_MAGMA/GWAS/MDD/MDD.phs001672.pha005122.txt'
mdd2019_path = '/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/13_MAGMA/GWAS/mdd2019edinburgh/PGC_UKB_depression_genome-wide.txt'
aud_path = here('processed-data', '10_MAGMA', 'AUD', 'AUDSumstats.txt.gz')
cud_path = here('processed-data', '10_MAGMA', 'CUD', 'CUDEASumstats.txt.gz')
ext_cannabis_path = here(
    'processed-data', '10_MAGMA', 'ext_cannabis', 'ExtSumstats.txt.gz'
)
lifetime_cannabis_path = here(
    'processed-data', '10_MAGMA', 'lifetime_cannabis',
    'LifetimeCannabisICC.txt.gz'
)
oud_path = here('processed-data', '10_MAGMA', 'OUD', 'OUDEASumstats.txt.gz')
sud_path = here('processed-data', '10_MAGMA', 'SUD2', 'SUDEA.txt.gz')
mdd_out_dir = here('processed-data', '10_MAGMA', 'MDD')
mdd2019_out_dir = here('processed-data', '10_MAGMA', 'MDD2019')
aud_out_dir = here('processed-data', '10_MAGMA', 'AUD')
cud_out_dir = here('processed-data', '10_MAGMA', 'CUD')
ext_cannabis_out_dir = here('processed-data', '10_MAGMA', 'ext_cannabis')
lifetime_cannabis_out_dir = here(
    'processed-data', '10_MAGMA', 'lifetime_cannabis'
)
oud_out_dir = here('processed-data', '10_MAGMA', 'OUD')
sud_out_dir = here('processed-data', '10_MAGMA', 'SUD2')
mdd_N = 1154267
mdd2019_N = 807553

################################################################################
#   Functions
################################################################################

#   Given an input tibble with columns SNP, CHR, BP, P, and N, lift the SNPs
#   from hg38 to hg19 using the provided chain file and return a lifted tibble
#   with the same columns
lift_df = function(input_df, chain) {
    this_gr = GRanges(
        seqnames = paste0('chr', input_df$CHR),
        ranges = IRanges(start = input_df$BP, width = 1),
        SNP = input_df$SNP
    )

    lifted = unlist(liftOver(this_gr, chain))
    message(
        sprintf(
            'Lifted over %d of %d SNPs to hg19',
            length(lifted), length(this_gr)
        )
    )

    out_df = tibble(
            SNP = lifted$SNP,
            CHR = gsub('chr', '', as.character(seqnames(lifted))),
            BP = start(lifted)
        ) |>
        left_join(input_df |> select(SNP, P, N), by = 'SNP')
    stopifnot(!any(is.na(out_df)))

    return(out_df)
}

#   Given an input tibble with columns SNP, CHR, BP, P, and N, write the SNP and
#   P-value input files for MAGMA to the specified output directory
to_input_files = function(input_df, out_dir) {
    #   SNP input
    input_df |>
        select(SNP, CHR, BP) |>
        write_tsv(file.path(out_dir, 'SNPs.tsv'))

    #   P-value input
    input_df |>
        select(SNP, P, N) |>
        write_tsv(file.path(out_dir, 'p_values.tsv'))
}

################################################################################
#   Main
################################################################################

#-------------------------------------------------------------------------------
#   MDD
#-------------------------------------------------------------------------------

chain = import.chain(chain_path)

fread(mdd_path, skip = 21) |>
    as_tibble() |>
    dplyr::rename(
        SNP = `SNP ID`, P = `P-value`, CHR = `Chr ID`, BP = `Chr Position`
    ) |>
    mutate(N = mdd_N) |>
    lift_df(chain) |>
    to_input_files(mdd_out_dir)

#-------------------------------------------------------------------------------
#   MDD2019
#-------------------------------------------------------------------------------

fread(mdd2019_path) |>
    as_tibble() |>
    dplyr::rename(SNP = MarkerName) |>
    mutate(N = mdd2019_N) |>
    to_input_files(mdd2019_out_dir)

#-------------------------------------------------------------------------------
#   AUD
#-------------------------------------------------------------------------------

read_tsv(aud_path, show_col_types = FALSE) |>
    dplyr::rename(
        SNP = SNP_ID, CHR = Chromsome, BP = Position, P = PValue, N = SampleSize
    ) |>
    to_input_files(aud_out_dir)

#-------------------------------------------------------------------------------
#   CUD
#-------------------------------------------------------------------------------

cud_df = read_tsv(cud_path, show_col_types = FALSE) |>
    to_input_files(cud_out_dir)

#-------------------------------------------------------------------------------
#   ext_cannabis
#-------------------------------------------------------------------------------

read_tsv(ext_cannabis_path, show_col_types = FALSE) |>
    dplyr::rename(BP = POS, P = PVAL) |>
    lift_df(chain) |>
    to_input_files(ext_cannabis_out_dir)

#-------------------------------------------------------------------------------
#   Lifetime Cannabis
#-------------------------------------------------------------------------------

read_tsv(
        lifetime_cannabis_path, col_types = c('ccidi')
    ) |>
    #   There are some missing chromosomes and X should be a number
    filter(!is.na(Chr)) |>
    mutate(Chr = as.integer(ifelse(Chr == 'X', 23, Chr))) |>
    dplyr::rename(CHR = Chr, BP = Bp) |>
    to_input_files(lifetime_cannabis_out_dir)

#-------------------------------------------------------------------------------
#   OUD
#-------------------------------------------------------------------------------

read_tsv(oud_path, show_col_types = FALSE) |>
    dplyr::rename(
        SNP = SNP_ID, CHR = Chrosome, BP = Position, P = PValue, N = Effective_N
    ) |>
    to_input_files(oud_out_dir)

#-------------------------------------------------------------------------------
#   SUD
#-------------------------------------------------------------------------------

read_tsv(sud_path, show_col_types = FALSE) |>
    dplyr::rename(CHR = Chr) |>
    #   A small number of rows had parsing issues; just drop them
    filter(!is.na(CHR)) |>
    to_input_files(sud_out_dir)

session_info()
