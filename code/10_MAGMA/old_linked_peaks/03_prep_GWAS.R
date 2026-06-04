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
input_paths = c(
    MDD = '/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/13_MAGMA/GWAS/MDD/MDD.phs001672.pha005122.txt',
    MDD2019 = '/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/13_MAGMA/GWAS/mdd2019edinburgh/PGC_UKB_depression_genome-wide.txt',
    AUD = here('processed-data', '10_MAGMA', 'AUD', 'AUDSumstats.txt.gz'),
    compulsive = here(
        'processed-data', '10_MAGMA', 'compulsive',
        'F1_CompulsiveDisorders_2025.txt.gz'
    ),
    CUD = here('processed-data', '10_MAGMA', 'CUD', 'CUDEASumstats.txt.gz'),
    ext_cannabis = here(
        'processed-data', '10_MAGMA', 'ext_cannabis', 'ExtSumstats.txt.gz'
    ),
    internalizing = here(
        'processed-data', '10_MAGMA', 'internalizing',
        'F4_Internalizing_2025.txt.gz'
    ),
    lifetime_cannabis = here(
        'processed-data', '10_MAGMA', 'lifetime_cannabis',
        'LifetimeCannabisICC.txt.gz'
    ),
    neurodev = here(
        'processed-data', '10_MAGMA', 'neurodev',
        'F3_Neurodevelopmental_2025.txt.gz'
    ),
    OUD = here('processed-data', '10_MAGMA', 'OUD', 'OUDEASumstats.txt.gz'),
    p_factor = here(
        'processed-data', '10_MAGMA', 'p_factor', 'PFactor_2025.txt.gz'
    ),
    SCZ_BPD = here(
        'processed-data', '10_MAGMA', 'SCZ_BPD',
        'F2_SchizophreniaBipolar_2025.txt.gz'
    ),
    SUD2 = here('processed-data', '10_MAGMA', 'SUD2', 'SUDEA.txt.gz'),
    SUD3 = here(
        'processed-data', '10_MAGMA', 'SUD3', 'F5_SubstanceUse_2025.txt.gz'
    )
)
out_dirs = here('processed-data', '10_MAGMA', names(input_paths))
names(out_dirs) = names(input_paths)
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

chain = import.chain(chain_path)

#-------------------------------------------------------------------------------
#   MDD
#-------------------------------------------------------------------------------

fread(input_paths['MDD'], skip = 21) |>
    as_tibble() |>
    dplyr::rename(
        SNP = `SNP ID`, P = `P-value`, CHR = `Chr ID`, BP = `Chr Position`
    ) |>
    mutate(N = mdd_N) |>
    lift_df(chain) |>
    to_input_files(out_dirs['MDD'])

#-------------------------------------------------------------------------------
#   MDD2019
#-------------------------------------------------------------------------------

#   P-value file only
fread(input_paths['MDD2019']) |>
    as_tibble() |>
    dplyr::rename(SNP = MarkerName) |>
    mutate(N = mdd2019_N) |>
    select(SNP, P, N) |>
    write_tsv(file.path(out_dirs['MDD2019'], 'p_values.tsv'))

#-------------------------------------------------------------------------------
#   AUD
#-------------------------------------------------------------------------------

read_tsv(input_paths['AUD'], show_col_types = FALSE) |>
    dplyr::rename(
        SNP = SNP_ID, CHR = Chromsome, BP = Position, P = PValue, N = SampleSize
    ) |>
    to_input_files(out_dirs['AUD'])

#-------------------------------------------------------------------------------
#   compulsive
#-------------------------------------------------------------------------------

read_tsv(input_paths['compulsive'], show_col_types = FALSE) |>
    to_input_files(out_dirs['compulsive'])

#-------------------------------------------------------------------------------
#   CUD
#-------------------------------------------------------------------------------

read_tsv(input_paths['CUD'], show_col_types = FALSE) |>
    to_input_files(out_dirs['CUD'])

#-------------------------------------------------------------------------------
#   ext_cannabis
#-------------------------------------------------------------------------------

read_tsv(input_paths['ext_cannabis'], show_col_types = FALSE) |>
    dplyr::rename(BP = POS, P = PVAL) |>
    lift_df(chain) |>
    to_input_files(out_dirs['ext_cannabis'])

#-------------------------------------------------------------------------------
#   internalizing
#-------------------------------------------------------------------------------

read_tsv(input_paths['internalizing'], show_col_types = FALSE) |>
    to_input_files(out_dirs['internalizing'])

#-------------------------------------------------------------------------------
#   Lifetime Cannabis
#-------------------------------------------------------------------------------

read_tsv(
        input_paths['lifetime_cannabis'], col_types = c('ccidi')
    ) |>
    #   There are some missing chromosomes and X should be a number
    filter(!is.na(Chr)) |>
    mutate(Chr = as.integer(ifelse(Chr == 'X', 23, Chr))) |>
    dplyr::rename(CHR = Chr, BP = Bp) |>
    to_input_files(out_dirs['lifetime_cannabis'])

#-------------------------------------------------------------------------------
#   neurodev
#-------------------------------------------------------------------------------

read_tsv(input_paths['neurodev'], show_col_types = FALSE) |>
    to_input_files(out_dirs['neurodev'])

#-------------------------------------------------------------------------------
#   OUD
#-------------------------------------------------------------------------------

read_tsv(input_paths['OUD'], show_col_types = FALSE) |>
    dplyr::rename(
        SNP = SNP_ID, CHR = Chrosome, BP = Position, P = PValue, N = Effective_N
    ) |>
    to_input_files(out_dirs['OUD'])

#-------------------------------------------------------------------------------
#   p_factor
#-------------------------------------------------------------------------------

read_tsv(input_paths['p_factor'], show_col_types = FALSE) |>
    to_input_files(out_dirs['p_factor'])

#-------------------------------------------------------------------------------
#   SCZ_BPD
#-------------------------------------------------------------------------------

read_tsv(input_paths['SCZ_BPD'], show_col_types = FALSE) |>
    to_input_files(out_dirs['SCZ_BPD'])

#-------------------------------------------------------------------------------
#   SUD2
#-------------------------------------------------------------------------------

read_tsv(input_paths['SUD2'], show_col_types = FALSE) |>
    dplyr::rename(CHR = Chr) |>
    #   A small number of rows had parsing issues; just drop them
    filter(!is.na(CHR)) |>
    to_input_files(out_dirs['SUD2'])

#-------------------------------------------------------------------------------
#   SUD3
#-------------------------------------------------------------------------------

read_tsv(input_paths['SUD3'], show_col_types = FALSE) |>
    to_input_files(out_dirs['SUD3'])

session_info()
