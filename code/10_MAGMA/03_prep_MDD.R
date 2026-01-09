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
mdd_out_dir = here('processed-data', '10_MAGMA', 'MDD')
mdd2019_out_dir = here('processed-data', '10_MAGMA', 'MDD2019')
aud_out_dir = here('processed-data', '10_MAGMA', 'AUD')
cud_out_dir = here('processed-data', '10_MAGMA', 'CUD')
ext_cannabis_out_dir = here('processed-data', '10_MAGMA', 'ext_cannabis')
lifetime_cannabis_out_dir = here(
    'processed-data', '10_MAGMA', 'lifetime_cannabis'
)
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
#   MDD
################################################################################

#-------------------------------------------------------------------------------
#   SNP input
#-------------------------------------------------------------------------------

mdd_df = fread(mdd_path, skip = 21) |>
    as_tibble()

mdd_gr = GRanges(
    seqnames = paste0('chr', mdd_df$`Chr ID`),
    ranges = IRanges(start = mdd_df$`Chr Position`, width = 1),
    SNP = mdd_df$`SNP ID`
)

chain = import.chain(chain_path)
lifted = unlist(liftOver(mdd_gr, chain))
message(
    sprintf(
        'Lifted over %d of %d SNPs to hg19',
        length(lifted), length(mdd_gr)
    )
)

tibble(
        SNP = lifted$SNP,
        CHR = gsub('chr', '', as.character(seqnames(lifted))),
        BP = start(lifted)
    ) |>
    write_tsv(file.path(mdd_out_dir, 'SNPs.tsv'))

#-------------------------------------------------------------------------------
#   P-value input
#-------------------------------------------------------------------------------

mdd_df |>
    dplyr::rename(SNP = `SNP ID`, P = `P-value`) |>
    mutate(N = mdd_N) |>
    select(SNP, P, N) |>
    filter(SNP %in% lifted$SNP) |>
    write_tsv(file.path(mdd_out_dir, 'p_values.tsv'))

################################################################################
#   MDD2019
################################################################################

#   P-value input only
mdd_df = fread(mdd2019_path) |>
    as_tibble() |>
    dplyr::rename(SNP = MarkerName) |>
    mutate(N = mdd2019_N) |>
    select(SNP, P, N) |>
    write_tsv(file.path(mdd2019_out_dir, 'p_values.tsv'))

################################################################################
#   AUD
################################################################################

aud_df = read_tsv(aud_path, show_col_types = FALSE) |>
    dplyr::rename(
        SNP = SNP_ID, CHR = Chromsome, BP = Position, P = PValue, N = SampleSize
    )

#   SNP input
aud_df |>
    select(SNP, CHR, BP) |>
    write_tsv(file.path(aud_out_dir, 'SNPs.tsv'))

#   P-value input
aud_df |>
    select(SNP, P, N) |>
    write_tsv(file.path(aud_out_dir, 'p_values.tsv'))

################################################################################
#   CUD
################################################################################

cud_df = read_tsv(cud_path, show_col_types = FALSE)

#   SNP input
cud_df |>
    select(SNP, CHR, BP) |>
    write_tsv(file.path(cud_out_dir, 'SNPs.tsv'))

#   P-value input
cud_df |>
    select(SNP, P, N) |>
    write_tsv(file.path(cud_out_dir, 'p_values.tsv'))

################################################################################
#   ext_cannabis
################################################################################

ext_cannabis_df = read_tsv(ext_cannabis_path, show_col_types = FALSE) |>
    dplyr::rename(BP = POS, P = PVAL)

#   It's hg38, so we lift to hg19
ext_cannabis_gr = GRanges(
    seqnames = paste0('chr', ext_cannabis_df$CHR),
    ranges = IRanges(start = ext_cannabis_df$BP, width = 1),
    SNP = ext_cannabis_df$SNP
)

#   Lift over
lifted = unlist(liftOver(ext_cannabis_gr, chain))
message(
    sprintf(
        'Lifted over %d of %d SNPs to hg19',
        length(lifted), length(ext_cannabis_gr)
    )
)

#   Convert back into a tibble with all columns
ext_cannabis_df = tibble(
        SNP = lifted$SNP,
        CHR = gsub('chr', '', as.character(seqnames(lifted))),
        BP = start(lifted)
    ) |>
    left_join(
        ext_cannabis_df |>
            select(SNP, P, N),
        by = 'SNP'
    )
stopifnot(!any(is.na(ext_cannabis_df)))

#   SNP input
ext_cannabis_df |>
    select(SNP, CHR, BP) |>
    write_tsv(file.path(ext_cannabis_out_dir, 'SNPs.tsv'))

#   P-value input
ext_cannabis_df |>
    select(SNP, P, N) |>
    write_tsv(file.path(ext_cannabis_out_dir, 'p_values.tsv'))

################################################################################
#   Lifetime Cannabis
################################################################################

lifetime_cannabis_df = read_tsv(
        lifetime_cannabis_path, col_types = c('ccidi')
    ) |>
    #   There are some missing chromosomes and X should be a number
    filter(!is.na(Chr)) |>
    mutate(Chr = as.integer(ifelse(Chr == 'X', 23, Chr))) |>
    dplyr::rename(CHR = Chr, BP = Bp)

#   SNP input
lifetime_cannabis_df |>
    select(SNP, CHR, BP) |>
    write_tsv(file.path(lifetime_cannabis_out_dir, 'SNPs.tsv'))

#   P-value input
lifetime_cannabis_df |>
    select(SNP, P, N) |>
    write_tsv(file.path(lifetime_cannabis_out_dir, 'p_values.tsv'))

session_info()
