#   The MDD GWAS summary statistics are in hg38. Lift them over to hg19. Also,
#   to ensure all inputs are consistent across GWAS datasets, ensure that
#   p-value input files all have columns SNP, P, N (affects MDD and MDD2019)

library(here)
library(rtracklayer)
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(sessioninfo)

chain_path = here('processed-data', '10_MAGMA', 'hg38ToHg19.over.chain')
mdd_path = '/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/13_MAGMA/GWAS/MDD/MDD.phs001672.pha005122.txt'
mdd2019_path = '/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/13_MAGMA/GWAS/mdd2019edinburgh/PGC_UKB_depression_genome-wide.txt'
mdd_out_dir = here('processed-data', '10_MAGMA', 'MDD')
mdd2019_out_dir = here('processed-data', '10_MAGMA', 'MDD2019')
mdd_N = 1154267
mdd2019_N = 807553

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

session_info()
