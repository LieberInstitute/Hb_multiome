#   The MDD GWAS summary statistics are in hg38. Lift them over to hg19

library(here)
library(rtracklayer)
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(sessioninfo)

chain_path = here('processed-data', '10_MAGMA', 'hg38ToHg19.over.chain')
mdd_path = '/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/13_MAGMA/GWAS/MDD/MDD.phs001672.pha005122.txt'
out_path = here('processed-data', '10_MAGMA', 'MDD', 'SNPs.tsv')

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
    write_tsv(out_path)

session_info()
