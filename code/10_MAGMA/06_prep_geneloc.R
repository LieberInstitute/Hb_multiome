#   Take the GTF used for the multiome data and lift it over to hg19. Everything
#   must be in hg19 for MAGMA (hg38 seems better, but it seems extremely
#   complicated to get the 1000 Genomes European plink files needed for MAGMA in
#   hg38)
library(rtracklayer)
library(here)
library(sessioninfo)

reference_gtf = '/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz'
chain_path = here('processed-data', '10_MAGMA', 'hg38ToHg19.over.chain')
out_path = here('processed-data', '10_MAGMA', 'hg19_gene_loc.tsv')

gtf = import(reference_gtf)
gtf = gtf[gtf$type == "gene"]

chain = import.chain(chain_path)
liftOver(gtf, chain) |>
    unlist() |>
    as.data.frame() |>
    as_tibble() |>
    mutate(chr = gsub("chr", "", seqnames)) |>
    select(gene_id, chr, start, end) |>
    write_tsv(out_path, col_names = FALSE)

session_info()
