#   Annotate DAR peaks with nearest genes and do GO by cell type on those genes.
#   Not technically TF-related

library(tidyverse)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # or species-appropriate
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg38)
library(clusterProfiler)
library(here)
library(sessioninfo)
library(GenomicRanges)
library(Signac)
library(rtracklayer)

peak_path = here(
    'processed-data', '06_peak_calling', '21_overlaping_FDRscores_TopHeatmap',
    'all_peaks_categorized.csv.gz'
)
plot_dir = here('plots', '11_TF', '03_GO_nearest_gene')
go_num_terms = 5
gtf_path = '/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz'
promoter_window = 2000  # +/- around TSS

dir.create(plot_dir, showWarnings = FALSE)

#   Load in DARs
peak_df = read_csv(peak_path, show_col_types = FALSE) |>
    filter(grepl('DAR', category)) |>
    distinct(peak_id, cell_type, .keep_all = TRUE)

peak_gr = peak_df |>
    pull(peak_id) |>
    unique() |>
    StringToGRanges(sep = c("-", "-"))

gtf = import(gtf_path) |>
    as.data.frame() |>
    as_tibble() |>
    filter(type == "gene") |>
    select(gene_id, gene_name)

#   Add nearest gene + basic annotation
ann = annotatePeak(
        peak_gr, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
        tssRegion = c(-1 * promoter_window, promoter_window),
        annoDb = "org.Hs.eg.db"
    ) |>
    as.data.frame() |>
    as_tibble() |>
    dplyr::rename(nearest_gene_name = SYMBOL) |>
    mutate(
        peak_id = paste(seqnames, start, end, sep = "-"),
        nearest_gene_id = gtf$gene_id[match(nearest_gene_name, gtf$gene_name)]
    ) |>
    select(peak_id, annotation, nearest_gene_id, nearest_gene_name)

peak_df = peak_df |>
    left_join(ann, by = "peak_id") |>
    #   For GO, an empirically linked gene is stronger evidence than using the
    #   nearest gene. Use whichever is available though
    mutate(gene_for_go = coalesce(link_gene_id, nearest_gene_id)) |>
    #   It's still a bit unclear why a considerable fraction of nearest genes
    #   don't have Ensembl IDs or symbols in the GTF. We'll only consider
    #   genes in the GTF for GO
    filter(!is.na(gene_for_go), gene_for_go %in% gtf$gene_id)

ego_df_list = list()
for (this_cell_type in unique(peak_df$cell_type)) {
    gene_set = peak_df |>
        filter(cell_type == this_cell_type) |>
        pull(gene_for_go) |>
        unique()

    #   Note the universe here-- we're constraining nearest genes to those in
    #   the GTF. This gives us a larger set actually than those genes expressed
    #   in the RNA assay (from which linked peaks were derived), which seems a
    #   favorable in this particular case since there are very few linked peaks
    #   relative to nearest genes
    ego = enrichGO(
        gene = gene_set, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP",
        universe = gtf$gene_id, pAdjustMethod= "BH", pvalueCutoff = 1,
        qvalueCutoff = 0.05
    )
    ego_df_list[[this_cell_type]] = ego@result |>
        as_tibble() |>
        mutate(cell_type = this_cell_type)
}

#   Custom dot plots by cell type
p = bind_rows(ego_df_list) |>
    group_by(cell_type) |>
    slice_min(p.adjust, n = go_num_terms) |>
    ungroup() |>
    mutate(
        gene_ratio = Count / as.integer(str_extract(GeneRatio, "(?<=/)[0-9]+")),
        log_fdr = -log10(p.adjust)
    ) |>
    ggplot(
            aes(
                x = cell_type, y = Description, color = log_fdr,
                size = gene_ratio
            )
        ) +
        geom_point() +
        scale_color_gradient(low = "red", high = "blue") +
        facet_wrap(~cell_type, scales = "free", ncol = 3) +
        theme_bw(base_size = 5)
pdf(file.path(plot_dir, 'GO_nearest_gene.pdf'))
print(p)
dev.off()

session_info()
