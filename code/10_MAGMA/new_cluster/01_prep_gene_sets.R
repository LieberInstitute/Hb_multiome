library(tidyverse)
library(here)
library(sessioninfo)

mean_ratio_threshold <- 1.05
max_num_genes <- 200

a <- readRDS("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/11_link_prep/04_registration_wrapper/model_results.rds")

out_dir= here('processed-data', '10_MAGMA', 'gene_sets')
dir.create(out_dir, showWarnings = FALSE)

# a$enrichment
enrich_df <- a$enrichment %>%
  as.data.frame()

# find logFC columns
logfc_cols <- grep("^logFC_", names(enrich_df), value = TRUE)

# change logFC columns to long format and calculate MeanRatio, then filter by MeanRatio threshold and select top genes for each cell type
marker_stats <- enrich_df %>%
  select(ensembl, gene, all_of(logfc_cols)) %>%
  pivot_longer(
    cols = all_of(logfc_cols),
    names_to = "cell_type",
    values_to = "logFC"
  ) %>%
  mutate(
    cell_type = sub("^logFC_", "", cell_type),
    MeanRatio = 2^logFC
  ) %>%
  filter(!is.na(MeanRatio), is.finite(MeanRatio)) %>%
  filter(MeanRatio > mean_ratio_threshold) %>%
  group_by(cell_type) %>%
  arrange(desc(MeanRatio), desc(logFC), .by_group = TRUE) %>%
  slice_head(n = max_num_genes) %>%
  ungroup() %>%
  transmute(
    set_id = cell_type,
    gene_id = ensembl,
    gene_name = gene,
    logFC = logFC,
    MeanRatio = MeanRatio
  )

# the number of genes in each set
table(marker_stats$set_id)

# save
readr::write_tsv(marker_stats, file.path(out_dir, "celltype_mean_ratio_top200.tsv"))

session_info()
