########################################################################
## explore and filter strong peak-gene links for evaluation 
## - Make several visualization to evaluate Peak scores
## - Make table with several confidence Peak scores
##
## Authors. CSC
## Date. Aug 11, 2025
## Recommended resources on interactive mode: srun --pty --mem=30GB --x11 bash
## Note. Seurat objects were created with module load conda_R/4.4.x
########################################################################

## I use EnsDb.Hsapiens.v86 for extracting gene names, positions, TSSs, chr locations, etc.
library("EnsDb.Hsapiens.v86")           # Gene annotation (GTF-style)
library("ggplot2")
library("patchwork")
library("tidyverse")
library("dplyr")
library("here")


## read input arguments
args = commandArgs(trailingOnly = TRUE)
p_met <- args[2]
w_size <- args[4]
# 0: pearson, 1e5
# 1: pearson, 5e4
# 2: spearman, 1e5
# 3: spearman, 2.5e4

# for testing
# p_met = "spearman"
# w_size = 2.5e4

## p_met:
# pearson -> peak-scores<0.2: likely due scATAC counts are ultra‑sparse; scRNA is zero‑inflated. Pearson r’s of 0.05–0.2 are common even for real links
# spearman -> as enhancer → gene relationships aren’t strictly linear; Pearson seems to underestimates. I will try spearman, more robust to nonlinearity/zeros

if (length(p_met) && length(w_size)) {
    message(
        "Processing job for peak-method:\n",
        p_met,
        "\nWindow-size\n",
        w_size
    )
    ## Use numeric comparison first, then assign string labels
    w_size_label <- case_when(
        isTRUE(all.equal(w_size, 25000))  ~ "2.5e4",
        isTRUE(all.equal(w_size, 50000))  ~ "5e4",
        isTRUE(all.equal(w_size, 100000)) ~ "1e5",
        TRUE                              ~ "00"
    )
    f_sufix <- paste0(".", p_met, ".", format(w_size, scientific = TRUE), ".cells_filtered_2perc")
    message("Processing: ", f_sufix)
} else {
    message("Input arguments missed")
    stop()
}


# Check/create directories

## clusters renamed for Spatial-Registration on Visium project
input_cvsDir <- here(
    "processed-data",
    "06_peak_calling",
    "00_link_peaks"
)
plotDir <- here(
  "plots",
  "06_peak_calling",
  "06_exploratory_peak_scores"
)
csvDir <- here(
    "processed-data",
    "06_peak_calling",
    "06_exploratory_peak_scores"
)


## Check directories
if (!dir.exists(plotDir)) {
  dir.create(plotDir)
}
if (!dir.exists(csvDir)) {
    dir.create(csvDir)
}

# load link peak-gene csv
gene_peaks_csv <- here(input_cvsDir, paste0("all_peak_gene_links", f_sufix, ".csv"))
if (file.exists(gene_peaks_csv)) {
    link_df <- read.csv(gene_peaks_csv)
    message("File loaded!")
} else {
    stop(paste("File not found:", gene_peaks_csv))
}

## for testing: ================================================================
#link_df <- read.csv(file = here(input_cvsDir, "all_peak_gene_links.spearman.5e4_test.csv"))
## for testing: ================================================================

colnames(link_df)
link_df <- link_df |>
    mutate(
        gene     = trimws(as.character(gene)),
        seqnames = as.character(seqnames),
        start    = as.numeric(start),
        end      = as.numeric(end)
    )
head(link_df)

message("Computing distance between peaks and TSS ...")

## Build TSS (strand-aware) GRanges table and 
#  compute the distance between each peak and its linked gene's TSS
gene_coords <- genes(EnsDb.Hsapiens.v86)
# tss_coords <- resize(gene_coords, width = 1, fix = "start")
# tss_coords <- keepStandardChromosomes(tss_coords, pruning.mode = "coarse")
# head(tss_coords)

# Build TSS table (strand-aware)
tss_coords  <- promoters(gene_coords, upstream = 0, downstream = 1) %>%   # 1bp TSS, respects strand
    keepStandardChromosomes(pruning.mode = "coarse")

# match UCSC-style peaks (chr1, chr2, etc.)
seqlevelsStyle(tss_coords) <- "UCSC"

tss_raw <- as.data.frame(tss_coords)
has_biotype <- "gene_biotype" %in% colnames(tss_raw)

# Prefer protein_coding, then take first per (gene_name, chr) - this avoid 1:many associations
tss_df <- tss_raw %>%
    mutate(gene_biotype = if (has_biotype) gene_biotype else NA_character_) %>%
    arrange(desc(gene_biotype == "protein_coding")) %>%  # prefer protein-coding where available
    group_by(gene_name, seqnames) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    transmute(
        gene_name,
        seqnames = as.character(seqnames),
        tss      = as.numeric(start),   # rename for clarity
        gene_strand = as.character(strand),
        gene_id
    )
head(tss_df)
#     gene_name seqnames       tss gene_strand gene_id        
#     <chr>     <chr>        <dbl> <chr>       <chr>          
# 1 5S_rRNA   chr1     143439605 +           ENSG00000252830
# 2 5S_rRNA   chr11    102057854 +           ENSG00000274097
# 3 5S_rRNA   chr17     37940790 -           ENSG00000277488
nrow(tss_df)
# [1] 56747

## check duplicates
dup_pairs <- tss_df %>%
    count(gene_name, seqnames, name = "n") %>%
    filter(n > 1)

if (nrow(dup_pairs) > 0) {
    print(head(dup_pairs, 10))
    stop("Non-unique (gene_name, seqnames) in tss_df: ", nrow(dup_pairs), " duplicates.")
}

# Join by gene + chromosome to avoid many-to-many 
colnames(link_df)
colnames(tss_df)
link_df2 <- link_df %>%
    left_join(tss_df, by = c("gene" = "gene_name", "seqnames" = "seqnames"))
head(link_df2, n = 3)
#     seqnames  start     end width strand      score  gene               peak
# 1     chr1 921198 1001138 79941      * 0.06070581 ISG15 chr1-920766-921629
# 2     chr1 960688 1000172 39485      * 0.07694616  HES4 chr1-960318-961058
# 3     chr1 960688 1001138 40451      * 0.06605429 ISG15 chr1-960318-961058
#     zscore     pvalue     tss gene_strand         gene_id
# 1 2.052715 0.02005013 1001138           + ENSG00000187608
# 2 1.758599 0.03932286 1000172           - ENSG00000188290
# 3 2.016967 0.02184947 1001138           + ENSG00000187608

# drop rows with no TSS match
n_before <- nrow(link_df2)
link_df2 <- link_df2 %>% filter(!is.na(tss))
message("Dropped ", n_before - nrow(link_df2), " rows with no TSS match.")
# Dropped 0 rows with no TSS match.
nrow(link_df2)

message("Distance between peaks and TSS completed!")

message("Computing Peak center and distance to TSS ...")

link_df2 <- link_df2 %>%
    mutate(
        peak_center       = (start + end) / 2,
        distance          = abs(peak_center - tss),
        signed_distance   = peak_center - tss,                       # genomic sign
        signed_by_strand  = ifelse(gene_strand == "-", -signed_distance, signed_distance),
        distance_kb       = distance / 1000
    )
head(link_df2, n=2)
# seqnames  start     end width strand      score  gene               peak
# 1     chr1 921198 1001138 79941      * 0.06070581 ISG15 chr1-920766-921629
# 2     chr1 960688 1000172 39485      * 0.07694616  HES4 chr1-960318-961058
#   zscore     pvalue     tss gene_strand         gene_id peak_center distance
# 1 2.052715 0.02005013 1001138           + ENSG00000187608      961168    39970
# 2 1.758599 0.03932286 1000172           - ENSG00000188290      980430    19742
# signed_distance signed_by_strand distance_kb
# 1          -39970           -39970      39.970
# 2          -19742            19742      19.742
message("Link gene-peak scores with TSS:")
summary(link_df2)
#sum(link_df2$distance > 1e5)  # should be ~0 if you used LinkPeaks(..., distance=1e5)
table(link_df2$gene_strand, useNA = "ifany")
    # -    + 
    # 2714 2786 

message("Peak center and distance to TSS completed!")


#===============================================================================

message("Building plots ...")

## Histogram TSS Scores
pdf(file = here(plotDir, 
                paste0("histogram_scores", f_sufix, ".pdf")), 
    width = 7, height = 5)

hist(link_df2$distance / 1000, breaks = 100,
     main = "Distance from Peaks to TSS",
     xlab = "Distance (kb)",
     col = "lightblue")
dev.off()

## Correlation vs Distance with smoothing
g1 <- ggplot(link_df2, aes(distance_kb, score)) +
    geom_point(alpha = 0.3, color = "steelblue") +
    geom_hline(yintercept = 0.2, linetype = "dashed", color = "red") +
    labs(
        x = "Distance from TSS (kb)",
        y = paste("Correlation Score", p_met),
        title = "Peak-Gene Correlation vs. Distance"
    ) + geom_smooth(method = "loess", se = FALSE, color = "darkred") +
    theme_minimal()

ggsave(here(plotDir, 
            paste0("distribution_scores", f_sufix, ".pdf")),
            g1, width = 8, height = 5)


#===============================================================================
# adding exploratory scores
# define high-confidence
# High: score ≥ 0.30 & FDR < 0.05
# Moderate: 0.20 ≤ score < 0.30 & FDR < 0.10
# Exploratory: 0.10 ≤ score < 0.20 & FDR < 0.10 (treat as hypotheses)

colnames(link_df2)
## add adjusted p-value using the Benjamini–Hochberg correction
link_df2 <- link_df2 %>%
    mutate(FDR = p.adjust(pvalue, method = "BH"))

# set tiers due we have high confidente peaks < 0.2 
link_df2 <- link_df2 %>%
    mutate(tier = case_when(
        score >= 0.30 & FDR < 0.05 ~ "High (>=0.30, FDR<0.05)",
        score >= 0.20 & FDR < 0.10 ~ "Moderate (0.20–0.30, FDR<0.10)",
        score >= 0.10 & FDR < 0.10 ~ "Exploratory (0.10–0.20, FDR<0.10)",
        TRUE ~ "Discarded"
    ))
# use plain ASCII hyphens
link_df$tier <- gsub("\u2013", "-", link_df2$tier)
head(link_df2)
table(link_df2$tier)

# quick view by distance (kb)
# Keep all data, no filtering of "Discarded" on the plot for visualization purposes
df_plot <- link_df2  

# Count total peaks and how many are below 0.1 to plot on discarted zone
count_below_01 <- sum(df_plot$score < 0.1, na.rm = TRUE)
count_below_02 <- sum((df_plot$score < 0.2 & df_plot$score > 0.1), na.rm = TRUE)
count_below_03 <- sum((df_plot$score < 0.3 & df_plot$score > 0.2), na.rm = TRUE)

g1 <- ggplot(df_plot, aes(x = distance/1000, y = score, color = tier)) +
    geom_point(alpha = 0.5, size = 0.8) +
    # trend over ALL tested links
    geom_smooth(
        data = df_plot,
        aes(x = distance_kb, y = score),
        method = "loess", se = FALSE, span = 0.8, color = "black", linewidth = 0.9) +
    # Threshold lines
    geom_hline(yintercept = 0.3, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 0.2, linetype = "dashed", color = "orange") +
    geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey50") +
    # Labels for thresholds
    annotate("text", x = max(df_plot$distance/1000)*1.02, y = 0.3, 
             label = paste("0.3 (", count_below_03, " peaks)"), hjust = 0.8, vjust = -0.5, color = "red") +
    annotate("text", x = max(df_plot$distance/1000)*1.02, y = 0.2, 
             label = paste("0.2 (", count_below_02, " peaks)"), hjust = 0.8, vjust = -0.5, color = "orange") +
    annotate("text", x = max(df_plot$distance/1000)*1.02, y = 0.1, 
             label = paste("0.1 (", count_below_01, " peaks)"), hjust = 0.8, vjust = -0.5, color = "grey50") +
    # Custom legend with count
    scale_color_manual(
        values = c(
            "High (>=0.30, FDR<0.05)"          = "#b2182b",
            "Moderate (0.20–0.30, FDR<0.10)"   = "#ef8a62",
            "Exploratory (0.10–0.20, FDR<0.10)" = "#67a9cf",
            "Discarded"                        = "grey80"
        ),
        name = paste0("Tier (Count < 0.1: ", count_below_01, ")")
    )  +
    labs(
        x = "Distance from TSS (kb)",
        y = "Correlation score",
        title = "Peak–gene links by tier",
        subtitle = paste(p_met, 
                         format(w_size, scientific = TRUE), "filtered genes < 2%")
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

ggsave(here(plotDir, 
            paste0("exploratory_scores_high_confidence", f_sufix, ".pdf")),
       g1, width = 8, height = 5,
       device = cairo_pdf)

#===============================================================================


## filtered peaks
write.csv(link_df2,
          file = here(csvDir, paste0("peak_gene_links_with_TSS_and_CC", f_sufix, ".csv")),
          row.names = FALSE)

message("Link Gene-Peak table with TSS distances and CC scores saved!")


#===============================================================================

## Check number of linked peaks per gene and viceverce
# Number of linked peaks per gene
peaks_per_gene <- link_df2 %>%
    count(gene, name = "n_peaks") %>%
    arrange(desc(n_peaks))

g1 <- ggplot(peaks_per_gene, aes(x = n_peaks)) +
    geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    labs(
        title = paste(p_met, format(w_size, scientific = TRUE), "filtered genes < 2%"),
        subtitle = "Peaks per gene",
        x = "Number of linked peaks per gene",
        y = "Number of genes"
    ) +
    theme_minimal()

# Number of linked genes per peak
genes_per_peak <- link_df2 %>%
    count(peak, name = "n_genes") %>%
    arrange(desc(n_genes))

g2 <- ggplot(genes_per_peak, aes(x = n_genes)) +
    geom_histogram(binwidth = 1, fill = "firebrick", color = "white") +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    labs(
        subtitle = "Genes per peak",
        x = "Number of linked genes per peak",
        y = "Number of peaks"
    ) +
    theme_minimal()

combined_plot <- g1 + g2
combined_plot

ggsave(here(plotDir, 
            paste0("histograms_link_peak_genes_peaks", f_sufix, ".pdf")),
       combined_plot, width = 8, height = 5)



# #===============================================================================
# ## Prepare Seurat for Heatmap
# 
# # Use Seurat with clusters renamed for Spatial-Registration on Visium project
# Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_HD.rds"
# seurat_name <- here(inputRDS_Dir, Seurat_base_name)
# SeuratOBJ <- readRDS(here(inputRDS_Dir, Seurat_base_name))
# levels(SeuratOBJ)
# 
# ## filter genes to those expressed in 3% of cells
# DefaultAssay(SeuratOBJ) <- "RNA"
# rna_counts <- GetAssayData(SeuratOBJ, assay="RNA", layer="data")
# length(rownames(rna_counts)) # [1] 36601
# #length(rownames(rna_counts)[Matrix::rowSums(rna_counts > 0)]) # 34738
# #length(rownames(rna_counts)[Matrix::rowSums(rna_counts > 0) > 0.02]) # 34738
# keep_genes <- rownames(rna_counts)[Matrix::rowSums(rna_counts > 0) > 0.02 * ncol(rna_counts)]
# length(keep_genes) # in count: [1] 14526
# 
# ## Set atac
# DefaultAssay(SeuratOBJ) <- "ATAC"
# class(SeuratOBJ[["ATAC"]])
# # make a readable base-name for plots
# Seurat_base_name <- str_extract(seurat_name, regex("C\\.\\w+"))
# Seurat_base_name <- sub("_renamed_visium$", "", Seurat_base_name)
# Seurat_base_name
# # C.leiden_lsi_r2
# 
# ##==============================================================================
# ## Make complex Heatmap
# 
# ## Prep: assign each peak to a cluster (by accessibility)
# ## average peaks × clusters
# avg_atac <- AggregateExpression(
#     SeuratOBJ, 
#     assays = "ATAC", 
#     group.by = "cluster_ann", 
#     layers = "counts"
# )$ATAC  # matrix: peaks x clusters
# 
# # for each peak, which cluster has highest average accessibility
# peak_cluster <- data.frame(
#     peak = rownames(avg_atac),
#     cluster_ann = colnames(avg_atac)[max.col(avg_atac, ties.method = "first")]
# )
# 
# # keep only links whose peaks are known in the assay
# links_cl <- link_df %>%
#     inner_join(peak_cluster, by = "peak") %>%
#     filter(!is.na(cluster_ann))
# 
# ## filter links of desired
# ## links_cl <- links_cl %>% filter(score > 0.3, pvalue < 0.05)
# 
# # bin definition (adjust if you like)
# bin_breaks <- c(0,1,2,3,5,10,Inf)
# bin_labels <- c("1","2","3","4–5","6–10",">10")
# 
# ## Distribution of number of linked peaks per gene, by cluster
# # Count linked peaks per (gene, cluster)
# peaks_per_gene <- links_cl %>%
#     group_by(cluster_ann, gene) %>%
#     summarise(n_peaks = n_distinct(peak), .groups = "drop") %>%
#     mutate(bin = cut(
#         n_peaks, breaks = bin_breaks,
#         labels = bin_labels, right = TRUE
#     ))
# 
# # Make a cluster × bin table
# dist_pg <- peaks_per_gene %>%
#     count(cluster_ann, bin, name = "n_genes") %>%
#     complete(cluster_ann, bin, fill = list(n_genes = 0)) %>%
#     group_by(cluster_ann) %>%
#     mutate(freq = n_genes / sum(n_genes)) %>%   # optional normalize by row
#     ungroup()
# 
# g_dp <- ggplot(dist_pg, aes(x = bin, y = cluster_ann, fill = freq)) +
#     geom_tile(color = "grey85") +
#     scale_fill_viridis_c(name = "Fraction of genes", option = "C") +
#     labs(
#         title = "Distribution of linked peaks per gene",
#         x = "Peaks per gene (bins)", y = "cluster_ann"
#     ) +
#     theme_minimal()
# 
# ## Distribution of number of linked genes per peak, by cluster
# # Count linked genes per (peak, cluster)
# genes_per_peak <- links_cl %>%
#     group_by(cluster_ann, peak) %>%
#     summarise(n_genes = n_distinct(gene), .groups = "drop") %>%
#     mutate(bin = cut(
#         n_genes, breaks = c(0,1,2,3,5,10,Inf),
#         labels = c("1","2","3","4-5","6-10",">10"), right = TRUE
#     ))
# 
# dist_gp <- genes_per_peak %>%
#     count(cluster_ann, bin, name = "n_peaks") %>%
#     complete(cluster_ann, bin, fill = list(n_peaks = 0)) %>%
#     group_by(cluster_ann) %>%
#     mutate(freq = n_peaks / sum(n_peaks)) %>%
#     ungroup()
# 
# g_dg <- ggplot(dist_gp, aes(x = bin, y = cluster_ann, fill = freq)) +
#     geom_tile(color = "grey85") +
#     scale_fill_viridis_c(name = "Fraction of peaks", option = "C") +
#     labs(
#         title = "Distribution of linked genes per peak",
#         x = "Genes per peak (bins)", y = "cluster_ann"
#     ) +
#     theme_minimal()
# 
# combined_plot <- g_dp + g_dg
# combined_plot
# 
# ggsave(here(plotDir, 
#             paste0("heatmaps_link_peak_genes_peaks", f_sufix, ".pdf")),
#        combined_plot, width = 12, height = 5)
# 
# 
# 
# ##  Top genes by total linked peaks, across clusters
# # Pick top N genes by total linked peaks (across all clusters)
# 
# N <- 50
# 
# top_genes <- links_cl %>%
#     group_by(gene) %>%
#     summarise(total_peaks = n_distinct(peak), .groups = "drop") %>%
#     slice_max(total_peaks, n = N) %>%
#     pull(gene)
# 
# mat_pg <- links_cl %>%
#     filter(gene %in% top_genes) %>%
#     group_by(cluster_ann, gene) %>%
#     summarise(n_peaks = n_distinct(peak), .groups = "drop") %>%
#     mutate(gene = factor(gene, levels = top_genes)) %>%
#     complete(cluster_ann, gene, fill = list(n_peaks = 0))
# 
# # Heatmap: genes (rows) × clusters (cols)
# g1 <- ggplot(mat_pg, aes(x = cluster_ann, y = gene, fill = n_peaks)) +
#     geom_tile() +
#     scale_fill_viridis_c(name = "# linked peaks") +
#     labs(
#         title = paste("Top", N, "genes: linked peaks per gene by cluster"),
#         x = "WNN cluster", y = "Gene"
#     ) +
#     theme_minimal() +
#     theme(axis.text.y = element_text(size = 7))
# 
# ggsave(here(plotDir, 
#             paste0("heatmap_Top50_genes_peaks", f_sufix, ".pdf")),
#        g1, width = 8, height = 5)
# 
# 
# 
# ## ====
# 
# library(ComplexHeatmap)
# library(circlize)
# 
# bin_breaks  <- c(0,1,2,3,5,10,Inf)
# bin_labels  <- c("1","2","3","4-5","6-10",">10")   # <-- ASCII hyphen
# 
# #1) Heatmap: Distribution of peaks per gene (binned), by cluster
# #Counts how many distinct peaks are linked to each gene within each cluster, bins those counts, and shows the fraction per cluster.
# 
# peaks_per_gene <- links_cl %>%
#     group_by(cluster_ann, gene) %>%
#     summarise(n_peaks = n_distinct(peak), .groups = "drop") %>%
#     mutate(bin = cut(n_peaks, breaks = bin_breaks, labels = bin_labels, right = TRUE))
# 
# dist_pg <- peaks_per_gene %>%
#     count(cluster_ann, bin, name = "n_genes") %>%
#     complete(cluster_ann, bin = factor(bin_labels, levels = bin_labels), fill = list(n_genes = 0)) %>%
#     group_by(cluster_ann) %>%
#     mutate(freq = n_genes / sum(n_genes)) %>%
#     ungroup()
# 
# # matrix: rows = clusters, cols = bins
# mat_pg <- dist_pg %>%
#     select(cluster_ann, bin, freq) %>%
#     pivot_wider(names_from = bin, values_from = freq) %>%
#     as.data.frame()
# 
# rownames(mat_pg) <- mat_pg$cluster_ann
# mat_pg$cluster_ann <- NULL
# mat_pg <- as.matrix(mat_pg)  # numeric matrix
# 
# # color scale
# col_fun_pg <- circlize::colorRamp2(c(0, max(mat_pg, na.rm = TRUE)), c("#F7FBFF", "#08306B"))
# 
# ht_pg <- Heatmap(
#     mat_pg,
#     name = "Frac genes",
#     col = col_fun_pg,
#     cluster_rows = TRUE,
#     cluster_columns = FALSE,
#     row_title = "WNN clusters",
#     column_title = "Peaks per gene (binned)",
#     heatmap_legend_param = list(direction = "horizontal"),
#     column_names_rot = 0
# )
# 
# # Save as PDF
# pdf(here(plotDir, "heatmap_frac_genes", f_sufix, ".pdf"),
#     width = 8, height = 6)
# draw(ht_pg)
# dev.off()
# 
# 
# ##2) Heatmap: Distribution of genes per peak (binned), by cluster
# ##Counts how many distinct genes each peak links to within each cluster, bins those counts, shows the fraction per cluster.
# 
# genes_per_peak <- links_cl %>%
#     group_by(cluster_ann, peak) %>%
#     summarise(n_genes = n_distinct(gene), .groups = "drop") %>%
#     mutate(bin = cut(n_genes, breaks = bin_breaks, labels = bin_labels, right = TRUE))
# 
# dist_gp <- genes_per_peak %>%
#     count(cluster_ann, bin, name = "n_peaks") %>%
#     complete(cluster_ann, bin = factor(bin_labels, levels = bin_labels), fill = list(n_peaks = 0)) %>%
#     group_by(cluster_ann) %>%
#     mutate(freq = n_peaks / sum(n_peaks)) %>%
#     ungroup()
# 
# mat_gp <- dist_gp %>%
#     select(cluster_ann, bin, freq) %>%
#     pivot_wider(names_from = bin, values_from = freq) %>%
#     as.data.frame()
# 
# rownames(mat_gp) <- mat_gp$cluster_ann
# mat_gp$cluster_ann <- NULL
# mat_gp <- as.matrix(mat_gp)
# 
# col_fun_gp <- circlize::colorRamp2(c(0, max(mat_gp, na.rm = TRUE)), c("#FFF5F0", "#7F0000"))
# 
# ht_gp <- Heatmap(
#     mat_gp,
#     name = "Frac peaks",
#     col = col_fun_gp,
#     cluster_rows = TRUE,
#     cluster_columns = FALSE,
#     row_title = "WNN clusters",
#     column_title = "Genes per peak (binned)",
#     heatmap_legend_param = list(direction = "horizontal"),
#     column_names_rot = 0
# )
# 
# # Save as PDF
# pdf(here(plotDir, "heatmap_frac_peaks", f_sufix, ".pdf"),
#     width = 8, height = 6)
# draw(ht_gp)
# dev.off()
# 
# #3) Heatmap: Top N genes by linked peaks across clusters
# #Rows = genes (top N by total linked peaks), columns = clusters, values = # linked peaks.
# 
# top_genes <- links_cl %>%
#     group_by(gene) %>%
#     summarise(total_peaks = n_distinct(peak), .groups = "drop") %>%
#     slice_max(total_peaks, n = N) %>%
#     pull(gene)
# 
# mat_top <- links_cl %>%
#     filter(gene %in% top_genes) %>%
#     group_by(cluster_ann, gene) %>%
#     summarise(n_peaks = n_distinct(peak), .groups = "drop") %>%
#     complete(cluster_ann, gene, fill = list(n_peaks = 0)) %>%
#     pivot_wider(names_from = cluster_ann, values_from = n_peaks) %>%
#     as.data.frame()
# 
# rownames(mat_top) <- mat_top$gene
# mat_top$gene <- NULL
# mat_top <- as.matrix(mat_top)
# 
# # Optional: row-wise z-score to show relative enrichment per gene
# mat_top_z <- t(scale(t(mat_top)))  # center/scale per gene
# 
# col_fun_top <- circlize::colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
# 
# ## rows = the top N genes (most variable or most linked, depending on your earlier filtering)
# ht_top <- Heatmap(
#     mat_top_z,
#     name = "Z (gene)",
#     col = col_fun_top,
#     cluster_rows = TRUE,
#     cluster_columns = TRUE,
#     show_row_dend = TRUE,
#     show_column_dend = TRUE,
#     column_title = "Cluster annotation",           # X-axis label
#     column_title_side = "bottom",                  # put at bottom if desired
#     row_title = paste0("Top ", N, " genes"),       # Y-axis label
#     show_row_names = TRUE,
#     row_names_gp = grid::gpar(fontsize = 7)
# )
# 
# ht_distributions <- ht_pg %v% ht_gp
# ht_distributions
# 
# # draw in one page
# pdf(file.path(plotDir, paste0("heatmaps_distributions_complexheatmap", f_sufix, ".pdf")),
#               width = 8, height = 10)
# draw(ht_distributions, heatmap_legend_side = "right", merge_legend = TRUE)
# dev.off()
# 
# # top-genes heatmap on its own
# pdf(file.path(plotDir, paste0("heatmap_top_genes_complexheatmap", f_sufix, ".pdf")),
#     width = 8, height = 10)
# draw(ht_top, heatmap_legend_side = "right")
# dev.off()

message("Plots done!!!")


# library("slurmjobs")
# job_single(
#   "06_exploratory_peak_scores",
#   create_shell = TRUE,
#   partition = "katun",
#   memory = "30G",
#   cores = 2,
#   logdir = "logs",
#   command = "Rscript 06_exploratory_peak_scores.R",
#   create_logdir = TRUE
# )

## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
