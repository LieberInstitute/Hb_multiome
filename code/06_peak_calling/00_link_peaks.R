########################################################################
## EDA: Compute LinkPeaks() and filter High-Confident Peaks 
##
## plots and tables to asses links at different thresholds 
##
## Authors. CSC
## Date. March 24, 2025
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
## Note. Seurat objects were created with module load conda_R/4.3.x
########################################################################

library("Seurat")
library("Signac")
library("ggplot2")
## I use BSgenome.Hsapiens.UCSC.hg38 for extracting DNA motifs, k-mers, sequence-based features and compute Tn5 bias correction
library("BSgenome.Hsapiens.UCSC.hg38")  # full reference genome sequence / actual DNA bases (A/T/C/G) for each chromosome
## I use EnsDb.Hsapiens.v86 for extracting gene names, positions, TSSs, chr locations, etc.
library("EnsDb.Hsapiens.v86")           # Gene annotation (GTF-style)
library("tidyverse")
library("dplyr")
library("purrr")
library("here")


## read input arguments
args = commandArgs(trailingOnly = TRUE)
p_met <- args[2]
w_size <- args[4]
# 0: pearson, 1e5
# 1: pearson, 5e4
# 2: spearman, 1e5
# 3: spearman, 5e4

# for testing
# p_met = "spearman"
# w_size = 5e4

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
    f_sufix <- paste0(".", p_met, ".", w_size, ".cells_filtered_3p")
} else {
    message("Input arguments missed")
    stop()
}


# Check/create directories

## clusters renamed for Spatial-Registration on Visium project
inputRDS_Dir <- here(
  "processed-data",
  "05_Clustering_ARCr",
  "17_wnn_clustering_final_ct"
)
inputCVS_Dir <- here(
  "processed-data",
  "05_Clustering_ARCr",
  "02_Hb_celltypes_from_seurat_reanalyze_v3",
  "cvs_files_markers"
)
plotDir <- here(
  "plots",
  "06_peak_calling",
  "00_link_peaks"
)
cvsDir <- here(
    "processed-data",
    "06_peak_calling",
    "00_link_peaks"
)

## Check directories
if (!dir.exists(plotDir)) {
  dir.create(plotDir)
}
if (!dir.exists(cvsDir)) {
    dir.create(cvsDir)
}

## Load Seurat
# Use Seurat with clusters renamed for Spatial-Registration on Visium project
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_HD.rds"
seurat_name <- here(inputRDS_Dir, Seurat_base_name)
SeuratOBJ <- readRDS(here(inputRDS_Dir, Seurat_base_name))
levels(SeuratOBJ)

## filter genes to those expressed in 3% of cells
DefaultAssay(SeuratOBJ) <- "RNA"
rna_counts <- GetAssayData(SeuratOBJ, assay="RNA", slot="data")
length(rownames(rna_counts)) # [1] 36601
keep_genes <- rownames(rna_counts)[Matrix::rowSums(rna_counts > 0) > 0.03 * ncol(rna_counts)]
length(keep_genes) # in count: [1] 14526

## Set atac
DefaultAssay(SeuratOBJ) <- "ATAC"
class(SeuratOBJ[["ATAC"]])
# make a readable base-name for plots
Seurat_base_name <- str_extract(seurat_name, regex("C\\.\\w+"))
Seurat_base_name <- sub("_renamed_visium$", "", Seurat_base_name)
Seurat_base_name
# C.leiden_lsi_r2

##==============================================================================

## Identifies cis-regulatory elements by linking chromatin-accessible peaks to gene expression using correlation (and optionally accounting for covariates).

## preprocessed peaks

## GC content correction
genome <- BSgenome.Hsapiens.UCSC.hg38

SeuratOBJ <- RegionStats(
    object = SeuratOBJ,
    genome = genome,
    assay = "ATAC"
)
SeuratOBJ

## see all chromosomes
table(seqnames(granges(SeuratOBJ)))
## see what are considered standard chromosomes
standardChromosomes(granges(SeuratOBJ))

## Even though this filtering doesn’t change anything in this dataset, it ensures reproducibility
## remove the features that correspond to chromosome scaffolds or other sequences instead of the (22+2) standard chromosomes
peaks.keep <- seqnames(granges(SeuratOBJ)) %in% standardChromosomes(granges(SeuratOBJ))
tryCatch(
    {
        SeuratOBJ <- SeuratOBJ[as.vector(peaks.keep), ]
    }, error = function(e) {
        message(e)
    })


## set a subset of genes to test
#hb_cannonical_genes <- c("GPR151",  "POU4F1", "TAC3")

## find peaks that are correlated with the expression of nearby genes
atac <- LinkPeaks(
    object = SeuratOBJ,
    peak.assay = "ATAC",
    expression.assay = "RNA",
    genes.use = keep_genes,
    method = p_met,
    distance = as.numeric(w_size)             # Only consider peaks within ±100 kb of gene TSS (cis-window)
)


##==============================================================================
## explore and filter strong peak-gene links

## inspect data
head(Links(atac))
# GRanges object with 5 ranges and 5 metadata columns:
#     seqnames              ranges strand |     score        gene
#        <Rle>           <IRanges>  <Rle> | <numeric> <character>
# [1]     chr5 146497025-146516190      * | 0.0515208      GPR151
# [2]     chr5 146516043-146516190      * | 0.0979604      GPR151
# [3]    chr13   78596294-78603560      * | 0.0569620      POU4F1
# [4]    chr13   78597457-78603560      * | 0.0579676      POU4F1
# [5]    chr13   78603450-78603560      * | 0.0594102      POU4F1
# peak    zscore      pvalue
# <character> <numeric>   <numeric>
# [1] chr5-146496517-14649..   2.71746 3.28927e-03
# [2] chr5-146515546-14651..   5.29336 6.00436e-08
# [3] chr13-78595767-78596..   4.35176 6.75248e-06
# [4] chr13-78597000-78597..   4.33686 7.22667e-06
# [5] chr13-78602870-78604..   5.15894 1.24175e-07

link_df <- as.data.frame(Links(atac))
summary(link_df$score)

write.csv(
    link_df,
    file = here(plotDir, paste0("all_peak_gene_links", f_sufix, ".csv")),
    row.names = FALSE
)
## for testing:
link_df <- read.csv(file = here(cvsDir, "all_peak_gene_links.spearman.1e5.csv"))
colnames(link_df)
link_df <- link_df |> 
    mutate(
        gene     = trimws(as.character(gene)),
        seqnames = as.character(seqnames),
        start    = as.numeric(start),
        end      = as.numeric(end)
    )
head(link_df)

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

# should be unique per (gene_name, chr)
stopifnot(!any(duplicated(tss_df[,c("gene_name","seqnames")])))

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
# [1] 5500

## Compute Peak center and distance to TSS
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
    geom_hline(yintercept = 0.3, linetype = "dashed", color = "red") +
    labs(
        x = "Distance from TSS (kb)",
        y = "Correlation Score",
        title = "Peak-Gene Correlation vs. Distance"
    ) + geom_smooth(method = "loess", se = FALSE, color = "darkred") +
    theme_minimal()

ggsave(here(plotDir, 
            paste0("distribution_scores", f_sufix, ".pdf")),
            g1, width = 8, height = 5)

## Highlighting high-confidence links
g2 <- ggplot(link_df2, aes(distance_kb, score, color = score > 0.3)) +
    geom_point(alpha = 0.4) +
    scale_color_manual(
        values = c("FALSE" = "steelblue", "TRUE" = "firebrick"),
        labels = c("FALSE" = "Score ≤ 0.3", "TRUE" = "Score > 0.3"),
        name = "Correlation Threshold"
    ) +
    labs(
        x = "Distance from TSS (kb)",
        y = "Correlation Score",
        title = "Peak-Gene Correlation vs. Distance",
        subtitle = "Red points indicate high-confidence peak-gene links (score > 0.3)"
    ) +
    theme_minimal()

ggsave(here(plotDir, 
            paste0("distribution_scores_high_confidence", f_sufix, ".pdf")),
       g2, width = 8, height = 5)

## filtered peaks
write.csv(link_df,
          file = here(cvsDir, paste0("filtered_peak_gene_links_with_distance", f_sufix, ".csv")),
          row.names = FALSE)

message("TSS Correlation scores completed!")




# library("slurmjobs")
# job_single(
#   "00_link_peaks",
#   create_shell = TRUE,
#   partition = "katun",
#   memory = "30G",
#   cores = 2,
#   logdir = "logs",
#   command = "Rscript -e \"options(width = 120); sessioninfo::session_info()\"",
#   create_logdir = TRUE
# )

## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
