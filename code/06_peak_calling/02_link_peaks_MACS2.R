########################################################################
## EDA: Compute LinkPeaks() and filter High-Confident Peaks 
## INPUT: Peaks generated with CallPeaks() - MACS2
##
## CVS tables with links peaks "global" and "local" with
## - Pearson and Spearman
## - 3 different open-windows sized (check below details) 
##
## Authors. CSC
## Date. August 18, 2025
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
## Note. Seurat objects were created with module load conda_R/4.3.x
########################################################################

library("Seurat")
library("Signac")
## I use BSgenome.Hsapiens.UCSC.hg38 for extracting DNA motifs, k-mers, sequence-based features and compute Tn5 bias correction
library("BSgenome.Hsapiens.UCSC.hg38")  # full reference genome sequence / actual DNA bases (A/T/C/G) for each chromosome
library("tidyverse")
library("tidyr")
library("stringr")
library("here")

#===============================================================================
# resolution_level = "Fine"     # 42 clusters
# resolution_level = "Broad"    # 8 cell-types
# resolution_level = "Mid"      # 18 cell-types
#===============================================================================

## read input arguments
args = commandArgs(trailingOnly = TRUE)
p_met <- args[2]
w_size <- args[4]
# 0: pearson, 1e5
# 1: pearson, 5e4
# 2: spearman, 1e5
# 3: spearman, 5e4

## for testing:
# p_met = "spearman"
# w_size = "5e4"
# resolution_level = "Mid" 

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
    f_sufix <- paste0(".", p_met, ".", w_size, ".cells_filtered_2perc")
} else {
    message("Input arguments missed")
    stop()
}


# Check/create directories

## clusters renamed for Spatial-Registration on Visium project
inputRDS_Dir <- here(
  "processed-data",
  "05_Clustering_ARCr",
  #"17_wnn_clustering_final_ct"
  "22_add_mid_level_clustering" # recent version with final wnn cell-types
)
cvsDir <- here(
    "processed-data",
    "06_peak_calling",
    #"00_link_peaks"
    "01_call_peaks_MACS2"
)

if (!dir.exists(cvsDir)) {
    stop("Peaks file should exist!")
}
   

##==============================================================================
## Set desired meta-data as current level

# Set Seurat identities from a metadata column
set_idents_from_meta <- function(seurat_obj, meta_col, level_order = NULL, na_fill = "Unknown") {
    ## double check level exist on meta-data
    if (!meta_col %in% colnames(seurat_obj@meta.data)) {
        stop("Meta column '", meta_col, "' not found in SeuratOBJ@meta.data")
    }
    # extract target vector
    target_vec <- as.character(seurat_obj[[meta_col]][, 1])
    
    # decide levels and keep appearance order
    if (is.null(level_order)) {
        level_order <- sort(unique(target_vec))
    }
    
    # only update if different from current Idents
    current_idents <- as.character(Idents(seurat_obj))
    if (!identical(current_idents, target_vec)) {
        seurat_obj <- SetIdent(seurat_obj, value = factor(target_vec, levels = level_order))
        message("Idents set from meta column '", meta_col, "'.")
    } else {
        message("Idents already match '", meta_col, "', nothing to do.")
    }
    
    return(seurat_obj)
}


##==============================================================================

## Load Seurat and make verification

# Use Seurat with clusters renamed for Spatial-Registration on Visium project
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_HD.rds"
seurat_name <- here(inputRDS_Dir, Seurat_base_name)

SeuratOBJ <- readRDS(here(inputRDS_Dir, Seurat_base_name))
DefaultAssay(SeuratOBJ) <- "RNA"

# colnames(SeuratOBJ@meta.data)
msg <- switch(resolution_level,
    "Fine" = paste0("Fine-level Cell-Types:\n ", paste(sort(unique(SeuratOBJ$cluster_ann)), collapse = "\n")),
    "Broad" = paste0("Broad-level Cell-Types:\n", paste(sort(unique(SeuratOBJ$merged_cluster)), collapse = "\n")),
    "Mid" = paste0("Mid-level Cell-Types:\n", paste(sort(unique(SeuratOBJ$mid_cluster)), collapse = "\n"))
)

message(msg)
message("Processing ", length(Cells(SeuratOBJ)), " cells")

## set resolution_level
meta_col <- case_when(
    resolution_level=="Fine" ~ "cluster_ann",
    resolution_level=="Broad" ~ "merged_cluster",
    resolution_level=="Mid" ~ "mid_cluster"
)

## set desired idents as current level
SeuratOBJ <- set_idents_from_meta(SeuratOBJ, meta_col = meta_col)
levels(SeuratOBJ)

message("Seurat loaded and ready!")


##==============================================================================
## filter genes to those expressed in 3% of cells

rna_counts <- GetAssayData(SeuratOBJ, assay="RNA", layer="data")
length(rownames(rna_counts)) # [1] 36601
#length(rownames(rna_counts)[Matrix::rowSums(rna_counts > 0)]) # 34738
#length(rownames(rna_counts)[Matrix::rowSums(rna_counts > 0) > 0.02]) # 34738

# filter rna count expressed in at least 2% of the cells
keep_genes <- rownames(rna_counts)[Matrix::rowSums(rna_counts > 0) > 0.02 * ncol(rna_counts)]
length(keep_genes) # in count: [1] 14526

## Set ATAC assay 
DefaultAssay(SeuratOBJ) <- "ATAC"
class(SeuratOBJ[["ATAC"]])

message("Chromatin loaded!")


##==============================================================================

## Identifies cis-regulatory elements by linking chromatin-accessible peaks to gene expression using correlation (and optionally accounting for covariates).

## pre-processed peaks / QC

message("Starting GC content correction ... ")

## GC content correction
genome <- BSgenome.Hsapiens.UCSC.hg38

SeuratOBJ <- RegionStats(
    object = SeuratOBJ,
    genome = genome,
    assay = "ATAC"
)
SeuratOBJ

message("GC content correction done!")

message("Searching LinkPeaks ... ")

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

message("Pre-processing ready ...")


##==============================================================================
## Process "global" Link peak-genes

message("Computing global link-peaks correlations ...")

## find peaks that are correlated with the expression of nearby genes 
atac <- LinkPeaks(
    object = SeuratOBJ,
    peak.assay = "ATAC",
    expression.assay = "RNA",
    genes.use = keep_genes,
    method = p_met,
    distance = as.numeric(w_size)             # Only consider peaks within ±100 kb of gene TSS (cis-window)
)

message("Global link-peaks correlations completed!")

## inspect data
head(Links(atac), n=3)
# GRanges object with 5 ranges and 5 metadata columns:
#     seqnames              ranges strand |     score        gene
#        <Rle>           <IRanges>  <Rle> | <numeric> <character>

link_df <- as.data.frame(Links(atac))
summary(link_df$score)

write.csv(
    link_df,
    file = here(cvsDir, paste0("global_link_peak_genes", f_sufix, ".csv")),
    row.names = FALSE
)

message("Global peaks saved!")



##==============================================================================
## Process "local" Link peak-genes (by cluster)

## find peaks by cluster correlated with the expression of nearby genes 

clusters <- levels(SeuratOBJ)

# make a list to store all subsets
seurat_subsets <- list()

for (clust in clusters) {
    message("Subsetting cluster: ", clust)
    
    seurat_subsets[[clust]] <- subset(
        SeuratOBJ,
        idents = clust
    )
}

# now access each subset by name

for (seurat_cluster in names(seurat_subsets)) {
    # seurat_cluster = "Endo"
    
    message("Processing seurat cluster: ", seurat_cluster)
    
    seurat_subset <- seurat_subsets[[seurat_cluster]]
    seurat_subset
    
    message("Total cells in cluster ", seurat_cluster, ": ", length(Cells(seurat_subset)))
    
    atac <- LinkPeaks(
        object = seurat_subset,
        peak.assay = "ATAC",
        expression.assay = "RNA",
        genes.use = keep_genes,
        method = p_met,
        distance = as.numeric(w_size)             # Only consider peaks within ±100 kb of gene TSS (cis-window)
    )
    
    message("Local link-peaks correlations completed!")
    ## inspect data
    head(Links(atac), n=3)

    ## prepare data to save cvs
    link_df <- as.data.frame(Links(atac))
    summary(link_df$score)
    
    f_name <- paste0(seurat_cluster, "_local_link_peak_genes", f_sufix, ".csv")
    write.csv(
        link_df,
        file = here(cvsDir, f_name),
        row.names = FALSE
    )
    
}


message("Local link-peaks correlations completed!")


message("All done!!!")


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
