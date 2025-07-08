s########################################################################
## Refine Annotation on WNN Clusters post Mean-Ratio Scores and HClustering Patterns
## INPUT:
##      (1) CSV with mean-ratio scores
##      (2) CSV with HCLust scores / Or likely manual curation
## OUPUT:
##      (1) Seurat with updated idents: clusters_ann (meta-data)
##      (2) Likely some UMAPS
## Authors. CSC 
## Date. Jul 30th, 2025
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
########################################################################

library("Seurat")
library("purrr")
library("tibble")
library("ggplot2")
library("dplyr")
library("tidyr")
library("stringr")
library("here")

## Directories

procData_Dir <- here("processed-data", "05_Clustering_ARCr")

# Seurat with Ident names 
inputSeurat_RDS <- here(procData_Dir, "05_rename_idents", 
                        "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium.rds")
# Meta-data summary 
input_ct_summary_CSV <- here(procData_Dir, "05_rename_idents", 
                                    "WNN_full_annotation_meta_data.csv")
# Mean-ratio data
inputMeanRatio_RDS <- here(procData_Dir, "13_wnn_geneExp_plt_mean_ratio_annotated", 
                            "marker_ranks_6_8_10_12_29k.RData") # marker_ranks_6k_12k_allk.RData

outputRDS_Dir <- here(procData_Dir, "16_rename_idents_post_mean_ratio_HClust")
outputCSV_Dir <- here(procData_Dir, "16_rename_idents_post_mean_ratio_HClust")
plotDir <- here("plots", "05_Clustering_ARCr", "16_rename_idents_post_mean_ratio_HClust")

## Check directories
if (!dir.exists(outputRDS_Dir)) {dir.create(outputRDS_Dir)}
if (!dir.exists(plotDir)) {dir.create(plotDir)}


## ==============================================================================
## load summary with cell_types and percentages by clusters 

message("=========== Loading Full Summary:\n ", basename(input_ct_summary_CSV))

summary_ct_df <- read.csv(input_ct_summary_CSV)

## verification
colnames(summary_ct_df)
summary_ct_df |>
    select(cluster, ct_frequency_in2ct, ct_CRegistration_fine) |>
    head()

## remove clusters<10 cells from full_annotation_df
summary_ct_df <- summary_ct_df[!summary_ct_df$number_cells<10, ]
summary_ct_df$number_cells
if (anyDuplicated(summary_ct_df$cluster)) { stop("There are duplicated clusters") }


## ==============================================================================
## load mean-ratio dataset 29k

message("=========== Loading Mean-Ratio dataset:\n ", basename(inputMeanRatio_RDS))

meanRatio_sets <- load(inputMeanRatio_RDS)
ls()
unique(marker_ranks_all$cellType.target)

## filter top 10 mean-ratio genes by cluster
marker_ranks_all_top10 <- marker_ranks_all |>
    group_by(cellType.target) |>
    slice_max(order_by = MeanRatio, n = 10, with_ties = FALSE) |>
    ungroup()

unique(marker_ranks_all_top10$cellType.target)
head(marker_ranks_all_top10)
nrow(marker_ranks_all_top10)



## ==============================================================================
## summarizes new meta-data for mean-ratio: mean_ratio_detail	mean_ratio_top2	    mean_ratio_repeated_ct
##                                          mean_ratio_range     mean_ratio_range

message("=========== Summarizing Mean-Ratio dataset")

## takes top 10 MeanRatio genes per cluster and collapse `cellType.2nd + MeanRatio` per cluster
collapsed_df <- marker_ranks_all_top10 |>
    group_by(cellType.target) |>
    summarise(
        comparisons = str_c(cellType.2nd, collapse = ", "),
        .groups = "drop"
    )
head(collapsed_df)
nrow(collapsed_df)

## count repeated clusters based on mean-ratio 'comparisons'
collapsed_summary_mean_ratio <- collapsed_df |>
    mutate(
        mean_ratio_detail = map_chr(comparisons, function(comp_str) {
            clusters <- str_split(comp_str, ",\\s*")[[1]] # # split by comma, trim spaces
            cluster_counts <- table(clusters) # count repeated clusters
            # format like "C.38.Inhib.Thal x2"
            summary_str <- paste0(names(cluster_counts), " x", cluster_counts, collapse = ", ")
            # insert \n every ~70 chars
            str_wrap(summary_str, width = 70)
        })
    )
#nchar(collapsed_summary_mean_ratio$mean_ratio_detail[7])
colnames(collapsed_summary_mean_ratio)
collapsed_summary_mean_ratio$comparisons <- NULL
head(collapsed_summary_mean_ratio)

## extract the top 2 most repeated patterns from each mean_ratio_detail
collapsed_summary_mean_ratio <- collapsed_summary_mean_ratio |>
    mutate(
        mean_ratio_top2 = map_chr(mean_ratio_detail, function(summary_str) {
            items <- str_split(summary_str, ",\\s*")[[1]] # # split by comma
            # extract the counts (xN) as numbers
            counts <- as.numeric(str_extract(items, "(?<= x)\\d+"))
            # order descending by counts
            top_items <- items[order(-counts)]
            # take top 2 and collapse back
            paste(top_items[1:min(2, length(top_items))], collapse = ", ")
        })
    )
head(collapsed_summary_mean_ratio)

## extract the top 1 most repeated patterns from each mean_ratio_detail
collapsed_summary_mean_ratio <- collapsed_summary_mean_ratio |>
    mutate(
        mean_ratio_repeated_ct = map_chr(mean_ratio_detail, function(summary_str) {
            items <- str_split(summary_str, ",\\s*")[[1]] # # split by comma
            # extract the counts (xN) as numbers
            counts <- as.numeric(str_extract(items, "(?<= x)\\d+"))
            # order descending by counts
            top_items <- items[order(-counts)]
            # take top 1
            top_items[1]
        })
    )
head(collapsed_summary_mean_ratio)
# ## clean mean_ratio_repeated_ct output, remove x4 characters
# collapsed_summary_mean_ratio <- collapsed_summary_mean_ratio |>
#     mutate(
#         #mean_ratio_repeated_ct = str_replace(mean_ratio_repeated_ct, "^C\\.\\d+\\.", "") |> str_remove("\\s.*")
#         mean_ratio_repeated_ct = str_remove(mean_ratio_repeated_ct, "\\s.*")
#     )
# head(collapsed_summary_mean_ratio$mean_ratio_repeated_ct)


message("=========== Summarizing Mean-Ratio ranges")

## takes top 10 MeanRatio genes per cluster and extract Range of MeanRatio
collapsed_df <- marker_ranks_all_top10 |>
    group_by(cellType.target) |>
    summarise(
        mean_ratio_ranges = str_c(
            round(MeanRatio, 2),
            collapse = ", "
        ),
        .groups = "drop"
    )
head(collapsed_df)

## Add mean_ratio_range per cluster
collapsed_summary_range <- collapsed_df |>     
    mutate(
        mean_ratio_range = map_chr(mean_ratio_ranges, function(x) {
            vals <- as.numeric(str_split(x, ",\\s*")[[1]]) # splits string into a vector by commas, removing spaces
            rng <- range(vals, na.rm = TRUE) 
            paste0("[", round(rng[2], 2), " - ", round(rng[1], 2), "]")
        })
    )
collapsed_summary_range$mean_ratio_ranges <- NULL
head(collapsed_summary_range)


message("=========== Joining Mean-Ratio summary and range summary")

mean_ratio_new_columns <- inner_join(collapsed_summary_mean_ratio, collapsed_summary_range, by = "cellType.target")
head(mean_ratio_new_columns)


## ==============================================================================

message("Integrating new meta-data into previous WNN_full_annotation_meta_data.csv") 

collapsed_summary_range <- mean_ratio_new_columns |>
    mutate(cluster = str_sub(cellType.target, 1, 4))
head(collapsed_summary_range)

## df with previous ct data
head(summary_ct_df[1:4])

## join both df
if (!identical(collapsed_summary_range$cluster, summary_ct_df$cluster)) { stop("Not equal cluster ID information on dataframes!") }

## check column names and rename if necessary
colnames(collapsed_summary_range)
colnames(summary_ct_df)
summary_ct_df <- summary_ct_df |> 
    rename(frequency_repeated_ct = all_cell_types_by_frequency)

WNN_full_annotation_df <- inner_join(
    summary_ct_df, 
    collapsed_summary_range, 
    by = "cluster"
    )
colnames(WNN_full_annotation_df)
head(WNN_full_annotation_df)


## ==============================================================================

message("Adding new meta-data for ambiguous clusters") 

## Create column with ambiguous clusters
colnames(WNN_full_annotation_df)

update_ct_info <- tribble(
    ~cluster, ~ct_ambiguous,        ~ct_final,     ~description_support,
    "C.03",   "Endo vs Excit.T",    "*Excit.Thal",  "SReg for Thal + mean-ratio for Thal + \nHClust in Thal clade", 
    "C.04",   "LHb4",               "*LHb.4",       "SReg for LHb + mean-ratio for Thal and LHb + \nVisiumHD for LHb + HClust in LHb clade",
    "C.06",   "LHb4",               "*LHb.4",       "mean-ratio for Thal and LHb + VisiumHD for LHb",
    "C.09",   "LHb4",               "*LHb.4",       "SReg for LHb + mean-ratio for Thal and LHb + \nVisiumHD for LHb + HClust in LHb clade",
    "C.13",   "LHb4",               "*LHb.4",       "SReg for LHb + mean-ratio for Thal + \nVisiumHD for LHb + HClust in LHb clade",
    "C.16",   "MHb vs LHb",         "*MHb.1.2",     "SReg for LHb + mean-ratio for MHb and LHb + \nVisiumHD for Hb + HClust support for MHb",
    "C.24",   "MHb vs LHb",         "*LHb.1.3.4",   "mean-ratio for LHb + HClust Hb clade",
    "C.30",   "MHb vs LHb",         "*MHb.LHb",     "mean-ratio for MHb and LHb + VisiumHD for Hb (low) + \nHClust in Hb clade",
    "C.31",   "LHb4",               "*LHb.4",       "mean-ratio for Thal and LHb + VisiumHD for Hb (low) + \nHClust in Hb clade",
    "C.34",   "Oligo",              "*MHb.1.2",     "mean-ratio for MHb + VisiumHD for Hb + \nHClust in MHb clade",
    "C.35",   "Excit.Thal",         "*Excit.Thal",  "SReg for Thal + HClust in Thal clade"
)

## ??
# C.25 still Thal ? -- In LHb clade
# C.31 still Thal ? -- In LHb clade
# C.24 

## Join to main table
WNN_full_annotation_with_ambiguous_df <- WNN_full_annotation_df |>
    left_join(update_ct_info, by = "cluster") |>
    mutate(
        ct_ambiguous = coalesce(ct_ambiguous, ""),
        ct_final = coalesce(ct_final, ""),
        description_support = coalesce(description_support, ""),
        HClust_support = ""
    )
colnames(WNN_full_annotation_with_ambiguous_df)
head(WNN_full_annotation_with_ambiguous_df)

## Update HClust_support column from HClustering results, only branches with two leaves from the same ct
cluster_hclust_map <- tribble(
    ~cluster, ~HClust_support,
    "C.02", "2-22",
    "C.22", "2-22",
    "C.07", "7-36",
    "C.36", "7-36",
    "C.11", "11-14",
    "C.14", "11-14",
    "C.13", "13-39",
    "C.39", "13-39",
    "C.28", "28-38",
    "C.01", "01-19",
    "C.19", "01-19",
    "C.12", "12-32",
    "C.32", "12-32",
    "C.03", "03-37",
    "C.37", "03-37",
    "C.15", "15-35",
    "C.35", "15-35"
)
cluster_hclust_map

## Joint to main table
WNN_full_annotation_with_ambiguous_df <- WNN_full_annotation_with_ambiguous_df |>
    select(-HClust_support) |>
    left_join(cluster_hclust_map, by = "cluster") |>
    mutate(
        HClust_support = coalesce(HClust_support, "")
    )
WNN_full_annotation_with_ambiguous_df$HClust_support


## copy ct_MeanRatio_support on ct_final only if ct_final empty or blank
WNN_full_annotation_with_ambiguous_df$ct_final
WNN_full_annotation_with_ambiguous_df <- WNN_full_annotation_with_ambiguous_df |>
    mutate(
        ct_final = if_else(
            ct_final == "", 
            ct_MeanRatio_support,
            ct_final
        )
    )
WNN_full_annotation_with_ambiguous_df$ct_final

## rearrange columns again
WNN_full_annotation_df <- 
    WNN_full_annotation_with_ambiguous_df |>
    select(
        cluster,
        ct_frequency_in2ct, 
        frequency_repeated_ct,
        ct_CRegistration_fine,
        ct_MeanRatio_support,
        mean_ratio_repeated_ct,
        mean_ratio_top2,
        mean_ratio_detail,
        mean_ratio_range,
        HClust_support,
        ct_final,
        description_support,
        number_cells, 
        cluster_percentage
        )
head(WNN_full_annotation_df)

## save summary WNN cluster annotations
f_name <- here(outputCSV_Dir, "WNN_full_annotation_meta_data.csv")
write.csv(WNN_full_annotation_df, f_name, row.names = FALSE)

# library("slurmjobs")
# job_single("16_rename_idents_post_mean_ratio_HClust", cores = 2, memory = 80, partition = "katun", create_shell = TRUE)


library("sessioninfo")
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()


