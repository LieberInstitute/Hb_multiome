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

## input directories

#inputRDS_Dir <- here("processed-data", "05_Clustering_ARCr", "01_clustering_std_method")
procData_Dir <- here("processed-data", "05_Clustering_ARCr")

inputSeurat_RDS <- here(procData_Dir, "05_rename_idents", 
                        "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium.rds")

input_ct_summary_CSV <- here(procData_Dir, "05_rename_idents", 
                                    "WNN_full_annotation_meta_data.csv")

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

summary_ct_df <- read.csv(input_ct_summary_CSV)

## verification
colnames(summary_ct_df)
summary_ct_df |>
    select(cluster, ct_frequency_in2ct, ct_CRegistration_fine) |>
    head()

# now remove from full_annotation_df
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
## summarizes new meta-data for mean-ratio: mean_ratio_detail	mean_ratio_detail	mean_ratio_range

message("=========== Summarizing Mean-Ratio dataset")

## takes top 10 MeanRatio genes per cluster and collapse `cellType.2nd + MeanRatio` per cluster
collapsed_df <- marker_ranks_all_top10 |>
    group_by(cellType.target) |>
    summarise(
        comparisons = str_c(cellType.2nd, collapse = ", "),
            #paste0(cellType.2nd, " (", round(MeanRatio, 2), ")"),
            #collapse = ", "
        #),
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
            summary_str
        })
    )
collapsed_summary_mean_ratio$mean_ratio_detail
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
head(collapsed_summary_mean_ratio$mean_ratio_repeated_ct)
## clean mean_ratio_repeated_ct output, remove x4 characters
collapsed_summary_mean_ratio <- collapsed_summary_mean_ratio |>
    mutate(
        #mean_ratio_repeated_ct = str_replace(mean_ratio_repeated_ct, "^C\\.\\d+\\.", "") |> str_remove("\\s.*")
        mean_ratio_repeated_ct = str_remove(mean_ratio_repeated_ct, "\\s.*")
    )
head(collapsed_summary_mean_ratio$mean_ratio_repeated_ct)


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
            vals <- as.numeric(str_split(x, ",\\s*")[[1]]) # splits your string into a vector by commas, removing spaces
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

message("Itegrating new meta-data into previous WNN_full_annotation_meta_data.csv") 

collapsed_summary_range <- mean_ratio_new_columns |>
    mutate(cluster = str_sub(cellType.target, 1, 4))
head(collapsed_summary_range)

## df with previous ct data
head(summary_ct_df[1:4])

## join both df
if (!identical(collapsed_summary_range$cluster, summary_ct_df$cluster)) { stop("Not equal cluster ID information on dataframes!") }

colnames(collapsed_summary_range)
colnames(summary_ct_df)

WNN_full_annotation_df <- inner_join(
    summary_ct_df, 
    collapsed_summary_range, 
    by = "cluster"
    )
colnames(WNN_full_annotation_df)
head(WNN_full_annotation_df)


WNN_full_annotation_df <- 
    WNN_full_annotation_df |>
    select(
        cluster,
        ct_frequency_in2ct, 
        ct_CRegistration_fine,
        mean_ratio_repeated_ct,
        mean_ratio_range,
        ct_MeanRatio_support,
        number_cells, 
        cluster_percentage,
        all_cell_types_by_frequency,
        mean_ratio_top2, 
        mean_ratio_detail)

head(WNN_full_annotation_df[1:7])

## save summary WNN cluster annotations
f_name <- here(outputCSV_Dir, "WNN_full_annotation_meta_data.csv")
write.csv(WNN_full_annotation_df, f_name, row.names = FALSE)



## ==============================================================================

message("Adding new meta-data for ambiguous clusters") 

# ## Create column with ambiguous clusters
# lookup_ct <- tribble(
#     ~cluster, ~ct_ambiguous,
#     "C.03", "Endo vs Excit.Thal",
#     "C.04", "LHb4",
#     "C.06", "LHb4",
#     "C.09", "LHb4",
#     "C.13", "LHb4",
#     "C.16", "MHb vs LHb",
#     "C.24", "MHb vs LHb",
#     "C.30", "MHb vs LHb",
#     "C.31", "LHb4",
#     "C.34", "Oligo",
#     "C.35", "Excit.Thal"
# )

# ## Join to main table 
# WNN_full_annotation_with_ambiguous_df <- WNN_full_annotation_df |>
#     left_join(lookup_ct, by = "cluster") |>
#     mutate(
#         ct_ambiguous = coalesce(ct_ambiguous, "")
#     )
# head(WNN_full_annotation_with_ambiguous_df)

colnames(WNN_full_annotation_df)

update_ct_info <- tribble(
    ~cluster, ~ct_ambiguous,        ~ct_final,     ~description_support,
    "C.03",   "Endo vs Excit.T",    "*Excit.Thal",  "SReg support for Thal + mean-ratio support too + HClust support too",
    "C.04",   "LHb4",               "*ExcitT.LHb.4","SReg support for LHb + similar mean-ratio for Thal and LHb + HClust support for both too",
    "C.06",   "LHb4",               "*ExcitT.LHb.4","mean-ratio support for Thal and LHb + HClust support for both too",
    "C.09",   "LHb4",               "*ExcitT.LHb.4","SReg support for LHb + similar mean-ratio for Thal and LHb + HClust support for both too",
    "C.13",   "LHb4",               "*Thal",        "SReg support for LHb + strong mean-ratio support too + HClust support too",
    "C.16",   "MHb vs LHb",         "*MHb.1.2",     "SReg support for LHb + mean-ratio support for MHb and LHb + HClust support for MHb",
    "C.24",   "MHb vs LHb",         "*LHb.1.3.4",   "strong mean-ratio support for LHb + HClust support too",
    "C.30",   "MHb vs LHb",         "*MHb.LHb",     "similar mean-ratio for MHb and LHb + HClust support for LHb",
    "C.31",   "LHb4",               "*ExcitT.LHb.4","strong mean-ratio support for Thal + HClust support for ExcitT.LHb.4",
    "C.34",   "Oligo",              "*MHb.1.2",     "mean-ratio support for MHb and LHb + HClust support too",
    "C.35",   "Excit.Thal",         "*Excit.Thal",  "SR support for Thal + HClust support too"
)
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

# ## Add additional columns with ct_final and description_support
# WNN_full_annotation_with_ambiguous_df <- 
#     WNN_full_annotation_with_ambiguous_df |>
#     mutate(HClust_support = "",
#            ct_final = "",
#            description_support=""
#            )
# WNN_full_annotation_with_ambiguous_df$HClust_support

## Update HClust_support from HClustering results, only branches with two leaves from the same ct
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

# # `ct_ambiguous` is assumed to be a vector of cluster names, e.g.:
# # ct_ambiguous <- c("C.02", "C.22", "C.07", "C.36", "C.11", "C.14", ...)
# WNN_full_annotation_with_ambiguous_df <- WNN_full_annotation_with_ambiguous_df %>%
#     mutate(
#         ct_final = if_else(
#             ct_ambiguous != "", 
#             "",                  # leave blank if ct_ambiguous has data
#             ct_MeanRatio_support # otherwise copy from support
#         )
#     )

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
        ct_CRegistration_fine,
        ct_ambiguous,
        mean_ratio_repeated_ct,
        mean_ratio_range,
        #ct_MeanRatio_support,
        HClust_support,
        ct_final,
        description_support,
        number_cells, 
        cluster_percentage,
        mean_ratio_top2, 
        mean_ratio_detail,
        all_cell_types_by_frequency)
head(WNN_full_annotation_df)

## save summary WNN cluster annotations
f_name <- here(outputCSV_Dir, "WNN_full_annotation_meta_data.csv")
write.csv(WNN_full_annotation_df, f_name, row.names = FALSE)



library("slurmjobs")
job_single("16_rename_idents_post_mean_ratio_HClust", cores = 2, memory = 80, partition = "katun", create_shell = TRUE)


library("sessioninfo")
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()


