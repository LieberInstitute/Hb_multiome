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
        ct_post_mean_ratio = map_chr(mean_ratio_detail, function(summary_str) {
            items <- str_split(summary_str, ",\\s*")[[1]] # # split by comma
            # extract the counts (xN) as numbers
            counts <- as.numeric(str_extract(items, "(?<= x)\\d+"))
            # order descending by counts
            top_items <- items[order(-counts)]
            # take top 1
            top_items[1]
        })
    )
head(collapsed_summary_mean_ratio$ct_post_mean_ratio)
## clean ct_post_mean_ratio output, remove x4 characters
collapsed_summary_mean_ratio <- collapsed_summary_mean_ratio |>
    mutate(
        #ct_post_mean_ratio = str_replace(ct_post_mean_ratio, "^C\\.\\d+\\.", "") |> str_remove("\\s.*")
        ct_post_mean_ratio = str_remove(ct_post_mean_ratio, "\\s.*")
    )
head(collapsed_summary_mean_ratio$ct_post_mean_ratio)


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

WNN_full_annotation_df <- inner_join(summary_ct_df, collapsed_summary_range, by = "cluster")
head(WNN_full_annotation_df)

colnames(WNN_full_annotation_df)

WNN_full_annotation_df <- 
    WNN_full_annotation_df |>
    select(cluster,
           ct_frequency_in2ct, 
           ct_CRegistration_fine,
           ct_post_mean_ratio,
           mean_ratio_range,
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

# ## Add extra columns o use later 
# WNN_full_annotation_with_ambiguous_df <- WNN_full_annotation_df |>
#     mutate(
#         ct_ambiguous = "",
#         ct_post_HCLust_manual = "",
#         ct_final = "",
#         annotation_reason = ""
#     )
# 
# colnames(WNN_full_annotation_with_ambiguous_df)

## Add ambiguous clusters
lookup_ct <- tibble(
    cluster = c("C.03", "C.04", "C.06", "C.09", "C.13", 
                "C.16", "C.24", "C.30", "C.31", "C.34", "C.35"),
    ct_ambiguous = c("Endo vs Excit.Thal", "LHb4", "LHb4", "LHb4", "LHb4", 
                     "MHb vs LHb", "MHb vs LHb", "MHb vs LHb", "LHb4", "Oligo", "Excit.Thal")
)
lookup_ct

WNN_full_annotation_with_ambiguous_df <- WNN_full_annotation_df |>
    left_join(lookup_ct, by = "cluster")
head(WNN_full_annotation_with_ambiguous_df)

# WNN_full_annotation_with_ambiguous_df <- WNN_full_annotation_with_ambiguous_df |>
#     mutate(ct_ambiguous = replace_na(ct_ambiguous, "Solved.on.review"))
# head(WNN_full_annotation_with_ambiguous_df)

## rearrange again
WNN_full_annotation_with_ambiguous_df <- 
    WNN_full_annotation_with_ambiguous_df |>
    select(cluster,
           ct_frequency_in2ct, 
           ct_CRegistration_fine,
           ct_ambiguous,
           mean_ratio_top2,
           ct_post_mean_ratio,
           mean_ratio_range,
           number_cells, 
           cluster_percentage,
           all_cell_types_by_frequency,
           mean_ratio_detail)
head(WNN_full_annotation_with_ambiguous_df)

## save summary WNN cluster annotations
f_name <- here(outputCSV_Dir, "WNN_full_annotation_meta_data.csv")
write.csv(WNN_full_annotation_df, f_name, row.names = FALSE)



## ==============================================================================
## load Seurat with WNN clusters
# 
# message("Loading Seurat:\n ", basename(inputSeurat_RDS))
# 
# SeuratOBJ <- readRDS(inputSeurat_RDS)
# SeuratOBJ
# 
# total_cells <- length(Cells(x = SeuratOBJ))
# message("Total cells: ", total_cells)
# message("Processing ", nrow(unique(SeuratOBJ[["seurat_clusters"]])), " clusters")


# # ==============================================================================
# ## Curated Manual Annotation based on: (1) gene marker match frequency (DD and LB) and (2) Clustering-Registration `Broad` and `Fine`
# # Details on: https://github.com/LieberInstitute/Hb_multiome/tree/8b614669db8f7b833e5a793fa01326f839c57cd8/data 
# cell_types_curated <- tribble(
#     ~cluster, ~ct_frequency_in2ct, ~ct_CRegistration_broad, ~ct_CRegistration_fine, ~ct_MeanRatio_support,
#     "C.01", "undetermined",        "Inhib.Thal",          "Inhib.Thal",             "Inhib.Thal",
#     "C.02", "DD_Oligo",            "Oligo",                "Oligo",                  "Oligo",
#     "C.03", "undetermined",        "Excit.Thal",           "Excit.Thal",             "Excit.Thal", # checking
#     "C.04", "undetermined",        "LHb",                  "LHb.4",                  "Excit.Thal", # checking
#     "C.05", "DD_LHb",              "LHb",                  "LHb.2.7",                "LHb.2.7",
#     "C.06", "DD_Excit.Thal",       "undetermined",         "undetermined",           "ExcitT.LHb.4",  # checking
#     "C.07", "DD_MHb",              "MHb",                  "MHb.2",                  "MHb.2",
#     "C.08", "undetermined",        "LHb",                  "LHb.4",                  "LHb.4",
#     "C.09", "undetermined",        "LHb",                  "LHb.4",                  "ExcitT.LHb.4",  # checking ?
#     "C.10", "DD_MHb",              "MHb",                  "MHb.1",                  "MHb.1",
#     "C.11", "DD_MHb",              "MHb",                  "MHb.1.2",                "MHb.1.2",
#     "C.12", "undetermined",        "Excit.Thal",           "Excit.Thal",             "Excit.Thal",
#     "C.13", "no-match",            "LHb",                  "LHb.4",                  "Thal",  # checking - HD
#     "C.14", "DD_MHb",              "MHb",                  "MHb.1",                  "MHb.1",
#     "C.15", "DD_Excit.Thal",       "Excit.Thal",           "Excit.Thal",             "Excit.Thal",
#     "C.16", "DD_MHb",              "LHb",                  "LHb.6",                  "MHb.1.2",  # checking
#     "C.17", "DD_Excit.Thal",       "Excit.Thal",           "Excit.Thal",             "Excit.Thal",
#     "C.18", "DD_LHb",              "LHb",                  "LHb.1.3.4",              "LHb.1.3.4",
#     "C.19", "DD_Inhib.Thal",       "Inhib.Thal",           "Inhib.Thal",             "Inhib.Thal",
#     "C.20", "DD_Astrocyte",        "Astrocyte",            "Astrocyte",              "Astrocyte",
#     "C.21", "DD_Astrocyte",        "Astrocyte",            "Astrocyte",              "Astrocyte",
#     "C.22", "undetermined",        "Oligo",                "Oligo",                  "Oligo",
#     "C.23", "DD_LHb",              "LHb",                  "LHb.1",                  "LHb.1",
#     "C.24", "DD_LHb",              "undetermined",         "undetermined",           "LHb", # checking
#     "C.25", "undetermined",        "Excit.Thal",           "Excit.Thal",             "Excit.Thal",
#     "C.26", "DD_OPC",              "OPC",                  "OPC",                    "OPC",
#     "C.27", "DD_Microglia",        "Microglia",            "Microglia",              "Microglia",
#     "C.28", "DD_Inhib.Thal",       "Inhib.Thal",           "Inhib.Thal",             "Inhib.Thal",
#     "C.29", "DD_Endo",             "Endo",                 "Endo",                   "Endo",
#     "C.30", "DD_LHb",              "undetermined",         "undetermined",           "MHb.LHb", # checking HD
#     "C.31", "DD_Excit.Thal",       "undetermined",         "undetermined",           "ExcitT.LHb.4", # checking HD
#     "C.32", "undetermined",        "Excit.Thal",           "Excit.Thal",             "Excit.Thal",
#     "C.33", "DD_LHb",              "LHb",                  "LHb.1.3",                "LHb.1.3",
#     "C.34", "DD_Oligo",            "undetermined",         "undetermined",           "Oligo", # checking
#     "C.35", "undetermined",        " Excit.Thal",          " Excit.Thal",            "Excit.Thal",
#     "C.36", "DD_MHb",              "MHb",                  "MHb.3",                  "MHb.3",
#     "C.37", "undetermined",        "Thal",                 "Thal",                   "Thal",
#     "C.38", "DD_Inhib.Thal",       "Inhib.Thal",           "Inhib.Thal",             "Inhib.Thal",
#     "C.39", "DD_Inhib.Thal",       "Inhib.Thal",           "Inhib.Thal",             "Inhib.Thal",
#     "C.40", "DD_LHb",              "LHb",                  "LHb.4",                  "LHb.4",
#     "C.41", "DD_Microglia",        "Microglia",            "Microglia",              "Microglia",
#     "C.42", "no-match",            "no-match",             "no-match",               "no-match"
# )
# 
# # double check spaces 
# cell_types_curated <- cell_types_curated |>
#     dplyr::mutate(across(everything(), ~ trimws(.)))
# 
# if (anyDuplicated(cell_types_curated$cluster)) { stop("There are duplicated clusters") }
# 
# message("Cell typers manually curated:")
# cell_types_curated <- as.data.frame(cell_types_curated)
# ## ct_MeanRatio_support was used to give support to cell_type_final, but this could be adjusted as the results support
# ## So, I duplicated the column to keept code consistent
# cell_types_curated$cell_type_final <- cell_types_curated$ct_MeanRatio_support
# cell_types_curated
# 
# ## inner join pseudo annotation (EDA) with curated annotation for supplementary material
# full_annotation_df <- left_join(cell_types_curated, summary_ct_df, by = "cluster")
# colnames(full_annotation_df)
# 
# # filter and sort columns
# full_annotation_df <- full_annotation_df[, c("cluster", "ct_frequency_in2ct", "ct_CRegistration_fine", "ct_MeanRatio_support", "cell_type_final", "number_cells", "cluster_percentage", "all_cell_types_by_frequency")]
# head(full_annotation_df)
# # cluster ct_frequency_in2ct ct_CRegistration_broad ct_CRegistration_fine
# # 1    C.01       undetermined             Inhib.Thal            Inhib.Thal
# # 2    C.02           DD_Oligo                  Oligo                 Oligo
# # 3    C.03       undetermined             Excit.Thal            Excit.Thal
# # 4    C.04       undetermined                    LHb                 LHb.4
# # 5    C.05             DD_LHb                    LHb               LHb.2.7
# # 6    C.06      DD_Excit.Thal           undetermined          undetermined
# # number_cells cluster_percentage
# # 1         3906             7.012%
# # 2         3371             6.052%
# # 3         3271             5.872%
# # 4         2911             5.226%
# # 5         2771             4.975%
# # 6         2667             4.788%
# # all_cell_types_by_frequency
# # 1 DD_Inhib.Thal (10), DD_LHb (1), LB_excitatory_neuron (1), LB_inhibitory_neuron (2), LB_Thalamus/MDm (1)
# # 2                                                     DD_Oligo (5), LB_neuron (2), LB_oligodendrocyte (1)
# # 3                                                        DD_Astrocyte (1), DD_Endo (4), DD_Excit.Thal (1)
# # 4                DD_Excit.Thal (6), DD_Inhib.Thal (1), DD_LHb (2), DD_OPC (1), LB_LHB neuron specific (1)
# # 5                   DD_Excit.Thal (1), DD_LHb (10), LB_Hb neuron specific (2), LB_LHB neuron specific (1)
# # 6                                                                           DD_Excit.Thal (6), DD_OPC (1)
# 
# ## We do not need this chunck any more, as we already have a polished annotation 
# ## Keep ct_frequency_in2ct only when ct_CRegistration_fine is "undetermined", otherwise use ct_CRegistration_fine
# # full_annotation_df <- full_annotation_df |>
# #     mutate(
# #         cell_type_final = if_else(
# #             ct_CRegistration_fine == "undetermined",
# #             ct_frequency_in2ct,
# #             ct_CRegistration_fine
# #         )
# #     )
# 
# # filter and sort columns
# full_annotation_df <- full_annotation_df[, c("cluster", "ct_frequency_in2ct", "ct_CRegistration_fine", "ct_MeanRatio_support", "cell_type_final", 
#                                              "number_cells", "cluster_percentage", "all_cell_types_by_frequency")]
# message("Updating full annotation summary ...")
# head(full_annotation_df)
# sort(unique(full_annotation_df$cell_type_final))
# 
# ## save summary WNN cluster annotations
# f_name <- here(outputCSV_Dir, "WNN_full_annotation_meta_data.csv")
# write.csv(full_annotation_df, f_name, row.names = FALSE)
# 
# ## create symlink path to data/ -- where I am actually storing summaries for speed searching
# dir_target <- here("data")
# 
# ## full symlink path
# symlink_path <- file.path(dir_target, basename(f_name))
# 
# # Create the symbolic link (if it doesn't already exist)
# if (!file.exists(symlink_path)) {
#     success <- tryCatch({
#         file.symlink(from = f_name, to = symlink_path)
#     }, warning = function(w) {
#         message("Warning: ", conditionMessage(w))
#         FALSE
#     }, error = function(e) {
#         message("Error: ", conditionMessage(e))
#         FALSE
#     })
#     
#     if (isTRUE(success)) {
#         message("Symlink created: ", symlink_path)
#     } else {
#         message("Failed to create symlink.")
#     }
# } else {
#     message("Symlink already exists: ", symlink_path)
# }
# 
# # ==============================================================================
# 
# ## Keep clusters with at least 10 cells
# 
# cluster_counts <- table(Idents(SeuratOBJ))
# valid_clusters <- names(cluster_counts[cluster_counts >= 10])
# SeuratOBJ <- subset(SeuratOBJ, idents = valid_clusters)
# # confirm clusters removed
# removed_clusters <- names(cluster_counts[cluster_counts < 10])
# message("Removed ID clusters with less than 10 cells: ", removed_clusters, "\n")
# levels(SeuratOBJ)
# 
# # now remove from full_annotation_df
# full_annotation_df <- full_annotation_df[!full_annotation_df$number_cells<10, ]
# head(full_annotation_df)
# sort(unique(full_annotation_df$cell_type_final))
# 
# ## extract original idents (clusters) and prepare the new ident names to rename Seurat clusters
# 
# # sanity check if any Seurat clusters was excluded
# Seurat_clusterIDS <- as.integer(levels(SeuratOBJ$seurat_clusters))
# Seurat_clusterIDS <- paste0("C.", sprintf('%02d',Seurat_clusterIDS))
# if (!identical(Seurat_clusterIDS, full_annotation_df$cluster)) {
#     stop("Cluster IDs in Seurat object do not match those in the annotation data frame.")
# }
# 
# # Extract numeric part from "C.XX" and convert to character
# full_annotation_df$cluster_id <- as.character(as.numeric(sub("C\\.", "", full_annotation_df$cluster)))
# # Use polished cell-types after inspect clustering-registration at fine res, plus mean-ratio scores
# new_names <- setNames(paste0(full_annotation_df$cluster, ".", full_annotation_df$cell_type_final), full_annotation_df$cluster_id)
# new_names
# 
# ## rename idents 
# SeuratOBJ <- RenameIdents(object = SeuratOBJ, new_names)
# head(levels(SeuratOBJ))
# head(Idents(SeuratOBJ))
# 
# ## Set new factor levels for identities to first plot Hb clusters
# all_clusters <- levels(SeuratOBJ)
# hb_clusters <- grep("MHb|LHb", all_clusters, value = TRUE)
# hb_clusters
# # [1] "C.04.LHb.4"     "C.05.LHb.2.7"   "C.07.MHb.2"     "C.08.LHb.4"    
# # [5] "C.09.LHb.4"     "C.10.MHb.1"     "C.11.MHb.1.2"   "C.13.LHb.4"    
# # [9] "C.14.MHb.1"     "C.16.LHb.6"     "C.18.LHb.1.3.4" "C.23.LHb.1"    
# # [13] "C.24.DD_LHb"    "C.30.DD_LHb"    "C.33.LHb.1.3"   "C.36.MHb.3"    
# # [17] "C.40.LHb.4"  
# no_hb_clust <- grep("MHb|LHb", all_clusters, value = TRUE, invert = TRUE)
# no_hb_clust
# 
# # ensure all clusters are included
# new_levels <- c(hb_clusters, setdiff(all_clusters, hb_clusters))
# new_levels
# # Apply the new order to Seurat object identities
# SeuratOBJ <- SetIdent(SeuratOBJ, value = factor(Idents(SeuratOBJ), levels = new_levels))
# 
# message("Renamed clusters:")
# levels(SeuratOBJ)
# 
# 
# message("Clusters renamed and sorted done!")
# 
# 
# ## =============================================================================
# ## Assign new "cluster_ann" column to Seurat meta-data for visualizations
# colnames(SeuratOBJ@meta.data)
# cluster_ann <- as.vector(Idents(SeuratOBJ))
# # assign new identities to the Seurat object
# SeuratOBJ$cluster_ann <- cluster_ann
# table(SeuratOBJ[["seurat_clusters"]])
# table(SeuratOBJ[["cluster_ann"]])
# 
# ## =============================================================================
# ## Add 3 meta-cluster as column: MHb, LHb and No-Habenula
# 
# # extract the ident ID for the 3 meta-groups
# LHb_clusters_to_merge <- grep("LHb", hb_clusters, value = TRUE)
# MHb_clusters_to_merge <- grep("MHb", hb_clusters, value = TRUE)
# LHb_clusters_to_merge
# # [1] "C.05.DD_LHb" "C.18.DD_LHb" "C.23.DD_LHb" "C.24.DD_LHb" "C.30.DD_LHb"
# # [6] "C.33.DD_LHb" "C.40.DD_LHb"
# MHb_clusters_to_merge
# # [1] "C.07.DD_MHb" "C.10.DD_MHb" "C.11.DD_MHb" "C.14.DD_MHb" "C.16.DD_MHb"
# # [6] "C.36.DD_MHb"
# no_hb_clust
# 
# # Add meta-data "merged_cluster" with 3 merged clusters classes: LHb, MHb and No-Hb clusters
# current_idents <- as.character(Idents(SeuratOBJ))
# # Assign merged labels
# merged_cluster <- ifelse(current_idents %in% LHb_clusters_to_merge, "LHb_merged",
#                          ifelse(current_idents %in% MHb_clusters_to_merge, "MHb_merged",
#                                 ifelse(current_idents %in% no_hb_clust, "No-Hb_merged", current_idents)))
# unique(merged_cluster)
# 
# # add to new metadata MERGED ident labels for further analysis
# SeuratOBJ$merged_cluster <- merged_cluster
# unique(SeuratOBJ$merged_cluster)
# #[1] "No-Hb_merged" "MHb_merged"   "LHb_merged" 
# 
# ## save RDS
# rds_file_name <- here(outputRDS_Dir, paste0(Seurat_base_name, "_renamed_visium.rds"))
# saveRDS(SeuratOBJ, rds_file_name)
# 
# message("Seurat with clusters renamed and merged saved!")
# 
# 
# ## =============================================================================
# ## Some visualizations: VPlots, DimPlot, Feature, DotPlot ...
# 
# message("Building some plots ...")
# 
# DefaultAssay(SeuratOBJ) <- "RNA"
# ## extract suffix name to give unique name to plots
# seurat_name <- str_extract(Seurat_base_name, pattern = "k[3:4]0\\_C\\.\\w*")
# # seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k40_C.leiden_lsi_r2
# 
# plt1 <- DimPlot(SeuratOBJ, 
#                 label = TRUE, 
#                 reduction = "wnn.umap",
#                 #group.by = "cluster_ann", 
#                 label.size = 3) + 
#     NoLegend() +
#     labs(title = paste0("**WNN Clusters: ", seurat_name))
# 
# tmp_name <- paste0(seurat_name, "_DimPlot_renamed.pdf")
# ggsave(plt1, filename = here(plotDir, tmp_name), height = 6, width = 6)
# message("WNN UMAP done!")
# 
# plt1 <- DimPlot(SeuratOBJ, 
#                 label = TRUE, 
#                 reduction = "wnn.umap",
#                 group.by = "merged_cluster", 
#                 label.size = 3) + 
#     NoLegend() +
#     labs(title = paste0("**WNN Merged Hb Clusters: ", seurat_name))
# 
# tmp_name <- paste0(seurat_name, "_Hb_merged_DimPlot_renamed.pdf")
# ggsave(plt1, filename = here(plotDir, tmp_name), height = 6, width = 6)
# message("WNN UMAP done!")
# 
# 
# ## Violin plots for the features selected
# 
# features <- c("POU4F1", "GPR151") # "TAC3"
#   
# plt1 <- VlnPlot(object = SeuratOBJ,
#                 features = features,
#                 layer = "data", 
#                 group.by = "cluster_ann",
#                 pt.size = 0) +
#   labs(x = paste0("**Clusters from WNN: ", seurat_name)) &
#   theme(text = element_text(size = 8), 
#         axis.text.x= element_text(size = 7), axis.text.y= element_text(size = 7)) 
# 
# tmp_name <- paste0(seurat_name, "_POU4F1_GPR151_Violin_plot.pdf")
# ggsave(plt1, filename = here(plotDir, tmp_name), height = 4, width = 17)
# message("Violin Plots done!")  
# 
# 
# # Visualize co-expression of two features simultaneously
# 
# plt1 <- FeaturePlot(SeuratOBJ, 
#                     features = features, 
#                     reduction = "wnn.umap",
#                     blend = TRUE) +
#   labs(title = paste0("**Clusters from WNN: ", seurat_name)) &
#   theme(text = element_text(size = 8), 
#         axis.text.x= element_text(size = 7), axis.text.y= element_text(size = 7),
#         plot.title=element_text(hjust=0.5)) 
# 
# tmp_name <- paste0(seurat_name, "_POU4F1_GPR151_FeaturePlot.pdf")
# ggsave(plt1, filename = here(plotDir, tmp_name), height = 3, width = 10)
# message("Feature Plots done!")  
# 
# 
# ## Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# plt1 <- DotPlot(SeuratOBJ, 
#                 group.by = "cluster_ann",
#                 features = c(features, "TAC3")) + 
#     RotatedAxis() +
#     labs(title = paste0("**Clusters from WNN: ", seurat_name)) &
#     theme(text = element_text(size = 8), 
#         axis.text.x= element_text(size = 7), axis.text.y= element_text(size = 7),
#         plot.title=element_text(hjust=0.5)) 
# 
# tmp_name <- paste0(seurat_name, "_POU4F1_GPR151_DotPlot.pdf")
# ggsave(plt1, filename = here(plotDir, tmp_name), height = 6, width = 6)
# message("DotPlot done!")  
# 
# message("Process completed!")
# 
# 
# # library("slurmjobs")
# # job_loop(
# #   loops = list(clustering_name = c("x1", "x2", "x3", "x4")),
# #   name = "05_rename_idents",
# #   cores = 2,
# #   create_shell = TRUE,
# #   partition = "katun"
# # )
# 
# library("sessioninfo")
# print('Reproducibility information:')
# Sys.time()
# proc.time()
# options(width = 120)
# session_info()


