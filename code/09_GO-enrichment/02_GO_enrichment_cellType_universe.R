########################################################################
## GO enrichment in cCRE and classified Linked OCR and DARs unlinked
##
## Authors. CSC
## Date. Oct 24, 2025
## Recommended resources on interactive mode: srun --pty --mem=20GB --x11 bash
########################################################################

library("dplyr")
library("purrr")
library("ggplot2")
library("tidyverse")
library("tidyr")
library("stringr")
library("here")
library("getopt")
library("org.Hs.eg.db")
library("clusterProfiler")

#===============================================================================
# resolution_level = "Mid"      # 18 cell-types
#===============================================================================


## Set directory names
input_rawLinks <- here(
    "processed-data",
    "06_peak_calling",
    "14_exploratory_pb_peak_scores_MACS2", # has ensembl id names
    "links_ct_merged",
    "all_links" 
)
inputCSV_merged_classified <- here(
    "processed-data",
    "06_peak_calling",
    "21_overlaping_FDRscores_TopHeatmap"
)
processedDir <- here(
    "processed-data",
    "09_GO-enrichment",
    "02_GO_enrichment"
)
plotDir <- here(
    "plots",
    "09_GO-enrichment",
    "02_GO_enrichment"
)

if (!dir.exists(processedDir)) { dir.create(processedDir, recursive = TRUE) }
if (!dir.exists(plotDir)) { dir.create(plotDir, recursive = TRUE) }

timestamp <- format(Sys.time(), "%Y%m%d_%H%M")


##==============================================================================
## Define universe genes for the different categories
##==============================================================================

message("Loading raw links ...")

## List all files matching the specific clustering resolution level
pattern = paste0("^Mid.*\\.csv$")
lst_peak_files <- list.files(
    path = input_rawLinks,
    pattern = pattern,
    full.names = TRUE
)

message("Found ", length(lst_peak_files), " raw links files:")
basename(lst_peak_files)


message("Making MASTER Entrez Map (Used for initial data conversion)...")

# Combine ALL gene_ids from ALL files to create the master map
master_link_df <- lst_peak_files |>
    purrr::map(~ read.csv(.x) |> select(cluster, gene_id)) |>
    list_rbind()
master_link_genes <- master_link_df |> 
    distinct(gene_id) |> 
    pull(gene_id)
# table(master_link_df$cluster)
# Astrocyte       Endo Excit.Thal Inhib.Thal      LHb.1    LHb.1.3  LHb.1.3.4 
# 334319      48030    1231642     982150     308624      28244     242985 
# LHb.2.7      LHb.4      LHb.7      MHb.1    MHb.1.2      MHb.2      MHb.3 
# 484856     914230      45796     266637     187704     394229      27774 
# Microglia      Oligo        OPC       Thal 
# 126002     419582     167940      76426 
# Perform the comprehensive mapping once
master_entrez_map <- clusterProfiler::bitr(
    master_link_genes,
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db,
    drop = TRUE
)


## Declare functions ===========================================================

load_entrez_cellType_universe <- function(
    clusterRes,  # Could be any of c("broad", "semi_broad", "mid")
    ct,
    lst_peak_paths) 
{
    # clusterRes = "Broad"
    # lst_peak_paths = lst_peak_files
    # ct = "LHb.2"
    
    message("Processing universe for clusterig resolution [", clusterRes, "]")
    #pattern_hb <- str_extract(ct, "(M|L)Hb")
    # Extract the target pattern from ALL file names (e.g., "Astrocyte", "MHb.2")
    target_patterns <- sub("^[^_]+_([^_]+)_.*$", "\\1", basename(lst_peak_paths))

    ## Determine which files to merge based on the cell type (ct) 
 
    if (clusterRes == "broad") { 
        # Identify ALL Habenula files (MHb and LHb) to merge for the 'broad' universe
        is_target_file <- grepl("^(M|L)Hb.*$", target_patterns)
        
    } else if (clusterRes == "semi_broad" || clusterRes == "mid") { 
        # For semi_broad and mid, the universe is ALL linked genes (no filtering by ct)
        is_target_file <- rep(TRUE, length(target_patterns))
        
    } else {
        stop("Invalid clusterRes value provided.")
    }
        
    ## load, combine and extract unique genes
    filtered_files <- lst_peak_paths[is_target_file]
    link_df <- filtered_files |>
        map(~ read.csv(.x) |> select(cluster, gene_id)) |>
        list_rbind()
    # table(link_df$cluster)
    link_genes <- link_df |> 
        distinct(gene_id) |> 
        pull(gene_id)
    # head(link_genes)
    #message("Total raw-links for [", ct , "] - ", clusterRes," level\n", length(link_genes))
    
    # === Entrez ID Mapping Section ===
    entrez_map <- clusterProfiler::bitr(
        link_genes,
        fromType = "ENSEMBL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db,
        drop = TRUE # Only keep mapped IDs
    )
    entrez_universe <- unique(entrez_map$ENTREZID)

    message("Mapped to ", length(entrez_universe), " Entrez IDs (DARs universe).")
    
    return(entrez_map)
    
}


## =============================================================================

## Function to Set go-universe depending on go-test category dataset

select_specific_go_universe <- function(test_name, go_universes) {
    
    if (test_name %in% c("Direct_Regulation_all", "Direct_Regulation_Enriched",
                         "Direct_Regulation_Depleted", "LinkPeaks_OCRs")) {
        selected <- go_universes[[test_name]]
        
    # } else if (test_name == "High_Interest_DARs") {
    #     selected <- go_universes$High_Interest_DARs
    #     
    # } else {
    #     stop(paste("Unknown test category:", test_name))
    }
    
    return(selected)

}



## =============================================================================
## Define GO Classification and Categories to test
## =============================================================================

message(" \n==== GO enrichment test for Classification-3 =====\n")
message("Loading cCRE and ORC file ...")

f_name = "overlaps_linkPeak_DARs_classified_thr_CC0.3_thr_DAR0.1.csv"
## load raw overlaps and 2-shared ct overlaps 
cCRE_df <- read.csv(here(inputCSV_merged_classified, f_name))
nrow(cCRE_df) # 6410
head(cCRE_df)
table(cCRE_df$type_classification)

cCRE_df <- cCRE_df |> dplyr::filter(type_classification=="classification-3")
nrow(cCRE_df) # 3205

## validation
unique(cCRE_df$type_classification)
as.data.frame(table(cCRE_df$type_classification, cCRE_df$category))
# Var1                    Var2 Freq
# 1 classification-3              Linked OCR  420
# 2 classification-3 Linked_DAR (-) depleted  657
# 3 classification-3 Linked_DAR (+) enriched  706
# 4 classification-3         Non-significant  354
# 5 classification-3            Unlinked DAR 1068

go_directReg = c("Linked_DAR (+) enriched", "Linked_DAR (-) depleted")
go_directReg_enriched = c("Linked_DAR (+) enriched")
go_directReg_depleted = c("Linked_DAR (-) depleted")
go_all_Linked = c("Linked_DAR (+) enriched", "Linked_DAR (-) depleted", "Linked OCR")
go_linked_OCRs = "Linked OCR"
# go_only_DARs = "Unlinked DAR"

lst_go_tests <- list(
    "Direct_Regulation_Enriched" = go_directReg_enriched,
    "Direct_Regulation_Depleted" = go_directReg_depleted,  
    "Direct_Regulation_all" = go_directReg,
    "LinkPeaks_OCRs" = go_all_Linked #,
    # "High_Interest_DARs" = go_only_DARs
)

# Initialize lst_go_universes with the master map (used by select_specific_go_universe)
lst_go_universes <- purrr::map(names(lst_go_tests), ~list(entrez_map = master_entrez_map))
names(lst_go_universes) <- names(lst_go_tests)
names(lst_go_universes)
# [1] "Direct_Regulation_Enriched" "Direct_Regulation_Depleted"
# [3] "Direct_Regulation_all"      "LinkPeaks_OCRs"

message("Master Map generated. Total unique Ensembl IDs mapped: ", nrow(master_entrez_map))


## =============================================================================
## Summarize stats for Hb clusters (validations)
## =============================================================================

## Generate Summary Tables for cCRE Categories Across Habenula (MHb/LHb) Clusters
purrr::map(lst_go_tests, function(.x) {
    # Filter the main data frame (cCRE_df) based on the current universe (.x)
    #    and restrict to Habenula (MHb/LHb) cell types
    cCRE_universe_df <- cCRE_df |>
        dplyr::filter(category %in% .x) |>
        dplyr::filter(grepl("MHb|LHb", cell_type))
    
    cat("\n--- Total Counts by Category ---\n")
    print(table(cCRE_universe_df$category))
    
    # Return the count data frame
    cCRE_universe_df |> dplyr::count(cell_type)
    
})


## =============================================================================
## Combine all filtered data frames into a named list for iteration

## iterate over the values of lst_go_tests (the vectors of categories)
list_df_cCRE <- lst_go_tests |>
    purrr::map(\(filter_vector) {
        cCRE_df |>
            dplyr::filter(category %in% filter_vector)
    })

## Check summary of list content
map_int(list_df_cCRE, nrow)
# Direct_Regulation_Enriched Direct_Regulation_Depleted 
# 706                        657 
# Direct_Regulation_all LinkPeaks_OCRs 
# 1363                       1783 


##==============================================================================
## Make DE_entrez data frame accordingly with the go-test to compute


make_DE_entrez_df <- function(
        test_name, # ge: "Direct_Regulation_all"
        cCRE_df,
        entrez_map)
    {
    
    if (test_name=="Direct_Regulation_all" || test_name=="Direct_Regulation_Enriched" || test_name=="Direct_Regulation_Depleted") {
     
        ########## Test "Direct_Regulation_all (s) " ##########
        
        DE_entrez <- cCRE_df |>
            # use this if duplicates are expected (because multiple peaks link to the same gene)
            left_join(entrez_map,
                      by = c("gene_id" = "ENSEMBL"), 
                      relationship = "many-to-many") |> 
            filter(!is.na(ENTREZID)) |>
            mutate(
                ## keep the original category direction
                category_direction = case_when(
                    category == "Linked_DAR (+) enriched" ~ "up",
                    category == "Linked_DAR (-) depleted" ~ "down"
            )) |>
            distinct(ENTREZID, .keep_all = TRUE)
        
    } else if ((test_name=="LinkPeaks_OCRs")) {
        
        ########## Test "LinkPeaks_OCRs" ########## 
        
        DE_entrez <- cCRE_df |>
            left_join(entrez_map,
                      by = c("gene_id" = "ENSEMBL"),
                      relationship = "many-to-many") |>
            filter(!is.na(ENTREZID)) |>
            mutate(
                category_direction = case_when(
                    category == "Linked_DAR (+) enriched" ~ "up",
                    category == "Linked_DAR (-) depleted" ~ "down",
                    category == "Linked OCR" & CCscore > 0 ~ "up",
                    category == "Linked OCR" & CCscore < 0 ~ "down"
                ),
                cell_type_clean = gsub("\\.", "-", cell_type)# ,
            ) |>
            distinct(ENTREZID, .keep_all = TRUE)
    
    } else if ((test_name=="High_Interest_DARs")) {
    
        ########## Test "DARs" ########## 
        
        DE_entrez <- cCRE_df |>
            left_join(entrez_map,
                      by = c("gene_id" = "ENSEMBL"),
                      relationship = "many-to-many") |>
            filter(!is.na(ENTREZID)) |>
            mutate(
                category_direction = case_when(
                    category == "Unlinked DAR" & logFC > 0  ~ "up",
                    category == "Unlinked DAR" & logFC < 0  ~ "down",
                    TRUE ~ "neutral"
                ),
                cell_type_clean = gsub("\\.", "-", cell_type)# ,
            ) |>
            distinct(ENTREZID, .keep_all = TRUE)
        
    } else {
        
        stop("Invalid test_name. Must be one of exixting tests!")
        
    }
    
    return(DE_entrez)
    
}
           

## =============================================================================
## Loop over each GO test by cell-type definition (UPDATED to use enrichGO per group)
## =============================================================================
## Run enrichment separately for each distinct value of DE_group
## Each GO enrichment is computed per cell-type (or per merged cell-type group) - not as one pooled dataset.

ont_list <- c("CC","BP","MF")
names(ont_list) <- ont_list
clustering_levels <- c("broad", "semi_broad", "mid")
# clustering_levels <- c("broad", "semi_broad")

# Tests: names(lst_go_tests)
# [1] "Direct_Regulation_Enriched" "Direct_Regulation_Depleted"
# [3] "Direct_Regulation_all"      "LinkPeaks_OCRs"            
# [5] "High_Interest_DARs"   

# Initialize list to store all results across all levels and tests
all_go_results <- list()

for (test_go in names(lst_go_tests)) {
    # testing: test_go = names(lst_go_tests)[[1]]
    
    message("\n==============================")
    message("Running GO enrichment for: ", test_go)
    message("==============================")
    
    ## Prepare category-specific subset
    cCRE_subset <- list_df_cCRE[[test_go]]
    if (is.null(cCRE_subset) || nrow(cCRE_subset) == 0) {
        message("No data found for ", test_go, " — skipping.")
        next
    }
    message("Subset size: ", nrow(cCRE_subset))
    # Subset size: 706
    
    ## Use the MASTER Entrez Map for initial gene ID conversion regardless of universe
    go_univ <- select_specific_go_universe(test_name = test_go, go_universes = lst_go_universes)
    ## Prepare entrez df
    DE_entrez <- make_DE_entrez_df(
        test_name = test_go,
        cCRE_df = cCRE_subset,
        entrez_map = go_univ$entrez_map
    )
    
    ## Loop over clustering resolution levels  =================================
    
    for (res_level in clustering_levels) {
        
        message("\n--- Running Clustering Resolution: ", res_level, " ---")
        
        if (hb_merged_ct == "mid") message("Using full mid-resolution clustering - no merging applied.")
        
        # # Load the universe based on the current clustering resolution
        # universe_map <- load_entrez_cellType_universe(
        #     clusterRes = hb_merged_ct, 
        #     lst_peak_paths = lst_peak_files
        # )
        # entrez_universe <- unique(universe_map$ENTREZID) # Dynamic universe
        
        ## Define Cell-Type Groups based on resolution
        DE_entrez_grouped <- DE_entrez |>
            mutate(
                cell_type_broad = case_when(
                    res_level == "broad" & grepl("^(LHb|MHb)", cell_type) ~ "Hb",
                    res_level == "semi_broad" & grepl("^LHb", cell_type) ~ "LHb",
                    res_level == "semi_broad" & grepl("^MHb", cell_type) ~ "MHb",
                    # # For mid, keep the original cell_type
                    res_level == "mid" ~ cell_type,
                    TRUE ~ cell_type # Fallback for non-Hb clusters if they were included
                ),
                DE_group = cell_type_broad
            ) |>
            distinct(DE_group, ENTREZID, .keep_all = TRUE) |>
            group_by(DE_group) |> 
            filter(n() >= 10) |> 
            ungroup()
        
        ## verifications
        ## table(DE_entrez_grouped$DE_group)
        # Astrocyte Excit.Thal Inhib.Thal      LHb.1  LHb.1.3.4    LHb.2.7      LHb.4 
        # 18        102        192         25         10         69         77 
        # MHb.1    MHb.1.2      MHb.2      Oligo 
        # 27         11        102         17
        # Filter the tibble to only show rows where DE_group is "LHb.1"
        # DE_entrez_grouped |>
        #     dplyr::filter(DE_group == "LHb.4") |>
        #     head()
        
        if (nrow(DE_entrez_grouped) == 0 || n_distinct(DE_entrez_grouped$DE_group) == 0) {
            message("No groups with >=10 genes for ", test_go, " at ", res_level, " — skipping.")
            next
        }
        
        # Extract list of unique cell groups/types to test
        cell_groups_to_test <- unique(DE_entrez_grouped$DE_group)
        
        # **Dynamically load the GO Universe for the current resolution level**
        # Since 'load_entrez_cellType_universe' function defines the universe 
        # based on the clustering resolution ('broad' is Habenula-only, 'mid/semi_broad' is all), 
        # I will use this resolution-specific universe for the enrichment of all groups
        # defined at that resolution. This correctly implements the 'cell-type has its own go-universe'
        
        universe_map <- load_entrez_cellType_universe(
            clusterRes = res_level, 
            lst_peak_paths = lst_peak_files 
        )
        entrez_universe <- unique(universe_map$ENTREZID)
        
        message("Using resolution-specific universe of ", length(entrez_universe), " genes for **ALL** groups at [", res_level, "].")
        
        ## checks: diagnostics vs universe
        n_unique <- length(unique(DE_entrez_tmp$ENTREZID))
        shared <- length(intersect(entrez_universe, DE_entrez_tmp$ENTREZID))
        message("Unique ENTREZ IDs: ", n_unique)
        message("Overlap with universe: ", shared, 
                " (", round(100 * shared / length(entrez_universe), 1), "%)")
        
        message("Running GO enrichment using background universe of ", length(entrez_universe), " genes...")
        
        ## GO terms summarized by broad cell types (LHb, MHb, Astrocyte, etc.) rather than by fine subclusters.
        
        suppressWarnings({
            go_result <- map(ont_list, ~compareCluster(
                ENTREZID ~ DE_group,
                data = DE_entrez_tmp,
                OrgDb = org.Hs.eg.db,
                fun = enrichGO,
                universe = entrez_universe,
                ont = .x,
                pAdjustMethod = "BH",
                pvalueCutoff = 0.1,
                qvalueCutoff = 0.2,
                readable = TRUE
            ))
        })
        
        # go_result for diagnostics: how many GO terms per ontology
        message("\nSummary of enriched terms:")
        map2(names(go_result), go_result, function(nm, x) {
            n_terms <- if (is.null(x)) 0 else nrow(x@compareClusterResult)
            message(sprintf("%s: %s terms", nm, n_terms))
        })
        
        
        ## testing: Preview & visualize top results for Biological Process
        if (!is.null(go_result$BP) && nrow(go_result$BP@compareClusterResult) > 0) {
            message("Top enriched BP terms:")
            print(head(go_result$BP@compareClusterResult, 5))
            
            p <- dotplot(go_result$BP, showCategory = 15) +
                ggtitle(paste(res_level, " GO BP enrichment"), 
                        subtitle = paste("Category: ", test_go))
                # ggtitle(paste("GO BP enrichment (", test_go, " — ", res_level, ")", sep = ""))
            
            plot_base <- sprintf("%s_GO_BP_%s_%s", res_level, test_go, timestamp)
            plot_pdf <- here(plotDir, paste0(plot_base, ".pdf"))
            ggsave(plot_pdf, plot = p, width = 8, height = 6)
            message("Saved plots: ", plot_pdf)
            
        } else {
            message("No significant BP terms for ", test_go)
        }
        
        ## Save results
        go_rds_name <- here(processedDir, sprintf("%s_GO_results_%s_%s.rds", res_level, test_go, timestamp))
        saveRDS(go_result, go_rds_name)
        message("Saved RDS: ", go_rds_name)
        
        ## Save summary table (flattened BP results if available)
        if (!is.null(go_result$BP) && nrow(go_result$BP@compareClusterResult) > 0) {
            bp_df <- go_result$BP@compareClusterResult
            go_csv_name <- here(processedDir, sprintf("%s_GO_results_BP_%s_%s.csv", res_level, test_go, timestamp))
            write.csv(bp_df, go_csv_name, row.names = FALSE)
            message("Saved CSV: ", go_csv_name)
        }
    
        message("Completed GO enrichment for ", res_level, " [test: ", test_go, "]")

    }
    
}


message("GO-Enrichment completed!")



## =========
# 
# # Remove NULL entries (e.g., CC = NULL)
# go_result_valid <- discard(go_result, is.null)
# 
# ## convert to table & extract compareClusterResult from each valid ontology and tag with ontology name
# compare_clus <- map2_dfr(go_result_valid, names(go_result_valid), function(x, nm) {
#     x@compareClusterResult |> mutate(ONTOLOGY = nm)
# })
# 
# # verify classes that survive
# map(go_result, function(x) {
#     if (is.null(x)) "NULL" else nrow(x@compareClusterResult)
# })
# 
# compare_clus |> count(DE_class_cluster, ONTOLOGY)
# # DE_class_cluster ONTOLOGY n
# # 1   LHb-1-3-4_down       MF 8
# # 2       LHb-2-7_up       BP 2
# # 3       MHb-2_down       MF 2
# 
# ## Save 
# f_name <- here(processedDir, sprintf("GO_compare_clus_%s.rds", peaks_classification_name))
# saveRDS(compare_clus, f_name)
# f_name <- here(processedDir, sprintf("GO_results_%s.csv", peaks_classification_name))
# write.csv(compare_clus, f_name, row.names = FALSE)
# 
# 
# #### dot plots ####
# 
# # Iterate over valid GO results (skip NULL or empty ones)
# go_result_valid <- discard(go_result, function(x) {
#     is.null(x) || nrow(x@compareClusterResult) == 0
# })
# 
# length(go_result_valid)
# map(go_result_valid, ~nrow(.@compareClusterResult))
# # pdf(f_name, width = 10, height = 10)
# # print(dotplot(go_result_valid[[2]])) # pass!
# # dev.off()
# 
# f_name = here(processedDir, sprintf("GO_dotplot_%s.pdf", peaks_classification_name))
# pdf(f_name, width = 10, height = 10)
# 
# ## one page per ontology 
# if (length(go_result_valid) == 0) {
#     message("No valid GO results to plot.")
#     dev.off()
# } else {
#     walk2(go_result_valid, names(go_result_valid), function(x, nm) {
#         # Check if the result for this ontology contains enough data to attempt plotting
#         if(nrow(x@compareClusterResult) > 0) {
#             
#             message("Plotting ontology: ", nm, " (", nrow(x@compareClusterResult), " terms)")
#             
#             # Use tryCatch to prevent a single failing plot from crashing the entire PDF
#             tryCatch({
#                 
#                 p <- dotplot(x,
#                              x = "DE_class_cluster", 
#                              showCategory = 5, 
#                              label_format = 60) +
#                     ggtitle(paste("GO Enrichment:", nm)) +
#                     theme_bw(base_size = 12) +
#                     theme(
#                         axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8),
#                         axis.text.y = element_text(size = 8),
#                         plot.title = element_text(hjust = 0.5, face = "bold")
#                     )
#                 print(p) 
#             }, error = function(e) {
#                 message("Skipping plot for ", nm, " due to error: ", conditionMessage(e))
#             })
#         } else {
#             message("Skipping plot for ", nm, ": No terms remaining after filtering.")
#         }
#     })
#     
#     dev.off()
#     
# }




## Reproducibility information
library(sessioninfo)
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()




