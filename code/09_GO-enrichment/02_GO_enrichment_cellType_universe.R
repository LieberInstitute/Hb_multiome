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
# inputCSV_merged_overlaps <- here(
#     "processed-data",
#     "06_peak_calling",
#     "19_Linkage_DARs_analysis"
# )
input_rawLinks <- here(
    "processed-data",
    "06_peak_calling",
    "13_pseudobulk_LinkPeaks_MACS2_split_ct",
    "links_ct_merged" # refers only to merged peaks - done to have idential genomic regions among datasets 
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

## For Direct-Regulation Categories (all measurable genes prefiltered) =========

# universe_linkPeaks <- read.csv(here(inputCSV_merged_overlaps, "ALL_LinkPeaks_signif_FDR02.csv")) |>
#     filter(!is.na(gene_id)) |>
#     distinct()

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


message("====================================================================\n")
message("Making Gene Universe for LinkPeaks Categories ...")

# message("Universe: LinkPeaks-based regulatory genes = ", nrow(universe_linkPeaks))
# head(universe_linkPeaks)
# 
# entrez_map_linkPeaks <- bitr(
#     unique(universe_linkPeaks$gene_id),
#     fromType = "ENSEMBL",
#     toType = "ENTREZID",
#     OrgDb = org.Hs.eg.db
# )
# entrez_universe_linkPeaks <- unique(entrez_map_linkPeaks$ENTREZID)
# message("Mapped to ", length(entrez_universe_linkPeaks), " Entrez IDs (LinkPeaks universe).")

load_cellType_universe <- function(
    clusterRes,  # Could be any of c("broad", "semi_broad", "mid")
    ct, 
    lst_peak_paths
    ) 
{
    # clusterRes = "Mid"
    # lst_peak_paths = lst_peak_files
    # ct="MHb.2"
    # ct="MHb.1.2"
    # ct="LHb.2"
    # ct="Oligo"

    message("Processing universe for ct: ", ct)
    pattern_hb <- str_extract(ct, "(M|L)Hb")
    # Extract the target pattern from ALL file names (e.g., "Astrocyte", "MHb.2")
    target_patterns <- sub("^[^_]+_([^_]+)_.*$", "\\1", basename(lst_peak_paths))

    ## Determine which files to merge based on the cell type (ct) ==============
    
    ## If pattern_hb is NA, it's an OTHER cell type (e.g., "Oligo")
    if (is.na(pattern_hb)) { 
        # Identify files matching the exact non-Hb cell type
        is_target_file <- grepl(paste0("^", ct, "$"), target_patterns)
    
    } else { ## # It is Habenula cell type (e.g., "MHb.2", "LHb.1.3")
        
        if (clusterRes=="broad") { ## we only have one Hb cell-type
            # Identify all "MHb" and "LHb" 
            is_target_file <- grepl("^(M|L)Hb.*$", target_patterns)
                        
        } else if (clusterRes=="semi-broad") {
            # Identify "MHb" or "LHb"
            if (pattern_hb=="MHb") {
                is_target_file <- grepl("^MHb.*$", target_patterns)
                
            } else (pattern_hb=="LHb") {
                is_target_file <- grepl("^LHb.*$", target_patterns)
                
            }
        } else if  (clusterRes=="mid") {
            # Identify files matching the exact Habenula sub-cluster
            is_target_file <- grepl(paste0("^", ct, "$"), target_patterns)
        }

    }
    ## =========================================================================
        
    # load, combine and extract unique genes
    filtered_files <- lst_peak_paths[is_target_file]
    link_df <- filtered_files |>
        map(~ read.csv(.x) |> select(cluster, gene)) |>
        list_rbind()
    table(link_df$cluster)
    link_genes <- link_df |> 
        distinct(gene) |> 
        pull(gene)
    head(link_genes)
    
    message("Total raw-links for [", ct , "] - ", clusterRes," level\n", length(link_genes))
    return(link_genes)
    
}



## For DARs Category (all measurable genes prefiltered) ========================
## all DARs with fdr < 0.1 (first filter)

# universe_dars <- read.csv(here(inputCSV_merged_overlaps, "ALL_DARs_signif_FDR01.csv")) |>
#     filter(!is.na(gene_id)) |>
#     distinct()
# 
# message("====================================================================\n")
# message("Universe: DARs-based regulatory genes = ", nrow(universe_dars))
#
# ## raw links reference to pull gene names
# raw_links_for_dars_path <- here("processed-data", "06_peak_calling", "14_exploratory_pb_peak_scores_MACS2", "links_ct_merged", "all_links")
# pattern <- paste0("^Mid.*\\.csv$")
# lst_peak_files <- list.files(
#     path = raw_links_for_dars_path,
#     pattern = pattern,
#     full.names = TRUE
# )
# #lst_peak_files = list.files(path = input_cvsDir)
# message("Link peak-genes raw files found:")
# lst_peak_files
# 
# raw_links_df <- lst_peak_files |>
#     map_dfr(read_csv, .id = "source_file",
#             show_col_types = FALSE)
# 
# universe_dars_cat <- raw_links_df # need to merge with universe_dars_cat to pull gene names --- in progress

# entrez_map_dars <- bitr(
#     unique(universe_dars$gene_id),
#     fromType = "ENSEMBL",
#     toType = "ENTREZID",
#     OrgDb = org.Hs.eg.db
# )
# entrez_universe_dars <- unique(entrez_map_dars$ENTREZID)
# message("Mapped to ", length(entrez_universe_dars), " Entrez IDs (DARs universe).")

## Combine all universes into one structured list

lst_go_universes <- list(
    Direct_Regulation = list(
        df = universe_linkPeaks,
        entrez_map = entrez_map_linkPeaks,
        entrez_universe = entrez_universe_linkPeaks
    ),
    LinkPeaks_OCRs = list(
        df = universe_linkPeaks,
        entrez_map = entrez_map_linkPeaks,
        entrez_universe = entrez_universe_linkPeaks
    ),
    High_Interest_DARs = list(
        df = universe_linkPeaks, # temporal assignation (universe_dars)
        entrez_map = entrez_map_linkPeaks, # entrez_map_dars
        entrez_universe = entrez_universe_linkPeaks #entrez_universe_dars
    )
)

names(lst_go_universes)
# [1] "Direct_Regulation"  "LinkPeaks_OCRs"     "High_Interest_DARs"


## =============================================================================

## Function to Set go-universe depending on go-test category dataset

select_specific_go_universe <- function(test_name, go_universes) {
    
    if (test_name %in% c("Direct_Regulation_all", "Direct_Regulation_Enriched",
                         "Direct_Regulation_Depleted", "LinkPeaks_OCRs")) {
        selected <- go_universes$Direct_Regulation
        
    } else if (test_name == "High_Interest_DARs") {
        selected <- go_universes$High_Interest_DARs
        
    } else {
        stop(paste("Unknown test category:", test_name))
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
go_only_DARs = "Unlinked DAR"

lst_go_tests <- list(
    "Direct_Regulation_Enriched" = go_directReg_enriched,
    "Direct_Regulation_Depleted" = go_directReg_depleted,  
    "Direct_Regulation_all" = go_directReg,
    "LinkPeaks_OCRs" = go_all_Linked,
    "High_Interest_DARs" = go_only_DARs
)



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
    
    # # Generate summary tables
    # cat("\n--- Counts by Cell Type and Category ---\n")
    # print(table(cCRE_universe_df$cell_type, cCRE_universe_df$category))
    # 
    # cat("\n--- Total Counts by Cell Type (Cluster) ---\n")
    # print(table(cCRE_universe_df$cell_type))
    
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
# High_Interest_DARs 
# 1068 


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
                #),cell_type_clean = gsub("\\.", "-", cell_type)# ,
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
## Loop over each GO test definition
## =============================================================================
## Run enrichment separately for each distinct value of DE_group
## Each GO enrichment is computed per cell-type (or per merged cell-type group) - not as one pooled dataset.

ont_list <- c("CC","BP","MF")
names(ont_list) <- ont_list
clustering_levels <- c("broad", "semi_broad", "mid")

# Tests: names(lst_go_tests)
# [1] "Direct_Regulation_Enriched" "Direct_Regulation_Depleted"
# [3] "Direct_Regulation_all"      "LinkPeaks_OCRs"            
# [5] "High_Interest_DARs"   

for (test_go in names(lst_go_tests)) {
    # testing: test_go = names(lst_go_tests)[[1]]
    
    
    message("\n==============================")
    message("Running GO enrichment for: ", test_go)
    message("==============================")
    
    
    ## Select GO universe and mappings for this test
    go_univ <- select_specific_go_universe(test_name = test_go, go_universes = lst_go_universes)
    entrez_universe <- go_univ$entrez_universe
    entrez_map <- go_univ$entrez_map

    ## Prepare category-specific subset (use the correct one per test)
    cCRE_subset <- list_df_cCRE[[test_go]]
    if (is.null(cCRE_subset) || nrow(cCRE_subset) == 0) {
        message("No data found for ", test_go, " — skipping.")
        next
    }
    message("Subset size: ", nrow(cCRE_subset))
    
    ## prepare entrez df
    DE_entrez <- make_DE_entrez_df(
        test_name = test_go,
        cCRE_df = cCRE_subset,
        entrez_map = entrez_map
    )
    
    ## Merge Hb sub-types into broad groups (LHb, MHb) ==========================
    for (hb_merged_ct in clustering_levels) {
        
        message("\n--- Running Hb merge mode: ", hb_merged_ct, " ---")
        if (hb_merged_ct == "mid") message("Using full mid-resolution clustering - no merging applied.")
        
        DE_entrez_tmp <- DE_entrez |>
            mutate(
                cell_type_broad = case_when(
                    hb_merged_ct == "broad" & grepl("^(LHb|MHb)", cell_type) ~ "Hb",
                    hb_merged_ct == "semi_broad" & grepl("^LHb", cell_type) ~ "LHb",
                    hb_merged_ct == "semi_broad" & grepl("^MHb", cell_type) ~ "MHb",
                    hb_merged_ct == "mid" ~ cell_type,
                    TRUE ~ cell_type
                ),
                DE_group = cell_type_broad
            ) |>
            distinct(DE_group, ENTREZID, .keep_all = TRUE) |>
            group_by(DE_group) |> 
            filter(n() >= 10) |> 
            ungroup()
        
        message("Final DE_entrez dimensions: ", nrow(DE_entrez_tmp), " rows, ", 
                length(unique(DE_entrez_tmp$DE_group)), " cell groups")
        
        ## test after filtering
        if (nrow(DE_entrez_tmp) == 0 || n_distinct(DE_entrez_tmp$DE_group) == 0) {
            message("No groups with >=10 genes for ", test_go, " — skipping enrichment.")
            next
        }
        
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
                ggtitle(paste(hb_merged_ct, " GO BP enrichment"), 
                        subtitle = paste("Category: ", test_go))
                # ggtitle(paste("GO BP enrichment (", test_go, " — ", hb_merged_ct, ")", sep = ""))
            
            plot_base <- sprintf("%s_GO_BP_%s_%s", hb_merged_ct, test_go, timestamp)
            plot_pdf <- here(plotDir, paste0(plot_base, ".pdf"))
            ggsave(plot_pdf, plot = p, width = 8, height = 6)
            message("Saved plots: ", plot_pdf)
            
        } else {
            message("No significant BP terms for ", test_go)
        }
        
        ## Save results
        go_rds_name <- here(processedDir, sprintf("%s_GO_results_%s_%s.rds", hb_merged_ct, test_go, timestamp))
        saveRDS(go_result, go_rds_name)
        message("Saved RDS: ", go_rds_name)
        
        ## Save summary table (flattened BP results if available)
        if (!is.null(go_result$BP) && nrow(go_result$BP@compareClusterResult) > 0) {
            bp_df <- go_result$BP@compareClusterResult
            go_csv_name <- here(processedDir, sprintf("%s_GO_results_BP_%s_%s.csv", hb_merged_ct, test_go, timestamp))
            write.csv(bp_df, go_csv_name, row.names = FALSE)
            message("Saved CSV: ", go_csv_name)
        }
    
        message("Completed GO enrichment for ", hb_merged_ct, " [test: ", test_go, "]")

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




