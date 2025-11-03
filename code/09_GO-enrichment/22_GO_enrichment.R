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
library("patchwork")
library("tidyverse")
library("tidyr")
library("stringr")
library("here")
library("getopt")
library("org.Hs.eg.db")
library("clusterProfiler")
library("rrvgo")
library("ComplexHeatmap")

#===============================================================================
# resolution_level = "Mid"      # 18 cell-types
#===============================================================================


## Set directory names
inputCSV_merged_overlaps <- here(
    "processed-data",
    "06_peak_calling",
    "19_Linkage_DARs_analysis"
)
inputCSV_merged_classified <- here(
    "processed-data",
    "06_peak_calling",
    "21_overlaping_FDRscores_TopHeatmap"
)
processedDir <- here(
    "processed-data",
    "09_GO-enrichment",
    "22_GO_enrichment"
)
plotDir <- here(
    "plots",
    "09_GO-enrichment",
    "22_GO_enrichment"
)

if (!dir.exists(processedDir)) { dir.create(processedDir, recursive = TRUE) }
if (!dir.exists(plotDir)) { dir.create(plotDir, recursive = TRUE) }


##==============================================================================
## defining universe from merged overlaps

# table containing all merged overlaps from 19 directory: all measurable genes

##  Linkable genes universe
universe_linkPeaks <- read.csv(here(inputCSV_merged_overlaps, "ALL_LinkPeaks_signif_FDR02.csv"))
nrow(linkPeaks_results_all) # 10948 ok

universe_overlaps_df <- read.csv(here(inputCSV_merged_overlaps, "Overlaps_LinkPeak_DARs_FDR0.2.csv"))
colnames(universe_overlaps_df)

head(universe_overlaps_df)

nrow(universe_overlaps_df)
# 23691

universe_df1 <- universe_overlaps_df |>
    filter(!is.na(gene_id)) |>
    distinct()
nrow(universe_df1)
# 23666

## Maybe Try a category-matched background:
# Global measurable background (current approach): all genes with any overlap link or DAR overlap (current) / for high signif vs all possible

# Linkable genes universe (more focused): all genes that appear in cCRE_df (genes with any LinkPeak correlation, regardless of DAR) / Linked OCR, Direct_Regulation, and High_Interest_DARs

# Cell-type–specific universe: each cell type, use only genes expressed/linked in that cell type as the background / increases biological relevance and specificity 



## ======

## Testing universe.  Map these to ENTREZ IDs
entrez_map <- bitr(universe_df1$gene_id,
                        fromType = "ENSEMBL",
                        toType = "ENTREZID",
                        OrgDb = org.Hs.eg.db) 

entrez_universe <- unique(entrez_map$ENTREZID)
message("Universe size: ", length(entrez_universe))
# 5817
# Warning message:
#     In bitr(universe_df1$gene_id, fromType = "ENSEMBL", toType = "ENTREZID",  :
#                 1.33% of input gene IDs are fail to map...


## =====

message("Loading cCRE and ORC file ...")

f_name = "overlaps_linkPeak_DARs_classified_thr_CC0.3_thr_DAR0.1.csv"
## load raw overlaps and 2-shared ct overlaps 
cCRE_df <- read.csv(here(inputCSV_merged_classified, f_name))
nrow(cCRE_df) # 6410
head(cCRE_df)
table(cCRE_df$type_classification)

cCRE_df <- cCRE_df |> filter(type_classification=="classification-3")

## validation
unique(cCRE_df$type_classification)
as.data.frame(table(cCRE_df$type_classification, cCRE_df$category))
# Var1                    Var2 Freq
# 1 classification-3              Linked OCR  420
# 2 classification-3 Linked_DAR (-) depleted  657
# 3 classification-3 Linked_DAR (+) enriched  706
# 4 classification-3         Non-significant  354
# 5 classification-3            Unlinked DAR 1068

peaks_classification_name = "peaks_classification3_hb"


## =============================================================================
## define GO universe 

## Primary Interest.** Canonical enhancer or promoter activity
go_directReg = c(
    "Linked_DAR (+) enriched",
    "Linked_DAR (-) depleted"
)

# add to the test
# go_directReg_enriched = c(
#     "Linked_DAR (+) enriched"
# )
# go_directReg_depleted = c(
#     "Linked_DAR (-) depleted"
# )

## Primary + Secondary Interest.** Regulatory link exists, but the element is not a strong DAR
go_Primary_Secondary_Interest = c(
    "Linked_DAR (+) enriched",
    "Linked_DAR (-) depleted",
    "Linked OCR"
)

## **High Interest.** Potential for distal regulation or non-coding targets
go_High_Interest_DARs = c(
    "Unlinked DAR"
)

# ## Secondary Interest.** Regulatory link exists, but not a strong DAR
# go_Secondary_Linked_OCR = c(
#     "Linked OCR"
# )

lst_go_tests <- list(
    "Direct_Regulation" = go_directReg,
    "Primary_Secondary_Interest" = go_Primary_Secondary_Interest,
    "High_Interest_DARs" = go_High_Interest_DARs
)

# Display the resulting list
lst_go_tests

## =========/


## Summarize stats for Hb clusters
## Generate Summary Tables for cCRE Categories Across Habenula (MHb/LHb) Clusters

purrr::map(lst_go_tests, function(.x) {
    # Filter the main data frame (cCRE_df) based on the current universe (.x)
    #    and restrict to Habenula (MHb/LHb) cell types
    cCRE_universe_df <- cCRE_df |>
        dplyr::filter(category %in% .x) |>
        dplyr::filter(grepl("MHb|LHb", cell_type))
    
    # Generate summary tables
    cat("\n--- Counts by Cell Type and Category ---\n")
    print(table(cCRE_universe_df$cell_type, cCRE_universe_df$category))
    
    cat("\n--- Total Counts by Cell Type (Cluster) ---\n")
    print(table(cCRE_universe_df$cell_type))
    
    cat("\n--- Total Counts by Category ---\n")
    print(table(cCRE_universe_df$category))
    
    # Return the count data frame
    cCRE_universe_df |> dplyr::count(cell_type)
    
})


## checks
colnames(cCRE_df)
table(cCRE_df$category)

## subset the specific go_universe

names(lst_go_tests)

cCRE_universe_Direct_Regulation <- cCRE_df |> 
    filter(category %in% lst_go_tests[["Direct_Regulation"]]) 
message("Direct_Regulation n = ", nrow(cCRE_universe_Direct_Regulation))

cCRE_universe_Primary_Secondary_Interest <- cCRE_df |> 
    filter(category %in% lst_go_tests[["Primary_Secondary_Interest"]]) 
message("Primary_Secondary_Interest n = ", nrow(cCRE_universe_Primary_Secondary_Interest))

cCRE_universe_df_High_Interest_DARs <- cCRE_df |> 
    filter(category %in% lst_go_tests[["High_Interest_DARs"]]) 
message("High_Interest_DARs n = ", nrow(cCRE_universe_df_High_Interest_DARs))


## Combine all filtered data frames into a named list for iteration
list_df_cCRE <- list(
    Direct_Regulation = cCRE_universe_Direct_Regulation,
    Primary_Secondary_Interest = cCRE_universe_Primary_Secondary_Interest,
    High_Interest_DARs = cCRE_universe_df_High_Interest_DARs
)
## Check summary of list content
map_int(list_df_cCRE, nrow)
# Direct_Regulation Primary_Secondary_Interest 
# 1363                       1783 
# High_Interest_DARs 
# 1068 


##==============================================================================
## Merge back and define DE_class + DE_class_cluster


make_DE_entrez_df <- function(
        test_name="Direct_Regulation",
        cCRE_df,
        entrez_map)
    {
    
    if (test_name=="Direct_Regulation") {
     
        ########## Test "Direct_Regulation" ##########
        
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
                ),
                cell_type_clean = gsub("\\.", "-", cell_type)# ,
            ) |>
            distinct(ENTREZID, .keep_all = TRUE)
        
        } else if ((test_name=="Primary_Secondary_Interest")) {
            
            ########## Test "Primary_Secondary_Interest" ########## 
            
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
        
                stop("Invalid test_name. Must be one of: 'Direct_Regulation', 'Primary_Secondary_Interest'")
        
        }
    
    return(DE_entrez)
    
}
           

## =============================================================================
## Loop over each GO test definition
## =============================================================================

ont_list <- c("CC","BP","MF")
names(ont_list) <- ont_list

clustering_levels <- c("broad", "semi_broad", "mid")
## testing
# clustering_levels = "broad"

for (test_go in names(lst_go_tests)) {
   
    message("\n==============================")
    message("Running GO enrichment for: ", test_go)
    message("==============================")
    
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
    
    ## Merge Hb subtypes into broad groups (LHb, MHb) ==========================
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
                ggtitle(paste("GO BP enrichment (", test_go, " — ", hb_merged_ct, ")", sep = ""))
            
            plot_base <- sprintf("%s_GO_BP_%s", hb_merged_ct, test_go)
            plot_pdf <- here(plotDir, paste0(plot_base, ".pdf"))
            ggsave(plot_pdf, plot = p, width = 8, height = 6)
            message("Saved plots: ", plot_pdf)
            
        } else {
            message("No significant BP terms for ", test_go)
        }
        
        ## Save results
        go_rds_name <- here(processedDir, sprintf("%s_GO_results_%s.rds", hb_merged_ct, test_go))
        saveRDS(go_result, go_rds_name)
        message("Saved RDS: ", go_rds_name)
        
        ## Save summary table (flattened BP results if available)
        if (!is.null(go_result$BP) && nrow(go_result$BP@compareClusterResult) > 0) {
            bp_df <- go_result$BP@compareClusterResult
            go_csv_name <- here(processedDir, sprintf("%s_GO_results_BP_%s.csv", hb_merged_ct, test_go))
            write.csv(bp_df, go_csv_name, row.names = FALSE)
            message("Saved CSV: ", go_csv_name)
        }
    
        message("Completed GO enrichment for ", hb_merged_ct, " [test: ", test_go, "]")

    }
    
}



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




