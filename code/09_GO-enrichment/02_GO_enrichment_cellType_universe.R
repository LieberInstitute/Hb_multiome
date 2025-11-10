########################################################################
## Performs cell-type-specific GO enrichment by dynamically defining both the foreground gene set and the background universe based on the clustering resolution.
##
## Authors. CSC
## Date. Oct 24, 2025
## Recommended resources on interactive mode: srun --pty --mem=20GB --x11 bash
########################################################################

# GO Enrichment Pipeline Summary (Setup and Data Flow)
# 
# I. Initialization and Data Setup
# Input Control: The script accepts the desired clustering resolution (broad, semi_broad, or mid) as a command-line argument.
# Master Mapping: All unique Ensembl IDs from all raw linkage files are mapped to ENTREZID to create a single master conversion map.
# Foreground Data: Classified cCRE data is loaded, and multiple test categories (e.g., enriched, depleted) are defined as the foreground gene sets.
# II. Dynamic Analysis Flow
# The script runs in nested loops:
#     Outer Loop: Processes each foreground test category defined (e.g., Direct Regulation All, Link Peaks OCRs).
#     Inner Loop: Processes the clustering resolution specified via the command line.
# 
# Dynamic Universe Definition: For each remaining cell group, calculates a unique background universe. This universe consists only of the Entrez IDs linked to raw peaks (pre-filtering) found within the cell types/files relevant to that specific group and resolution.

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
# Define resolution_level
#===============================================================================
# Specification for argument parsing: -r or --res_level <value>
spec = matrix(c(
    'res_level', 'r', 1, "character", "Clustering resolution level (broad, semi_broad, mid)."
), byrow=TRUE, ncol=5)
opt = getopt(spec)

if (is.null(opt$res_level)) { stop(getopt(spec, usage = TRUE)) }

res_level_arg <- opt$res_level
message("Processing job array for resolution level: ", res_level_arg)

# Define the single-element vector for the inner loop
res_level <- c(res_level_arg)

# ==============================================================================

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

message("Runnig go-enrichment at [", res_level, "] resolution")


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

message("Total raw links by cell-type ....")
table(master_link_df$cluster)

# Perform the comprehensive mapping once
master_entrez_map <- clusterProfiler::bitr(
    master_link_genes,
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db,
    drop = TRUE
)


## Declare functions ===========================================================

## Dynamic Universe Definition ***********
load_entrez_cellType_universe <- function(
    clusterRes,  # Could be any of c("broad", "semi_broad", "mid")
    ct,          # for Broad=Hb / for Semi-Broad= MHb or LHb
    lst_peak_paths) 
{
    # clusterRes = "Broad"
    # lst_peak_paths = lst_peak_files
    # ct = "LHb.2"
    
    message("Processing universe for clusterig resolution [", clusterRes, "] and cell type [", ct, "]")
    
    # Extract the target pattern from ALL file names (e.g., "Astrocyte", "MHb.2")
    target_patterns <- sub("^[^_]+_([^_]+)_.*$", "\\1", basename(lst_peak_paths))

    ## Determine which files to merge based on the cell type (ct) 
 
    if (clusterRes == "broad") { 
        if (ct == "Hb") {
            # Combine ALL Habenula files for the 'Hb' group's universe
            is_target_file <- grepl("^(M|L)Hb.*$", target_patterns)
        } else {
            # Match files exactly to the specific cell type for non-habenula groups
            is_target_file <- target_patterns == ct
        }
        
    } else if (clusterRes == "semi_broad") { 
        # Semi-Broad: Merge files that match the broad LHb or MHb group (e.g., 'LHb' or 'MHb')
        # Here, 'ct' is expected to be 'LHb' or 'MHb'
        pattern <- paste0("^", ct, ".*$")
        is_target_file <- grepl(pattern, target_patterns)
        
    } else if (clusterRes == "mid") {
        # Mid: Use ONLY the file corresponding to the specific cell type (ct)
        is_target_file <- target_patterns == ct
        
    } else {
        stop("Invalid clusterRes value provided.")
        
    }
        
    ## load, combine and extract unique genes
    filtered_files <- lst_peak_paths[is_target_file]
    
    link_df <- filtered_files |>
        map(~ read.csv(.x) |> select(cluster, gene_id)) |>
        list_rbind()
    
    ## Filter by cluster if necessary (ensuring only linked genes relevant to 'ct' are used)
    link_genes <- link_df |> 
        distinct(gene_id) |> 
        pull(gene_id)

    # === Entrez ID Mapping Section ===
    entrez_map <- clusterProfiler::bitr(
        link_genes,
        fromType = "ENSEMBL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db,
        drop = TRUE # Only keep mapped IDs
    )
    entrez_universe <- unique(entrez_map$ENTREZID)

    message("Mapped to ", length(entrez_universe), " Entrez IDs (", ct, " specific universe)")
    
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
#nrow(cCRE_df) # 6410
#head(cCRE_df)
#table(cCRE_df$type_classification)

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
            left_join(entrez_map, by = c("gene_id" = "ENSEMBL")) |> #,
                      #relationship = "many-to-many") |> it is ignored in newer versions
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
  
         
check_go_results_validity <- function(go_result_list) {
    # Check if any ontology result is an S4 object and has rows in compareClusterResult
    has_valid_data <- purrr::map_lgl(go_result_list, function(res_obj) {
        if (is.null(res_obj)) {
            return(FALSE)
        }
        # Check if it has the required slot and if that slot has more than 0 rows
        return(methods::is(res_obj, "compareClusterResult") && 
                   methods::slot(res_obj, "compareClusterResult") |> nrow() > 0)
    })
    
    # Return TRUE if at least one ontology has valid data
    return(any(has_valid_data))
}

##################################################################
# Loop (GO Test Category)
##################################################################

ont_list <- c("CC","BP","MF")
names(ont_list) <- ont_list

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
    
    message("\n--- Running Clustering Resolution: ", res_level, " ---")
    
    if (res_level == "mid") message("Using full mid-resolution clustering - no merging applied.")
    
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
    
    ## verification
    ## table(DE_entrez_grouped$DE_group)
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
    
    ## starts enrichment test ====>
    
    # **Dynamically load the GO Universe for the current resolution level and cell-type**
    # Run individual enrichGO tests per DE_group with its unique universe
    
    go_result <- map(ont_list, function(ont_type) {
        
        # Generate a named list of cell groups to ensure proper naming of results
        named_groups <- set_names(cell_groups_to_test)
        
        # Run enrichGO for each group and ontology type
        list_of_enrichments <- map(named_groups, function(group_name) {
            # group_name = "Excit.Thal"
            # Critical Step!!! - Load the universe specific to the current 'group_name' and 'res_level'
            universe_map <- load_entrez_cellType_universe(
                clusterRes = res_level, 
                ct = group_name, # specific cell group name. ge. LHb.4
                lst_peak_paths = lst_peak_files 
            )
            entrez_universe <- unique(universe_map$ENTREZID) 
            
            # Subset genes for the current group
            genes_to_test <- DE_entrez_grouped |>
                filter(DE_group == group_name) |>
                pull(ENTREZID)
            
            # Check for empty universe or insufficient test genes
            if(length(entrez_universe) == 0 || length(genes_to_test) < 5) {
                message(sprintf("  -> Skipping %s (%s genes, Universe size: %s)...", 
                                group_name, length(genes_to_test), length(entrez_universe)))
                return(NULL) 
            }
            
            # Ensure the genes to test are present in the universe (required by enrichGO)
            genes_to_test <- intersect(genes_to_test, entrez_universe)
            
            if(length(genes_to_test) < 5) return(NULL) # Skip if too few genes after filtering
            
            message(sprintf("Testing %s (%s genes/Universe %s) for %s...", 
                            group_name, length(genes_to_test), 
                            length(entrez_universe), ont_type))
            
            suppressWarnings({
                enrich_res <- enrichGO(
                    gene = genes_to_test,
                    OrgDb = org.Hs.eg.db,
                    keyType = "ENTREZID",
                    universe = entrez_universe, # Use the now-specific universe
                    ont = ont_type,
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.1,
                    qvalueCutoff = 0.2,
                    readable = TRUE
                )
            })
            
            # Tag the result with the group name (return full S4 to before combining)
            if (!is.null(enrich_res) && nrow(enrich_res) > 0) {
                enrich_res@result$Cluster <- group_name
                # Only return the data frame result, merge_result can combine data frames
                return(enrich_res)
            } else {
                return(NULL)
            }
        }) |> compact() # Remove NULL results (groups with no enrichment)
        
        # Combine all enrichments for this ontology into a single compareClusterResult
        if (length(list_of_enrichments) > 0) {
            final_result <- clusterProfiler::merge_result(list_of_enrichments)
            # Rename the 'Cluster' column to 'DE_group' to match expectations for plotting
            colnames(final_result@compareClusterResult)[
                colnames(final_result@compareClusterResult) == "Cluster"
            ] <- "DE_group"
            
            return(final_result)
        } else {
            return(NULL)
        }
        
    })
    ## ends enrichment test ====/
    
    names(go_result) <- ont_list
    
    ## Save and plot results
    is_valid_go_result <- check_go_results_validity(go_result)
    
    if (is_valid_go_result) {
        message("\nValid GO enrichment results found. Proceeding with saving and plotting.")
        
        ## Save Full Results (All Ontologies: RDS)
        go_rds_name <- here::here(processedDir, 
                                  sprintf("%s_GO_results_ALL_ONTOLOGIES_%s_%s.rds", res_level, test_go, timestamp))
        saveRDS(go_result, go_rds_name)
        message("Saved full GO results object (ALL ontologies): ", go_rds_name)
        
        ## Plot Biological Process (BP) Results (PDF)
        # Check specifically for BP before plotting
        bp_result <- go_result$BP # Define the BP result S4 object
        
        if (!is.null(bp_result) && nrow(bp_result@compareClusterResult) > 0) {
            message("Generating and saving GO BP dot plot...")
            
            bp_df <- bp_result@compareClusterResult # access df
            bp_df$GeneRatio <- as.character(bp_df$GeneRatio)
            bp_df_clean <- bp_df[!is.na(bp_df$GeneRatio) & bp_df$GeneRatio != "", ]
            
            if (nrow(bp_df_clean) == 0) {
                message("Warning: All GO BP results were removed during cleanup; cannot plot.")
                next # Skip plotting and move to next resolution/test if inside a loop
            } 
            p <- ggplot(bp_df_clean, 
                        aes(x = DE_group, 
                            y = reorder(Description, Count), 
                            size = Count, 
                            color = p.adjust)) +
                geom_point() +
                scale_color_gradient(low = "red", high = "blue", name = "Adj. p-value") +
                labs(
                    title = paste(res_level, "GO BP enrichment"),
                    subtitle = paste("Category:", test_go),
                    x = "Cell Type",
                    y = "GO Biological Process"
                ) +
                theme_minimal(base_size = 11) +
                theme(plot.title = element_text(face = "bold"))
                
            plot_base <- sprintf("%s_GO_BP_%s_%s", res_level, test_go, timestamp)
            plot_pdf <- here(plotDir, paste0(plot_base, ".pdf"))
            ggsave(plot_pdf, plot = p, width = 8, height = 6)
            message("Saved plot: ", plot_pdf)
        }

    } else {
        message("Skipping saving and plotting: No significant GO terms found across any ontology.")
    }
    
}


message("GO-Enrichment completed!")



## Reproducibility information
library(sessioninfo)
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()




