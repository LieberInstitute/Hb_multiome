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

universe_overlaps_df <- read.csv(here(inputCSV_merged_overlaps, "Overlaps_LinkPeak_DARs_FDR0.2.csv"))
colnames(universe_overlaps_df)
# > colnames(universe_overlaps_df)
# [1] "peak_id_links"    "CCscore"          "gene_name"        "gene_id"         
# [5] "FDR_CC"           "cluster"          "tss"              "gene_strand"     
# [9] "distance"         "distance_kb"      "signed_distance"  "signed_by_strand"
# [13] "peak_id"          "cell_type"        "FDR_threshold"    "logFC"           
# [17] "fdr_dars"  
head(universe_overlaps_df)
# peak_id_links     CCscore  gene_name         gene_id     FDR_CC
# 1 chr1-1745993-1746550 -0.40853659 AL645728.1 ENSG00000279244 0.02898562
# 2 chr1-1745993-1746550 -0.40853659 AL645728.1 ENSG00000279244 0.02898562
# 3 chr1-3910780-3911122  0.03048837       DFFB ENSG00000169598 0.16370976
# 4 chr1-3910780-3911122  0.03048837       DFFB ENSG00000169598 0.16370976
# 5 chr1-3910780-3911122  0.03048837       DFFB ENSG00000169598 0.16370976
# 6 chr1-3910780-3911122  0.03048837       DFFB ENSG00000169598 0.16370976
# cluster     tss gene_strand distance distance_kb signed_distance
# 1 Astrocyte 1579756           +    83258      83.258           83258
# 2 Astrocyte 1579756           +    83258      83.258           83258
# 3 Astrocyte 3857267           +    26842      26.842           26842
# 4 Astrocyte 3857267           +    26842      26.842           26842
# 5 Astrocyte 3857267           +    26842      26.842           26842
# 6 Astrocyte 3857267           +    26842      26.842           26842
# signed_by_strand              peak_id  cell_type FDR_threshold      logFC
# 1            83258 chr1-1745993-1746550  Astrocyte        FDR0.2  0.4569419
# 2            83258 chr1-1745993-1746550  Microglia        FDR0.2  1.1139565
# 3            26842 chr1-3910780-3911122  Astrocyte        FDR0.2 -0.9013231
# 4            26842 chr1-3910780-3911122 Inhib.Thal        FDR0.2 -0.3233352
# 5            26842 chr1-3910780-3911122      LHb.4        FDR0.2  0.2851914
# 6            26842 chr1-3910780-3911122      Oligo        FDR0.2 -1.4320263
nrow(universe_overlaps_df)
# 23691

universe_df1 <- universe_overlaps_df |>
    filter(!is.na(gene_id)) |>
    distinct()
nrow(universe_df1)
# 23666


## ======

## Testing universe.  Map these to ENTREZ IDs
entrez_map <- bitr(universe_df1$gene_id,
                        fromType = "ENSEMBL",
                        toType = "ENTREZID",
                        OrgDb = org.Hs.eg.db) 
# |> pull(ENTREZID) |> unique()

entrez_universe <- unique(entrez_map$ENTREZID)
message("Universe size: ", length(entrez_universe))
# 5817
# Warning message:
#     In bitr(universe_df2$gene_id, fromType = "ENSEMBL", toType = "ENTREZID",  :
#                 1.33% of input gene IDs are fail to map...





## Build a background universe before subsetting for Direct_Regulation

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
    ## need to add non-signifcance to the entrez universe 
)

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

# tests
cCRE_df

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


## testing
colnames(cCRE_df)
table(cCRE_df$category)

## subset the specific go_universe
cCRE_universe_df <- cCRE_df |> 
    filter(category %in% lst_go_tests[["Direct_Regulation"]]) 
cCRE_universe_df |> nrow()
# 1363

# ## subset the specific go_universe
# cCRE_universe_df <- cCRE_df |> 
#     filter(category %in% lst_go_tests[["Primary_Secondary_Interest"]]) 
# cCRE_universe_df |> nrow()
# # 1783


## =============================================================================
## Merge Hb sub-types into broad groups (LHb, MHb)

cCRE_universe_df <- cCRE_universe_df |>
    mutate(
        cell_type_broad = case_when(
            grepl("^LHb", cell_type) ~ "LHb",
            grepl("^MHb", cell_type) ~ "MHb",
            TRUE ~ cell_type
        )
    )

table(cCRE_universe_df$cell_type)
## Check how many entries per broad group
table(cCRE_universe_df$cell_type_broad)


## =============================================================================
## Build DE_entrez dataframe for GO enrichment
## =============================================================================

# ## ENTREZID look up
# names(cCRE_universe_df)
# 
# ## gene symbol and gene entrez name
# cCRE_universe_df[c("gene_name", "gene_id")]
# 
# # Convert ENSEMBL -> ENTREZID
# entrez_search <- bitr(cCRE_universe_df$gene_id, 
#                       fromType = "ENSEMBL", 
#                       toType = "ENTREZID", 
#                       OrgDb = "org.Hs.eg.db")
# # Warning message:
# #     In bitr(cCRE_universe_df$gene_id, fromType = "ENSEMBL", toType = "ENTREZID",  :
# #                 1.43% of input gene IDs are fail to map.



##==============================================================================
## Merge back and define DE_class + DE_class_cluster

########## Test "Direct_Regulation" ##########

DE_entrez <- cCRE_universe_df |>
    # use this if duplicates are expected (because multiple peaks link to the same gene)
    left_join(entrez_map,
              by = c("gene_id" = "ENSEMBL"),
              relationship = "many-to-many") |>
    filter(!is.na(ENTREZID)) |>
    mutate(
        ## keep the original category direction
        category_direction = case_when(
            category == "Linked_DAR (+) enriched" ~ "up",
            category == "Linked_DAR (-) depleted" ~ "down"#,
            #TRUE ~ "neutral"
        ),
        cell_type_clean = gsub("\\.", "-", cell_type)# ,
        # DE_class_cluster = paste0(cell_type_clean, "_", DE_class)
    ) %>%
    distinct(ENTREZID, .keep_all = TRUE)
           

message("Final DE_entrez dimensions: ", nrow(DE_entrez), " rows, ", length(unique(DE_entrez$DE_group)), " cell groups")
message("Background universe size: ", length(entrez_universe))


## group by cell_type for enrichment
DE_entrez <- DE_entrez |>
    # mutate(DE_group = cell_type_clean) |>
    mutate(DE_group = cell_type_broad) |>
    distinct(DE_group, ENTREZID, .keep_all = TRUE)
head(DE_entrez)
table(DE_entrez$cell_type_broad, DE_entrez$DE_group)
# down  up
# Astrocyte    68  18
# Endo          0   1
# Excit.Thal   66 101
# Inhib.Thal  407 183
# LHb          10 169
# MHb          23 135
# Microglia     0   5
# Oligo        41  15
# OPC           0   6


## keep groups with at least 10 genes
DE_entrez <- DE_entrez |>
    group_by(DE_group) |>
    filter(n() >= 10) |>
    ungroup()

table(DE_entrez$DE_group)
unique(DE_entrez$cell_type_broad)

## checks
# Check how many unique genes per DE_group
DE_entrez |> count(DE_group)
# DE_group       n
# <chr>      <int>
# 1 Astrocyte     86
# 2 Excit.Thal   167
# 3 Inhib.Thal   590
# 4 LHb          179
# 5 MHb          158
# 6 Oligo         56


# test_ids <- DE_entrez |> filter(DE_group == "LHb") |> pull(ENTREZID)
# ego_test <- enrichGO(gene = test_ids, OrgDb = org.Hs.eg.db,
#                      universe   = entrez_universe,   # correct background
#                      ont = "BP", 
#                      pvalueCutoff = 1, 
#                      qvalueCutoff = 0.2,
#                      readable = TRUE)
# head(ego_test)



#### Run GO ####

ont_list <- c("CC","BP","MF")
names(ont_list) <- ont_list

## GO terms summarized by broad cell types (LHb, MHb, Astrocyte, etc.) rather than by fine subclusters.

go_result <- map(ont_list, ~compareCluster(
               ENTREZID ~ DE_group,
               data = DE_entrez,
               OrgDb = org.Hs.eg.db,
               fun = enrichGO,
               universe = entrez_universe,
               ont = .x,
               pAdjustMethod = "BH",
               pvalueCutoff = 0.1,
               qvalueCutoff = 0.2,
               readable = TRUE
               ))

# go_result for diagnostics: how many GO terms per ontology
map2(names(go_result), go_result, function(nm, x) {
    n_terms <- if (is.null(x)) 0 else nrow(x@compareClusterResult)
    message(sprintf("%s: %s terms", nm, n_terms))
})



# Remove NULL entries (e.g., CC = NULL)
go_result_valid <- discard(go_result, is.null)

## convert to table & extract compareClusterResult from each valid ontology and tag with ontology name
compare_clus <- map2_dfr(go_result_valid, names(go_result_valid), function(x, nm) {
    x@compareClusterResult |> mutate(ONTOLOGY = nm)
})

# verify classes that survive
map(go_result, function(x) {
    if (is.null(x)) "NULL" else nrow(x@compareClusterResult)
})

compare_clus |> count(DE_class_cluster, ONTOLOGY)
# DE_class_cluster ONTOLOGY n
# 1   LHb-1-3-4_down       MF 8
# 2       LHb-2-7_up       BP 2
# 3       MHb-2_down       MF 2

## Save 
f_name <- here(processedDir, sprintf("GO_compare_clus_%s.rds", peaks_classification_name))
saveRDS(compare_clus, f_name)
f_name <- here(processedDir, sprintf("GO_results_%s.csv", peaks_classification_name))
write.csv(compare_clus, f_name, row.names = FALSE)


#### dot plots ####

# Iterate over valid GO results (skip NULL or empty ones)
go_result_valid <- discard(go_result, function(x) {
    is.null(x) || nrow(x@compareClusterResult) == 0
})

length(go_result_valid)
map(go_result_valid, ~nrow(.@compareClusterResult))
# pdf(f_name, width = 10, height = 10)
# print(dotplot(go_result_valid[[2]])) # pass!
# dev.off()

f_name = here(processedDir, sprintf("GO_dotplot_%s.pdf", peaks_classification_name))
pdf(f_name, width = 10, height = 10)

## one page per ontology 
if (length(go_result_valid) == 0) {
    message("No valid GO results to plot.")
    dev.off()
} else {
    walk2(go_result_valid, names(go_result_valid), function(x, nm) {
        # Check if the result for this ontology contains enough data to attempt plotting
        if(nrow(x@compareClusterResult) > 0) {
            
            message("Plotting ontology: ", nm, " (", nrow(x@compareClusterResult), " terms)")
            
            # Use tryCatch to prevent a single failing plot from crashing the entire PDF
            tryCatch({
                
                p <- dotplot(x,
                             x = "DE_class_cluster", 
                             showCategory = 5, 
                             label_format = 60) +
                    ggtitle(paste("GO Enrichment:", nm)) +
                    theme_bw(base_size = 12) +
                    theme(
                        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8),
                        axis.text.y = element_text(size = 8),
                        plot.title = element_text(hjust = 0.5, face = "bold")
                    )
                print(p) 
            }, error = function(e) {
                message("Skipping plot for ", nm, " due to error: ", conditionMessage(e))
            })
        } else {
            message("Skipping plot for ", nm, ": No terms remaining after filtering.")
        }
    })
    
    dev.off()
    
}

# library("slurmjobs")
# job_single(
#   "22_GO_enrichment",
#   create_shell = TRUE,
#   partition = "katun",
#   memory = "30G",
#   cores = 2,
#   logdir = "logs",
#   command = "Rscript 22_GO_enrichment.R",
#   create_logdir = FALSE
# )

# 
# pdf(file = here(plot_dir, sprintf("GO_dotplot_%s_2plus.pdf", opt$datatype)), width = 10, height = 10)
# walk2(go_result, names(go_result), function(gr, ont){
#     
#     gr@compareClusterResult <- gr@compareClusterResult |> filter(!grepl("0", DE_class_cluster), Count >= 2)
#     
#     if(nrow(gr@compareClusterResult) == 0) return(NULL)
#     
#     print(
#         dotplot(gr, 
#                 x = "DE_class_cluster", 
#                 showCategory = 3, 
#                 label_format = 60)  +
#             ggtitle(paste("GO Enrichment:", ont, " (2+ genes)")) +
#             theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
#     )
# })
# dev.off()
# 
# pdf(file = here(plot_dir, sprintf("GO_dotplot_%s_5plus.pdf", opt$datatype)), width = 10, height = 10)
# walk2(go_result, names(go_result), function(gr, ont){
#     
#     gr@compareClusterResult <- gr@compareClusterResult |> filter(!grepl("0", DE_class_cluster), Count >= 5)
#     
#     if(nrow(gr@compareClusterResult) == 0) return(NULL)
#     
#     print(
#         dotplot(gr, 
#                 x = "DE_class_cluster", 
#                 showCategory = 3, 
#                 label_format = 60)  +
#             ggtitle(paste("GO Enrichment:", ont)) +
#             theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
#     )
# })
# dev.off()
# 
# # for sn fine plot by cell type
# if(opt$datatype == "sn_fine"){
#     
#     ## by broad cell types
#     broad_cell_types <- unique(jaffelab::ss(cluster_levels, "\\."))
#     
#     pdf(file = here(plot_dir, sprintf("GO_dotplot_%s_cell_type.pdf", opt$datatype)), width = 10, height = 10)
#     
#     map(broad_cell_types, function(ct){
#         go_result_ct <- map2(go_result, names(go_result), function(gr, ont){
#             # subset
#             gr@compareClusterResult <- gr@compareClusterResult |> filter(grepl(ct, Cluster))
#             
#             if(nrow(gr@compareClusterResult ) > 1){
#                 # dotplot
#                 print(
#                     dotplot(gr,
#                             x = "DE_class_cluster",
#                             showCategory = 5,
#                             label_format = 60)  +
#                         ggtitle(sprintf("GO Enrichment:%s - %s", ont, ct)) +
#                         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
#                 )}
#             
#             # return(gr@compareClusterResult |> nrow())
#         })
#         # return(go_result_ct)
#     })
#     
#     dev.off()
#     
#     pdf(file = here(plot_dir, sprintf("GO_dotplot_%s_cell_type_2plus.pdf", opt$datatype)), width = 10, height = 10)
#     map(broad_cell_types, function(ct){
#         go_result_ct <- map2(go_result, names(go_result), function(gr, ont){
#             # subset
#             gr@compareClusterResult <- gr@compareClusterResult |> filter(grepl(ct, Cluster), Count >2)
#             
#             if(nrow(gr@compareClusterResult ) > 1){
#                 # dotplot
#                 print(
#                     dotplot(gr,
#                             x = "DE_class_cluster",
#                             showCategory = 5,
#                             label_format = 60)  +
#                         ggtitle(sprintf("GO Enrichment:%s - %s (2+ genes)", ont, ct)) +
#                         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
#                 )}
#             
#             # return(gr@compareClusterResult |> nrow())
#         })
#         # return(go_result_ct)
#     })
#     dev.off()    
#     
#     pdf(file = here(plot_dir, sprintf("GO_dotplot_%s_cell_type_5plus.pdf", opt$datatype)), width = 10, height = 10)
#     map(broad_cell_types, function(ct){
#         go_result_ct <- map2(go_result, names(go_result), function(gr, ont){
#             # subset
#             gr@compareClusterResult <- gr@compareClusterResult |> filter(grepl(ct, Cluster), Count >5)
#             
#             if(nrow(gr@compareClusterResult ) > 1){
#                 # dotplot
#                 print(
#                     dotplot(gr,
#                             x = "DE_class_cluster",
#                             showCategory = 5,
#                             label_format = 60)  +
#                         ggtitle(sprintf("GO Enrichment:%s - %s (5+ genes)", ont, ct)) +
#                         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
#                 )}
#             
#             # return(gr@compareClusterResult |> nrow())
#         })
#         # return(go_result_ct)
#     })
#     dev.off()
#     
#     
#     
#     
#     go_result_Oligo3 <- map(go_result, function(gr){
#         gr@compareClusterResult <- gr@compareClusterResult |> filter(grepl("Oligo.3", Cluster) | grepl("Astro.3", Cluster))
#         return(gr)
#     })
#     
#     map(go_result_Oligo3, ~.x@compareClusterResult |> count(Cluster))
#     
#     pdf(file = here(plot_dir, sprintf("GO_dotplot_%s_Oligo.3.pdf", opt$datatype)), width = 6, height = 6)
#     walk2(go_result_Oligo3, names(go_result_Oligo3), 
#           ~print(
#               dotplot(.x, 
#                       x = "DE_class_cluster", 
#                       showCategory = 5, 
#                       label_format = 60)  +
#                   ggtitle(paste("GO Enrichment:", .y)) +
#                   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
#           )
#     )
#     dev.off()    
#     
#     ## select genes
#     go_result_select_genes <- map(go_result, function(gr){
#         gr@compareClusterResult <- gr@compareClusterResult  |> filter(grepl("FOS", geneID) |
#                                                                           grepl("TLR2", geneID) | 
#                                                                           grepl("STAT1", geneID) | 
#                                                                           grepl("STAT4", geneID))
#         
#         return(gr)
#     })
#     
#     map(go_result_select_genes, ~.x@compareClusterResult |> count(Cluster))
#     
#     pdf(file = here(plot_dir, sprintf("GO_dotplot_%s_select_genes.pdf", opt$datatype)), width = 8, height = 8)
#     walk2(go_result_select_genes, names(go_result_select_genes), 
#           ~print(
#               dotplot(.x, 
#                       x = "DE_class_cluster", 
#                       showCategory = 10, 
#                       label_format = 60)  +
#                   ggtitle(paste("GO Enrichment:", .y, "FOS|TLR2|STAT1|STAT4")) +
#                   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
#           )
#     )
#     dev.off()
#     
#     
# }




## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()




