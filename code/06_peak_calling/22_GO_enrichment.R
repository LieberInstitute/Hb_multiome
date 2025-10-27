########################################################################
## GO enrichment in cCRE and classified Linked OCR and DARs unlinked
##
## Authors. CSC
## Date. Oct 24, 2025
## Recommended resources on interactive mode: srun --pty --mem=20GB --x11 bash
########################################################################

#library("pheatmap")
#library("reshape2")
library("dplyr")
library("purrr")
library("ggplot2")
library("patchwork")
#library("ggrepel")
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
inputCSV_cCRE_ORC_Dir <- here(
    "processed-data",
    "06_peak_calling",
    "21_overlaping_FDRscores_TopHeatmap"
)
processedDir <- here(
    "processed-data",
    "06_peak_calling",
    "22_GO_enrichment"
)
plotDir <- here(
    "plots",
    "06_peak_calling",
    "22_GO_enrichment"
)

if (!dir.exists(processedDir)) { dir.create(processedDir) }
if (!dir.exists(plotDir)) { dir.create(plotDir) }


##==============================================================================

message("Loading cCRE and ORC file ...")

f_name = "overlaps_linkPeak_DARs_classified_thr_CC0.3_thr_DAR0.1.csv"

## load raw overlaps and 2-shared ct overlaps 
cCRE_df <- read.csv(here(inputCSV_cCRE_ORC_Dir, f_name))
nrow(cCRE_df) # 6410
head(cCRE_df)

## validation
unique(cCRE_df$type_classification)
# [1] "classification-2" "classification-3"
# as.data.frame(table(cCRE_df$type_classification, cCRE_df$category))

peaks_classification3 = c("Linked_DAR (+) enriched", "Linked_DAR (-) depleted", "Linked OCR", "Unlinked DAR", "Non-significant")
peaks_classification3_universe = c("Linked_DAR (+) enriched", "Linked_DAR (-) depleted")

cCRE_universe_df <- cCRE_df |> 
    filter(type_classification == "classification-3"
           & category %in% peaks_classification3_universe)

table(cCRE_universe_df$type_classification)
table(cCRE_universe_df$category)
head(cCRE_universe_df)

cluster_levels <- cCRE_universe_df$cell_type |> unique()

cCRE_universe_df |> count(cluster)
# cluster   n
# 1   Astrocyte  19
# 2        Endo   2
# 3  Excit.Thal 213
# 4  Inhib.Thal 278
# 5       LHb.1  29
# 6     LHb.1.3   6
# 7   LHb.1.3.4  42
# 8     LHb.2.7 161
# 9       LHb.4 236
# 10      MHb.1  30
# 11    MHb.1.2   4
# 12      MHb.2 323
# 13      MHb.3   1
# 14  Microglia   1
# 15        OPC   2
# 16      Oligo  11
# 17       Thal   5

summary(cCRE_universe_df$fdr_dars)


## =============================================================================
## Build DE_entrez dataframe for GO enrichment
## =============================================================================

## ENTREZID look up
names(cCRE_universe_df)
## gene symbol and gene entrez name
cCRE_universe_df[c("gene_name", "gene_id")]

## 
# Convert ENSEMBL -> ENTREZID
entrez_search <- bitr(cCRE_universe_df$gene_id, 
                      fromType = "ENSEMBL", 
                      toType = "ENTREZID", 
                      OrgDb = "org.Hs.eg.db")
# Warning message:
#     In bitr(cCRE_universe_df$gene_id, fromType = "ENSEMBL", toType = "ENTREZID",  :
#                 1.43% of input gene IDs are fail to map.
head(entrez_search)
message("Mapped ", nrow(entrez_search), " gene IDs to ENTREZID")
# Mapped 1248 gene IDs to ENTREZID

# DE_entrez <- DE_data |> 
#     left_join(entrez_search, by = c("gene_id" = "ENSEMBL"), relationship = "many-to-many") |>
#     filter(!is.na(ENTREZID)) |>
#     mutate(DE_class = case_when(vlmf_logFC > 0 & vlmf_adj.P.Val < 0.05 ~ "up",
#                                 vlmf_logFC < 0 & vlmf_adj.P.Val < 0.05 ~ "down",
#                                 TRUE ~ "None"),
#            DE_class_cluster = paste0(gsub("\\.", "-", cluster), "_",DE_class)) ## doesn't like .  in cluster names

## Merge back and define DE_class + DE_class_cluster

DE_entrez <- cCRE_universe_df |>
    # use this if duplicates are expected (because multiple peaks link to the same gene)
    left_join(entrez_search,
              by = c("gene_id" = "ENSEMBL"),
              relationship = "many-to-many") |>
    filter(!is.na(ENTREZID)) |>
    mutate(
        ## Use logFC sign only — since FDR filtering was done upstream
        DE_class = case_when(
            logFC > 0 ~ "up",
            logFC < 0 ~ "down",
            TRUE ~ "None"
        ),
        ## Optionally keep the original category direction if meaningful
        category_direction = case_when(
            category == "Linked_DAR (+) enriched" ~ "up",
            category == "Linked_DAR (-) depleted" ~ "down",
            TRUE ~ "neutral"
        ),
        cell_type_clean = gsub("\\.", "-", cell_type),
        DE_class_cluster = paste0(cell_type_clean, "_", DE_class)
    ) %>%
    distinct(ENTREZID, .keep_all = TRUE) %>%
    select(cell_type, cell_type_clean, gene_name, gene_id, ENTREZID,
           logFC, fdr_dars, category, DE_class, category_direction,
           DE_class_cluster)

## keep all mappings but avoid inflating the enrichment gene universe
DE_entrez <- DE_entrez |>
    distinct(cell_type, ENTREZID, .keep_all = TRUE)


DE_entrez |> count(cluster)
# DE_entrez |> filter(DE_class != "None") |> count(DE_class_cluster) |> head()


#### Run GO ####
universe <- unique(DE_entrez$ENTREZID)
length(universe)

ont_list <- c("CC","BP","MF")
names(ont_list) <- ont_list

go_result <- map(ont_list, ~compareCluster(ENTREZID ~ DE_class_cluster,
                                           data = DE_entrez |> filter(DE_class != "None"), 
                                           OrgDb = org.Hs.eg.db,
                                           fun = enrichGO,
                                           universe = universe,
                                           ont = .x, ##ALL,CC,BP,MF
                                           pAdjustMethod = "BH",
                                           pvalueCutoff = 0.05,
                                           # qvalueCutoff = 0.05,
                                           readable = TRUE))

## control point 
saveRDS(go_result, file = here(processedDir, sprintf("GO_result_%s.rds", "cCRE")))



go_result <- readRDS(here(data_dir, sprintf("GO_result_%s.rds", opt$datatype)))

## convert to table
compare_clus <- map2_dfr(go_result, names(go_result), ~.x@compareClusterResult |> mutate(ONTOLOGY = .y))
compare_clus |> count(DE_class_cluster, ONTOLOGY)

# Save 
saveRDS(compare_clus, file = here(data_dir, sprintf("GO_compare_clus_%s.rds", opt$datatype)))
write.csv(compare_clus, file = here(data_dir, sprintf("GO_results_%s.csv", opt$datatype)), row.names = FALSE)


#### dot plots ####

pdf(file = here(plot_dir, sprintf("GO_dotplot_%s.pdf", opt$datatype)), width = 10, height = 10)
walk2(go_result, names(go_result), 
      ~print(
          dotplot(.x, 
                  x = "DE_class_cluster", 
                  showCategory = 3, 
                  label_format = 60)  +
              ggtitle(paste("GO Enrichment:", .y, " (5+ genes)")) +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
      )
)

dev.off()

pdf(file = here(plot_dir, sprintf("GO_dotplot_%s_2plus.pdf", opt$datatype)), width = 10, height = 10)
walk2(go_result, names(go_result), function(gr, ont){
    
    gr@compareClusterResult <- gr@compareClusterResult |> filter(!grepl("0", DE_class_cluster), Count >= 2)
    
    if(nrow(gr@compareClusterResult) == 0) return(NULL)
    
    print(
        dotplot(gr, 
                x = "DE_class_cluster", 
                showCategory = 3, 
                label_format = 60)  +
            ggtitle(paste("GO Enrichment:", ont, " (2+ genes)")) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
    )
})
dev.off()

pdf(file = here(plot_dir, sprintf("GO_dotplot_%s_5plus.pdf", opt$datatype)), width = 10, height = 10)
walk2(go_result, names(go_result), function(gr, ont){
    
    gr@compareClusterResult <- gr@compareClusterResult |> filter(!grepl("0", DE_class_cluster), Count >= 5)
    
    if(nrow(gr@compareClusterResult) == 0) return(NULL)
    
    print(
        dotplot(gr, 
                x = "DE_class_cluster", 
                showCategory = 3, 
                label_format = 60)  +
            ggtitle(paste("GO Enrichment:", ont)) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
    )
})
dev.off()

# for sn fine plot by cell type
if(opt$datatype == "sn_fine"){
    
    ## by broad cell types
    broad_cell_types <- unique(jaffelab::ss(cluster_levels, "\\."))
    
    pdf(file = here(plot_dir, sprintf("GO_dotplot_%s_cell_type.pdf", opt$datatype)), width = 10, height = 10)
    
    map(broad_cell_types, function(ct){
        go_result_ct <- map2(go_result, names(go_result), function(gr, ont){
            # subset
            gr@compareClusterResult <- gr@compareClusterResult |> filter(grepl(ct, Cluster))
            
            if(nrow(gr@compareClusterResult ) > 1){
                # dotplot
                print(
                    dotplot(gr,
                            x = "DE_class_cluster",
                            showCategory = 5,
                            label_format = 60)  +
                        ggtitle(sprintf("GO Enrichment:%s - %s", ont, ct)) +
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
                )}
            
            # return(gr@compareClusterResult |> nrow())
        })
        # return(go_result_ct)
    })
    
    dev.off()
    
    pdf(file = here(plot_dir, sprintf("GO_dotplot_%s_cell_type_2plus.pdf", opt$datatype)), width = 10, height = 10)
    map(broad_cell_types, function(ct){
        go_result_ct <- map2(go_result, names(go_result), function(gr, ont){
            # subset
            gr@compareClusterResult <- gr@compareClusterResult |> filter(grepl(ct, Cluster), Count >2)
            
            if(nrow(gr@compareClusterResult ) > 1){
                # dotplot
                print(
                    dotplot(gr,
                            x = "DE_class_cluster",
                            showCategory = 5,
                            label_format = 60)  +
                        ggtitle(sprintf("GO Enrichment:%s - %s (2+ genes)", ont, ct)) +
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
                )}
            
            # return(gr@compareClusterResult |> nrow())
        })
        # return(go_result_ct)
    })
    dev.off()    
    
    pdf(file = here(plot_dir, sprintf("GO_dotplot_%s_cell_type_5plus.pdf", opt$datatype)), width = 10, height = 10)
    map(broad_cell_types, function(ct){
        go_result_ct <- map2(go_result, names(go_result), function(gr, ont){
            # subset
            gr@compareClusterResult <- gr@compareClusterResult |> filter(grepl(ct, Cluster), Count >5)
            
            if(nrow(gr@compareClusterResult ) > 1){
                # dotplot
                print(
                    dotplot(gr,
                            x = "DE_class_cluster",
                            showCategory = 5,
                            label_format = 60)  +
                        ggtitle(sprintf("GO Enrichment:%s - %s (5+ genes)", ont, ct)) +
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
                )}
            
            # return(gr@compareClusterResult |> nrow())
        })
        # return(go_result_ct)
    })
    dev.off()
    
    
    
    
    go_result_Oligo3 <- map(go_result, function(gr){
        gr@compareClusterResult <- gr@compareClusterResult |> filter(grepl("Oligo.3", Cluster) | grepl("Astro.3", Cluster))
        return(gr)
    })
    
    map(go_result_Oligo3, ~.x@compareClusterResult |> count(Cluster))
    
    pdf(file = here(plot_dir, sprintf("GO_dotplot_%s_Oligo.3.pdf", opt$datatype)), width = 6, height = 6)
    walk2(go_result_Oligo3, names(go_result_Oligo3), 
          ~print(
              dotplot(.x, 
                      x = "DE_class_cluster", 
                      showCategory = 5, 
                      label_format = 60)  +
                  ggtitle(paste("GO Enrichment:", .y)) +
                  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
          )
    )
    dev.off()    
    
    ## select genes
    go_result_select_genes <- map(go_result, function(gr){
        gr@compareClusterResult <- gr@compareClusterResult  |> filter(grepl("FOS", geneID) |
                                                                          grepl("TLR2", geneID) | 
                                                                          grepl("STAT1", geneID) | 
                                                                          grepl("STAT4", geneID))
        
        return(gr)
    })
    
    map(go_result_select_genes, ~.x@compareClusterResult |> count(Cluster))
    
    pdf(file = here(plot_dir, sprintf("GO_dotplot_%s_select_genes.pdf", opt$datatype)), width = 8, height = 8)
    walk2(go_result_select_genes, names(go_result_select_genes), 
          ~print(
              dotplot(.x, 
                      x = "DE_class_cluster", 
                      showCategory = 10, 
                      label_format = 60)  +
                  ggtitle(paste("GO Enrichment:", .y, "FOS|TLR2|STAT1|STAT4")) +
                  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
          )
    )
    dev.off()
    
    
}


#### rrvgo ####
# 
# reducedTerms_list <- map(ont_list, function(o){
#     ## loop ontologies
#     go_analysis <- go_result[[o]]@compareClusterResult
#     go_analysis$Cluster <- droplevels(go_analysis$Cluster)
#     ## get clusters
#     clusters <- levels(go_analysis$Cluster)
#     names(clusters) <- clusters
#     
#     message(Sys.time(), " - ", o)
#     
#     ont_reducedTerms <- map(clusters, function(clus){
#         ## calc similarity scores
#         go_analysis_c <- go_analysis |>
#             filter(Cluster == clus,
#                    !is.na(qvalue)) ## why are some NA?
#         
#         if(nrow(go_analysis_c) == 0) return(NULL)
#         
#         message(sprintf("--%s (%i)--", clus, nrow(go_analysis_c)))
#         
#         
#         simMatrix <- calculateSimMatrix(go_analysis_c$ID,
#                                         orgdb="org.Hs.eg.db",
#                                         ont=o,
#                                         semdata = GOSemSim::godata(annoDb = "org.Hs.eg.db", ont = o),
#                                         method="Rel")
#         
#         ## group terms based on similarity
#         scores <- setNames(-log10(go_analysis_c$qvalue), go_analysis_c$ID)
#         reducedTerms <- NA
#         
#         try(reducedTerms <- reduceSimMatrix(simMatrix,
#                                             scores,
#                                             threshold=0.7,
#                                             orgdb="org.Hs.eg.db"))
#         
#         return(reducedTerms)
#     }
#     )
#     
#     return(ont_reducedTerms)
# })
# 
# 
# # map_depth(reducedTerms_list, 2, length)
# 
# ## save
# saveRDS(reducedTerms_list, file = here(data_dir, sprintf("GO_reduced_terms_%s.rds", opt$datatype)))
# 
# reducedTerms_list2 <- list_transpose(reducedTerms_list)
# reducedTerms_list2 <- reducedTerms_list2[order(names(reducedTerms_list2))]
# 
# # map_depth(reducedTerms_list2, 2, length)
# # map_depth(reducedTerms_list2, 2, is.null)
# 
# # ## plot treemap
# # pdf(here(plot_dir, "treemap_test.pdf"))
# # treemapPlot(reducedTerms, title = "Oligo")
# # dev.off()
# 
# ## rm empty results
# reducedTerms_list2 <- reducedTerms_list2[map_lgl(reducedTerms_list2, function(rt) !all(map_lgl(rt, ~all(is.null(.x)))))]
# 
# # reducedTerms_list2[["Inhib-Vip_down"]]["CC"]
# 
# pdf(here(plot_dir, sprintf("GO_treemap_%s.pdf", opt$datatype)))
# walk2(reducedTerms_list2, names(reducedTerms_list2), function(rt, clus_name){
#     
#     rt <- rt[!map_lgl(rt, is.null)]
#     map2(rt, names(rt), ~try(treemapPlot(.x, title = paste(clus_name, .y))))
#     
# })
# dev.off()
# 
# 
# #### GO Heatplots ####
# 
# reducedTerms_list_long <- map(reducedTerms_list, ~do.call("rbind", .x))
# 
# ## check most common parent terms
# map(reducedTerms_list_long, ~.x |> count(parentTerm) |> arrange(-n) |> head())
# 
# ## check most common terms
# map(reducedTerms_list_long, ~.x |> count(term) |> arrange(-n) |> head())
# 
# ## check GO terms w/ most genes
# compare_clus |> arrange(-Count) |> head()
# 
# ## source functions
# source(here("code", "13_compile_DGE", "GO_logFC_heatmap.R"))
# 
# # compare_clus |> filter(grepl("MAPT", geneID))
# 
# #### GO heatmap by datatype ####
# if(opt$datatype == "Visium"){
#     
#     go_terms1 <- c("myelin assembly", 
#                    "cell-cell recognition",
#                    "myelination", 
#                    "ensheathment of neurons", 
#                    "axon ensheathment")
#     
#     go_terms1 %in% go_terms_myelination
#     
#     go_stats_multi1 <- get_go_DE_stats_multi(go_list = go_terms1)
#     
#     pdf(here(plot_dir, sprintf("GO_logFC_heatmap_%s.pdf", opt$datatype)))
#     
#     # GO_logfc_Heatmap(get_go_DE_stats("sperm−egg recognition"))
#     GO_logfc_Heatmap(get_go_DE_stats("learning"))
#     GO_logfc_Heatmap(get_go_DE_stats("myelin assembly"))
#     GO_logfc_Heatmap(get_go_DE_stats("central nervous system myelination"))
#     GO_logfc_Heatmap(get_go_DE_stats("oligodendrocyte differentiation"))
#     GO_logfc_Heatmap(get_go_DE_stats("oligodendrocyte differentiation"))
#     GO_logfc_Heatmap(go_stats_multi1)
#     
#     dev.off()
#     
#     parent_term_heatmap(search_term = c("myelination",
#                                         "oligodendrocyte differentiation",
#                                         "calcium ion transmembrane import into cytosol",
#                                         "negative regulation of developmental growth",
#                                         "dendrite terminus",
#                                         "mitotic spindle midzone",
#                                         "dipeptidase activity"), 
#                         pdf_suffix = "ParentTerms")
#     
#     
# } else if(opt$datatype == "sn_broad"){
#     
#     pdf(here(plot_dir, sprintf("GO_logFC_heatmap_MF_%s.pdf", opt$datatype)))
#     # Astro Up
#     GO_logfc_Heatmap(get_go_DE_stats_multi(c("channel activity", "passive transmembrane transporter activity","GABA receptor activity")))
#     # Astro Down
#     GO_logfc_Heatmap(get_go_DE_stats_multi(c("exopeptidase activity", "dipeptidase activity","oligopeptide binding")))
#     # Inhib UP
#     GO_logfc_Heatmap(get_go_DE_stats_multi(c("heparan sulfate sulfotransferase activity", "proteoglycan sulfotransferase activity","microfilament motor activity")))
#     # Macro down
#     GO_logfc_Heatmap(get_go_DE_stats_multi(c("G protein-coupled purinergic nucleotide receptor activity","semaphorin receptor activity","purinergic nucleotide receptor activity")))
#     # Micro down
#     GO_logfc_Heatmap(get_go_DE_stats_multi(c("protein folding chaperone", "hydrolase activity, hydrolyzing O-glycosyl compounds","unfolded protein binding")))
#     # Oligo up
#     GO_logfc_Heatmap(get_go_DE_stats_multi(c("collagen binding", "high voltage-gated calcium channel activity")))
#     
#     dev.off()
#     
#     
#     
#     ## select GO terms
#     pdf(here(plot_dir, sprintf("GO_logFC_heatmap_%s.pdf", opt$datatype)))
#     
#     # Astro Up
#     GO_logfc_Heatmap(get_go_DE_stats(go_term = "modulation of chemical synaptic transmission"))
#     
#     GO_logfc_Heatmap(get_go_DE_stats(go_term = "oligodendrocyte differentiation"))
#     GO_logfc_Heatmap(get_go_DE_stats(go_term = "heparin proteoglycan"))
#     
#     GO_logfc_Heatmap(get_go_DE_stats(go_term = "cell fate commitment"))
#     GO_logfc_Heatmap(get_go_DE_stats(go_term = "passive transmembrane transporter activity"))
#     
#     # inhib
#     GO_logfc_Heatmap(get_go_DE_stats_multi(c("central nervous system myelination", "oligodendrocyte differentiation")))
#     
#     dev.off()
#     
#     # go_lookup("dipeptidase activity")
#     # parent_term_lookup("peptid")
#     # compare_clus |> filter(grepl("pepti", Description))
#     compare_clus |> filter(grepl("hydrolase activity", Description))
#     
#     parent_term_heatmap(search_term = c("synaptic membrane",
#                                         "oligodendrocyte differentiation",
#                                         "neurotransmitter transport",
#                                         "cellular response to calcium ion",
#                                         "learning or memory",
#                                         "early endosome membrane",
#                                         "external encapsulating structure"), 
#                         pdf_suffix = "ParentTerms")
#     
#     
#     pdf(here(plot_dir, sprintf("GO_logFC_heatmap_%s-TOP.pdf", opt$datatype)), height = 10, width = 10)
#     GO_logfc_Heatmap(get_go_DE_stats_multi(c("synaptic membrane", 
#                                              "external encapsulating structure",
#                                              
#                                              "regulation of trans-synaptic signaling",
#                                              "extracellular structure organization",
#                                              "regulation of membrane potential")),
#                      title = "cellular response to calcium ion")
#     dev.off()
# } else if(opt$datatype == "sn_fine"){
#     
#     # compare_clus |> filter(grepl("MAPT", geneID))
#     # compare_clus |> filter(grepl("FOS", geneID))
#     
#     ## select GO terms
#     pdf(here(plot_dir, sprintf("GO_logFC_heatmap_%s.pdf", opt$datatype)), width = 12, height = 12)
#     
#     GO_logfc_Heatmap(get_go_DE_stats(go_term = "neuronal cell body"))
#     GO_logfc_Heatmap(get_go_DE_stats(go_term = "synaptic membrane"))
#     GO_logfc_Heatmap(get_go_DE_stats(go_term = "metal ion transmembrane transporter activity"))
#     GO_logfc_Heatmap(get_go_DE_stats(go_term = "response to calcium ion"))
#     GO_logfc_Heatmap(get_go_DE_stats(go_term = "oligodendrocyte differentiation"))
#     GO_logfc_Heatmap(get_go_DE_stats(go_term = "main axon"))
#     GO_logfc_Heatmap(get_go_DE_stats_multi(c("main axon", "neuronal cell body")))
#     GO_logfc_Heatmap(get_go_DE_stats_multi(c("memory", "cognition")))
#     GO_logfc_Heatmap(get_go_DE_stats_multi(c("central nervous system myelination", "oligodendrocyte differentiation")))
#     
#     dev.off()
#     
#     # parent_term_lookup("calcium")
#     
#     parent_term_heatmap(search_term = c("myelin sheath",  # down in Astro.3 + Oligo.3
#                                         "myelination",
#                                         "cell-cell junction",
#                                         "synaptic membrane",  # Oligo.3 Up
#                                         "metal ion transmembrane transporter activity",
#                                         "regulation of metal ion transport",
#                                         "response to calcium ion",
#                                         "actin filament-based movement"), 
#                         pdf_suffix = "ParentTerms",
#                         height = 14, width = 12)
#     
# }


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()




