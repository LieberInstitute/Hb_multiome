########################################################################
## Read cellranger-arc clusters to look for cell-types (Fast cluster exploration)
## INPUT:
## Gene-marker list
## CSV file with clusters produced by cellranger-arc
## 
## Authors. CSC
## Date. Jan 26th, 2024 / last md. XXX 
########################################################################

# load libraries
library(tidyverse)
library(dplyr)
library(here)

here::here()




############  Initials ############

## commandArgs scans the arguments which have been supplied when the current R script was invoked (from shell sh)
sample_tmp <- commandArgs(trailingOnly = TRUE)
#sample_tmp <- args[1]
# testing
#sample_tmp <- 'S1_Hb_KDM,human'  # testing HUMAN tissue
#sample_tmp <- 'S2_Hb_KDM,human'  # testing HUMAN tissue
#sample_tmp <- '2_HPC_KDM,human'  # testing HUMAN tissue
#sample_tmp <- '2_HPC_KDM,human'  # testing HUMAN tissue
sample_data = unlist(strsplit(sample_tmp,","))

s_sample <- sample_data[[1]]
s_tissue <- sample_data[[2]]
message('Processing sample: ',s_sample, ' from ', s_tissue, ' tissue.')


############################# Marker gene list ################################
# [MUST HAVE AT LEAST TWO GENES PER MARKER FOR FUNCTION TO WORK CORRECTLY]
Erick.markers.custom = list(
    'neuron' = c('SYT1', 'SNAP25'), # 'GRIN1','MAP2'),    # NOTE: MT-ATP8, Fake assignations for testing purposes. CSC
    #'neuron' = c('SYT1', 'SNAP25'), # 'GRIN1','MAP2'),
    'excitatory_neuron' = c('SLC17A6', 'SLC17A7'), # 'SLC17A8'),
    'inhibitory_neuron' = c('GAD1', 'GAD2'), #'SLC32A1'),
    'mediodorsal thalamus'= c('EPHA4','PDYN', 'LYPD6B', 'LYPD6', 'S1PR1', 'GBX2', 
                              'RAMP3', 'COX6A2', 'SLITRK6', 'DGAT2', "ADARB2"), # LYPD6B and ADARB2*, EPHA4 may not be best
    'Pf/PVT' = c("INHBA", "NPTXR"), 
    'Hb neuron specific'= c('POU4F1','GPR151'), #'STMN2','NR4A2','VAV2','LPAR1'), #CALB2 and POU2F2 were thrown out because not specific enough]
    'MHB neuron specific' = c('TAC1','CHAT','CHRNB4', "TAC3"),#'SLC17A7'
    'LHB neuron specific' = c('HTR2C','MMRN1', "ANO3"),#'RFTN1'
    'oligodendrocyte' = c('MOBP', 'MBP'), # 'PLP1'),
    'oligodendrocyte_precursor' = c('PDGFRA', 'VCAN'), # 'CSPG4', 'GPR17'),
    'microglia' = c('C3', 'CSF1R'), #'C3'),
    'astrocyte' = c('GFAP', 'AQP4'),
    "Endo/CP" = c("TTR", "FOLR1", "FLT1", "CLDN5")
)
Erick.markers.custom 



#############################  DGE gene lists from cellranger-arc ################################


# Extract number of clusters for the given sample
path_cellranger_clusters_df <- here(paste0('processed-data/cellrangerARC/', s_sample, '/outs/analysis/clustering/gex/graphclust'),
                                 'clusters.csv')
cellranger_clusters <- as.data.frame(read.csv(path_cellranger_clusters_df, header = TRUE))
cellr_clusters <- unique(cellranger_clusters['Cluster'])
clusters <- cellr_clusters[['Cluster']]

message('Looking gene markers for ', length(clusters), ' clusters for sample ', s_sample)


# Extract DGE genes for all clusters for the given sample
# Read path to cellranger-arc DGE clusters
path_cellranger_DGE_clust_df <- here(paste0('processed-data/cellrangerARC/', s_sample, '/outs/analysis/clustering/gex/graphclust'),
                                    'differential_expression.csv')
path_cellranger_DGE_clust_df
# [1] "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/cellrangerARC/S1_Hb_KDM/outs/analysis/clustering/graphclust/differential_expression.csv"
cellr_clusters <- as.data.frame(read.csv(path_cellranger_DGE_clust_df, header = TRUE))
head(cellr_clusters, n=3)



############################# Marker gene list ################################

# testing
Erick.markers.custom$neuron
# [1] "SYT1"   "SNAP25"

length(Erick.markers.custom)     #levels/categories of cell-types

#testing
# numb <- 9 
# Cluster_number <- paste0('Cluster.',numb,'.')
# Cluster_number
# f <-  paste0(Cluster_number,'Adjusted.p.value')
# 
# top10_DGE_clust <- cellr_clusters %>% 
#     dplyr::arrange(get(f)) %>%
#     select(c(Feature.ID, Feature.Name), 1, starts_with(Cluster_number)) %>%
#     dplyr::filter(get(f) < 0.05) %>%
#     slice_head(n = 10)
# head(top10_DGE_clust)
# 
# 
gm_lst <- as.list(Erick.markers.custom)
gm_lst <- as.vector(Erick.markers.custom)
names(gm_lst[1])

#testing
# i_pos <- 0
# for ( x in gm_lst ) {
#     i_pos <- i_pos+1
#     cell_type <- names(gm_lst[i_pos])
#     print(cell_type)
#     print(x)
# }



####### testing chunk    ####### 

clusters <- sort(clusters)
# empty df to save cell-types found
all_gene_match <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("Feature.ID", "Feature.Name", "Cluster.Adjusted.p.value", "cell-type", "cluster"))

for (clust in clusters) {
    
    #test. clust <- 9
    Cluster_number <- paste0('Cluster.',clust,'.')
    f <-  paste0(Cluster_number,'Adjusted.p.value')
    #print(f)    # Cluster.1.Adjusted.p.value
    top10_DGE_clust <- cellr_clusters %>% 
        dplyr::arrange(get(f)) %>%
        #select(c(Feature.ID, Feature.Name, f)) %>%
        select(c(Feature.ID, Feature.Name), 1, starts_with(f)) %>%
        dplyr::filter(get(f) < 0.05) %>%
        slice_head(n = 10)
    
    # Parse cell-types on the given cluster
    if ( nrow(top10_DGE_clust) > 0 ) {
        
        #print(head(top10_DGE_clust, n=3))
        #map(top10_DGE_clust, map_celltypes)
        #bnew <- TRUE    # first row in all_gene_match df 
        
        i_pos <- 0      # reset gene-marker list
        
        for ( x in gm_lst ) {
            # testing
            #x <- 1 # neuron
            #x <- 10 # excitatory_neuron
            
            i_pos <- i_pos+1    # extract celltype name 
            cell_type <- names(gm_lst[i_pos])
            #message(' Looking gene markers for ', cell_type)
            gene_markers <- x   
            # testing
            #gene_markers <- c('SYT1', 'SNAP25', 'MT-ATP8')
            #gene_markers = c('PDGFRA', 'VCAN')
            
            # match cell-types in the filtered DGE df
            gene_match <- top10_DGE_clust %>% filter_all(any_vars(. %in% gene_markers))
            
            if ( nrow(gene_match) > 0 ) {
                # rename column to allow rbind
                names(gene_match)[names(gene_match) == f ] <- "Cluster.Adjusted.p.value"

                # PENDIENTE PARSEAR POR RENGLON CUANDO HAY MAS DE UNA COINCIDENCIA. CSC
                #for (r in 1:nrow(gene_match)) {
                #    gene_match[r, ][c('cell-type', 'cluster')]  <- c(cell_type, clust)
                #    print(gene_match[r, ])

                    gene_match[c('cell-type', 'cluster')]  <- c(cell_type, clust)
                    all_gene_match <- rbind(all_gene_match, gene_match)
                    
                #}
                
            }
            
        }
    }
    
} 

all_gene_match
path_cellranger_clusters_markers <- here('processed-data/05_DiffExpr_Clustering', paste0(s_sample,'_cell_types.csv'))
write.csv(all_gene_match, path_cellranger_clusters_markers, row.names=FALSE)


##################### Load gene-markers from Erik's habenula list ##################### 



# 
# 
# my_markers <- as.data.frame(read.csv(s_file_name, header = TRUE))
# 
# ## Build another data frame for the Nelson gene-markers list to compare gene-markers
# ## Notice Nelson's markers are probably comming from mouse donors because these are in lower case
# nlst_Astrocytes <- toupper(list("Aqp4", "Gfap"))
# nlst_Oligodendrocytes <- toupper(list("Mbp", "Mobp"))
# nlst_Oligodendrocyte_precursor <- toupper(list("Pdgfra", "Vcan"))
# 
# nlst_Microglia <- list("Csf1r", "C3", "Tmem119")
# nlst_Choroid_plexus <- list("Ttr", "Folr1")
# nlst_Vascular_endothelial <- list("Cldn5", "Flt1")
# nlst_Neurons_general <- list("Syt1", "Snap25")
# nlst_IN_neurons <- list("Gad1", "Gad2", "Sst", "Pvalb", "Htr3a")
# nlst_EX_neurons <- list("Slc17a7", "Slc17a6", "Nrgn")
# n_lst_EXsub_GranCells <- list("Prox1", "Glis3")
# n_lst_EXsub_CA4 <- list("Calb2", "Csf2rb2")
# n_lst_EXsub_CA3 <- list("Nectin3")
# n_lst_EXsub_CA2 <- list("Ntsr2", "Rgs14")
# n_lst_EXsub_CA1 <- list("Man1a", "Fibcd1", "Tbc1d1")
# n_lst_EXsub_SubProS <- list("Fn1", "Nts", "Ndst4", "Dcn")
# 
# nelson_gene_markers <- list(nlst_Astrocytes, nlst_Oligodendrocytes, nlst_Oligodendrocyte_precursor)
# nelson_gene_markers$gene
# # View cluster levels
# View(Idents(SeuratOBJ))
# clusters <- 0:cluster_idx
# 
# map_celltypes <- function(gene_marker_lst) {
#     ## Parse known gene markers to annotate cluster cell-types 
#     
#     show(gene_marker_lst)
#     for (gene in (gene_marker_lst)) {
#         my_markers %>% 
#             rowwise %>% 
#             mutate(firstZero = names(.)[match(0, c_across(everything()), 
#                 nomatch = ncol(.))]) %>% 
#                 ungroup
#     }
#     
#     #new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
#     #                     "NK", "DC", "Platelet")
#     names(new.cluster.ids) <- levels(pbmc)
#     pbmc <- RenameIdents(pbmc, new.cluster.ids)
#     DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
#     
# }


#map(top10All_DGE_df, map_celltypes)


message(' Cell types assingned')




################################################################################


library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
#"2023-04-04 12:42:26 EDT"
proc.time()
options(width = 120)
session_info()

