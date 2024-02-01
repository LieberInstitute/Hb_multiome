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

source(here("code/functions_custom", "remote_DGE_marker_gene_lists.R"))       # Call functions to read paths


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

# We have access to get_erik_markers_genes_HPC() and get_bukola_markers_genes_Hb()
markers.custom <- get_bukola_markers_genes_Hb()
markers.custom 
markers.custom$neuron
# [1] "SYT1"   "SNAP25"

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
path_cellranger_DGE_clust_df      # ~/Hb_multiome/processed-data/cellrangerARC/S1_Hb_KDM/outs/analysis/clustering/graphclust/differential_expression.csv"
cellr_clusters <- as.data.frame(read.csv(path_cellranger_DGE_clust_df, header = TRUE))
head(cellr_clusters[1:5], n=3)
# Feature.ID Feature.Name Cluster.1.Mean.Counts Cluster.1.Log2.fold.change
# 1 ENSG00000243485  MIR1302-2HG                     0                   7.110087
# 2 ENSG00000237613      FAM138A                     0                   7.110087
# 3 ENSG00000186092        OR4F5                     0                   7.110087
# Cluster.1.Adjusted.p.value
# 1                          1
# 2                          1
# 3                          1


############################# Marker gene list ################################

    #levels/categories of cell-types
message('Parsing ', length(markers.custom), ' cell-types levels')

# testing clusters
numb <- 9
Cluster_number <- paste0('Cluster.',numb,'.')
Cluster_number
f <-  paste0(Cluster_number,'Adjusted.p.value')

top10_DGE_clust <- cellr_clusters %>%
    dplyr::arrange(get(f)) %>%
    select(c(Feature.ID, Feature.Name), 1, starts_with(Cluster_number)) %>%
    dplyr::filter(get(f) < 0.05) %>%
    slice_head(n = 10)
head(top10_DGE_clust)

# get a vector wit all marker genes 
gm_lst <- as.vector(as.list(markers.custom))
#gm_lst <- as.vector(markers.custom)
names(gm_lst[1])

# # testing vector levels
i_pos <- 0
for ( x in gm_lst ) {
    i_pos <- i_pos+1
    print(paste('i_post', i_pos, ' x ', x))
    cell_type <- names(gm_lst[i_pos])
    print(cell_type)
    #print(x)
}



####### testing chunk    ####### 

clusters <- sort(clusters)
# empty df to save cell-types matched 
all_gene_match <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("Feature.ID", "Feature.Name", "Cluster.Adjusted.p.value", "cell-type", "cluster"))

for (clust in clusters) {

    # Read cluster x and extract the 10 ten most relevant genes
    Cluster_number <- paste0('Cluster.',clust,'.')
    f <-  paste0(Cluster_number,'Adjusted.p.value')

    top10_DGE_clust <- cellr_clusters %>% 
        dplyr::arrange(get(f)) %>% #select(c(Feature.ID, Feature.Name, f)) %>%
        select(c(Feature.ID, Feature.Name), 1, starts_with(f)) %>%
        dplyr::filter(get(f) < 0.05) %>%
        slice_head(n = 10)
    
    # Parse each gene in the top10 list against the marker gene list provided 
    if ( nrow(top10_DGE_clust) > 0 ) {
        
        i_pos <- 0      # reset gene-marker list position
        
        for ( x in gm_lst ) {

            i_pos <- i_pos+1                        # to extract cell type position
            cell_type <- names(gm_lst[i_pos])       # to extract cell type name
            gene_markers <- x                       # marker genes string vector
            
            # Match top10genes with the marker genes for the cell-type x 
            gene_match <- top10_DGE_clust %>% filter_all(any_vars(. %in% gene_markers))
            print(gene_match)
            
            # add matched genes to a dataframe
            if ( nrow(gene_match) > 0 ) {
                # rename column to allow rbind
                names(gene_match)[names(gene_match) == f ] <- "Cluster.Adjusted.p.value"
                gene_match['cell-type']  <- cell_type
                gene_match['cluster']  <- clust
                all_gene_match <- rbind(all_gene_match, gene_match)
            }
            
        }
        
    }
    
} 

# save the matched genes for the corresponding sample 
all_gene_match
path_cellranger_clusters_markers <- here('processed-data/05_DiffExpr_Clustering', paste0(s_sample,'_cell_types.csv'))
write.csv(all_gene_match, path_cellranger_clusters_markers, row.names=FALSE)



message(' Cell type identification in clusters done!')



################################################################################


library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
#"2023-04-04 12:42:26 EDT"
proc.time()
options(width = 120)
session_info()

