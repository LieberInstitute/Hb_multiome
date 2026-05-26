#Using MetaNeighbors GO term testing to identify the gene sets that characterize cell-type variability
#Can also use KEGG terms

#Do we want to focus on just the habenula cell-types?



library(SingleCellExperiment)
library(MetaNeighbor)
library(MetaMarkers)
library(qs2)
library(here)
library(sessioninfo)

library(org.Hs.eg.db)
library(GO.db)


message('Loading data...')

multiome_path = here('processed-data', '05_03_annotation_adjustments','06_refined_annotations')
multiome_sce = qs_read(paste0(multiome_path, '/refined_annotation_multiomeHab_SCE.qs2'))

assay(multiome_sce, 'cpm') = MetaMarkers::convert_to_cpm(assay(multiome_sce, 'counts'))

#Drop ATAC assay
altExps(multiome_sce) <- NULL
gc()

#Filter for just the habenula cell-types
hab_celltypes = c('Inhib_LHb_4.1','Inhib_LHb_4.2','LHb.4','LHb.2.7','LHb.1.3.4','MHb.1','MHb.1.2','MHb.2','MHb.3')
multiome_hab_sce = multiome_sce[ ,multiome_sce$refined_mid_cluster %in% hab_celltypes]

#And filter for the non-neurons
multiome_nonHab_sce = multiome_sce[ ,!multiome_sce$refined_mid_cluster %in% hab_celltypes]



plot_path = here('plots','05_03_annotation_adjustments', '10_celltype_func_annot')
data_path = here('processed-data','05_03_annotation_adjustments', '10_celltype_func_annot')

if (!dir.exists(data_path)) dir.create(data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)


message('prepping go terms...')
#GOSOURCEDATE: 2025-02-06
#First set up the GO annotations, adapting from the MetaNeighbor code https://github.com/gillislab/MetaNeighbor-Protocol/blob/ae5f25eae996881518b1f2c6e72110ff8c1a1cfd/data/make_go_sets.R
full_go_name = function(go_id) {
    go_term = GOTERM[[go_id]]
    return(paste(go_id, Term(go_term), Ontology(go_term), sep = "|"))
}

#Check if the GO sets are already saved, if not make them
if (file.exists(file.path(data_path, "go_human.rds"))) {
    go_sets = readRDS(file.path(data_path, "go_human.rds"))
} else {

  go_sets = as.list(org.Hs.egGO2ALLEGS)
  all_genes = unique(unlist(go_sets))
  entrez_to_symbols = unlist(mget(all_genes, org.Hs.egSYMBOL))
  go_sets = lapply(go_sets, function(gs) { unique(entrez_to_symbols[gs]) })
  names(go_sets) = sapply(names(go_sets), full_go_name)

  saveRDS(go_sets, here(data_path, "go_human.rds"))
  
}


#Prep the GO terms, filtering down to the present genes and then filtering gene sets to be between 10 and 100 genes
known_genes = rownames(multiome_sce)

go_sets = lapply(go_sets, function(gene_set) {gene_set[gene_set %in% known_genes]})
min_size = 10
max_size = 100
go_set_size = sapply(go_sets, length)
go_sets = go_sets[go_set_size >= min_size & go_set_size <= max_size]
length(go_sets)

# 6719 gene sets to test


#Run MetaNeighbor on the GO gene sets

message('Starting MetaNeighbor functional annotation...')
aurocs = MetaNeighbor(dat = multiome_sce,
  experiment_labels = multiome_sce$orig.ident,
  celltype_labels = multiome_sce$refined_mid_cluster,
  genesets = go_sets, 
  fast_version = TRUE, bplot = FALSE, batch_size = 50)

write.table(aurocs, here(data_path, 'multiomeHab_functional_aurocs.txt'))


message('Starting MetaNeighbor functional annotation on the non-Hab neurons...')
aurocs = MetaNeighbor(dat = multiome_nonHab_sce,
  experiment_labels = multiome_nonHab_sce$orig.ident,
  celltype_labels = multiome_nonHab_sce$refined_mid_cluster,
  genesets = go_sets, 
  fast_version = TRUE, bplot = FALSE, batch_size = 50)

write.table(aurocs, here(data_path, 'multiomeNonHab_GO_functional_aurocs.txt'))



print("Reproducibility information:")
options(width = 120)
session_info()







  