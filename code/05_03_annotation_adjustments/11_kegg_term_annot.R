#Similar to the GO term annotation of cell-types, do the same with the kegg gene sets
#using the clusterprofiler R package to get the Kegg terms
 



library(SingleCellExperiment)
library(MetaNeighbor)
library(MetaMarkers)
library(qs2)
library(here)
library(sessioninfo)

library(clusterProfiler)
library(org.Hs.eg.db)


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


data_path = here('processed-data','05_03_annotation_adjustments', '11_kegg_term_annot')

if (!dir.exists(data_path)) dir.create(data_path)



#Check out the KEGG gene sets
if(file.exists(file.path(data_path, "kegg_human.rds"))) {
    kegg_sets = readRDS(file.path(data_path, "kegg_human.rds"))
} else {
  kegg_db = download_KEGG(species = 'hsa', keggType = "KEGG", keyType = "kegg")

  kegg_terms = kegg_db$KEGGPATHID2NAME$to

  kegg_sets = lapply(kegg_terms, function(gs) {kegg_id = kegg_db$KEGGPATHID2NAME$from[which(kegg_db$KEGGPATHID2NAME$to == gs)]
      kegg_genes = kegg_db$KEGGPATHID2EXTID$to[which(kegg_db$KEGGPATHID2EXTID$from == kegg_id)]
      gene_symbols = bitr(kegg_genes, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
      return(gene_symbols$SYMBOL)}
  )

  names(kegg_sets) = kegg_terms

  saveRDS(kegg_sets, here(data_path, "kegg_human.rds"))
}

#Prep the KEGG terms, filtering down to the present genes and then filtering gene sets to be between 10 and 100 genes
known_genes = rownames(multiome_sce)

kegg_sets = lapply(kegg_sets, function(gene_set) {gene_set[gene_set %in% known_genes]})
min_size = 10
max_size = 100
kegg_set_size = sapply(kegg_sets, length)
kegg_sets = kegg_sets[kegg_set_size >= min_size & kegg_set_size <= max_size]
length(kegg_sets)


# 216 gene sets to test


#Run MetaNeighbor on the KEGG gene sets

message('Starting MetaNeighbor functional annotation...')
aurocs = MetaNeighbor(dat = multiome_hab_sce,
  experiment_labels = multiome_hab_sce$orig.ident,
  celltype_labels = multiome_hab_sce$refined_mid_cluster,
  genesets = kegg_sets, 
  fast_version = TRUE, bplot = FALSE, batch_size = 50)

write.table(aurocs, here(data_path, 'multiomeHab_KEGG_functional_aurocs.txt'))


message('Starting MetaNeighbor functional annotation on the non-Hab neurons...')
aurocs = MetaNeighbor(dat = multiome_nonHab_sce,
  experiment_labels = multiome_nonHab_sce$orig.ident,
  celltype_labels = multiome_nonHab_sce$refined_mid_cluster,
  genesets = kegg_sets, 
  fast_version = TRUE, bplot = FALSE, batch_size = 50)

write.table(aurocs, here(data_path, 'multiomeNonHab_KEGG_functional_aurocs.txt'))



print("Reproducibility information:")
options(width = 120)
session_info()

  