#Using MetaNeighbors GO term testing to identify the gene sets that characterize cell-type variability
#Can also use KEGG terms

#Do we want to focus on just the habenula cell-types?



library(SingleCellExperiment)
library(MetaNeighbor)
library(here)

library(org.Hs.eg.db)
library(GO.db)


plot_path = here('plots','05_03_annotation_adjustments', '10_celltype_func_annot')
data_path = here('processed-data','05_03_annotation_adjustments', '10_celltype_func_annot')

if (!dir.exists(data_path)) dir.create(data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)




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




  