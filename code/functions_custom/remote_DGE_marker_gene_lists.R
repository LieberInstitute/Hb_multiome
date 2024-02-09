########################################################################
##
## FUNCTIONS TO BUILD MARKER GENE LISTS FOR HABENULA
##
## Authors. CSC
## Date. Feb 1st, 2024
## 
########################################################################

# load libraries
library(tidyverse)
library(dplyr)
library(here)
library(readxl)

here::here()



get_erik_markers_genes_HPC <- function() {
    
    # Load the Erik's marker gene list based on HPC from Human Brain

    markers.custom = list(
      'neuron' = c('SYT1', 'SNAP25'), #'SNAP25', 'GRIN1','MAP2'),
      'excitatory_neuron' = c('SLC17A6', 'SLC17A7'), # 'SLC17A8'),
      'inhibitory_neuron' = c('GAD1', 'GAD2'), #'SLC32A1'),
      'mediodorsal thalamus'= c('EPHA4','PDYN', 'LYPD6B', 'LYPD6', 'S1PR1', 'GBX2', 'RAMP3', 'COX6A2', 'SLITRK6', 'DGAT2'),
      'Hb neuron specific'= c('POU2F2','POU4F1','GPR151','CALB2'),#,'GPR151','POU4F1','STMN2','CALB2','NR4A2','VAV2','LPAR1'),
      'MHB neuron specific' = c('TAC1','CHAT','CHRNB4'),#'TAC3','SLC17A7'
      'LHB neuron specific' = c('HTR2C','MMRN1'),#'RFTN1'
      'oligodendrocyte' = c('MOBP', 'MBP'), # 'PLP1'),
      'oligodendrocyte_precursor' = c('PDGFRA', 'VCAN'), # 'CSPG4', 'GPR17'),
      'microglia' = c('C3', 'CSF1R'), #'C3'),
      'astrocyte' = c('GFAP', 'AQP4')
    )
    
    return(markers.custom)
    
}

get_bukola_markers_genes_Hb <- function() {
    
    # Load the Bukola/Louise's marker gene list based on Habenula from Human Brain
    
    markers.custom = list(
        'neuron' = c('SYT1', 'SNAP25'), # 'GRIN1','MAP2'),
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
    
    
    return(markers.custom)
    
}


get_Top50r_markers_genes_Hb <- function() {
    
    ## Load the Top 50 marker gene for habenula from Human Brain
   
    ## read top50 by ratio gene marker list
    s_path_name <- here('data', 'sfigu_top_50_MarkerGenes_Table.xlsx')
    Hb_gene_markers <- as.data.frame(read_excel(s_path_name, na = "---")) #sheet = "data"
    #f <-  c('LHb, MHb') 
    ## Note: cellType.target is the column you want to use. That is the "target" cell type that the data corresponds to, 
    ## the second cellType is the second highest non-target cell type (so the cell type we are comparing the target cell type to)
    
    top50_LHb_genes_byratio <- Hb_gene_markers %>% 
        dplyr::arrange(cellType.target) %>% 
        select(c(cellType.target, rank_ratio, Symbol, log.p.value), 1) %>%   #, starts_with(f)
        dplyr::filter((cellType.target == 'LHb')) #| (cellType.target == 'MHb') %>%
        #slice_head(n = 50)
    
    top50_MHb_genes_byratio <- Hb_gene_markers %>% 
        dplyr::arrange(cellType.target) %>% 
        select(c(cellType.target, rank_ratio, Symbol, log.p.value), 1) %>%
        dplyr::filter((cellType.target == 'MHb')) #| (cellType.target == 'MHb') %>%
    
    #nrow(top50_LHb_genes_byratio)
    #tail(top50_LHb_genes_byratio, n=5)
    
    if ( (nrow(top50_LHb_genes_byratio) > 0) &  (nrow(top50_MHb_genes_byratio) > 0) ) {
        markers.custom = list(
            'LHb' = top50_LHb_genes_byratio$Symbol, # 'GRIN1','MAP2'),
            'MHb' = top50_MHb_genes_byratio$Symbol # 'SLC17A8'),
        )
    }
    
    return(markers.custom)
    
}
