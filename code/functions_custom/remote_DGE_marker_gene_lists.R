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
    top50_gene_markers <- as.data.frame(read_excel(s_path_name, na = "---")) #sheet = "data"
    ## Note: cellType.target is the column you want to use. That is the "target" cell type that the data corresponds to, 
    ## the second cellType is the second highest non-target cell type (so the cell type we are comparing the target cell type to)
    ## `mean-ratio` in the excel sheet is called `ratio`
    cell_types <- unlist(top50_gene_markers |> distinct(cellType.target))
    length(cell_types)
    
    library(purrr)
    names(cell_types) <- paste0("DD_", cell_types)
    markersDD.list <- map(cell_types, ~top50_gene_markers |>
                            dplyr::filter(cellType.target == .x & ratio > 1) |> 
                            select(c(cellType.target, Symbol)) |>
                           pull(Symbol)) 
                            
    # markersDD.list <- list()
    # 
    # for (cell_class in cell_types) {
    #   # for testing:  cell_class <- "Endo"
    #   #               cell_class <- "Astrocyte"
    # 
    #   top50_class_name <- paste0("DD_", cell_class)
    #   top50_DD_class_symbols <- top50_gene_markers |>
    #     dplyr::filter(cellType.target == cell_class & ratio > 1) |> 
    #     select(c(cellType.target, Symbol)) 
    #   #print(top50_class_name)
    #   #print(head(top50_DD_class_symbols))
    #   #if (length(top50_DD_class_symbols) > 0) {
    #     if (!length(markersDD.list) == 0)  {
    #         markersDD.list[[top50_class_name]] <- top50_DD_class_symbols$Symbol
    #     } else {
    #         markersDD.list <- list(top50_class_name = top50_DD_class_symbols$Symbol)
    #         names(markersDD.list) <- top50_class_name
    #     }
    #   #}
    #   
    # }
    # 
    message("Added ", names(markersDD.list), " to Data-Driven markers list")

    # top50_LHb_genes_byratio <- top50_gene_markers |> 
    #     dplyr::arrange(cellType.target) |>
    #     select(c(cellType.target, rank_ratio, Symbol, logFC, log.p.value), 1) %>%   #, starts_with(f)
    #     dplyr::filter(cellType.target == 'LHb')
    # 
    # top50_MHb_genes_byratio <- top50_gene_markers %>% 
    #     dplyr::arrange(cellType.target) %>% 
    #     select(c(cellType.target, rank_ratio, Symbol, log.p.value), 1) %>%
    #     dplyr::filter((cellType.target == 'MHb')) # (cellType.target == 'MHb') %>%
    # 
    # top50_ThalamusExcit_genes_byratio <- top50_gene_markers %>% 
    #   dplyr::arrange(cellType.target) %>% 
    #   select(c(cellType.target, rank_ratio, Symbol, log.p.value), 1) %>%
    #   dplyr::filter((cellType.target == 'Excit.Thal')) 
    # top50_ThalamusInhib_genes_byratio <- top50_gene_markers %>% 
    #   dplyr::arrange(cellType.target) %>% 
    #   select(c(cellType.target, rank_ratio, Symbol, log.p.value), 1) %>%
    #   dplyr::filter((cellType.target == 'Inhib.Thal')) 
    
    # if ( (nrow(top50_LHb_genes_byratio) > 0) &  (nrow(top50_MHb_genes_byratio) > 0) ) {
    #     markers.custom = list(
    #         'LHb_putative' = top50_LHb_genes_byratio$Symbol, 
    #         'MHb_putative' = top50_MHb_genes_byratio$Symbol, 
    #         'Thal_excit_putative' = top50_ThalamusExcit_genes_byratio$Symbol,
    #         'Thal_inhib_putative' = top50_ThalamusInhib_genes_byratio$Symbol 
    #     )
    # }
    
    return(markersDD.list)
    
}


get_erik_and_Hb_markers_genes <- function() {
    
    # ONLY literature marker gene list 
    
    # markersTop50 <- get_Top50r_markers_genes_Hb()
    # MHb_putative = c(markersTop50$MHb_putative)
    # LHb_putative = c(markersTop50$LHb_putative)
    # ThalE_putative = c(markersTop50$Thal_excit_putative)
    # ThalI_putative = c(markersTop50$Thal_inhib_putative)
    
    markers.custom = list(
        'neuron' = c('SYT1', 'SNAP25'), #'SNAP25', 'GRIN1','MAP2'),
        'excitatory_neuron' = c('SLC17A6', 'SLC17A7'), # 'SLC17A8'),
        'inhibitory_neuron' = c('GAD1', 'GAD2'), #'SLC32A1'),
        #'mediodorsal thalamus'= c('EPHA4','PDYN', 'LYPD6B', 'LYPD6', 'S1PR1', 'GBX2', 'RAMP3', 'COX6A2', 'SLITRK6', 'DGAT2'),
        'Hb neuron specific'= c('POU2F2','POU4F1','GPR151','CALB2'),#,'GPR151','POU4F1','STMN2','CALB2','NR4A2','VAV2','LPAR1'),
        'MHB neuron specific' = c('TAC1','CHAT','CHRNB4'),#'TAC3','SLC17A7'
        'LHB neuron specific' = c('HTR2C','MMRN1'),#'RFTN1'
        'oligodendrocyte' = c('MOBP', 'MBP'), # 'PLP1'),
        'oligodendrocyte_precursor' = c('PDGFRA', 'VCAN'), # 'CSPG4', 'GPR17'),
        'microglia' = c('C3', 'CSF1R'), #'C3'),
        'astrocyte' = c('GFAP', 'AQP4'),
        "Endo/CP" = c("TTR", "FOLR1", "FLT1", "CLDN5"),
        "Thalamus broad" =  "NECAB1",
        "Thalamus/MDm/Endo" =  "GBX2",
        "Thalamus/MDm/Endo-" =  "TNNT1",
        "Thalamus/MDm" =  "MEIS2"
        # 'MHb_putative' = MHb_putative,
        # 'LHb_putative' = LHb_putative,
        # 'ThalE_putative' = ThalE_putative,
        # 'ThalI_putative' = ThalI_putative
    )
    
    markers.custom
    
    return(markers.custom)
    
}
