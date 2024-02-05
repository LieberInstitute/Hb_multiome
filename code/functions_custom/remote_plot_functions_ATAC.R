########################################################################
##
## FUNCTIONS TO HANDLE ATAC PLOTS 
##
## Authors. CSC
## Date. Feb 5th, 2024
## LM: XXX
########################################################################

library(ggplot2)
library(patchwork)

get_Vplots_blackR_ATAC  <- function(seuratOBJ) {

    if ((!is.null(seuratOBJ@meta.data$pct_reads_in_peaks)) && (!is.null(seuratOBJ@meta.data$blacklist_ratio))) { 
        p1 <- VlnPlot(object = seuratOBJ, features = c("pct_reads_in_peaks","blacklist_ratio"), group.by = "orig.ident", ncol = 2)  
    } else {
        p1 <- VlnPlot(object = seuratOBJ, features = c("pct_reads_in_peaks","blacklist_ratio"), group.by = "orig.ident", ncol = 2) 
    }
    p1 <- p1 &
        theme(#legend.position = 'none',
            axis.text.x = element_text(angle=0, hjust=1, size=8),  #10
            axis.text.y = element_text(size=8),  #10
            axis.title.x = element_blank(),
            axis.title.y = element_blank()) # &labs(title = "", x = 'Samples', y ="")
    
    return(p1)
    
} 
# Plot Peaks in black ratio and ATAC main feature scoreds
#p1_BlackR <- VlnPlot(SeuratOBJ, features = c("pct_reads_in_peaks","blacklist_ratio"), ncol = 2)


get_Vplots_main_ATAC  <- function(seuratOBJ) {
    # TSS.enrichment will fail if you do not have enough memory
    if (!is.null(seuratOBJ@meta.data$TSS.enrichment)) { 
        p1 <- VlnPlot(object = seuratOBJ, features = c("nCount_ATAC", "nFeature_ATAC", "TSS.enrichment"), 
                      group.by = "orig.ident") # , ncol = 4 
    } else {
        p1 <- VlnPlot(object = seuratOBJ, features = c("nCount_ATAC", "nFeature_ATAC", "nucleosome_signal"), 
                      group.by = "orig.ident") # , ncol = 3 
    }
    p1 <- p1 &
        theme(#legend.position = 'none',
            axis.text.x = element_text(angle=0, hjust=1, size=8),  #10
            axis.text.y = element_text(size=8),  #10
            axis.title.x = element_blank(),
            axis.title.y = element_blank()) # &labs(title = "", x = 'Samples', y ="")
    
    return(p1)
    
} 
#p1_ATAC <- VlnPlot(SeuratOBJ, features = c("nCount_ATAC", "nFeature_ATAC", "nucleosome_signal", "TSS.enrichment"), ncol = 4)