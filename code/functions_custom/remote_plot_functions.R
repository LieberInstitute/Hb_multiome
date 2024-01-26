########################################################################
##
## FUNCTIONS TO HANDLE SEURAT PLOTS 
##
## Authors. CSC
## Date. May 17th, 2023
## LM: Jan.2024
########################################################################

library(ggplot2)
library(patchwork)

get_Vplots_integrated_GEX  <- function(seuratOBJ, sfeature, stitle, imax) {
    # CSC - Sep 06, 2023 / Sep 28th, 2023
    # Generic function to build violin plot for one Seurat Object 
    # Arguments:
    #       @seuratOBJ: the seurat object
    #       @sfeature: any of 'nFeature_RNA', 'nCount_RNA', percent.mt', 'percent.ribo'
    #       @stitle: title to be used
    #       @imax: y-axis limit

    # Call Vln plot from seurat suite. (ggplot2 wrapped)    
    if (sfeature=='nFeature_RNA') {ylabel<-'Gene counts (Abs)'}
    if (sfeature=='nCount_RNA') {ylabel<-'UMI counts (Abs)'}
    if (sfeature=='percent.mt') {ylabel<-'Percentual mt levels'}
    if (sfeature=='percent.ribo') {ylabel<-'Percentual ribo levels'}
    #& geom_boxplot(width=0.1,fill="white")

    p1 <- VlnPlot(object = SeuratOBJ, features = c(sfeature), y.max = imax,
                  group.by = 'orig.ident', pt.size = 0) & geom_boxplot() &
        theme(legend.position = 'none',
              axis.text.x = element_text(angle=0, hjust=1, size=8),  #10
              axis.text.y = element_text(size=8),  #10
              axis.title.x = element_blank(),
              axis.title.y = element_blank()) &
              labs(title = "", x = 'Samples', y ="")
              #ylabs("") &
              #ylab("") &
        ggtitle(stitle)

    return(p1)
} 


get_Vplots_main_GEX  <- function(seuratOBJ, largest_genes = FALSE) {
#Old function name: f_plt_Vplots_main_GEX
    # Violin plot for a Seurat Object unique or composed for several assays  
    if ((!is.null(seuratOBJ@meta.data$percent.ribo)) && (!is.null(seuratOBJ@meta.data$percent.Largest.Gene))) { 
        p1 <- VlnPlot(object = seuratOBJ, features = c('nCount_RNA','nFeature_RNA',
                                                 'percent.mt', 'percent.ribo', 
                                                 'percent.Largest.Gene'), group.by = "orig.ident")  
    } else {
        p1 <- VlnPlot(object = seuratOBJ, features = c('nCount_RNA','nFeature_RNA','percent.mt'), group.by = "orig.ident") 
    }
    p1 <- p1 &
        theme(#legend.position = 'none',
              axis.text.x = element_text(angle=0, hjust=1, size=8),  #10
              axis.text.y = element_text(size=8),  #10
              axis.title.x = element_blank(),
              axis.title.y = element_blank()) # &labs(title = "", x = 'Samples', y ="")
    
    return(p1)
    
} 
        # # Customizing the output
    # # width	Width (defaults to 480 pixels)
    # # height	Height (defaults to 480 pixels)
    # # unit	Size unit (“px”, “in”, “cm” and “mm”)
    # # pointsize	Size of the plotted text
    # # bg	Initial background color
    # # res	Resolution in ppi
    # # type	“cairo”, “Xlib”, “quartz”
    # # (and “cairo-png” for PNG)    
    # 
    #     # Call the pdf command to start the plot
    # # Color model (cmyk is required for most publications)
    # pdf(fname, width = 8, height = 7, bg = "white", colormodel = "cmyk", paper = "A4") 
    # print(p1)
    # # Closing the graphical device
    # dev.off() 


# Visualize number of cells per sample 
get_plt_cells_by_sample  <- function(genes_per_cell) {
    genes_per_cell %>%
        ggplot(aes(x=orig.ident, fill=orig.ident)) +
        geom_bar(alpha = 0.7) +
        theme_classic() +
        #    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
        theme(plot.title = element_text(hjust=0.5)) +   #, face="bold"
        xlab("Sample ID") +
        ggtitle("Cells per sample") + 
        theme(plot.title = element_text(size = 10, face = "bold"))
              
} 

# Visualize the number UMIs/transcripts per cell
get_plt_UMIs_per_cell  <- function(genes_per_cell) {
    genes_per_cell %>%
        ggplot(aes(color=orig.ident, x=nCount_RNA, fill=orig.ident)) +
        geom_density(alpha = 0.2) +
        scale_x_log10() +
        theme_classic() +
        theme(plot.title = element_text(hjust=0.5)) +   #, face="bold"
        ylab("Log10 Cell Density") +
        xlab("UMIs") +
        ggtitle("UMIs per cell") +  
        theme(plot.title = element_text(size = 10, face = "bold"))
        geom_vline(xintercept = 500, linetype=2)     # should generally be above 500
}

# Plot the distribution of genes per cell (histogram)
get_plt_genes_per_cell_density <- function(genes_per_cell) {
    genes_per_cell %>%
        ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) +
        geom_density(alpha = 0.2) +
        scale_x_log10() +
        theme_classic() +
        theme(plot.title = element_text(hjust=0.5)) +
        geom_vline(xintercept = 300) +
        ylab("Log10(UMIs)") +
        xlab("Gene-counts") +
        ggtitle("Genes density by cell") 
}

# Plot the distribution of genes detected per cell via boxplot
get_plt_genes_per_cell_boxplot  <- function(genes_per_cell) {
    genes_per_cell %>%
        ggplot(aes(x=orig.ident, y=log10(nFeature_RNA), fill=orig.ident)) +
        geom_boxplot(alpha = 0.7) +
        theme_classic() +
        theme(axis.text.x = element_text(vjust = 1, hjust=1)) +
        theme(plot.title = element_text(hjust=0.5)) +
        ylab("Log10(gene-counts)") +
        xlab("") +
        ggtitle("Genes distribution by cell")
}

# Plot a linear model among the UMIs vs the Genes against the MT levels
get_plt_UMIS_genes_MT_geomlm  <- function(genes_per_cell, i_vline = 0, i_hline = 0) {

    p1 <- genes_per_cell %>%
            ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt, group.by = 'orig.ident')) + # MTRatio
            #    ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=MTRatio)) + # MTRatio
            geom_point() +
            scale_colour_gradient(low = "gray90", high = "black") +
            stat_smooth(method=lm) +
            scale_x_log10() +
            scale_y_log10() +
            theme_classic() +
            # if (i_vline > 0) {
            #     geom_vline(xintercept = i_vline, linetype=2)} +
            # if (i_hline > 0) {
            #     geom_hline(yintercept = i_hline, linetype=2)} +
            # facet_wrap(~orig.ident) +
            ylab("log10(genes-counts)") +
            xlab("log10(UMI-counts)") +
            ggtitle('UMIs/Genes by MT levels')

    return(p1)

}

