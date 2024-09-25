########################################################################
##
## Additional FUNCTIONS for SEURAT OBJECTS 
##
########################################################################

get_basic_stats_GEX <- function(seuratOBJ, s_tissue) {
## Calculate statistics over the GEX seurat assay: nCount, nFeature, MITO, Ribo, #cells, etc.
## INPUT: 
##      @seuratOBJ: Seurat obj 
##      @s_tissue : string used to distinguish btw human and mouse gene MITO phenotype.    
## OUTPUT:
##      @tab_stats: table with the statistics about the GEX Seurat Assay
##      Print a file with the same data for further analysis
##
    
tryCatch( {

    # Initialize metrics
    s_sample = toString(unique(seuratOBJ@meta.data$orig.ident))
    totalCells <- ncol(seuratOBJ)
    # Continuous sample quantile types 4: linear interpolation of the Empirical Cumulative Distribution Function 
    nCount = quantile(seuratOBJ$nCount_RNA, na.rm = TRUE, type = 4)       #  if na.rm = true, any NA and NaN's are removed
    nFeature = quantile(seuratOBJ$nFeature_RNA, na.rm = TRUE, type = 4)   
    #nMT = quantile(seuratOBJ$percent.mt, na.rm = TRUE, type = 4)
    #Calculate deciles and extract the 10th and 90th percentiles of the vector
    nMT = quantile(seuratOBJ$percent.mt, c(.01,.10,.25,.50,.75,.90,.99), na.rm = TRUE, type = 4)

    # Identify the number if mitochondrial genes by their names starting gene name pattern
    if (s_tissue=='human') {
        mt_genes = length(grep("^MT-",rownames(seuratOBJ@assays$RNA@counts),value = TRUE))
    } else { # it is mouse
        mt_genes = length(grep("^mt",rownames(seuratOBJ@assays$RNA@counts),value = TRUE))
    }
    
    # Identify the number of ribosomal genes, that usually tend to be very highly represented, 
    #           This can vary between cell types,check how prevalent they are in the data.
    # if (s_tissue=='human') {
    #     ribo_genes = length(grep("^RP[LS]",rownames(seuratOBJ@assays$RNA@counts),value = TRUE))
    # } else {    # it is mouse
    #     ribo_genes = length(grep("^Rp[ls]",rownames(seuratOBJ@assays$RNA@counts),value = TRUE))
    # }
    # print(ribo_genes)  
    # Build a table with the stats applied and the outputs gotten
    # Table with the cellranger-arc `gene expression statistics`
    #mtx_stats <- matrix(c(nCount, nFeature, nMT), ncol=15, byrow=TRUE)
    mtx_stats <- matrix(c(nCount, nFeature, nMT), ncol=17, byrow=TRUE)
    #class(mtx_stats)
    #dim(mtx_stats)
    colnames(mtx_stats) <- c('nCount0%','nCount25%','nCount50%','nCount75%','nCount100%',
                             'nGenes0%','nGenes25%','nGenes50%','nGenes75%','nGenes100%',
                             'nMito0%','nMito10p','nMito25%','nMito50%','nMito75%','nMito90p','nMito99%')
    print("x") 
    # add additional columns: sample name, total cells, #mito genes and #ribo genes.
    
    #message('nCount, nFeature and nMT statistics done.')
    mtx_stats<-cbind(mtx_stats,mt_genes)    # number of mito genes
    #mtx_stats<-cbind(mtx_stats,ribo_genes)  # number of ribo genes  
    sample <- s_sample
    mtx_stats<-cbind(mtx_stats,totalCells)
    mtx_stats<-cbind(mtx_stats,sample)
    # convert to table
    tab_stats <- as.table(mtx_stats)
    print(tab_stats)

    return(tab_stats) }
    , error = function(e) {print('An error ocurred. Verify you have an GEX object active') })
}


l_get_perc_largest_genes  <- function(seuratOBJ) {
    # Get the Percentage of Largest Gene    
    # Inspired by https://www.bioinformatics.babraham.ac.uk/training/10XRNASeq/seurat_workflow.html
    #          Having a high proportion of your data dominated by a single gene could either give biological context or indicate #          a technical problem. This fuctions is to monitor this variable.
    
    # It is find MALAT1 normally largest by some distance - it’s a non-coding nuclear gene expressed at very high levels
    seuratOBJ.nomalat <- seuratOBJ[rownames(seuratOBJ) != "MALAT1",]
    seuratOBJ.nomalat$largest_count <- apply(seuratOBJ.nomalat@assays$RNA@counts, 2, max)
    seuratOBJ.nomalat$largest_index <- apply(seuratOBJ.nomalat@assays$RNA@counts, 2, which.max)
    seuratOBJ.nomalat$largest_gene <- rownames(seuratOBJ.nomalat)[seuratOBJ.nomalat$largest_index]
    seuratOBJ.nomalat$percent.Largest.Gene <- 100 * seuratOBJ.nomalat$largest_count / seuratOBJ.nomalat$nCount_RNA
    seuratOBJ$largest_gene <- seuratOBJ.nomalat$largest_gene
    seuratOBJ$percent.Largest.Gene <- seuratOBJ.nomalat$percent.Largest.Gene
    #VlnPlot(seuratOBJ, features=c("nCount_RNA","percent.mt", "percent.ribo","percent.Largest.Gene"))
}

# Taken from https://github.com/satijalab/seurat/issues/5343
# function for inspecting a seurat object
f_InspectSeurat <- function(seuratObj) {
    library(Seurat)
    options(Seurat.object.assay.version = 'v5')
    
    print(paste0('Seurat object size: ', format(object.size(seuratObj), units = 'auto')))
    nCells <- ncol(seuratObj)
    if (nrow(seuratObj@meta.data) != nCells) {
        print(paste0('Metadata rows not equal to ncol: ', nCells, ' / ', nrow(seuratObj@meta.data)))
    }
    
    for (assayName in names(seuratObj@assays)) {
        print(paste0('Assay: ', assayName))
        for (slotName in c('counts', 'data')) { #, 'scale.data'
            dat <- Seurat::GetAssayData(seuratObj, assay = assayName, slot = slotName)  # Use 'layer' for SeuratObject v5 otherwise use  'GetAssayData' --deprecable 
            print(paste0('slot: ', slotName, ', size: ', format(object.size(x = dat), units = 'auto')))
            
            if (!is(dat, 'sparseMatrix')) {
                print(paste0('Non-sparse! Assay: ', assayName, ', slot: ', slotName))
            }
            
            if (ncol(dat) > 0 && ncol(dat) != nCells) {
                print(paste0('Assay cols not equal to ncol(object): ', nCells, ' / ', ncol(dat)))
            }
        }  
    }
    
    print('All slots:')
    for (slotName in slotNames(seuratObj)) {
        val <- object.size(slot(seuratObj, slotName))
        if (val > 1000) { 
            print(paste0(slotName, ': ', format(val, units = 'auto')))
        }
    }
    
}


get_basic_stats_GEX <- function(seuratOBJ, s_tissue) {
  ## Calculate statistics over the GEX seurat assay: nCount, nFeature, MITO, Ribo, #cells, etc.
  ## INPUT: 
  ##      @seuratOBJ: Seurat obj 
  ##      @s_tissue : string used to distinguish btw human and mouse gene MITO phenotype.    
  ## OUTPUT:
  ##      @tab_stats: table with the statistics about the GEX Seurat Assay
  ##      Print a file with the same data for further analysis
  ##
  
  tryCatch( {
    
    # Initialize metrics
    s_sample = toString(unique(seuratOBJ$orig.ident))
    totalCells <- ncol(seuratOBJ)

    # Continuous sample quantile types 4: linear interpolation of the Empirical Cumulative Distribution Function 
    nCount = round(quantile(seuratOBJ$nCount_RNA, na.rm = TRUE, type = 4), digits = 2)       #  if na.rm = true, any NA and NaN's are removed
    nFeature = round(quantile(seuratOBJ$nFeature_RNA, na.rm = TRUE, type = 4), digits = 2)   
    
    #nMT = quantile(seuratOBJ$percent.mt, na.rm = TRUE, type = 4)
    #Calculate deciles and extract the 10th and 90th percentiles of the vector
    nMT = round(quantile(seuratOBJ$percent.mt, c(.01,.10,.25,.50,.75,.90,.99), na.rm = TRUE, type = 4), digits = 2)
   
    # # Identify the number if mitochondrial genes by their names starting gene name pattern
    # if (s_tissue=='human') {
    #   mt_genes = length(grep("^MT-",rownames(seuratOBJ@assays$RNA@counts),value = TRUE))
    # } else { # it is mouse
    #   mt_genes = length(grep("^mt",rownames(seuratOBJ@assays$RNA@counts),value = TRUE))
    # }
    
    # # Identify the number of ribosomal genes, that usually tend to be very highly represented, 
    # #           This can vary between cell types,check how prevalent they are in the data.
    # if (s_tissue=='human') {
    #   ribo_genes = length(grep("^RP[LS]",rownames(seuratOBJ@assays$RNA@counts),value = TRUE))
    # } else {    # it is mouse
    #   ribo_genes = length(grep("^Rp[ls]",rownames(seuratOBJ@assays$RNA@counts),value = TRUE))
    # }
    
    # Build a table with the stats applied and the outputs gotten
    # Table with the cellranger-arc `gene expression statistics`
    #mtx_stats <- matrix(c(nCount, nFeature, nMT), ncol=15, byrow=TRUE)
    mtx_stats <- matrix(c(nCount, nFeature, nMT), ncol=17, byrow=TRUE)

    colnames(mtx_stats) <- c('nCount0%','nCount25%','nCount50%','nCount75%','nCount100%',
                             'nGenes0%','nGenes25%','nGenes50%','nGenes75%','nGenes100%',
                             'nMito0%','nMito10%','nMito25%','nMito50%','nMito75%','nMito90p','nMito99%')
    
    # add additional columns: sample name, total cells, #mito genes and #ribo genes.
    
    #message('nCount, nFeature and nMT statistics done.')
    mtx_stats<-cbind(mtx_stats,NULL)    # number of mito genes
    mtx_stats<-cbind(mtx_stats,NULL)  # number of ribo genes  
    sample <- s_sample
    mtx_stats<-cbind(mtx_stats,totalCells)
    mtx_stats<-cbind(mtx_stats,sample)
    # convert to table
    tab_stats <- as.table(mtx_stats)
    
    return(tab_stats) }
    , error = function(e) {print('An error ocurred. Verify you have an GEX object active') })
}