########################################################################
##
## FUNCTIONS TO HANDLE SEURAT OBJECTS 
##
## Authors. CSC
## Date. May 4th, 2023
##
## Functions added:
##      get_seurat_obj()

########################################################################

# Old version used for first implementation
# load_seurat_obj  <- function(SeuratOBJ, s_seurat_name) {
# 
#     # Bug found .... need extra validation to check if the object exists ... working on it ..
#     
#     
#     
#     # Validate the Seurat object exists, of does exit
#     tryCatch( {
#         if ("Seurat" %in% class(SeuratOBJ)) return(SeuratOBJ)
#     }, error = function(e) { message('Seurat object not found. It will be loaded from disk.') })
# 
#     # If not exist, load from disk
#     if (s_seurat_name == '') stop('Seurat name not provided') 
#     message('Loading Seurat from disk if exist.')
#     if (s_seurat_name == 'pbmc3k') { SeuratOBJ <- readRDS(here('processed-data','pbmc3k.GEX.RDS')) }
#     if (s_seurat_name == 'pbmc10k') { SeuratOBJ <- readRDS(here('processed-data','pbmc10k.GEX.RDS')) }
#     if (s_seurat_name == 'hippo-42_1') { SeuratOBJ <- readRDS(here('processed-data','hippo-42_1.GEX.RDS')) }
#     if (s_seurat_name == 'hippo-42_4') { SeuratOBJ <- readRDS(here('processed-data','hippo-42_4.GEX.RDS')) }
#     # Filtered objects by GEX
#     if ( s_seurat_name == 'pbmc3k.filtered.GEX.M1' ) { SeuratOBJ <- load(here('processed-data','pbmc3k.filtered.GEX.M1.Prob.RData')) }
#     if ( s_seurat_name == 'pbmc3k.filtered.GEX.M2' ) { SeuratOBJ <- load(here('processed-data','pbmc3k.filtered.GEX.M2.SD.RData')) }
#     if ( s_seurat_name == 'pbmc10k.filtered.GEX.M1' ) { SeuratOBJ <- load(here('processed-data','pbmc10k.filtered.GEX.M1.Prob.RData')) }
#     if ( s_seurat_name == 'pbmc10k.filtered.GEX.M2' ) { SeuratOBJ <- load(here('processed-data','pbmc10k.filtered.GEX.M2.SD.RData')) }    
#     # Integrated objects
#     if ( s_seurat_name == 'pbmc10k.combined' ) { SeuratOBJ <- readRDS(here('processed-data','pbmc10k.combined.GEX.RDS')) }
#     
#     # Ensure Seurat object is a Seurat object
#     tryCatch( {
#         if (!("Seurat" %in% class(SeuratOBJ))) stop("Loading from disk fails. Check the object exists and the sample name assigned.")
#         message('Seurat loaded successfully')
#         #head(SeuratOBJ, n = 3)
#     }, error = function(e) { message('An error ocurred. Verify you have assigned a valid sample name') })
#         
#     return(SeuratOBJ)
# 
# }


get_seurat_obj_QCed <- function(SeuratObj, b_std_filtering=TRUE) {
# CSC. 01.2024
# Filter out low quality cells
#        if b_std_filtering equal TRUE apply standard filtering
    
    if (b_std_filtering) {    
        SeuratObj <- subset(
            x = SeuratObj,
            subset = nCount_ATAC < 100000 &
                nCount_RNA < 25000 &
                nCount_ATAC > 1000 &
                nCount_RNA > 1000 &
                nucleosome_signal < 2 &
                TSS.enrichment > 1
        )
    } else {
        # some customize filtering
    }
    return(SeuratObj)
    
}

get_seurat_obj <- function(seuratName, s_bc_mtx, s_tissue, s_meta, b_additional_feat = FALSE) {

    # This function creates a Seurat object with rna counts and meta-data attached
    #       @seuratName         string with the name to be assigned to the Seurat object
    #       @s_bc_mtx           string with the filtered barcode mtx file path
    #       @s_tissue           string used to distinguish btw human and mouse gene MITO phenotype.
    #       @s_meta             string with the meta_data file path 
    #       @b_additional_feat  Boolean to required calculate additional features; RIBO and largest genes (******** OPT *********)
    
    # filtered bc mtx
    mtx <- Read10X_h5(s_bc_mtx)             # Returns a sparse matrix with rows and columns labeled
    rna_counts <- mtx$`Gene Expression`
    print(head(rna_counts, n = 3))
    message('GEX counts loaded successfully')

    metadata <- read.csv(file = s_meta, header = TRUE, row.names = 1)
    # subset specific fields in the meta data df
    meta_tmp = c('atac_peak_region_fragments','atac_fragments')
    meta = metadata[meta_tmp]
    print(head(meta, n = 3))
    
    seur_obj <- CreateSeuratObject(
        counts = rna_counts,
        assay = "RNA",
        project = seuratName,
        meta.data = meta
    )
    
    message('Meta-data attached successfully')

    # Ensure seurat_obj is a Seurat object
    if (!("Seurat" %in% class(seur_obj))) { stop("Seurat object not found") } 
    #Layers(SeuratOBJ[["RNA"]])
    #if (!("counts" %in% Layers(seur_obj[["RNA"]]))) { stop("RNA counts layer not found!") } 
    #if (!("data" %in% Layers(seur_obj[["RNA"]]))) { stop("Data counts not found!") }
    
    # Calculate MITO levels
    seur_obj$log10GenesPerUMI <- log10(seur_obj$nFeature_RNA) / log10(seur_obj$nCount_RNA)
    
    # Select the correct Symbol for each genome. By default is assumed MT, belonging to humans.
    # For mouse genome we use the Symbol: mt-. 
    # https://www.ncbi.nlm.nih.gov/nuccore/34538597 
    
    # MITO GENES for human and mouse genomes. GRCh38 and mm10, respectively
    if (s_tissue=='human') {
        seur_obj[["percent.mt"]] <- PercentageFeatureSet(seur_obj, pattern = "^MT-")
    } else { # it is mouse
        seur_obj[["percent.mt"]] <- PercentageFeatureSet(seur_obj, pattern = "^Mt")
    }

    # RIBO GENES for human and mouse genomes. GRCh38 and mm10, respectively    
    if (s_tissue=='human') {
        seur_obj[["percent.ribo"]] <- PercentageFeatureSet(seur_obj, pattern = "^RP[LS]")
    } else { # it is mouse
        seur_obj[["percent.ribo"]] <- PercentageFeatureSet(seur_obj, pattern = "^Rp[ls]")
    }
    seur_obj[["MTRatio"]] <- seur_obj$percent.mt / 100 
    
    message('^MITO and ^RIBO levels processed successfully')
    # Calculate percentages of largest genes by single cell
    # Use the "b_additional_feat" property carefully because of grows the size object significantly.  
    # if (b_additional_feat==TRUE) {        
    #     seur_obj <- l_get_perc_largest_genes(seur_obj)
    #     message('Additional largest genes by cell features attached successfully')
    # }  
    
    message('Seurat completed successfully!')
    return(seur_obj) 
}


f_create_seurat_RAW  <- function(seuratName, s_bc_mtx) {
    # This function creates a Seurat object from RAW to process with empty::droplets for objects interaction
    #       @seuratName         string with the name to be assigned to the Seurat object
    #       @s_bc_mtx           string with the barcode file path

    # filtered bc mtx
    mtx <- Read10X_h5(s_bc_mtx)             # Returns a sparse matrix with rows and columns labeled
    rna_counts <- mtx$`Gene Expression`
    message('GEX counts loaded successfully')
    
    seur_obj <- CreateSeuratObject(
        counts = rna_counts,
        assay = "RNA",
        project = seuratName,
        #meta.data = meta
    )
    #UpdateSeuratObject(seur_obj)
    print(seur_obj)
    message('Seurat assay completed successfully')
    
    # Ensure seurat_obj is a Seurat object
    if (!("Seurat" %in% class(seur_obj))) { stop("Seurat object not found") } 
    if (!("counts" %in% Layers(seur_obj[["RNA"]]))) { stop("RNA counts layer not found!") } 
    if (!("data" %in% Layers(seur_obj[["RNA"]]))) { stop("Data counts not found!") }
    
    message('Process completed successfully. Seurat object saved as ', here("processed-data", seuratName))
    
    return(seur_obj) 
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
    if (s_tissue=='human') {
        ribo_genes = length(grep("^RP[LS]",rownames(seuratOBJ@assays$RNA@counts),value = TRUE))
    } else {    # it is mouse
        ribo_genes = length(grep("^Rp[ls]",rownames(seuratOBJ@assays$RNA@counts),value = TRUE))
    }
    
    # Build a table with the stats applied and the outputs gotten
    # Table with the cellranger-arc `gene expression statistics`
    #mtx_stats <- matrix(c(nCount, nFeature, nMT), ncol=15, byrow=TRUE)
    mtx_stats <- matrix(c(nCount, nFeature, nMT), ncol=17, byrow=TRUE)
    #class(mtx_stats)
    #dim(mtx_stats)
    colnames(mtx_stats) <- c('nCount0%','nCount25%','nCount50%','nCount75%','nCount100%',
                             'nGenes0%','nGenes25%','nGenes50%','nGenes75%','nGenes100%',
                             'nMito0%','nMito10p','nMito25%','nMito50%','nMito75%','nMito90p','nMito99%')

    # add additional columns: sample name, total cells, #mito genes and #ribo genes.
    
    #message('nCount, nFeature and nMT statistics done.')
    mtx_stats<-cbind(mtx_stats,mt_genes)    # number of mito genes
    mtx_stats<-cbind(mtx_stats,ribo_genes)  # number of ribo genes  
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

# from https://github.com/satijalab/seurat/issues/5343
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
    
    # print('Reductions:')
    # for (reductionName in names(seuratObj@reductions)) {
    #     print(paste0(reductionName, ', size: ', format(object.size(x = seuratObj@reductions[reductionName]), units = 'auto')))
    # }
    
    print('All slots:')
    for (slotName in slotNames(seuratObj)) {
        val <- object.size(slot(seuratObj, slotName))
        if (val > 1000) { 
            print(paste0(slotName, ': ', format(val, units = 'auto')))
        }
    }
    
    # if (nrow(seuratObj@meta.data) != nCells) {
    #     print(paste0('Metadata rows not equal to ncol: ', nCells, ' / ', nrow(seuratObj@meta.data)))
    # }
}
