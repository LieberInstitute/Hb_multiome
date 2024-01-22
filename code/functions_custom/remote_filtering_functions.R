########################################################################
##
## FUNCTIONS TO HANDLE SEURAT OBJECTS 
##
## Authors. HT + CSC
## Date. May 4th, 2023
## Seurat filtering functions. HT 19/05/2023
## Pre-processed/subset functions. CS 08/06/2023
########################################################################

get_preprocessed_subset <- function(seurat_obj, nCountRNA_h, nCountRNA_l, mito_p, nCount_ATAC_h, nCount_ATAC_l, nucleosome_grp, tss_s) {
    # This function subsets a Seurat object with chromatin data attached 
    #       @seurat_obj         string with the name to be assigned to the Seurat object
    #       @nCountRNA_h        higher number of genes per cell to kept
    #       @nCountRNA_l        lower number of genes per cell to kept
    #       @mito_p             higher mitochondrial percentage of genes per cell to kept
    #       @nCount_ATAC_h      higher ATAC counts per cell to kept
    #       @nCount_ATAC_l      lower ATAC counts per cell to kept
    #       @nucleosome_grp     higher nucleosome group to kept, derived from the nucleosome signal pattern
    #       @tss_s              higher TSS (Transcription Start Sites) group to kept, derived from the TSS enrichment score
    
    # Ensure seurat_obj is a Seurat object
    if(!("Seurat" %in% class(seurat_obj))) stop("Seurat object should be a Seurat object")
    
    # Check if required slots exist
    required_slots <- c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'nucleosome_group', 'high.tss')
    if(any(!required_slots %in% names(seurat_obj@meta.data))) 
        stop("Seurat object should have the slots: nCount_RNA, nFeature_RNA, and percent.mt, nucleosome_group, high.tss")
    
    # Subset Seurat obj with chromatin data based on the given thresholds 
    # i.e: nCount_ATAC < 7e4 & nCount_ATAC > 5e3 & nCount_RNA < 25000 & nCount_RNA > 1000 & percent.mt < 20
    seurat_obj.subset <- subset(x = seurat_obj, subset = (nCount_RNA < nCountRNA_h & nCount_RNA > nCountRNA_l & percent.mt < mito_p) &
                 (nCount_ATAC < nCount_ATAC_h & nCount_ATAC > nCount_ATAC_l) &
                 (nucleosome_group < nucleosome_grp & TSS.enrichment > tss_s)
    )

    message('Process completed successfully')
    
    # Return the filtered Seurat object
    return(seurat_obj.subset)
    
    }

# Check specific slots: GetAssayData(object = pbmc3k[["RNA"]], slot = "data")[1:5,1:5]
# Subseting seurat object to smaller one as a matter to efficiently deal with the memory resources and implement processes quicker 
# seurat10k_subset <- subset(pbmc10k, cells = colnames(pbmc10k)[1:100])
# print(ncol(seurat10k_subset))
# head(pbmc10k@meta.data, n=3)
# colnames(pbmc3k[[]])

#### Method 1 ####
# GEX method 1 (Probabilities)
get_seurat_GEX_filteringM1p <- function(seurat_obj, lowRNA_prob = 0.01, highRNA_prob = 0.99, lowFeature_prob = 0.01, 
                                      highMT_prob = 0.90, save_combined=FALSE) {
  # Ensure seurat_obj is a Seurat object
  if(!("Seurat" %in% class(seurat_obj))) stop("seurat_obj should be a Seurat object")
  
  # Check if required slots exist
  required_slots <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
  if(any(!required_slots %in% names(seurat_obj@meta.data))) 
      stop("seurat_obj should have the slots nCount_RNA, nFeature_RNA, and percent.mt")
  
  # Get the name of seurat_obj
  seurat_name <- deparse(substitute(seurat_obj))
  
  # Calculate the thresholds
  countLOW = quantile(seurat_obj$nCount_RNA, probs=lowRNA_prob)
  countHIGH = quantile(seurat_obj$nCount_RNA, probs=highRNA_prob)
  featureLOW = quantile(seurat_obj$nFeature_RNA, probs=lowFeature_prob)
  mitHIGH = round(quantile(seurat_obj$percent.mt, probs=highMT_prob))
  
  # Print the thresholds
  message(paste("Seurat Object: ", seurat_name, 
                "\nThresholds: ", 
                "\ncountLOW_UMIs: ", countLOW,"%",
                "\ncountHIGH_UMIs: ", countHIGH, "%",
                "\nfeatureLOW_genes: ", featureLOW, "%",
                "\nmitoHIGH: ", mitHIGH,"%"))
  
  # Subset the Seurat object based on the thresholds
  seuratOBJ.filtered.M1 <- subset(x = seurat_obj, subset = (nFeature_RNA >= featureLOW) &
                                    (nCount_RNA >= countLOW)  &
                                    (nCount_RNA < countHIGH) &
                                    (percent.mt < mitHIGH))
  
    if (save_combined) {
        # Check if processed_data directory exists, if not create it
        if (!dir.exists(here("processed-data"))) {
        dir.create(here("processed-data"))
        }
        
        # Define the file name based on combined parameter
        file_name <- ifelse(save_combined, paste0(seurat_name, ".combined.filtered.GEX.M1.Prob.RDS"), 
                          paste0(seurat_name, ".filtered.GEX.M1.Prob.RDS"))
        
        # Save the filtered Seurat object using here package
        write_rds(seuratOBJ.filtered.M1, here("processed-data", file_name), compress = ('gz'))
        message('Filtered object saved as ', here("processed-data", file_name))
    }      
    
    message('Process completed successfully')
    # Return the filtered Seurat object
    return(seuratOBJ.filtered.M1)
}

#### Method 2 ####
# GEX method 2 (Standard Deviation)
get_seurat_GEX_filteringM2sd <- function(seurat_obj, iSD=2, save_combined=FALSE){
  
    # Ensure seurat_obj is a Seurat object
    if(!("Seurat" %in% class(seurat_obj))) stop("seurat_obj should be a Seurat object")
    
    # Check if required slots exist
    required_slots <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
    if(any(!required_slots %in% names(seurat_obj@meta.data))) stop("seurat_obj should have the slots nCount_RNA, nFeature_RNA, and percent.mt")
    
    # Get the name of seurat_obj
    seurat_name <- deparse(substitute(seurat_obj))
    
    # Calculate count, feature, and mt thresholds
    count.max <- round(mean(seurat_obj$nCount_RNA) + iSD * sd(seurat_obj$nCount_RNA), digits = -2)
    count.min <- max(round(mean(seurat_obj$nCount_RNA) - iSD * sd(seurat_obj$nCount_RNA), digits = -2), 0)
    feat.max <- round(mean(seurat_obj$nFeature_RNA) + iSD * sd(seurat_obj$nFeature_RNA), digits = -2)
    feat.min <- max(round(mean(seurat_obj$nFeature_RNA) - iSD * sd(seurat_obj$nFeature_RNA), digits = -2), 0)
    mt.min <- round(mean(seurat_obj$percent.mt) - iSD * sd(seurat_obj$percent.mt))
    mt.max <- round(mean(seurat_obj$percent.mt) + iSD * sd(seurat_obj$percent.mt))
    
    # Print the thresholds
    message(paste("Seurat Object: ", seurat_name, 
                "\nThresholds: ", 
                "\ncount_min_UMIs: ", count.min,
                "\ncount_max_UMIs: ", count.max, 
                "\nfeature_min_genes: ", feat.min, 
                "\nfeature_max_genes: ", feat.max,
                "\nmt_min: ", mt.min,
                "\nmt_max: ", mt.max))
  
    # Subset the Seurat object based on the thresholds
    seurat_obj.filtered <- subset(x = seurat_obj, subset = (
            ((nFeature_RNA > feat.min) & (nFeature_RNA < feat.max)) &
            ((nCount_RNA < count.max) & (nCount_RNA > count.min)) & 
            (percent.mt < mt.max))
            )
  
    if (save_combined) {
        # Check if processed_data directory exists, if not create it
        if (!dir.exists(here("processed-data"))) {
            dir.create(here("processed-data"))
        }
      
        # Define the file name based on combined parameter
        file_name <- ifelse(save_combined, paste0(seurat_name, ".combined.filtered.GEX.M2.SD.RData"), 
                          paste0(seurat_name, ".filtered.GEX.M2.SD.RData"))
        
        # Save the filtered Seurat object using here package
        save(seurat_obj.filtered, file = here("processed-data", file_name))
        message('Filtered object saved as ', here("processed-data", file_name))
    }

    message('Process completed successfully')
    # Return the filtered Seurat object
    return(seurat_obj.filtered)
}
