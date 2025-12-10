## Authors. HT / CS
## Date. May 30th, 2023

# RNA analysis by SCTransform workflow
run_DimRed_scRNA <- function(seuratOBJ, umap_dims = 1:50, b_mito_confounders = FALSE, verbose = FALSE) {
    # First set the RNA assay active, then perform SCTransform normalization, followed by PCA and UMAP dimensionality reduction
    # INPUT:
    #   @seuratOBJ          Seurat obj
    #   @umap_dims          string describing the object (optional / provide short-name)
    #   @b_mito_confounders Boolean to remove confounding sources of variation, in this case mitochondrial mapping percentage
    #   @verbose            Boolean to output verbose in the screen
    # OUTPUT:
    #   @seuratOBJ          Seurat obj with dimensionality reduction attached
    
    # NOTE:
    # SCTransform() is a Seurat alternative to the NormalizeData, FindVariableFeatures and ScaleData workflow, 
    #               both ways use regularized negative binomial regression to normalize UMI counts when using VST in the FindVariableFeatures()

    # Ensure Seurat object is a Seurat object
    tryCatch( {
        if (!("Seurat" %in% class(SeuratOBJ))) stop("Seurat object not provided. Check the object exists!")
        message('Seurat with chromatin loaded successfully')
        #head(SeuratOBJ, n = 3)
    }, error = function(e) { message('An error ocurred. Verify you have assigned a Seurat object.') })
    
    tryCatch( {            
        # Set the default assay to RNA
        DefaultAssay(seuratOBJ) <- "RNA"

        if (verbose) { message('Running SCTransform on assay: RNA') }
        
        if (b_mito_confounders) {
            
            # Perform SCTransform normalization
            seuratOBJ <- SCTransform(seuratOBJ, 
                                     vars.to.regress = "percent.mt", verbose = verbose)

        } else {
            
            # Perform SCTransform normalization
            seuratOBJ <- SCTransform(seuratOBJ, verbose = verbose)

        }
        if (verbose) { message('SCTransform normalization done.') }
        
        # Run PCA
        seuratOBJ <- RunPCA(seuratOBJ)
        if (verbose) { message('PCA done.') }
        
        # Run UMAP dimensionality reduction
        #seuratOBJ <- Seurat::RunUMAP(seuratOBJ, umap.method = "umap-learn", metric = 'correlation',
        seuratOBJ <- RunUMAP(seuratOBJ, #method = "glmGamPoi",
                                     dims = umap_dims, reduction.name = 'umap.rna', 
                                     reduction.key = 'rnaUMAP_')
        ## Warning: Default method for RunUMAP has changed from calling Python UMAP via reticulate.... needs correction CSC
        if (verbose) { message('UMAP done.') }
        
    }, error = function(e) { message('An error ocurred. Not dimensionality reduction done!') })
            
    return(seuratOBJ)
    
}

# ATAC analysis. 
run_DimRed_scATAC <- function(seuratOBJ, umap_dims = 1:10, quantile = 'q0', min_dist = NULL, num_epocs = NULL, verbose = TRUE) {
  # This function uses LSI/LSA or Latent Semantic Indexing/Analysis for analyse ATAC data, some parametrization is permitted to improve the result as follow:
  #     @seuratOBJ is the Seurat object with ATAC assay included
  #     @umpa_dim is the dimensions to use as input features, used only if features is NULL. [5 to 50]
  #     @quantile, cutoff for feature to be included in the VariableFeatures for the object. By default is set to 'q0'   
  #     @min_dist, controls how tightly the embedding is allowed compress points together. Larger values ensure embedded points are moreevenly
  #                  distributed, while smaller allow the algorithm to optimize more accurately. Sensible values are in [0.001 to 0.5]
  #     @num_epocs, larger values result in more accurate embeddings. If NULL is specified, a value will be selected based on the size of 
  #                     the input dataset (200 for large datasets, 500 for small).    
  #     @verbose is to silence the verbose TRUE default mode when processing data
    
  # First is set the ATAC assay active
  DefaultAssay(seuratOBJ) <- "ATAC"
  
  # Perform TF-IDF (Term Frequency Inverse Document Frequency on a matrix normalization. Scale factor used as Default, equal to 10000.
  seuratOBJ <- RunTFIDF(seuratOBJ)
  
  # Find top features. The total counts and percentile rank for each feature is stored in the feature metadata for the assay.
  # To only compute the feature metadata, set min.cutoff=NA
  seuratOBJ <- FindTopFeatures(seuratOBJ, min.cutoff = quantile)  # q0 to q5; Ex. q5 is the top 95% most common features as the VariableFeatures
  
  # Perform SVD (Singular Value Decomposition) using irlba (Implicitly Restarted Lanczos Bidiagonalization Algorithm), in other words uses Fast Truncated SVD and PC Analysis for Sparse Matrices)
  seuratOBJ <- RunSVD(seuratOBJ)
  
  # Perform UMAP dimensionality reduction. We first exclude the first dimension as this is typically correlated with sequencing depth
  seuratOBJ <- RunUMAP(seuratOBJ, reduction = 'lsi', 
                       dims =umap_dims,         # [5 to 50]
                       #min.dist = min_dist,     # [0.001 to 0.5]
                       n.epochs = num_epocs,    # [200-500]                                       
                       reduction.name = "umap.atac", reduction.key = "atacUMAP_")
  
  # Find neighbors. Construct a shared NN graph by calculating the overlap (Jaccard index) between every cell and its k.param nearest neighbors.
  # By default uses PCA, here is called LSI. k.param defines k for the k-nearest neighbor algorithm, k.param = 20
  seuratOBJ <- FindNeighbors(object = seuratOBJ, reduction = 'lsi', dims = 2:umap_dims)
  
  # This is an alternative approach for weighted nearest neighbor graph
  # seuratOBJ <- FindMultiModalNeighbors(seuratOBJ, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
  # seuratOBJ <- RunUMAP(seuratOBJ, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
   
  # Find clusters. Calculate k-nearest neighbors and construct the SNN graph, then optimize the modularity function to determine clusters
  seuratOBJ <- FindClusters(object = seuratOBJ, verbose = verbose, algorithm = 3)
  # Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 
  #           3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires the leiden alg python.
  
  message('Process completed successfully. ')
  
  return(seuratOBJ)
}

#### run_DimRed_Integration function ####
run_DimRed_Integration <- function(seuratOBJ, resolution = 0.5, 
                                   multimodal_dims = 30, verbose = TRUE) {
  
  # Find Multi-Modal Neighbors
  seuratOBJ <- FindMultiModalNeighbors(seuratOBJ, reduction.list = list("pca", "lsi"), 
                                  dims.list = list(1:multimodal_dims, 2:multimodal_dims))
  
  # Perform UMAP dimensionality reduction
  seuratOBJ <- RunUMAP(seuratOBJ, nn.name = "weighted.nn", reduction.name = "wnn.umap", 
                  reduction.key = "wnnUMAP_")
  
  # Find clusters
  # Algorithm   1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 
  #             3 = SLM algorithm). SLM (Subspace Learning Machine)
  seuratOBJ <- FindClusters(seuratOBJ, graph.name = "wsnn", 
                       algorithm = 3, resolution = resolution, verbose = verbose)
  
  return(seuratOBJ)
}

#### plot_ATAC_UMAP ####
plot_dimension_reduction <- function(object, label = TRUE) {
  obj_name <- deparse(substitute(object))
  name_parts <- unlist(strsplit(obj_name, "\\."))[1]
  title <- paste("ATAC UMAP dimensionality reduction -", name_parts, sep = " ")
  
  plot <- DimPlot(object = object, label = label) + NoLegend() + ggtitle(title)
  return(print(plot))
}


############ Reproducibility information ####################

library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
#"2023-04-04 12:42:26 EDT"
proc.time()
options(width = 120)
session_info()
 
# [1] "Reproducibility information:"
# > # Last modification
#     > Sys.time()
# [1] "2023-06-14 11:25:32 EDT"
# > #"2023-04-04 12:42:26 EDT"
#     > proc.time()
# user  system elapsed 
# 0.196   0.067   4.963 