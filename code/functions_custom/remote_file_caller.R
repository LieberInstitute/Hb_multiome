########################################################################
##
## FUNCTIONS TO HANDLE SEURAT OBJECTS 
##
## Authors. CSC
## Date. Aug 7th, 2023
## Files caller. 
## Goal: find and load a file for a particular use. Ex. Load a raw dataset, load a metadata file, etc. 
## 
## Functions added:
##      get_raw_barcode_mtx 
##      get_filtered_barcode_mtx
##      get_job_array_targets
########################################################################


get_raw_barcode_mtx <- function(s_experiment_name) {
    # This function find a barcode file and return the full path name 
    #       @s_experiment_name         experiment name to be used to look for the raw barcode mtx

    # experiment name must to exist
    # if(!(s_experiment_name %in% lst_experiments)) stop("Experiment name does not match any declared experiment.")
    print(paste('Recovering path for raw barcode matrix ',s_experiment_name))
    
    # Assign the file paths to load the raw barcode mtx
    
    
    # These are the new Habenula samples (LIEBER 2024)
    if (s_experiment_name == 'S1_Hb_KDM' || s_experiment_name == 'S2_Hb_KDM') {
        base_path <- '/users/csoto/Hb_multiome/processed-data/cellrangerARC'
        if ( s_experiment_name == 'S1_Hb_KDM') {
            s_featured_bc_mtx <- here(paste0(base_path,"/S1_Hb_KDM/outs"), "raw_feature_bc_matrix.h5") }
        else {  # S2_Hb_KDM
            s_featured_bc_mtx <- here(paste0(base_path,"/S2_Hb_KDM/outs"), "raw_feature_bc_matrix.h5") }
    }
    
    # These are the public PBMC datasets 
    if (s_experiment_name == 'pbmc3k' || s_experiment_name == 'pbmc10k') {
        base_path <- '/users/csoto/cellranger-arc-public/raw-data'   
        if ( s_experiment_name == 'pbmc3k') {
            s_featured_bc_mtx <- here(paste0(base_path,"/PBMC_CellSorted_ARC2_0_0"), "pbmc_granulocyte_sorted_3k_raw_feature_bc_matrix.h5") }
        else {
            s_featured_bc_mtx <- here(paste0(base_path,"/PBMC_CellSorted_ARC2_0_0"), "pbmc_granulocyte_sorted_10k_raw_feature_bc_matrix.h5") }
    }

    # These are the samples 42_* (2020)
    if ( s_experiment_name == 'hippo42_1' || s_experiment_name == 'hippo42_4' ) {
        base_path <- '/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/rafael_rotation/cellranger_rerun'
        if ( s_experiment_name == 'hippo42_1') {
            s_featured_bc_mtx = here(paste0(base_path,"/42_1/outs"), "raw_feature_bc_matrix.h5") }
        else {
            s_featured_bc_mtx = here(paste0(base_path,"/42_4/outs"), "raw_feature_bc_matrix.h5") }
    }    
    
    # These are the new HPC samples (2023)
    if ( s_experiment_name == '3_HPC_KDM' || s_experiment_name == '2_HPC_KDM' || s_experiment_name == '1_HPC_KDM' ) {
        base_path <- '/users/csoto/HPC_multiome_pilot/processed-data/cellranger_run_fast_version/'
        if ( s_experiment_name == '3_HPC_KDM') {
            s_featured_bc_mtx <- here(paste0(base_path,"/3_HPC_KDM/outs"), "raw_feature_bc_matrix.h5") }
        else if ( s_experiment_name == '2_HPC_KDM' ) {
            s_featured_bc_mtx <- here(paste0(base_path,"/2_HPC_KDM/outs"), "raw_feature_bc_matrix.h5") }
        else if ( s_experiment_name == '1_HPC_KDM' ) {
            s_featured_bc_mtx <- here(paste0(base_path,"/1_HPC_KDM/outs"), "raw_feature_bc_matrix.h5") }
    }    
    
    
    #message('Process completed successfully')
    #print(paste('Reading',s_featured_bc_mtx))
    return(s_featured_bc_mtx)
    
}



get_filtered_barcode_mtx <- function(s_experiment_name) {
    # Read the filtered barcode file and return the full path name 
    #       @s_experiment_name         experiment name to be used to look for the raw barcode mtx
    
    # experiment name must to exist
   #if(!(s_experiment_name %in% lst_experiments)) stop("Experiment name does not match any declared experiment.")

    # These are the new Habenula samples (LIEBER 2024)
    if (s_experiment_name == 'S1_Hb_KDM' || s_experiment_name == 'S2_Hb_KDM') {
        
        base_path <- 'processed-data/cellrangerARC/'

        
    } else {

        # These are the new HPC samples (LIEBER 2023)
        base_path <- paste0('processed-data/cellranger_run_fast_version/')

    }
    
    base_path <- paste0(base_path, s_experiment_name,"/outs")
    s_featured_bc_mtx <- here(base_path, "filtered_feature_bc_matrix.h5" )
    message('Filtered barcode matrix H5 found: ', s_featured_bc_mtx)

    return(s_featured_bc_mtx)
    
}

# new function. CSC.01.2024
get_ATAC_barcode_tsv <- function(s_experiment_name) {
    # Read fragment metadata and return the full path name
    #       @s_experiment_name         experiment name to be used
    # NOTE. Since v2.0.0 to v2.0.2 is used the same name standard. Be careful in future cellranger-arc releases. (Sep, 2023)
    
    if (s_experiment_name == 'S1_Hb_KDM' || s_experiment_name == 'S2_Hb_KDM') {
        
        base_path <- 'processed-data/cellrangerARC/'
    
    } else {
        
        # These are the new HPC samples (LIEBER 2023)
        base_path <- paste0('processed-data/cellranger_run_fast_version/')
        
    }    

    base_path <- paste0(base_path, s_experiment_name,"/outs")
    s_atac_path <- here(base_path,"atac_fragments.tsv.gz")
    message('ATAC fragments (TSV) file found: ', s_atac_path)

    return(s_atac_path)
    
}



# Function name changed: get_ATAC_metadata_paths by get_ATAC_barcode_paths
get_ATAC_barcode_paths <- function(s_experiment_name) {
    # This function call a fragment metadata and return the full path name
    #       @s_experiment_name         experiment name to be used
    # Return a list with 2 paths, one for the filtered barcode mtx (H5) and other for the fragments (TSV)  
    
    atac_paths <- list()

    # Assign the file paths to load the filtered barcode mtx
    # NOTE. Since v2.0.0 to v2.0.2 is used the same name standard. Be careful in future cellranger-arc releases. (Sep, 2023)
    
    base_path <- paste0("processed-data/cellranger_run_fast_version/", s_experiment_name,"/outs")

    # Read barcode filtered mtx path (*.h5)
    atac_paths[1] <- here(base_path,"filtered_feature_bc_matrix.h5")
    # Read fragments (*.tsv.gz) paths
    atac_paths[2] <- here(base_path,"atac_fragments.tsv.gz")
    message('ATAC fragments (TSV) and barcode (h5) path files found')

    return(atac_paths)

}

get_metadata_path <- function(s_experiment_name) {
    # This function build metadata obj with specific features
    #       @s_experiment_name         experiment name reference
    
    # NOTE. Since v2.0.0 to v2.0.2 is used the same name standard. Be careful in future cellranger-arc releases. (Sep, 2023)
    # CSC 09.19.2023
    
    # Read meta-data csv file
    # These are the new Habenula samples (LIEBER 2024)
    if (s_experiment_name == 'S1_Hb_KDM' || s_experiment_name == 'S2_Hb_KDM') {
        
        base_path <- paste0('processed-data/cellrangerARC/')

    } else {
        
        # These are the new HPC samples (LIEBER 2023)
        base_path <- paste0("processed-data/cellranger_run_fast_version/")

    }    
    
    base_path <- paste0(base_path, s_experiment_name,"/outs")
    meta_path_file <- here(base_path, "per_barcode_metrics.csv")
    message('Meta-data (csv) file found: ', meta_path_file)
    return(meta_path_file)
    
}

get_job_array_targets <- function(lst_experiments) {
    # This function create a txt file with the full paths to the raw barcode matrices from cellranger-arc data
    #       @lst_experiments: list of experiment names
    
    i = 1
    # Recover the full file paths
    for (s_experiment_name in lst_experiments) {
        
        message("Reading raw barcode path for sample ", s_experiment_name)
        # main function with the paths for all the experiments
        sample_path <- get_raw_barcode_mtx(s_experiment_name)    
        message("Raw barcode path: ", sample_path)
        
        # Build the array target file
        if (i == 1) {
            # Table with the paths to raw-data
            tab_paths <- matrix(c(sample_path), ncol=1, byrow=TRUE)
            #colnames(tab_paths) <- c('File_path')
            tab_paths <- as.table(tab_paths)
            i <- i+1
        } else {
            # Add Row to table with the paths to raw data using rbind()
            new_row = c(sample_path)
            tab_paths = rbind(tab_paths,new_row)
        }
        
    }

    rownames(tab_paths) <- NULL
    #message(tab_paths)
    
    # # Export the table in txt format
    # s_file_name <- here('code/', 'array_targets.txt')
    # write.table(tab_paths, s_file_name, append = FALSE, 
    #             row.names = FALSE, col.names = FALSE, quote = FALSE)
    return(tab_paths)
    
}

