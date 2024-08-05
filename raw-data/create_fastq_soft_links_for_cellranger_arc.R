

# Require symbolic access to fasq raw data in the directory to process
# path_atac_July11_2024 -> /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2024-07-09-Psomagen/AN00019739_10X_RawData_Outs

# Track fasta multiome data from directories

library(stringr)
library(here)

# specified if any of this libraies is not provided 
b_scATAC_lib = TRUE
b_scRNA_lib = FALSE

here::here()
data_path <- here("/Hb_multiome/raw-data/")
system('ls -l')

# create directories for symbolic links
main_dir <- here("raw-data/")
fasta_dir <- "FASTQ_2024/"
path_atac_name <- "path_atac_July11_2024/"
path_snrna_name <- "path_snrna_2024" # check this related chunk when rna samples arrived. csc

## check if sub directories exists, otherwise create it 
if (!file.exists(file.path(main_dir, fasta_dir))){ dir.create(file.path(main_dir, fasta_dir)) }

# Main FASTQ container
main_dir <- here(paste0("raw-data/", fasta_dir))  # raw-data/FASTQ_2024/

# Sub-folders
if (!file.exists(file.path(main_dir, "GEX"))){ dir.create(file.path(main_dir, "GEX")) }
if (!file.exists(file.path(main_dir, "ATAC"))){ dir.create(file.path(main_dir, "ATAC")) }


######### Read the snRNAseq fastq files corresponding to the project. #########

if (b_scRNA_lib) {
  
    # Relative paths to create the soft links
    snRNAseq_path <- here('raw-data/path_snrna/')
    file_list <- list.files(path=snRNAseq_path, pattern=".*Hb.*\\.fastq\\.gz", full.names=FALSE, recursive = FALSE)
    # equivalent command in linux: $ ls path_snrna/ | grep .*Hb.*\\.fastq\\.gz
    file_list
    # > file_list
    # [1] "37---1C-Hb-KDM-Hb_S17_L001_R1_001.fastq.gz"
    # [2] "37---1C-Hb-KDM-Hb_S17_L001_R2_001.fastq.gz"
    # [3] "37---1C-Hb-KDM-Hb_S17_L002_R1_001.fastq.gz"
    # ...

    # GEX fast to rename
    ln_GEX_files <- length(file_list)
    
    # Read first folder in the directory
    if (!ln_GEX_files>0) { stop('GEX fastq files NOT found!') }
    
    # Read first fastq file, extract base-name (sDS_name) and create sub-dir
    if (ln_GEX_files>0) {
        ffn <- file_list[1]
        ln_file_name <- nchar(ffn)
        # extract base name
        sDS_name <- substring(ffn, 6, ln_file_name)
        # first base name and sub-folder name
        new_subdir <- substring(sDS_name, 1, 9)     #* Need to be defined, used temporal nomenclature, CSC
        # Main FASTQ container
        main_dir <- here("raw-data/FASTQ/GEX/")
        fasta_dir <- new_subdir  # 1C-Hb-KDM
        if (!file.exists(file.path(main_dir, fasta_dir))){ dir.create(file.path(main_dir, fasta_dir)) }
    } 
    
    #sDS_name
    # [1] "1C-Hb-KDM-Hb_S17_L001_R1_001.fastq.gz"
    
    # GEX fastq files
    # Parse the directory and create symbolic links for cellranger-arc pipeline
    here(snRNAseq_path) #"/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/raw-data/path_snrna/"
    
    for (f in file_list) {
        
        # testing
        # f <- "37---1C-Hb-KDM-Hb_S17_L001_R1_001.fastq.gz"
        # ...
      
        ln_file_name <- nchar(f)
        new_name <- substring(f, 6, ln_file_name)
        print(new_name)
        
        # read file, assign base_name and create subdirectory --> * Need to be re-defined, used temporal nomenclature, CSC
        if (!new_name==sDS_name) { 
            sDS_name <- substring(f, 6, ln_file_name) 
            new_subdir <- substring(sDS_name, 1, 9)     
            if (!new_subdir==fasta_dir) { 
                fasta_dir <- new_subdir 
                if (!file.exists(file.path(main_dir, fasta_dir))) { dir.create(file.path(main_dir, fasta_dir)) }  
            }
        } else {    
            sDS_name <- substring(f, 6, ln_file_name) 
            new_subdir <- substring(sDS_name, 1, 9)    
            #main_dir <- here("raw-data/FASTQ/GEX/")
            fasta_dir <- new_subdir
            if (!file.exists(file.path(main_dir, fasta_dir))) { dir.create(file.path(main_dir, fasta_dir)) }
        }
        
        # assign the symbolic link to this directory
        raw_path <- paste0(snRNAseq_path,'/')
        symbolic_args <- paste0('ln -s ', raw_path, f, ' ', 'raw-data/FASTQ/GEX/',fasta_dir, '/', sDS_name)
        print(symbolic_args)
        # ln -s /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/raw-data/path_snrna//37---1C-Hb-KDM-Hb_S17_L001_R1_001.fastq.gz raw-data/FASTQ/GEX/1C-Hb-KDM/1C-Hb-KDM-Hb_S17_L001_R1_001.fastq.gz
        system(symbolic_args)
            
    }
    
    message(ln_GEX_files,' files renamed with symbolic links for GEX assay')    
    system('ls raw-data/FASTQ/GEX/ -l')
}



# ATAC fasta files are arranged in folders
# Parse by directory and create symbolic links for cellranger-arc pipeline

if (b_scATAC_lib) {
    
    ## point to ATAC
    main_dir <- here(main_dir, "ATAC")
    snATACseq_path <- here("raw-data", path_atac_name)   
    dir_list <- list.dirs(path=snATACseq_path, full.names=FALSE, recursive = FALSE)
    dir_list
    # [1] "3A_Hb_KDM" "4A_Hb_KDM" "5A_Hb_KDM" "6A_Hb_KDM"
    
    ## Nmber of ATAC fasta files to rename
    ln_ATAC_files <- length(dir_list)
    
    # Read first directory
    if (!ln_ATAC_files>0) { stop('ATAC fastq files NOT found!') }
    fasta_subdir <- dir_list[1]
    if (!file.exists(file.path(main_dir, fasta_subdir))) { dir.create(file.path(main_dir, fasta_subdir)) }  
    
    ## Parse the atac files into the directory
    for (dd in dir_list) {
        
        # values for testing
        # dd <- "3A_Hb_KDM"
        # dd <- "6A_Hb_KDM"
        
        ## create subdirectoty if does not exist
        new_subdir <- dd
        # create sub-directory
        if (!new_subdir==fasta_subdir) { 
            fasta_subdir <- new_subdir 
            if (!file.exists(file.path(main_dir, fasta_subdir))) { dir.create(file.path(main_dir, fasta_subdir)) }  
        }
        
        ## get all fastq files found recursively in the given directory
        file_list <- list.files(path=paste0(snATACseq_path,'/',dd), pattern = "*.fastq.gz", full.names=FALSE, recursive = TRUE)
        file_list
        # [1] "2252KWLT4/3A_Hb_KDM_S3_L008_I1_001.fastq.gz"
        # [2] "2252KWLT4/3A_Hb_KDM_S3_L008_R1_001.fastq.gz"

        ## directory name assigned to retrieve the files / Usually the `main soft link` manually created to access the sequenced data
        raw_path <- paste0(snATACseq_path, dd, '/')
        # [1] "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/raw-data/path_atac//1A_Hb_KDM-1/"
        
        # create soft link to ATAC fasta files
        for (f in file_list) {
            # Data for testing / It should be accordingly with your own data
            # f <- "1A_Hb_KDM_S61_R1_001.fastq.gz"
            # f <- "1A_Hb_KDM_S62_I1_001.fastq.gz"
            
            # in this case we do not require rename files, so we only assign the symbolic link to the file
            source_file <-  here(raw_path, f)
            target_file <-  here("raw-data", fasta_dir, "ATAC", dd, basename(f))
            sys_command <- paste('ln -s ', source_file, target_file)
            system(sys_command)
            
        }
    
    }
}
    
message(ln_ATAC_files,' symbolic/soft links created for ATAC assay')    
#system('ls raw-data/FASTQ_2024/ATAC/ -l')
# drwxrws---+ 2 csoto lieber_lcolladotor 6 Jul 30 16:12 3A_Hb_KDM
# drwxrws---+ 2 csoto lieber_lcolladotor 6 Jul 30 16:12 4A_Hb_KDM
# drwxrws---+ 2 csoto lieber_lcolladotor 6 Jul 30 16:12 5A_Hb_KDM
# drwxrws---+ 2 csoto lieber_lcolladotor 6 Jul 30 16:12 6A_Hb_KDM






