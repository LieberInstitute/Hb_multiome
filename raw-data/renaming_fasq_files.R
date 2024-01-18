

# Require symbolic access to fasq raw data in the directory to process
# Ex. path_atac -> /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/24-01-02_SPag110823_ATAC/
# Ex. path_snrna -> /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/24-01-02_SPag110823/

# Tracking multiome data from directories
library(stringr)
library(here)

here::here()

data_path <- here("/Hb_multiome/raw-data/")
#here(data_path)
system('ls -l')

# create directories for symbolic links
main_dir <- here("raw-data/")
sub_dir <- "FASTQ"
# check if sub directories exists 
if (!file.exists(file.path(main_dir, sub_dir))){ dir.create(file.path(main_dir, sub_dir)) }
# Main FASTQ container
main_dir <- here("raw-data/FASTQ/")
# Sub-folders
sub_dir <- "GEX"
if (!file.exists(file.path(main_dir, sub_dir))){ dir.create(file.path(main_dir, sub_dir)) }
sub_dir <- "ATAC"
if (!file.exists(file.path(main_dir, sub_dir))){ dir.create(file.path(main_dir, sub_dir)) }


######### Read the snRNAseq fastq files corresponding to the project. #########

# Fastq relative paths to create the soft links
snRNAseq_path <- here('raw-data/path_snrna/')
file_list <- list.files(path=snRNAseq_path, pattern=".*Hb.*\\.fastq\\.gz", full.names=FALSE, recursive = FALSE)
# equivalent command in linux: $ ls path_snrna/ | grep .*Hb.*\\.fastq\\.gz
file_list
# > file_list
# [1] "37---1C-Hb-KDM-Hb_S17_L001_R1_001.fastq.gz"
# [2] "37---1C-Hb-KDM-Hb_S17_L001_R2_001.fastq.gz"
# [3] "37---1C-Hb-KDM-Hb_S17_L002_R1_001.fastq.gz"
# ...
# [9] "38---2C-Hb-KDM-Hb_S18_L001_R1_001.fastq.gz"
# [10] "38---2C-Hb-KDM-Hb_S18_L001_R2_001.fastq.gz"
# [11] "38---2C-Hb-KDM-Hb_S18_L002_R1_001.fastq.gz"
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
    sub_dir <- new_subdir  # 1C-Hb-KDM
    if (!file.exists(file.path(main_dir, sub_dir))){ dir.create(file.path(main_dir, sub_dir)) }
} 

#sDS_name
# [1] "1C-Hb-KDM-Hb_S17_L001_R1_001.fastq.gz"

# GEX fastq files
# Parse the directory and create symbolic links for cellranger-arc pipeline
here(snRNAseq_path) #"/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/raw-data/path_snrna/"

for (f in file_list) {
    
    # testing
    # f <- "37---1C-Hb-KDM-Hb_S17_L001_R1_001.fastq.gz"
    # f <- "37---1C-Hb-KDM-Hb_S17_L001_R2_001.fastq.gz"
    # f <- "38---2C-Hb-KDM-Hb_S18_L001_R1_001.fastq.gz"
    # f <- "38---2C-Hb-KDM-Hb_S18_L001_R2_001.fastq.gz"
    
    ln_file_name <- nchar(f)
    new_name <- substring(f, 6, ln_file_name)
    print(new_name)
    
    # read file, assign base_name and create subdirectory --> * Need to be re-defined, used temporal nomenclature, CSC
    if (!new_name==sDS_name) { 
        sDS_name <- substring(f, 6, ln_file_name) 
        new_subdir <- substring(sDS_name, 1, 9)     
        if (!new_subdir==sub_dir) { 
            sub_dir <- new_subdir 
            if (!file.exists(file.path(main_dir, sub_dir))) { dir.create(file.path(main_dir, sub_dir)) }  
        }
    } else {    
        sDS_name <- substring(f, 6, ln_file_name) 
        new_subdir <- substring(sDS_name, 1, 9)    
        #main_dir <- here("raw-data/FASTQ/GEX/")
        sub_dir <- new_subdir
        if (!file.exists(file.path(main_dir, sub_dir))) { dir.create(file.path(main_dir, sub_dir)) }
    }
    
    # assign the symbolic link to this directory
    raw_path <- paste0(snRNAseq_path,'/')
    symbolic_args <- paste0('ln -s ', raw_path, f, ' ', 'raw-data/FASTQ/GEX/',sub_dir, '/', sDS_name)
    print(symbolic_args)
    # ln -s /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/raw-data/path_snrna//37---1C-Hb-KDM-Hb_S17_L001_R1_001.fastq.gz raw-data/FASTQ/GEX/1C-Hb-KDM/1C-Hb-KDM-Hb_S17_L001_R1_001.fastq.gz
    system(symbolic_args)
        
}

message(ln_GEX_files,' files renamed with symbolic links for GEX assay')    
system('ls raw-data/FASTQ/GEX/ -l')




# ATAC Fast files are arrenged in folders
# Parse by directory and create symbolic links for cellranger-arc pipeline

main_dir <- here("raw-data/FASTQ/ATAC/")
snATACseq_path <- here('raw-data/path_atac/')   #'/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/24-01-02_SPag110823_ATAC/'
#here(snATACseq_path) # "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/raw-data/path_atac/"

dir_list <- list.files(path=snATACseq_path, full.names=FALSE, recursive = FALSE)
dir_list
# [1] "1A_Hb_KDM-1" "1A_Hb_KDM-2" "1A_Hb_KDM-3" "1A_Hb_KDM-4" "2A_Hb_KDM-1"
# [6] "2A_Hb_KDM-2" "2A_Hb_KDM-3" "2A_Hb_KDM-4"

# ATAC fasta files to rename
ln_ATAC_files <- length(file_list)

# Read first folder in the directory
if (!ln_ATAC_files>0) { stop('ATAC fastq files NOT found!') }
sub_dir <- dir_list[1]
if (!file.exists(file.path(main_dir, sub_dir))) { dir.create(file.path(main_dir, sub_dir)) }  

# Parse files into each subfolder
for (dd in dir_list) {
    # values for testing
    # dd <- "1A_Hb_KDM-1"
    # dd <- "1A_Hb_KDM-2"
    
    new_subdir <- dd
    
    # create sub-directory
    if (!new_subdir==sub_dir) { 
        sub_dir <- new_subdir 
        if (!file.exists(file.path(main_dir, sub_dir))) { dir.create(file.path(main_dir, sub_dir)) }  
    }
    
    # list of files contained in the given directory
    file_list <- list.files(path=paste0(snATACseq_path,'/',dd), full.names=FALSE, recursive = FALSE)
    file_list
    # [1] "1A_Hb_KDM_S61_I1_001.fastq.gz" "1A_Hb_KDM_S61_R1_001.fastq.gz"
    # [3] "1A_Hb_KDM_S61_R2_001.fastq.gz" "1A_Hb_KDM_S61_R3_001.fastq.gz"
    
    raw_path <- paste0(snATACseq_path,'/', dd, '/')
    # [1] "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/raw-data/path_atac//1A_Hb_KDM-1/"
    
    # rename ATAC fasta files
    for (f in file_list) {
        # f <- "1A_Hb_KDM_S61_I1_001.fastq.gz"
        # f <- "1A_Hb_KDM_S61_R1_001.fastq.gz"
        # f <- "1A_Hb_KDM_S62_I1_001.fastq.gz"
        
        # in this case we do not require rename files, so we only assign the symbolic link to the file
        symbolic_args <- paste0('ln -s ', raw_path, f, ' ', 'raw-data/FASTQ/ATAC/', dd, '/', f)
        print(symbolic_args)
        # ln -s /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/raw-data/path_atac//1A_Hb_KDM_S61_I1_001.fastq.gz raw-data/FASTQ/ATAC/1A_Hb_KDM_S61_I1_001.fastq.gz
        system(symbolic_args)
        
    }

}

message(ln_ATAC_files,' symbolic links created for ATAC assay')    
#system('ls raw-data/FASTQ/ATAC/ -l')






