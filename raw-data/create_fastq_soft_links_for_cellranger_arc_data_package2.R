
## GEX and/or ATAC Fasta files are arranged in folders
##    Parse each sub-directory and create symbolic links for cellranger-arc libraries

## For Habenula multiome we have these soft links to access the rawData:
##    path_atac_July11_2024 -> /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2024-07-09-Psomagen/AN00019739_10X_RawData_Outs
##    path_rna_Aug05_2024 -> /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2024-08-02_Psomagen/AN00020424_10X_RawData_Outs

library(stringr)
library(here)

# specified if any of this libraies is not provided 
b_scATAC_lib = TRUE
b_scRNA_lib = TRUE

here::here()
## softlinks pointing to fastq files in disk 
path_atac_name <- "path_atac_July11_2024/"
path_rna_name <- "path_rna_Aug05_2024/" 

# Main FASTQ container
main_dir <- here("raw-data/FASTQ_2024/")  # raw-data/FASTQ_2024/
if (!file.exists(file.path(main_dir))){ dir.create(file.path(main_dir)) }


######### This chunk is to create the GEX Fastq files soft links #########

if (b_scRNA_lib) {
    
    ## Create sub directories for symbolic links
    if (!file.exists(file.path(main_dir, "GEX"))){ dir.create(file.path(main_dir, "GEX")) }
  
    ## point to GEX
    main_subdir <- here(main_dir, "GEX")  
    ## soft-link to fastq raw data
    snRNAseq_path <- here("raw-data", path_rna_name)
    dir_list <- list.files(path=snRNAseq_path, pattern="*C", full.names=FALSE, recursive = FALSE)
    dir_list
    # [1] "4C_Hb_KDM" "5C_Hb_KDM" "6C_Hb_KDM" "7C_Hb_KDM" "8C_Hb_KDM"
    
    ## Validate if dir(s) available and create sub-directories inside GEX/ subdir
    ln_RNA_files <- length(dir_list)
    if (!ln_RNA_files>0) { stop('RNA fastq files NOT found!') }
    fasta_subdir <- dir_list[1]
    if (!file.exists(file.path(main_subdir, fasta_subdir))) { dir.create(file.path(main_subdir, fasta_subdir)) }  
      
    ## Parse the GEX files into the sub-directory
    for (new_subdir in dir_list) {
      
      # values for testing
      # new_subdir <- "4C_Hb_KDM"
      
      # create sub-directory inside GEX/
      if (!new_subdir==fasta_subdir) { 
        fasta_subdir <- new_subdir 
        if (!file.exists(file.path(main_subdir, fasta_subdir))) { dir.create(file.path(main_subdir, fasta_subdir)) }  
      }
      
      ## get all Fastq files found recursively in the given directory
      file_list <- list.files(path=paste0(snRNAseq_path, '/', new_subdir), pattern = "*.fastq.gz", full.names=FALSE, recursive = TRUE)
      file_list
      # [1] "22F3JLLT4/4C_Hb_KDM_S4_L007_I1_001.fastq.gz"
      # [2] "22F3JLLT4/4C_Hb_KDM_S4_L007_I2_001.fastq.gz"
      # ...
      
      raw_path <- paste0(snRNAseq_path, new_subdir, '/')
      # "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/raw-data/path_rna_Aug05_2024/4C_Hb_KDM/"
      
      # create soft links to each GEX Fastq file
      for (f in file_list) {
        # Data for testing:
        # f <- "22F3JLLT4/4C_Hb_KDM_S4_L007_I1_001.fastq.gz"
        
        # in this case we do not require rename files, so we only assign the symbolic link to the file
        source_file <-  here(raw_path, f)
        target_file <-  here(main_subdir, new_subdir, basename(f))
        sys_command <- paste('ln -s ', source_file, target_file)
        system(sys_command)
        
      }
      
    }
    
}

message(' Tree of softlinks created in ', main_subdir)
system(paste0('tree -L 2 ', main_subdir))


######### This chunk is to create the ATAC Fastq files soft links #########

if (b_scATAC_lib) {
  
    ## Create sub directories for symbolic links
    if (!file.exists(file.path(main_dir, "ATAC"))){ dir.create(file.path(main_dir, "ATAC")) }
    
    ## point to ATAC
    main_subdir <- here(main_dir, "ATAC")
    snATACseq_path <- here("raw-data", path_atac_name)   
    dir_list <- list.dirs(path=snATACseq_path, full.names=FALSE, recursive = FALSE)
    dir_list
    # [1] "3A_Hb_KDM" "4A_Hb_KDM" "5A_Hb_KDM" "6A_Hb_KDM"
    
    ## Validate if dir(s) available and create sub-directories inside GEX/ subdir
    ln_ATAC_files <- length(dir_list)
    if (!ln_ATAC_files>0) { stop('ATAC fastq files NOT found!') }
    fasta_subdir <- dir_list[1]
    if (!file.exists(file.path(main_subdir, fasta_subdir))) { dir.create(file.path(main_subdir, fasta_subdir)) }  
    
    ## Parse the ATAC files into the sub directory
    for (new_subdir in dir_list) {
        
        # values for testing
        # new_subdir <- "3A_Hb_KDM"

        # create sub-directory inside ATAC/
        if (!new_subdir==fasta_subdir) { 
            fasta_subdir <- new_subdir 
            if (!file.exists(file.path(main_subdir, fasta_subdir))) { dir.create(file.path(main_subdir, fasta_subdir)) }  
        }
        
        ## get all fastq files found recursively in the given directory
        file_list <- list.files(path=paste0(snATACseq_path,'/',new_subdir), pattern = "*.fastq.gz", full.names=FALSE, recursive = TRUE)
        file_list
        # [1] "2252KWLT4/3A_Hb_KDM_S3_L008_I1_001.fastq.gz"
        # [2] "2252KWLT4/3A_Hb_KDM_S3_L008_R1_001.fastq.gz"
        # ...

        ## directory name assigned to retrieve the files / Usually the `main soft link` manually created to access the sequenced data
        raw_path <- paste0(snATACseq_path, new_subdir, '/')
        # [1] "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/raw-data/path_atac//1A_Hb_KDM-1/"
        
        # create soft link to each ATAC Fastq files
        for (f in file_list) {
            # Data for testing:
            # f <- "1A_Hb_KDM_S61_R1_001.fastq.gz"
            # f <- "1A_Hb_KDM_S62_I1_001.fastq.gz"
            
            # in this case we do not require rename files, so we only assign the symbolic link to the file
            source_file <-  here(raw_path, f)
            target_file <-  here(main_subdir, new_subdir, basename(f))
            sys_command <- paste('ln -s ', source_file, target_file)
            system(sys_command)
            
        }
    
    }
}
    
message(' Tree of softlinks created in ', main_subdir)
system(paste0('tree -L 2 ', main_subdir))

