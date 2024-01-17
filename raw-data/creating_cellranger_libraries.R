

# Tracking multiome data from directories
library(stringr)
library(here)

here::here()

data_path <- here("/Hb_multiome/raw-data/FASTQ/GEX")

# Read the GEX fastq files (symbolic links)

# Fastq relative paths to create the soft links
snRNAseq_path <- here("raw-data/FASTQ/GEX/")
file_list <- list.files(path=snRNAseq_path, full.names=FALSE, recursive = FALSE)
file_list
# [1] "1C-Hb-KDM-Hb_S17_L001_R1_001.fastq.gz"
# [2] "1C-Hb-KDM-Hb_S17_L001_R2_001.fastq.gz"
# [3] "1C-Hb-KDM-Hb_S17_L002_R1_001.fastq.gz"
# [4] "1C-Hb-KDM-Hb_S17_L002_R2_001.fastq.gz"
# [5] "1C-Hb-KDM-Hb_S17_L003_R1_001.fastq.gz"
# [6] "1C-Hb-KDM-Hb_S17_L003_R2_001.fastq.gz"
# [7] "1C-Hb-KDM-Hb_S17_L004_R1_001.fastq.gz"
# [8] "1C-Hb-KDM-Hb_S17_L004_R2_001.fastq.gz"
# ...

# GEX fastq(s) list
ln_GEX_files <- length(file_list)
message(ln_GEX_files, ' GEX files to add')

# cellranger library layout
# fastqs,sample,library_type
# /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/HPC_multiome_pilot/raw-data/FASTQ/GEX/1HPC_S1_L003, 1G_HPC_KDM, Gene Expression
# /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/HPC_multiome_pilot/raw-data/FASTQ/ATAC/1HPC_S1_L003, 1A_HPC_KDM, Chromatin Accessibility
# /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/HPC_multiome_pilot/raw-data/FASTQ/ATAC/1HPC_S2_L003, 1A_HPC_KDM, Chromatin Accessibility
# /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/HPC_multiome_pilot/raw-data/FASTQ/ATAC/1HPC_S3_L003, 1A_HPC_KDM, Chromatin Accessibility
# /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/HPC_multiome_pilot/raw-data/FASTQ/ATAC/1HPC_S4_L003, 1A_HPC_KDM, Chromatin Accessibility

## Library name csv file
s_sample <- 'libraries_S1_2024_habenula'

# extract base name
f <- file_list[1]
# [1] "1A_Hb_KDM_S61_I1_001.fastq.gz"
base_name <- substring(f, 1, 9)
print(base_name)

message('Integrating ... wait ')

cellranger_data <- data_path

i <- 1
for (filex in file_list) {
    # filex <- "1A_Hb_KDM_S61_I1_001.fastq.gz"
    # Build the tables
    if (i == 1) {
        # create header
        df_all <- 'fastqs,sample,library_type'
        i <- i+1
    } else {
        # Append rows
        new_row = snRNAseq_path
        df_all = rbind(df_all,new_row)
        df_all
    }
}




# Read the ATAC fastq files (symbolic links)

# Fastq relative paths to create the soft links
snATACseq_path <- here("raw-data/FASTQ/ATAC/")
file_list <- list.files(path=snATACseq_path, full.names=FALSE, recursive = FALSE)
file_list
# [1] "1A_Hb_KDM_S61_I1_001.fastq.gz" "1A_Hb_KDM_S61_R1_001.fastq.gz"
# [3] "1A_Hb_KDM_S61_R2_001.fastq.gz" "1A_Hb_KDM_S61_R3_001.fastq.gz"
# [5] "1A_Hb_KDM_S62_I1_001.fastq.gz" "1A_Hb_KDM_S62_R1_001.fastq.gz"
# [7] "1A_Hb_KDM_S62_R2_001.fastq.gz" "1A_Hb_KDM_S62_R3_001.fastq.gz"
# [9] "1A_Hb_KDM_S63_I1_001.fastq.gz" "1A_Hb_KDM_S63_R1_001.fastq.gz"
# [11] "1A_Hb_KDM_S63_R2_001.fastq.gz" "1A_Hb_KDM_S63_R3_001.fastq.gz"
# [13] "1A_Hb_KDM_S64_I1_001.fastq.gz" "1A_Hb_KDM_S64_R1_001.fastq.gz"
# [15] "1A_Hb_KDM_S64_R2_001.fastq.gz" "1A_Hb_KDM_S64_R3_001.fastq.gz"
# ...

# GEX fastq(s) list
ln_ATAC_files <- length(file_list)

message(ln_ATAC_files, ' ATAC files to add')
