# Hb_multiome

Human Habenula Multiome Project

### Overview

#### Local project name: Hb_multiome

Allocated in: `/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome`

#### Fastq GEX and ATAC Raw-data

First sequencing for `Human Habenula Pilot`:

snRNA-Seq fastq libraries in: `/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/24-01-02_SPag110823/`

snATAC-Seq fastq libraries in: `/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/24-01-02_SPag110823_ATAC/`

<br> Second sequencing for `Human Habenula Neurotypical Controls` (August 2024):

snRNA-Seq fastq libraries in: `/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2024-07-09-Psomagen/AN00019739_10X_RawData_Outs`

snATAC-Seq fastq libraries in: `/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2024-08-02_Psomagen/AN00020424_10X_RawData_Outs`

Notes: - Second sequencing have `Visium` sequencing files too. - `RNA multiome files` have a `C` letter; `Visum files` have a `v` lower case letter.

<br>

````{=html}
<!--

#### Soft links into raw-data/ subfolder:

These links are used to track and rename fastq files to run the `cellranger-arc` pipeline.

```
lrwxrwxrwx  1 csoto lieber_lcolladotor 80 Jan 11 15:02 path_atac -> /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/24-01-02_SPag110823_ATAC/
lrwxrwxrwx  1 csoto lieber_lcolladotor 75 Jan 11 15:03 path_snrna -> /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/24-01-02_SPag110823/
```
-->
````

Example of tree FASTQ directory structure (raw-data):

```         
15:04 raw-data $ tree FASTQ/
FASTQ/
├── ATAC
└── GEX
```

````{=html}
<!--
Identified files for ATACseq:
```
drwxrws---+ 2 rmiller lieber_lmh_storage 6 Jan  2 11:16 1A_Hb_KDM-1
drwxrws---+ 2 rmiller lieber_lmh_storage 6 Jan  2 11:16 1A_Hb_KDM-2
drwxrws---+ 2 rmiller lieber_lmh_storage 6 Jan  2 11:16 1A_Hb_KDM-3
drwxrws---+ 2 rmiller lieber_lmh_storage 6 Jan  2 11:15 1A_Hb_KDM-4
drwxrws---+ 2 rmiller lieber_lmh_storage 6 Jan  2 11:15 2A_Hb_KDM-1
drwxrws---+ 2 rmiller lieber_lmh_storage 6 Jan  2 11:12 2A_Hb_KDM-2
drwxrws---+ 2 rmiller lieber_lmh_storage 6 Jan  2 11:15 2A_Hb_KDM-3
drwxrws---+ 2 rmiller lieber_lmh_storage 6 Jan  2 11:13 2A_Hb_KDM-4
```

Identified files for snRNAseq:
```
14:05 Hb_2024 $ ls -l | grep Hb
-rw-r-----+ 1 rmiller lieber_lmh_storage 1889340602 Jan  2 16:36 37---1C-Hb-KDM-Hb_S17_L001_R1_001.fastq.gz
-rw-r-----+ 1 rmiller lieber_lmh_storage 4432968966 Jan  2 11:57 37---1C-Hb-KDM-Hb_S17_L001_R2_001.fastq.gz
-rw-r-----+ 1 rmiller lieber_lmh_storage 1908702363 Jan  2 16:29 37---1C-Hb-KDM-Hb_S17_L002_R1_001.fastq.gz
-rw-r-----+ 1 rmiller lieber_lmh_storage 4475186382 Jan  2 11:40 37---1C-Hb-KDM-Hb_S17_L002_R2_001.fastq.gz
-rw-r-----+ 1 rmiller lieber_lmh_storage 1908480908 Jan  2 16:29 37---1C-Hb-KDM-Hb_S17_L003_R1_001.fastq.gz
-rw-r-----+ 1 rmiller lieber_lmh_storage 4464333226 Jan  2 11:44 37---1C-Hb-KDM-Hb_S17_L003_R2_001.fastq.gz
-rw-r-----+ 1 rmiller lieber_lmh_storage 1908465304 Jan  2 16:30 37---1C-Hb-KDM-Hb_S17_L004_R1_001.fastq.gz
-rw-r-----+ 1 rmiller lieber_lmh_storage 4472526848 Jan  2 11:40 37---1C-Hb-KDM-Hb_S17_L004_R2_001.fastq.gz
-rw-r-----+ 1 rmiller lieber_lmh_storage 1697727336 Jan  2 16:58 38---2C-Hb-KDM-Hb_S18_L001_R1_001.fastq.gz
-rw-r-----+ 1 rmiller lieber_lmh_storage 3851384409 Jan  2 13:53 38---2C-Hb-KDM-Hb_S18_L001_R2_001.fastq.gz
-rw-r-----+ 1 rmiller lieber_lmh_storage 1725024267 Jan  2 16:55 38---2C-Hb-KDM-Hb_S18_L002_R1_001.fastq.gz
-rw-r-----+ 1 rmiller lieber_lmh_storage 3907813188 Jan  2 13:35 38---2C-Hb-KDM-Hb_S18_L002_R2_001.fastq.gz
-rw-r-----+ 1 rmiller lieber_lmh_storage 1726565648 Jan  2 16:55 38---2C-Hb-KDM-Hb_S18_L003_R1_001.fastq.gz
-rw-r-----+ 1 rmiller lieber_lmh_storage 3902949402 Jan  2 13:39 38---2C-Hb-KDM-Hb_S18_L003_R2_001.fastq.gz
-rw-r-----+ 1 rmiller lieber_lmh_storage 1722948027 Jan  2 16:56 38---2C-Hb-KDM-Hb_S18_L004_R1_001.fastq.gz
-rw-r-----+ 1 rmiller lieber_lmh_storage 3900694608 Jan  2 13:43 38---2C-Hb-KDM-Hb_S18_L004_R2_001.fastq.gz
```
-->
````

Script used to create soft links to Fastq RawData:

```         
create_cellranger_ARC_libraries.R
```

Script to create the `csv library file` to input in the `CellRanger-ARC` pipeline:

```         
create_fastq_soft_links_for_cellranger_arc.R
```

<br>

IMPORTANT NOTES

1.  These scripts are located into `raw-data/` directory.

2.  These is the only exception to set R Scripts in a different directory to `code/` dir.

<br>

CSC

Aug 06th, 2024
