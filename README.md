# Hb_multiome

Human Habenula Multiome Project

### Overview

#### Local project name: Hb_multiome

Allocated in: `/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome`


#### Fastq multiome snRNA-seq and snATAC-seq datasets from `Human Habenula Neurotypical Controls`

JHPCE Path to `first`  sequenced dataset (January 2024):

snRNA-Seq fastq libraries in: `/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/24-01-02_SPag110823/`

snATAC-Seq fastq libraries in: `/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/24-01-02_SPag110823_ATAC/`

Path to `second`  sequenced dataset (July- early August 2024):

snRNA-Seq fastq libraries in: `/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2024-07-09-Psomagen/AN00019739_10X_RawData_Outs/`

snATAC-Seq fastq libraries in: `/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2024-08-02_Psomagen/AN00020424_10X_RawData_Outs/`

Path to `third`  sequenced dataset (late August 2024):

snRNA-Seq fastq libraries in: `/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2024-08-26_Psomagen/`

snATAC-Seq fastq libraries in: `/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2024-08-27_Psomagen/`

<!-- 

Notes: - Second sequencing have `Visium` sequencing files too. - `RNA multiome files` have a `C` letter; `Visum files` have a `v` lower case letter.

<br>

````{=html}

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


To create soft links to the fastq multiome data use the following scripts:

```         
ls create_cellranger_ARC_libraries* -l
```

To create the `csv library naming convention cellrangerARC` files use the following scripts:

```         
ls create_fastq_soft_links_for_cellranger_arc* -l
```

<br>

Note these R scripts are located into the `raw-data/` directory. These is the only exception to set R Scripts in a different directory to `code/` dir.

<br>


#### Summary reports to cellranger pipelines:

[CellrangerARC web summary report](https://github.com/LieberInstitute/Hb_multiome/tree/732545d59dcbb1087621414eb39c1b093832aa5e/processed-data/cellrangerARC_summary_rpts)

[Cellranger-count web summary report](https://github.com/LieberInstitute/Hb_multiome/tree/732545d59dcbb1087621414eb39c1b093832aa5e/processed-data/cellrangerGEX_summary_rpts)

[Cellranger-atac web summary report](https://github.com/LieberInstitute/Hb_multiome/tree/732545d59dcbb1087621414eb39c1b093832aa5e/processed-data/cellrangerATAC_summary_rpts)



CSC
Sep. 2024
