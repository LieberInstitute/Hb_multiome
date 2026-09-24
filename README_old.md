# Hb_multiome

<br>

## Project Description

<br> 

------------------------------------------------------------------------

## Local Data in JHPCE

#### Github Project Name: Hb_multiome

Allocated in: `/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome`

<br>

#### Fastq multiome snRNA-seq and snATAC-seq datasets from `Human Habenula Neurotypical Controls`

**JHPCE Path to `first` sequenced dataset (January 2024):**

snRNA-Seq fastq libraries in: `/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/24-01-02_SPag110823/`

snATAC-Seq fastq libraries in: `/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/24-01-02_SPag110823_ATAC/`

<br>

**Path to `second` sequenced dataset (July- early August 2024):**

snRNA-Seq fastq libraries in: `/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2024-07-09-Psomagen/AN00019739_10X_RawData_Outs/`

snATAC-Seq fastq libraries in: `/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2024-08-02_Psomagen/AN00020424_10X_RawData_Outs/`

<br>

**Path to `third` sequenced dataset (late August 2024):**

snRNA-Seq fastq libraries in: `/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2024-08-26_Psomagen/`

snATAC-Seq fastq libraries in: `/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2024-08-27_Psomagen/`

<br>

Example of tree FASTQ directory structure (raw-data):

```         
15:04 raw-data $ tree FASTQ/
FASTQ/
├── ATAC
└── GEX
```

<br>

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


------------------------------------------------------------------------

## Core Results

## 10x Genomics Cell Ranger Summary Reports:

[CellrangerARC web summary report](https://github.com/LieberInstitute/Hb_multiome/tree/732545d59dcbb1087621414eb39c1b093832aa5e/processed-data/cellrangerARC_summary_rpts)

[Cellranger-count web summary report](https://github.com/LieberInstitute/Hb_multiome/tree/732545d59dcbb1087621414eb39c1b093832aa5e/processed-data/cellrangerGEX_summary_rpts)

[Cellranger-atac web summary report](https://github.com/LieberInstitute/Hb_multiome/tree/732545d59dcbb1087621414eb39c1b093832aa5e/processed-data/cellrangerATAC_summary_rpts)

<br> 

------------------------------------------------------------------------

## RNA-Multiome

*(Section content to be added later)*

📁 **Data Folder:**\
[RNA-Multiome Processed Data](https://github.com/LieberInstitute/Hb_multiome/tree/main/processed-data/05_RNA_multiome)

📊 **Related Plots:**\
[RNA-Multiome Plots](https://github.com/LieberInstitute/Hb_multiome/tree/main/plots/05_RNA_multiome)

<br> 

------------------------------------------------------------------------

## ATAC-Multiome

### LinkPeaks–DARs Overlap Analysis

This analysis identifies significant chromatin regions that serve as both cis-regulatory elements (LinkPeaks) and differentially accessible regions (DARs). Overlaps were assessed across cell types:

-   **Resolution:** 18 mid-resolution WNN clusters\
-   **Threshold:** False Discovery Rate (FDR) ≤ 0.2

<br> 

#### Key Files:

| File | Description |
|------------------------|------------------------------------------------|
| [`ALL_DARs_signif_FDR02.csv`](https://github.com/LieberInstitute/Hb_multiome/tree/32b27fc1015923f46d4c44fe742e776e9e0a67e6/processed-data/06_peak_calling/19_Linkage_DARs_analysis/ALL_DARs_signif_FDR02.csv) | Identifies DARs representing chromatin opening/closing in specific cell types |
| [`Overlaps_LinkPeak_DARs_FDR0.2.csv`](https://github.com/LieberInstitute/Hb_multiome/tree/32b27fc1015923f46d4c44fe742e776e9e0a67e6/processed-data/06_peak_calling/19_Linkage_DARs_analysis/Overlaps_LinkPeak_DARs_FDR0.2.csv) | Validates that correlated peaks are also cell-type-specific in accessibility |

📁 **Data Folder:**\
[19_Linkage_DARs_analysis (processed-data)](https://github.com/LieberInstitute/Hb_multiome/tree/ad1f1062c61f9f3d8909cc75f063e286c5273caa/processed-data/06_peak_calling/19_Linkage_DARs_analysis)

📊 **Related Plots:**\
[Linkage & DAR Overlap Visualizations (plots)](https://github.com/LieberInstitute/Hb_multiome/tree/ad1f1062c61f9f3d8909cc75f063e286c5273caa/plots/06_peak_calling/19_Linkage_DARs_analysis)

<br><br>

------------------------------------------------------------------------

### LinkPeaks–DARs Overlap Summaries (Unique vs Shared)

This section summarizes and visualizes the overlap structure between LinkPeaks and DARs at FDR = 0.2, distinguishing `unique` (cell-type–specific) and `shared` (multi-cell-type) regulatory relationships.

<br> 

#### Key Files:

| File | Description |
|------|-------------|
| [`overlaps_summary_linkPeak_DARs_unique_ct_FDR0.2.csv`](https://github.com/LieberInstitute/Hb_multiome/blob/040a7df1fa0007ba88fddf1508efe93de0602deb/processed-data/06_peak_calling/20_explore_overlapping_linkpeaks_DARs/overlaps_summary_linkPeak_DARs_unique_ct_FDR0.2.csv) | Summarizes all LinkPeaks that overlap DARs at FDR = 0.2, reporting for each peak the number of DARs and distinct cell types it overlaps. Used to distinguish unique (cell-type–specific) and shared (multi-cell-type) regulatory peaks. ≈ 5,603 peaks identified. |
| [`overlaps_unique_linkPeak_DARs_detail_FDR0.2.csv`](https://github.com/LieberInstitute/Hb_multiome/blob/040a7df1fa0007ba88fddf1508efe93de0602deb/processed-data/06_peak_calling/20_explore_overlapping_linkpeaks_DARs/overlaps_unique_linkPeak_DARs_detail_FDR0.2.csv) | Detailed table of unique LinkPeaks overlapping DARs (FDR = 0.2), including full accessibility direction (`logFC`) and per-cell-type annotation. Contains ≈ 2,900 rows representing LinkPeak–DAR pairs. |

📁 Data folder:  
[20_explore_overlapping_linkpeaks_DARs (processed-data)](https://github.com/LieberInstitute/Hb_multiome/tree/040a7df1fa0007ba88fddf1508efe93de0602deb/processed-data/06_peak_calling/20_explore_overlapping_linkpeaks_DARs)

📊 Related plots:  
[Unique/Shared Overlap Visualizations (plots)](https://github.com/LieberInstitute/Hb_multiome/tree/040a7df1fa0007ba88fddf1508efe93de0602deb/plots/06_peak_calling/20_explore_overlapping_linkpeaks_DARs)

------------------------------------------------------------------------


*Cynthia SC*

*October 2025*
