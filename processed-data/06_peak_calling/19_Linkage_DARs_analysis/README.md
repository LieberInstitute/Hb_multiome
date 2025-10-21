# LinkPeaks–DARs Overlap Analysis

Results for significant/selected chromatin regions that are both:

-   **LinkPeaks**: Putative cis-regulatory elements linked to gene expression through correlation\
-   **DARs**: Differentially Accessible Regions (DARs) based on cell-type-specific chromatin accessibility

Overlaps were computed to pinpoint cis-regulatory elements (CREs) of interesr across cell types.

-   Tested in: 18 mid-resolution WNN clusters\
-   FDR threshold: 0.2

------------------------------------------------------------------------

<br> 

#### Key Result Files

|  |  |
|-----------------------|------------------------------------------------|
| [`ALL_DARs_signif_FDR02.csv`](https://github.com/LieberInstitute/Hb_multiome/blob/32b27fc1015923f46d4c44fe742e776e9e0a67e6/processed-data/06_peak_calling/19_Linkage_DARs_analysis/ALL_DARs_signif_FDR02.csv) | All significant DARs across clusters (FDR \< 0.2), based on pseudobulk `voomLmFit` modeling. \~325,580 open chromatin regions (OCRs) identified |
| [`Overlaps_LinkPeak_DARs_FDR0.2.csv`](https://github.com/LieberInstitute/Hb_multiome/blob/32b27fc1015923f46d4c44fe742e776e9e0a67e6/processed-data/06_peak_calling/19_Linkage_DARs_analysis/Overlaps_LinkPeak_DARs_FDR0.2.csv) | Genomic overlaps between LinkPeaks and selected DARs (FDR \< 0.2). \~20,650 overlapping regions with merged metadata 
|  |  |


Related plots:\
[Linkage & DAR Overlap Visualizations](https://github.com/LieberInstitute/Hb_multiome/tree/32b27fc1015923f46d4c44fe742e776e9e0a67e6/plots/06_peak_calling/19_Linkage_DARs_analysis)

------------------------------------------------------------------------

<br>

#### Columns description: `ALL_DARs_signif_FDR02.csv`

<br>

|  |  |
|----------------------------|--------------------------------------------|
| `FDR_threshold` | FDR label (e.g., "FDR0.2") |
| `cell_type` | Cluster/cell type where the DAR was detected |
| `peak_id` | Genomic coordinates in `chr-start-end` format |
| `logFC` | Log₂ fold change in accessibility. Positive = more accessible ("open"); Negative = less accessible ("closed") |
| `fdr` | Adjusted p-value from differential accessibility testing |
| `is_significant` | Boolean flag for peaks below FDR threshold |

------------------------------------------------------------------------

<br>

#### Columns description: `Overlaps_LinkPeak_DARs_FDR0.2.csv`

Each row represents a shared peak that is both:

-   Correlated with gene expression (LinkPeaks)\

-   Differentially accessible in a specific cell type (DARs)

<br>

####  From LinkPeaks

| Column | Description |
|----------------------------|--------------------------------------------|
| `peak_id_links` | LinkPeaks identifier (`chr-start-end`) |
| `CCscore` | Spearman correlation between accessibility and gene expression |
| `gene_name` | Correlated gene symbol |
| `gene_id` | Ensembl gene ID |
| `FDR_CC` | FDR from correlation test |
| `cluster` | Cell type where the link was detected |
| `distance`, `tss`, `signed_distance` | Distances from peak to gene TSS |
| `signed_by_strand` | Signed distance adjusted by gene strand (+/–) |

<br>

#### From DARs Merged Columns

| Column          | Description                      |
|-----------------|----------------------------------|
| `peak_id`       | Genomic coordinate from DAR      |
| `cell_type`     | Cell type where DAR was detected |
| `logFC`         | Log₂ fold change from DARs       |
| `fdr_dars`      | FDR from DARs test               |
| `FDR_threshold` | FDR threshold used (0.2)         |

------------------------------------------------------------------------

<br>

#### Summary

|  |  |  |
|------------------|--------------------------|----------------------------|
| LinkPeaks | Spearman correlation (`CCscore`, `FDR_CC`) | Functional association between a peak and a gene |
| DARs | Differential accessibility (`logFC`, `fdr_dars`) | Structural chromatin change in a cell type |
| Overlap | Exact coordinate match | Indicates that the same region is both correlated and cell-type-specific |

These overlapping peaks represent candidate **cis-regulatory elements (CREs)** that are:

-   Correlated with transcriptional activity
-   Differentially accessible in a specific cell type

<br><br>

