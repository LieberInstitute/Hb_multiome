#### This table is a filtered set of **peak–gene links** that passed the **“Exploratory” tier** threshold (`score > 0.10`), with each row representing a unique association between one ATAC peak and one gene within a specific cluster.

------------------------------------------------------------------------

### Column breakdown

#### Genomic location of the peak

-   **seqnames**, **start**, **end**, **width**, **strand** → coordinates of the ATAC peak in the genome.

#### Link statistics

-   **score** → correlation strength between peak accessibility and gene expression (Pearson/Spearman).
-   **zscore**, **pvalue**, **FDR** → statistical measures of link significance.
-   **tier** → confidence category based on thresholds.

#### Gene information

-   **gene**, **gene_id** → linked gene name and Ensembl ID.
-   **tss** → transcription start site position.
-   **gene_strand** → gene strand orientation.

#### Distance metrics

-   **peak_center**, **distance** → distance from peak center to gene TSS (absolute and signed).
-   **signed_by_strand** → signed distance adjusted for gene strand direction.
-   **distance_kb** → same in kilobases.

#### Peak identity and cluster assignment

-   **peak** → unique peak identifier.
-   **cluster_ann** → WNN multiome cell cluster the peak was assigned to after pseudobulking.

------------------------------------------------------------------------

Essentially, this is a table summarizing working dataset for downstream **“Peaks per gene” / “Genes per peak” heatmaps**,\
each row says:\
*In cluster X, peak Y is linked to gene Z with this correlation and distance.*
