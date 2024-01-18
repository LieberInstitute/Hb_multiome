# HPC multiome pilot project 
### Hippocampus brain tissue from Human and Mouse

***

Here, I have compiled all the essential information regarding the R devel environment and a brief description of the scripts used in this directory.
Please refer to the scripts for a wide explanation of each one.<br>

<br>

### R scripts Empty-drops scripts

<!--
Table.1 R scripts
| Name                            | Description                                                         |
|---------------------------------|---------------------------------------------------------------------|
| 05_DimRed_GEX_ATAC.R            | This script runs *dimensionality reduction* and *clustering* in cellranger-ARC datasets. First, it loads Seurat object with ATAC assay and then it runs PCA for GEX (k-means) and LSA/LSI (Latent Semantic Indexing/Analysis) for ATAC. Lastly it integrates both outputs. | 
| 05b_DimRed_GEX_ATAC_plt_integration.R  | A combination of different plots is created for visualization purposes. | 
| 05_DimRed_GEX_ATAC_SLURM_loop.sh  | Bash script to run as slurm job-array. It requires *array_targets_names.txt* (samples-names) as input. Experiments are run as job-arrays: **$ sbatch 05_DimRed_GEX_ATAC_SLURM_loop.sh** | 
-->

**code $ ls -1  02_emptydrops_for_multiome/**

02b_empty_droplets_stats_integration.R

02c_empty_dropets_vs_cellrangerarc.R

02_empty_droplets_stats_loop.R

slurm_emptydrops_stats.sh 


**Execute in slurm server with the next command:**

$ sbatch slurm_emptydrops_stats.sh

<br>


### Integrating Plots for Comprehensive Data Visualization

This script built 5 different plots as follow:

<!-- *<sample_name>_DimRed_scATAC.png*         UMAP for scRNAseq, PC for gene expression.

*<sample_name>_DimRed_scRNA.png*          UMAP for scATACseq, LSI for chromatin accessibility.

-->

<br>


[Go to main page](../../README.md)
<br><br>


*Cynthia SC*

Oct 30th, 2023

<br><br>

***


