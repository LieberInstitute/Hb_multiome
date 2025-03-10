
## Working in a bash job to compute: (IN PROGRESS)

################################################################################
## Batch correction and analysis of multiome data for WNN clustering with Seurat
##
## This script gather all the slurm jobs to perform:
## Batch correction on multiome data for improved WNN clustering, identify cell 
## types and generate a diverse set of plots for exploratory data analysis (EDA)
## and, optionally, evaluate clustering performance using well-known metrics.
##
## Author: CSC
## Date: Feb, 2025
################################################################################

## Workflow:
## 1. prepare dataset with outliers removed based on customized barcodes detected
## 2. compute batch correction with harmony on rna and atac
## 3. Plot combinations of clustering for EDA
## 4. Cell-type identification with gene-markers reference from both `Literature-based` 
##.   and `Data-Driven` (Human Habenula Project) genes
## Optionally, evaluate clusters using Jaccard and Silhouette scores with labeled clusters. 
##    You may first run a quick annotation using 05_rename_idents.sh



MAINDIR="/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome"
CODEDIR="${MAINDIR}/code"
PROCESSEDIR="${MAINDIR}/processed-data"
PLOTDIR="${MAINDIR}/plots"


######################  Seruat prepared without outliers (1) ###################

# cd code/03_pseudobulking
## change directory
SUBDIR="03_pseudobulking"

cd ${CODEDIR}/${SUBDIR}
echo "Current code dir: ${CODEDIR}/${SUBDIR}/"

03_pseudobulking/07_prepare_seurat_reanalyze_without_outliers.R


######################  Compute batch correction (2) ###########################

## - run snn and lsi and run harmonize on both assay. plus plot unintgerated, integrated and tsne for futher reference)

rm logs/08_harmony_CR_ARCr.txt
rm -f ../../plots/03_pseudobulking/cellrangerARC_reanalyze/*.png
rm -f ../../processed-data/03_pseudobulking/cellrangerARC_reanalyze/seurat.norm_counts_ARCr_harmony_atac_rna_QCed.rds

sbatch 08_harmony_CR_ARCr.sh


###################### Plot combinations of clustering for EDA (3) #############

## change directory
SUBDIR="05_Clustering_ARCr"

cd ${CODEDIR}/${SUBDIR}
echo "Current code dir: ${CODEDIR}/${SUBDIR}/"

cd code/05_Clustering_ARCr
rm -f logs/01_clustering_std_method_*.txt
rm logs/01_clustering_wnn_method_*.txt
rm ../../plots/05_Clustering_ARCr/01_clustering_std_method/*.png
rm ../../processed-data/05_Clustering_ARCr/01_clustering_std_method/*.rds
rm ../../processed-data/05_Clustering_ARCr/01_clustering_std_method/cvs_files/*.csv

## previous versions: 
## 01_clustering_std_method.sh: has harmonyzed RNA data 

echo "Running `01_clustering_std_method_v2.sh` with harmonyzed data for both RNA and ATAC"

sbatch 01_clustering_std_method_v2.sh


###################### Cell-type identification with gene-markers reference (4) #############

## sub-subdir: 02_Hb_celltypes_from_seurat_reanalyze_v3

rm -f logs/05_rename_idents*.txt
rm -f ../../processed-data/05_Clustering_ARCr/02_Hb_celltypes_from_seurat_reanalyze_v3/cvs_files_markers/*cluster_info.csv
rm -f ../../processed-data/05_Clustering_ARCr/02_Hb_celltypes_from_seurat_reanalyze_v3/cvs_files_markers/*cellTypes_integrated_top*.csv

## The version v3 is an improved script to identify cell-types based in an integrated reference using `literature-based` and `data-driven` genes with no repetitions (remove redundant genes to avoid cell-type duplications in the reports)

sbatch 02_Hb_celltypes_from_seurat_reanalyze_v3.sh



################################################################################
## This chunk OPTIONAL to make a sort of fast clustering evaluation using **Jaccard** and **Silhouette**
################################################################################

## Quickly annotates (renames identities) WNN clusters in Seurat objects using a gene marker reference. Convenient for evaluating cluster concordance and proximity, providing a pseudo-annotation for easy identification of target clusters (Habenula)
## Note. Requires to run first `01_clustering_std_method_v2.sh` and `02_Hb_celltypes_from_seurat_reanalyze_v3.sh`

# cd code/05_Clustering_ARCr

outSUBDIR="05_rename_idents"

rm -f logs/05_rename_idents*.txt
rm -f ${PLOTDIR}/${SUBDIR}/${outSUBDIR}/*.pdf
rm -f ${PROCESSEDIR}/${SUBDIR}/${outSUBDIR}/*.rds

sbatch 05_rename_idents.sh
# squeue -u csoto

# cd code/05_Clustering_ARCr

## manage directories

outSUBDIR="06_jaccard_annotated"

rm -f logs/06_jaccard_annotated*.txt
rm -f ${PLOTDIR}/${SUBDIR}/${outSUBDIR}/*_Silhouette_WNN.png
rm -f ${PROCESSEDIR}/${SUBDIR}/${outSUBDIR}/*_Silhouette_WNN.cvs


sbatch 06_jaccard_annotated.sh




03_cell_types_perc_WNN_leidenR1_knn30_v2.sh


