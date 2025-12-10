#!/bin/bash
#SBATCH --partition=katun	
#SBATCH --job-name=run_all_peaks
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH --mem=2GB
#SBATCH -o /dev/null
#SBATCH -e /dev/null
# SBATCH --mail-type=ALL

## Run all for LinkPeaks and DARs in pseudobulk multiome assays
## - Peaks from macs2 --> reduce() 

echo "**** Job starts ****"
echo "Run run_all_peaks"

date

echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job name (script): ${SLURM_JOB_NAME}"
echo "Hostname (computer node): ${HOSTNAME}"
echo ""
echo "Job id: ${SLURM_JOBID}"

MAINDIR="/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome"
CODEDIR="${MAINDIR}/code"
PROCESSEDIR="${MAINDIR}/processed-data"
PLOTDIR="${MAINDIR}/plots"

## change dir
SUBDIR="06_peak_calling"
cd ${CODEDIR}/${SUBDIR}
pwd

echo "===== Call macs2  ................."

## rm previous log files and output files
rm -f ${CODEDIR}/${SUBDIR}/logs/01_call_peaks_MACS2*.txt
rm -f ${PROCESSEDIR}/${SUBDIR}/macs_*.csv
id1=$(sbatch --parsable 01_call_peaks_MACS2.sh)    
echo $id1

echo "===== Merge macs2 peaks  ................."

id2=$(sbatch --parsable --dependency=afterok:$id1 11_peaks_merge_MACS2.sh)
echo $id2

echo "===== Pseudobulk Merge macs2 peaks  ................."

id3=$(sbatch --parsable --dependency=afterok:$id2 12_pseudobulk_MACS2.sh)
echo $id3

echo "===== LinkPeaks in pb assays with merge macs2 peaks / by cell-type  ................."

rm -f ${CODEDIR}/${SUBDIR}/logs/13_pseudobulk_LinkPeaks_MACS2_split*.log
rm -r ${PROCESSEDIR}/${SUBDIR}/links_ct_merged
# rm -r ${PROCESSEDIR}/${SUBDIR}/links_ct_not_merged
# rm -r ${PROCESSEDIR}/${SUBDIR}/links_global_regular # this are comming from cellRanger-ARC (global measures not split by ct)

id4=$(sbatch --parsable --dependency=afterok:$id3 13_pseudobulk_LinkPeaks_MACS2_split_ct.sh)


## Make distribution plots to evaluate Peak scores
sbatch --dependency=afterok:$id4 14_exploratory_pb_peak_scores_MACS2.sh

echo "===== MACS2 peaks widths at different clustering resolutions  ................."

sbatch --parsable --dependency=afterok:$id3 15_pseudobulk_compare_link_peak_gene_distributions_MACS2.R


echo "===== DARs (search markers) in pb assays with merge macs2 peaks / by cell-type  ................."

id5=$(sbatch --parsable --dependency=afterok:$id3 16_pseudobulk_DARs_MACS2_reduced.sh)










echo "===== DARs (voomLmFit) in pb assays with merge macs2 peaks / by cell-type  ................."

id6=$(sbatch --parsable --dependency=afterok:$id3 17_pseudobulk_DARs_MACS2_reduced_voomLmFit.sh)


echo "===== Visualizations for DARs by cellType  ................."

rm -f ${CODEDIR}/${SUBDIR}/logs/18_pseudobulk_DARs_Volcano_Violin.txt
rm -f ${PLOTDIR}/${SUBDIR}/18_pseudobulk_DARs_Volcano_Violin/*.pdf
# sbatch 18_pseudobulk_DARs_Volcano_Violin.sh
id6=$(sbatch --parsable --dependency=afterok:$id6 18_pseudobulk_DARs_Volcano_Violin.sh)


echo "===== 19 xxxx  ................."

rm -f ${CODEDIR}/${SUBDIR}/logs/...
sbatch <- 19_Linkage_DARs_analysis.R


echo "===== 20 xxxx  ................."

rm -f ${CODEDIR}/${SUBDIR}/logs/20_explore_overlapping_linkpeaks_DARs.txt.txt
rm -f ${PROCESSEDIR}/${SUBDIR}/20_explore_overlapping_linkpeaks_DARs/overlaps*.csv
rm -f ${PLOTDIR}/${SUBDIR}/20_explore_overlapping_linkpeaks_DARs/overlaps*.pdf

# sbatch 20_explore_overlapping_linkpeaks_DARs.sh
id6=$(sbatch --parsable --dependency=afterok:$id6 20_explore_overlapping_linkpeaks_DARs.sh)


echo "Done!!-------------------------------------------------------------------"

echo "**** Job ends ****"

date



## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/