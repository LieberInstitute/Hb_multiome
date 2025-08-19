#!/bin/bash
#SBATCH --partition=katun	
#SBATCH --job-name=run_all_peaks
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH --mem=2GB
#SBATCH -o /dev/null
#SBATCH -e /dev/null
# SBATCH --mail-type=ALL

00_link_peaks.sh
06_exploratory_peak_scores.sh


set -e


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

echo "===== Make Coverage Plots for cannonical genes and top5 DEG ................."
## rm previous log files and output files
rm -f ${CODEDIR}/${SUBDIR}/logs/01_coverage_basic.txt
rm -f ${PLOTDIR}/${SUBDIR}/*.png
rm -f ${PLOTDIR}/${SUBDIR}/*.pdf
sbatch 01_coverage_basic.sh
echo "Done!!-------------------------------------------------------------------"

echo "===== Call Peaks - 01_coverage_basic.R .................................."

# id1=$(sbatch --parsable 01_clustering_std_method_v2.sh)
# echo $id1
echo "Done!!-------------------------------------------------------------------"




echo "**** Job ends ****"


date



## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/