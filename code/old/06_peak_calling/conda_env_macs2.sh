## CallPeaks() from Signac is a function to search peaks in ATAC data
## It internally call macs2 package, so you need to install it properly
## - Here, you will find a draft with suggestions to create an isolated env for macs2 to avoid dependency conflicts
## - Plus a pseudo slurm script to activate your conda env

 /users/csoto/.conda/envs/macs2_conda3_env/bin/macs2 --version

#=====================================
# if necessary (optional), update your conda
conda update -n base -c conda-forge conda

#=====================================
# JHPCE env preparation
# conda_R/4.4.x is too slow, hangs too long the installation
module unload conda_R/4.4.x 
which macs
## Unloading conda_R/4.4.x

# So, the alternative was to unload the module and kept conda/3-24.3.0, which is faster
module list
## Currently Loaded Modules:
##     1) JHPCE_ROCKY9_DEFAULT_ENV   2) JHPCE_tools/3.0   3) conda/3-24.3.0

#=====================================
# Create env and Install macs2

# 1. Create a clean environment with Python 3.8 (recommended)
conda create -n macs2_conda3_env python=3.8

## Why Python 3.8?
# - Stable version
# - Supported by Signac and most R workflows
# - Works reliably with Python 3.6–3.9, especially using conda install -c bioconda macs2

# 2. Activate it
conda activate macs2_conda3_env

# 3. (optional) prioritize channels over others when resolving packages 
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# 4. Install MACS2 from bioconda. Takes long time .... so consider use bioconda
conda install macs2

# or ensures Conda won’t waste time trying multiple combinations of packages across channels
# I used this ONE
conda install -c bioconda macs2


#=====================================
# If it still hangs too long — Try to priorize channel to strict mode

# simplifies dependency resolution
conda config --set channel_priority strict  
conda config --show channels
# channels:
#       - conda-forge
#       - bioconda
#       - defaults

#=====================================

# Alternatively use Mamba which is 10–100× faster at solving environments (you need to install it)

conda install -c conda-forge mamba
mamba install -c bioconda macs2 

#=====================================

# Confirm Installation
which macs2
macs2 --version

# seurat_atac <- CallPeaks(
#     object = seurat_atac,
#     group.by = "seurat_clusters",
#     macs2.path = "/home/your_user/miniconda3/envs/macs2_env/bin/macs2"
# )

#====================================
# for slurm

#!/bin/bash
#SBATCH --job-name=macs2_job
#SBATCH --output=macs2_out.log
#SBATCH --error=macs2_err.log
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4

# Load any necessary modules (optional)
# module load anaconda3

#  find your Conda environment path:
conda info --envs

# Conda setup (do NOT skip)
# Initialize Conda (adjust to your installation path)
source /home/your_user/miniconda3/etc/profile.d/conda.sh

# Activate your env
conda activate macs2_env

# Optional: confirm it's active
echo "Conda env: $CONDA_DEFAULT_ENV"
which macs2
macs2 --version

# call peaks
macs2 callpeak -t fragments.bed -f BED -g hs -n atac_peaks --nomodel --shift -100 --extsize 200


