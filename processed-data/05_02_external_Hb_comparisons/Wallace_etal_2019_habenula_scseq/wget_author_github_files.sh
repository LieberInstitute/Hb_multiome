#!/bin/bash
#SBATCH --partition=transfer
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out

wget https://raw.githubusercontent.com/mwall2017/habenula_indrops/master/Habenula_Seurat_meta.RData
wget https://raw.githubusercontent.com/mwall2017/habenula_indrops/master/lhb_Seurat_meta.RData
wget https://raw.githubusercontent.com/mwall2017/habenula_indrops/master/mhb_Seurat_meta.RData






