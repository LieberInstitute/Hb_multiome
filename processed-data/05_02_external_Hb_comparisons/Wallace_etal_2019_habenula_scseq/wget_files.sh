#!/bin/bash
#SBATCH --partition=transfer
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out

wget -O habData.csv "https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/2VFWF6/4MT4I4"
wget -O hab_batch1.rds "https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/2VFWF6/YA7RBW"
wget -O ori_author_read_ME.docx "https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/2VFWF6/WLPRKB"



