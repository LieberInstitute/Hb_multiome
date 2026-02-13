#!/bin/bash
#SBATCH --partition=transfer
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out


#10x files
wget -O 10X/GSE109158_RAW.tar 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE109158&format=file'
tar -xvf 10X/GSE109158_RAW.tar -C 10X


#Smartseq files
wget -O smartSeq/GSE109159_AH7GKNBGXY.RSEM.genes.counts.matrix.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE109159&format=file&file=GSE109159%5FAH7GKNBGXY%2ERSEM%2Egenes%2Ecounts%2Ematrix%2Egz'
wget -O smartSeq/GSE109159_AH7GKNBGXY.RSEM.genes.tpm.matrix.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE109159&format=file&file=GSE109159%5FAH7GKNBGXY%2ERSEM%2Egenes%2Etpm%2Ematrix%2Egz'
wget -O smartSeq/GSE109159_HGJMWBGXY.RSEM.genes.counts.matrix.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE109159&format=file&file=GSE109159%5FHGJMWBGXY%2ERSEM%2Egenes%2Ecounts%2Ematrix%2Egz'
wget -O smartSeq/GSE109159_HGJMWBGXY.RSEM.genes.tpm.matrix.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE109159&format=file&file=GSE109159%5FHGJMWBGXY%2ERSEM%2Egenes%2Etpm%2Ematrix%2Egz'
wget -O smartSeq/GSE109159_HHKCVBGXY.RSEM.genes.counts.matrix.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE109159&format=file&file=GSE109159%5FHHKCVBGXY%2ERSEM%2Egenes%2Ecounts%2Ematrix%2Egz'
wget -O smartSeq/GSE109159_HHKCVBGXY.RSEM.genes.tpm.matrix.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE109159&format=file&file=GSE109159%5FHHKCVBGXY%2ERSEM%2Egenes%2Etpm%2Ematrix%2Egz' 





