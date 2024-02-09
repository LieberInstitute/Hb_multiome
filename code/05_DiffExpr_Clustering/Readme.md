# Habenula multiome project 
### Identification of cell types 

***

<br>


## PART 1. Identification of cell types in cellranger-arc clusters


Clusters generated for *cellranger-arc* were extracted from each sample's `~/outs/analysis/clustering/gex/graphclust/` directory and used to search for cell types given a list of custom marker genes.


### R scripts


Table 1. R scripts
| Name                            | Description                                                         |
|---------------------------------|---------------------------------------------------------------------|
| 05_Hb_celltypes_from_cellranger_files.R   | Reads a Cellranger-arc DGE (cluster) CSV file to identify cell types (fast cluster scan) given a defined marker gene list.  |
| remote_DGE_marker_gene_lists.R    |   This is the marker gene list builder. It is located in *code/functions_custom*  |

<br>


### Output files


Table 2. Outputs
| Name                            | Description                                                         |
|---------------------------------|---------------------------------------------------------------------|
| *<sample_name>*_cell_types_Top50r_gm10.csv   | Table with cell-types found in each cluster using the *Top50r_gm* marker list. Considering only the Top 10 DGE genes by cluster   |
| *<sample_name>*_cell_types_Top50r_gm20.csv   | Table with cell-types found in each cluster using the *Top50r_gm* marker list. Considering only the Top 20 DGE genes by cluster   |
| *<sample_name>*_cell_types_erik_gm10.csv      | Table with cell-types found in each cluster using the *erik_gm* marker list. Considering only the Top 10 DGE genes by cluster   |
| *<sample_name>*_cell_types_erik_gm20.csv      | Table with cell-types found in each cluster using the *erik_gm* marker list. Considering only the Top 20 DGE genes by cluster   |

<br>
Outputs are located at: `~/Habenula_R01_LIBD4270/Hb_multiome/processed-data/05_DiffExpr_Clustering/`


<br>


### Marker Gene Lists


*We employ two marker gene lists:*
<br>


1) Erik's marker gene list. Named `erik_gm`

```
> markers.custom 
$neuron
[1] "SYT1"   "SNAP25"

$excitatory_neuron
[1] "SLC17A6" "SLC17A7"

$inhibitory_neuron
[1] "GAD1" "GAD2"

$`mediodorsal thalamus`
 [1] "EPHA4"   "PDYN"    "LYPD6B"  "LYPD6"   "S1PR1"   "GBX2"    "RAMP3"  
 [8] "COX6A2"  "SLITRK6" "DGAT2"   "ADARB2" 

$`Pf/PVT`
[1] "INHBA" "NPTXR"

$`Hb neuron specific`
[1] "POU4F1" "GPR151"

$`MHB neuron specific`
[1] "TAC1"   "CHAT"   "CHRNB4" "TAC3"  

$`LHB neuron specific`
[1] "HTR2C" "MMRN1" "ANO3" 

$oligodendrocyte
[1] "MOBP" "MBP" 

$oligodendrocyte_precursor
[1] "PDGFRA" "VCAN"  

$microglia
[1] "C3"    "CSF1R"

$astrocyte
[1] "GFAP" "AQP4"

$`Endo/CP`
[1] "TTR"   "FOLR1" "FLT1"  "CLDN5"
```
<br>


2) Top50 by ratio marker gene list (curated list focus on Habenula). Named `Top50r_gm`

```
> markers.custom                                                              
$LHb
 [1] "HTR4"       "BVES"       "NRP1"       "HTR2C"      "CDH4"      
 [6] "COL25A1"    "CBLN2"      "CACNA1I"    "WSCD2"      "SOX1-OT"   
[11] "AC092969.1" "MDGA1"      "GALR1"      "TRAM2"      "ADCYAP1"   
[16] "TEX2"       "AL078604.4" "AL391807.1" "TLL1"       "NTNG2"     
[21] "EPHA5-AS1"  "PCDH11X"    "LINC01876"  "KAZN-AS1"   "RFTN1"     
[26] "LGI2"       "OVCH1"      "ZNF235"     "ZNF208"     "LINC01572" 
[31] "STS"        "SLC12A8"    "LINC00886"  "AL357873.1" "CCDC85A"   
[36] "CLVS2"      "SLC7A14"    "ARHGAP28"   "TMEM196"    "FAM81A"    
[41] "AC007614.1" "KIAA1217"   "PXDNL"      "CLEC2L"     "EPHA5"     
[46] "PRKCB"      "ADCY6"      "DCHS2"      "GLIS1"      "GREB1L"    

$MHb
 [1] "CHAT"       "LINC01307"  "NEUROD1"    "CHRNB4"     "LINC02143" 
 [6] "AC114321.1" "AC104170.1" "AC079760.2" "AC024610.2" "AC022382.2"
[11] "LINC02240"  "AC103843.1" "GNG8"       "SLC5A7"     "GPR149"    
[16] "CERKL"      "ROBO3"      "LINC00616"  "AC008415.1" "AC109466.1"
[21] "KCNMA1-AS1" "CHRNB3"     "CYP1B1-AS1" "CPNE9"      "GPR139"    
[26] "ACOXL"      "CD24"       "CHRNA3"     "SCUBE1"     "RASGRP1"   
[31] "KIAA2012"   "CASZ1"      "PDE2A"      "LINC02505"  "GIPR"      
[36] "AC008667.1" "MPP7"       "AC003991.1" "VAV2"       "THSD7B"    
[41] "CACNG5"     "GSDME"      "AL512662.2" "AC004594.1" "GDPD5"     
[46] "AL137139.2" "PLCH2"      "CNGB1"      "GPR39"      "AC015845.2"
```
<br>


### Top DGE genes

<br>

We evaluated two criteria:

1) The first top 10 DGE genes by cluster.

2) The first top 20 DGE genes by cluster.

<br>

DGE genes were organized by *adjusted p-value* and posterior *values<0.05* were kept. A subset of the top 10 or top 20 genes were then used to search for markers. Below is provided an example:

<br>

```
# Feature.ID Feature.Name Cluster.1.Mean.Counts Cluster.1.Log2.fold.change
# 1 ENSG00000243485  MIR1302-2HG                     0                   7.110087
# 2 ENSG00000237613      FAM138A                     0                   7.110087
# 3 ENSG00000186092        OR4F5                     0                   7.110087
# Cluster.1.Adjusted.p.value
# 1                          1
# 2                          1
# 3                          1
```


<br>
<br>



[Go to main page](../../README.md)
<br><br>


*Cynthia SC*

Feb 9th, 2024

<br><br>



#### Useful resources consulted:



<br><br>


***


