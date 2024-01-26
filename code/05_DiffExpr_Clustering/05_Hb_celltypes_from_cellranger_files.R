########################################################################
## Read cellranger-arc clusters to look for cell-types (Fast cluster exploration)
## INPUT:
## Gene-marker list
## CSV file with clusters produced by cellranger-arc
## 
## Authors. CSC
## Date. Jan 26th, 2024 / last md. XXX 
########################################################################

# load libraries
library(here)

here::here()

# Samples for testing
#sample_tmp <- '1_HPC_KDM,human'  # testing HUMAN tissue
#sample_tmp <- '3_HPC_KDM,mouse'  # testing MOUSE tissue 
#sample_tmp <- 'hippo42_1,human'    # testing human tissue -- smaller 615 cells

s_sample <- sample_tmp %>% str_split(",") %>%flatten_chr() %>% .[[1]]
s_tissue <- sample_tmp %>% str_split(",") %>%flatten_chr() %>% .[[2]]
message("Processing sample: ", s_sample, ' / ', s_tissue)


############  Initials ############

## read our gene markers (those found in our clusters)
s_file_name <- here('processed-data/06_Diff_expr_genes/Gene_markers', paste0(s_sample,'_cluster.Top5markers_POS_minPC025_log2FC025.csv'))
my_markers <- as.data.frame(read.csv(s_file_name, header = TRUE))

## Build another data frame for the Nelson gene-markers list to compare gene-markers
## Notice Nelson's markers are probably comming from mouse donors because these are in lower case
nlst_Astrocytes <- toupper(list("Aqp4", "Gfap"))
nlst_Oligodendrocytes <- toupper(list("Mbp", "Mobp"))
nlst_Oligodendrocyte_precursor <- toupper(list("Pdgfra", "Vcan"))

nlst_Microglia <- list("Csf1r", "C3", "Tmem119")
nlst_Choroid_plexus <- list("Ttr", "Folr1")
nlst_Vascular_endothelial <- list("Cldn5", "Flt1")
nlst_Neurons_general <- list("Syt1", "Snap25")
nlst_IN_neurons <- list("Gad1", "Gad2", "Sst", "Pvalb", "Htr3a")
nlst_EX_neurons <- list("Slc17a7", "Slc17a6", "Nrgn")
n_lst_EXsub_GranCells <- list("Prox1", "Glis3")
n_lst_EXsub_CA4 <- list("Calb2", "Csf2rb2")
n_lst_EXsub_CA3 <- list("Nectin3")
n_lst_EXsub_CA2 <- list("Ntsr2", "Rgs14")
n_lst_EXsub_CA1 <- list("Man1a", "Fibcd1", "Tbc1d1")
n_lst_EXsub_SubProS <- list("Fn1", "Nts", "Ndst4", "Dcn")

nelson_gene_markers <- list(nlst_Astrocytes, nlst_Oligodendrocytes, nlst_Oligodendrocyte_precursor)

# Current markers used here comes from Erik Nelson project
# Note: this nomenclature could come from Mouse data, due usually the nomenclature is that gene names are all caps for human (SLC17A7) (Slc17a7)
# Non-neuronal:
# Astrocytes: Aqp4, Gfap
# Oligodendrocytes: Mbp, Mobp
# Oligodendrocyte precursor cells: Pdgfra, Vcan
# Microglia: Csf1r, C3, Tmem119
# Choroid plexus: Ttr, Folr1
# Vascular/endothelial: Cldn5, Flt1
# Neurons, general: Syt1, Snap25
# Inhibitory neurons: Gad1, Gad2, Sst, Pvalb, Htr3a
# Excitatory neurons: Slc17a7, Slc17a6, Nrgn
# Excitatory subtypes:
# Granule cells: Prox1, Glis3
# CA4: Calb2, Csf2rb2
# CA3: Nectin3
# CA2: Ntsr2, Rgs14
# CA1: Man1a, Fibcd1, Tbc1d1
# Sub/ProS: Fn1, Nts, Ndst4, Dcn

nelson_gene_markers$gene
# View cluster levels
View(Idents(SeuratOBJ))
clusters <- 0:cluster_idx

map_celltypes <- function(gene_marker_lst) {
    ## Parse known gene markers to annotate cluster cell-types 
    
    show(gene_marker_lst)
    for (gene in (gene_marker_lst)) {
        my_markers %>% 
            rowwise %>% 
            mutate(firstZero = names(.)[match(0, c_across(everything()), 
                nomatch = ncol(.))]) %>% 
                ungroup
    }
    
    #new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
    #                     "NK", "DC", "Platelet")
    names(new.cluster.ids) <- levels(pbmc)
    pbmc <- RenameIdents(pbmc, new.cluster.ids)
    DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
    
}


for (clust in (clusters)) {
    message('Integrating explorative plots for gene markers sample ', s_sample)
    message('Processing cluster: ', clust)
    # map and integrate the plots
    map(clust, map_celltypes)
}

message(' All gene markers parsed')







message("Tasks completed successfully!")


################################################################################


library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
#"2023-04-04 12:42:26 EDT"
proc.time()
options(width = 120)
session_info()

