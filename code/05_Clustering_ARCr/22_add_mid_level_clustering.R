########################################################################
## Add mid-level resolution clustering
##
## Authors. CSC
## Date. August 15, 2025
## Recommended resources on interactive mode: srun --pty --mem=80GB --x11 bash
## Note. Seurat objects were created with module load conda_R/4.3.x
########################################################################

library("Seurat")
library("Signac")
library("ggplot2")
library("tidyverse")
library("dbplyr")
library("colorspace") # make color gradients 
library("here")

# Check/create directories
## clusters renamed for Spatial-Registration on Visium project

inputRDS_Dir <- here(
    "processed-data",
    "05_Clustering_ARCr",
    "17_wnn_clustering_final_ct"
)
outputRDS_Dir <- here(
    "processed-data",
    "05_Clustering_ARCr",
    "22_add_mid_level_clustering"
)
plotDir <- here(
    "plots",
    "05_Clustering_ARCr",
    "22_add_mid_level_clustering"
)
outputCSV_Dir <- here(
    "processed-data", 
    "05_Clustering_ARCr",
    "22_add_mid_level_clustering"
)

## Check directories
if (!dir.exists(plotDir)) {
    dir.create(plotDir)
}
if (!dir.exists(outputCSV_Dir)) {
    dir.create(outputCSV_Dir)
}

message("Loading Seurat ... ")

# Use Seurat with final ct-annotations
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_HD.rds"
seurat_name <- here(inputRDS_Dir, Seurat_base_name)
SeuratOBJ <- readRDS(here(inputRDS_Dir, Seurat_base_name))
SeuratOBJ

##==============================================================================
## Adding mid level clustering resolution

# Define mid-level
clusters_to_merge <- c(SeuratOBJ$cluster_ann)

# Create a vector mapping old cluster names to MID-LEVEL clusters
cluster_map <- sapply(clusters_to_merge, function(x) {
    if (grepl("LHb.4", x)) { return("LHb.4")
    } else if (grepl("LHb.2.7", x)) { return("LHb.2.7")
    } else if (grepl("LHb.1$", x)) { return("LHb.1")
    } else if (grepl("LHb.1.3$", x)) { return("LHb.1.3")
    } else if (grepl("LHb.1.3.4", x)) { return("LHb.1.3.4")
    } else if (grepl("LHb.7", x)) { return("LHb.7")
    } else if (grepl("MHb.1$", x)) { return("MHb.1")
    } else if (grepl("MHb.1.2", x)) { return("MHb.1.2")
    } else if (grepl("MHb.2", x)) { return("MHb.2")
    } else if (grepl("MHb.3", x)) { return("MHb.3")
    } else if (grepl("Inhib.Thal", x)) { return("Inhib.Thal")
    } else if (grepl("Excit.Thal", x)) { return("Excit.Thal")
    } else if (grepl("Astrocyte", x)) { return("Astrocyte")
    } else if (grepl("OPC", x)) { return("OPC")
    } else if (grepl("Oligo", x)) { return("Oligo")
    } else if (grepl("Microglia", x)) { return("Microglia")
    } else if (grepl("Thal", x)) { return("Thal")
    } else if (grepl("Endo", x)) { return("Endo")
    } else { return(NA_character_)  # unmatched case
    }
})

names(cluster_map) <- clusters_to_merge
# map current identities to LHb/MHb
mid_idents <- plyr::revalue(as.character(Idents(SeuratOBJ)), cluster_map)
# assign new identities
Idents(SeuratOBJ) <- mid_idents
## add as meta=data
SeuratOBJ$mid_cluster <- mid_idents

# verify mid-levels
if (length(unique(SeuratOBJ$cluster_ann[is.na(SeuratOBJ$mid_cluster)])) > 0) {
    stop(
        "All clusters should be assigned to a mid-level resolution\n",
        "Levels not assigned\n",
        paste(unique(SeuratOBJ$cluster_ann[is.na(SeuratOBJ$mid_cluster)]), collapse = "\n")
    )
} else {
    message("Mid-level clusters added:")
    data.frame(mid_cluster = sort(unique(SeuratOBJ$mid_cluster)))
    nrow(data.frame(mid_cluster = sort(unique(SeuratOBJ$mid_cluster))))
}

# colnames(SeuratOBJ@meta.data)
msg <- paste0("Mid-level Cell-Types:\n", paste(sort(unique(SeuratOBJ$mid_cluster)), collapse = "\n"))
message(msg)

## =============================================================================

## save RDS
rds_file_name <- here(outputRDS_Dir, Seurat_base_name)
saveRDS(SeuratOBJ, rds_file_name)

message("Seurat with mid-level clusters meta-data saved!")


## =============================================================================
## Picked up Hex-color codes similar across cell-type

my_colors <- c(
    LHb = "#1f78b4",
    MHb = "#ad1d8c",
    Oligo = "#384a08",
    Astrocyte = "#532222", 
    OPC = "#829454",
    Microglia = "#141b02",
    Endo = "#d95f02",
    Inhib_Thal = "#9a9fe7",
    Excit_Thal = "#42467b",
    Thal = "#4d55b7"
)

## assign color gradients to mid resolution clusters based on Broad cell-types

# extract LHb and MHb clusters
cluster_levels <- levels(SeuratOBJ)
cluster_levels
# [1] "Excit.Thal" "LHb.4"      "Inhib.Thal" "Astrocyte"  "MHb.1.2"   
# [6] "LHb.1"      "OPC"        "Oligo"      "Microglia"  "LHb.2.7"   
# [11] "Endo"       "LHb.1.3.4"  "MHb.1"      "MHb.2"      "MHb.3"     
# [16] "Thal"       "LHb.1.3"    "LHb.7" 
LHb_clusters <- grep("LHb", cluster_levels, value = TRUE)
MHb_clusters <- grep("MHb", cluster_levels, value = TRUE)

# Create tonal gradients for LHb and MHb
LHb_colors <- sequential_hcl(length(LHb_clusters), h = 210, c = 80, l = c(30, 80))
MHb_colors <- sequential_hcl(length(MHb_clusters), h = 320, c = 80, l = c(30, 80))

# Build full cluster color map
my_colors_mid <- setNames(rep("#bdbdbd", length(cluster_levels)), cluster_levels)
# Excit.Thal      LHb.4 Inhib.Thal  Astrocyte    MHb.1.2      LHb.1        OPC 
# "#bdbdbd"  "#bdbdbd"  "#bdbdbd"  "#bdbdbd"  "#bdbdbd"  "#bdbdbd"  "#bdbdbd" 
# Oligo  Microglia    LHb.2.7       Endo  LHb.1.3.4      MHb.1      MHb.2 
# "#bdbdbd"  "#bdbdbd"  "#bdbdbd"  "#bdbdbd"  "#bdbdbd"  "#bdbdbd"  "#bdbdbd" 
# MHb.3       Thal    LHb.1.3      LHb.7 
# "#bdbdbd"  "#bdbdbd"  "#bdbdbd"  "#bdbdbd" 

my_colors_mid[LHb_clusters] <- LHb_colors
my_colors_mid[MHb_clusters] <- MHb_colors

# assign base color for other types from your existing palette
for (category in c("Oligo", "Astrocyte", "OPC", "Microglia", "Endo", "Inhib.Thal", "Excit.Thal", "Thal")) {
    matched <- grep(category, cluster_levels, value = TRUE)
    my_colors_mid[matched] <- my_colors[[gsub("\\.", "_", category)]]
}
#scales::show_col(my_colors_fine)

## =============================================================================

message("Processing UMAP ...")

## extract suffix name to give unique name to plots
seurat_name <- str_extract(Seurat_base_name, pattern = "k[3:4]0\\_C\\.\\w*")

Reductions(SeuratOBJ)

plt1 <- DimPlot(SeuratOBJ, 
                label = TRUE, 
                reduction = "wnn.umap",
                group.by = "mid_cluster",
                label.size = 3,
                cols = my_colors_mid) + 
    NoLegend() +
    labs(title = "WNN cell types (Mid-resolution)")

ggsave(here(plotDir, "WNN_umap_mid.pdf"), plt1, width = 7, height = 7)

# library("slurmjobs")
# job_single(
#     "22_add_mid_level_clustering",
#     cores = 2,
#     partition = "katun",
#     memory = "80G",
#     create_shell = TRUE
#     )

## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()








