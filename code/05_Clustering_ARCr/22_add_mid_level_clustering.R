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





