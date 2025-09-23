########################################################################
## Differential Accessibility (DA) using voomLmFit / pseudobulk multiome assays with merged peaks
## ## refs:
## adapted from https://github.com/LieberInstitute/spatialLIBD/blob/40da043d0235e01a12a7f52a0b367d3850bad9e8/R/registration_stats_enrichment.R#L40
##
## Authors. CSC
## Date. Sep 16, 2025
## Recommended resources mem=30GB
########################################################################
library("spatialLIBD")
library("ggplot2")
library("ggrepel") 
library("dplyr")
library("stringr")
library("here")


# Testing spearman at 5e4 on macs2 peaks 
resolution_level = "Mid"
lfc_thresh <- 0.25 # 2^0.25 ≈1.189 
FDR_thr = 0.20

# ## read input arguments
# args = commandArgs(trailingOnly = TRUE)
# clus <- args[1]
# if (is.na(clus) || !nzchar(clus)) stop("Missing cluster_name argument")

## Check/create directories
input_cvsDir <- here( 
    "processed-data",
    "06_peak_calling",
    "17_pseudobulk_DARs_MACS2_reduced_voomLmFit"
)
# these are mac2 peaks merged and normalized 
# Seurat_base_name <- "Mid_pseudobulk.spearman.5e5_merged_peaks.rds"

processed_Dir <- here(
    "processed-data",
    "06_peak_calling",
    "18_pseudobulk_DARs_Volcano_Violin"
)
plot_Dir <- here(
    "plots",
    "06_peak_calling",
    "18_pseudobulk_DARs_Volcano_Violin"
)

if (!dir.exists(plot_Dir)) {
    dir.create(plot_Dir)
}
if (!dir.exists(processed_Dir)) {
    dir.create(processed_Dir)
}

##==============================================================================

## load DARs

pattern = paste0("^voomlmFit_DAR_peaks_ALL_in_.*\\.csv$")
lst_DARs_cvs <- list.files(
    path = input_cvsDir,
    pattern = pattern
)

message("DAR files found:")
lst_DARs_cvs

##==============================================================================

## create volcano
create_volcano <- function(res_enrich, clus,FDR_thr, lfc_thresh, plot_Dir) {
    
    plot_data <- res_enrich |>
        select(
            peak_id,
            logFC = paste0("logFC_", clus),
            fdr = paste0("fdr_", clus)
        ) |>
        mutate(
            is_significant = case_when(
                fdr < FDR_thr ~ "Significant",
                #fdr < FDR_thr & abs(logFC) > lfc_thresh ~ "Significant",
                TRUE ~ "Not Significant"
            )
        )
    
    # Calculate the number of all and significant DARs
    n_all_dars <- nrow(plot_data)
    n_significant_dars <- sum(plot_data$is_significant == "Significant")
    
    plot_title <- paste0("Differentially Accessible Regions (DARs) in ", clus)
    plot_subtitle <- paste0("Total DARs: ", n_all_dars, " | Significant DARs: ", n_significant_dars)
    
    f_name <- here(plot_Dir, paste0("Volcano_", clus, "_voomLmFit.pdf"))
    pdf(f_name, width = 7, height = 6)  
    
    g1 <- ggplot(plot_data, aes(x = logFC, y = -log10(fdr), color = is_significant)) +
        geom_point(alpha = 0.5, size = 1) +
        scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
        geom_hline(yintercept = -log10(FDR_thr), linetype = "dashed", color = "blue") +
        geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", color = "blue") +
        labs(
            title = plot_title,
            subtitle = plot_subtitle,
            x = "Log2 Fold Change (logFC)",
            y = "-Log10(FDR)",
            color = "Significance"
        ) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0))
    
    print(g1)
    dev.off()
    message("Volcano done!")

}


## Parse csv DAR files and plot Volcano and ViolinPlot for the 5 "Up/Down" DARs by cellType

for (ct_DARs in lst_DARs_cvs) {
    # ct_DARs = lst_DARs_cvs[2] # "Endo"
    
    message("Processing:\n", ct_DARs)
    
    # load DARs
    DAR_peaks_cvs <- here(input_cvsDir, ct_DARs)
    if (file.exists(DAR_peaks_cvs)) {
        DARs_df <- read.csv(DAR_peaks_cvs)
        message("File loaded!")
    } else {
        stop(paste("File not found:", DAR_peaks_cvs))
    }
    #colnames(DARs_df)

    # get cluster name. ge. "Endo"
    clust_name <- sub("^voomlmFit_DAR_peaks_ALL_in_(.*)\\.csv$", "\\1", ct_DARs)
        
    create_volcano(DARs_df, clust_name, FDR_thr, lfc_thresh, plot_Dir)
    
    # # save enrichment stats
    # f_name <- here(output_Dir, paste0("voomlmFit_DAR_peaks_ALL_in_", clus, ".csv"))
    # write.csv(res_enrich, f_name, row.names = FALSE)
    # 
    # # Identify top up/down peaks based on FDR
    # fdr_cols <- grep("^fdr_", colnames(res_enrich), value = TRUE)
    # 
    # # cluster-specific significant peaks
    # res_sig_ct <- res_enrich |>
    #     # filter(.data[[paste0("fdr_", clus)]] < FDR_thr) |>
    #     filter(.data[[paste0("fdr_", clus)]] < FDR_thr,
    #            abs(.data[[paste0("logFC_", clus)]]) > lfc_thresh) |>
    #     mutate(
    #         logFC_ct = .data[[paste0("logFC_", clus)]],
    #         direction = case_when(
    #             logFC_ct >  0 ~ "Up",    # opening
    #             logFC_ct <  0 ~ "Down",  # closing
    #             TRUE ~ "NS"              # should not occur if you filtered
    #         )
    #     )
    # 
    # if (nrow(res_sig_ct) > 0) { 
    #     f_name <- here(output_Dir, paste0("voomlmFit_DAR_peaks_", clus, ".csv"))
    #     write.csv(res_sig_ct, f_name, row.names = FALSE)
    #     message("Enrichment statistics saved [", clus, "]")    
    # }

}

message("Plots done!")


# library("slurmjobs")
# job_single(
#   "18_pseudobulk_DARs_Volcano_Violin",
#   create_shell = TRUE,
#   partition = "katun",
#   memory = "30G",
#   cores = 2,
#   logdir = "logs",
#   command = "Rscript 18_pseudobulk_DARs_Volcano_Violin.R",
#   create_logdir = FALSE
# )


## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
