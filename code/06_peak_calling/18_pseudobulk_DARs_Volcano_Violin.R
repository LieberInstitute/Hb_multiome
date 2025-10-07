########################################################################
## Visualizations for DARs by cellType
## - data computer with modified version for registration_stats_enrichment.R using voomLmFit
## Refs:
## adapted from https://github.com/LieberInstitute/spatialLIBD/blob/40da043d0235e01a12a7f52a0b367d3850bad9e8/R/registration_stats_enrichment.R#L40
##
## Authors. CSC
## Date. Sep 23, 2025
## Recommended resources mem=30GB
########################################################################

library("spatialLIBD")
library("purrr")
library("ggplot2")
library("ggrepel") 
library("patchwork")
library("dplyr")
library("stringr")
library("here")


# Testing spearman at 5e4 on macs2 peaks 
resolution_level = "Mid"
lfc_thresh <- 0.25 # 2^0.25 ≈ 1.189 
FDR_thr = c(0.10, 0.20)

## Check/create directories
input_cvsDir <- here( 
    "processed-data",
    "06_peak_calling",
    "17_pseudobulk_DARs_MACS2_reduced_voomLmFit"
)

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

## load DARs enriched peaks 

pattern = paste0("^voomlmFit_DAR_peaks_ALL_in_.*\\.csv$")
lst_DARs_enrich_cvs <- list.files(
    path = input_cvsDir,
    pattern = pattern
)

message("DAR files found:")
lst_DARs_enrich_cvs

##==============================================================================

create_barPlot_significant_DARs <- function(res_enrich, clus, FDR_thr, lfc_thresh) {
    ## Barplot of Up vs Down DARs per cluster
    # testing: 
    # res_enrich = res_enrich
    # clust_name = "Endo"
    # FDR_thr = 0.2
    
    message("Starting Barplot for cluster ", clus)
    
    # Prepare counts of Up and Down DARs for the specified cluster
    up_down_counts <- res_enrich |>
        mutate(
            logFC = .data[[paste0("logFC_", clus)]], # Access column by dynamic name
            fdr = .data[[paste0("fdr_", clus)]]      # Access column by dynamic name
        ) |>
        filter(
            fdr < FDR_thr                           # Filter for significant DARs
            # abs(logFC) > lfc_thresh               # Filter by logFC threshold
        ) |>
        mutate(
            direction = case_when(
                # logFC > lfc_thresh ~ "Up",
                # logFC < -lfc_thresh ~ "Down",
                logFC > 0 ~ "Up", # Up/Down is now based on positive/negative logFC
                logFC < 0 ~ "Down",
                TRUE ~ "Not"
            )
        ) |>
        filter(direction != "Not") |>
        group_by(direction) |>
        summarise(n = n(), .groups = "drop")
    
    g1 <- ggplot(up_down_counts, aes(x = "", y = n, fill = direction)) +
        geom_col(position = "dodge") +
        scale_fill_manual(values = c("Up" = "red", "Down" = "blue")) +
        labs(
            title = clus,
            x = "",
            y = "Number of DARs",
            fill = "Direction"
        ) +
        theme_minimal() +
        geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.5)
    
    message("Up/Down Barplot done!")
    return(g1)
    
}


##==============================================================================

## create volcano
create_volcano <- function(res_enrich, clus, FDR_thr, lfc_thresh, plot_Dir) {
    # res_enrich=res_enrich
    # clus=clust_name
    # FDR_thr=0.2
        
    message("Starting Volcano for cluster ", clus)
    
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
    plot_subtitle <- paste0("Total DARs: ", n_all_dars, " | Significant DARs: ", n_significant_dars, " | FDR = ", FDR_thr)
    
    g1 <- ggplot(plot_data, aes(x = logFC, y = -log10(fdr), color = is_significant)) +
        geom_point(alpha = 0.5, size = 1) +
        scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
        geom_hline(yintercept = -log10(FDR_thr), linetype = "dashed", color = "blue") +
        geom_vline(xintercept = 0, linetype = "solid", color = "darkgrey") +
        labs(
            title = plot_title,
            subtitle = plot_subtitle,
            x = "Log2 Fold Change (logFC)",
            y = "-Log10(FDR)",
            color = "Significance"
        ) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0))
    
    f_name <- here(plot_Dir, paste0("Volcano_", clus, "_voomLmFit_FDR", FDR_thr, ".pdf"))
    pdf(f_name, width = 7, height = 6)  
    print(g1)
    dev.off()
    
    print("Volcano done!")

}


##==============================================================================


## Parse csv DAR files and plot Volcano and ViolinPlot for the 5 "Up/Down" DARs by cellType

# empty list to store the plots
barPlot_list <- list()

for (ct_DARs in lst_DARs_enrich_cvs) {
    # ct_DARs = lst_DARs_enrich_cvs[1] # "Endo"
    
    message("Processing:\n", ct_DARs)
    
    # load DARs
    DAR_peaks_enrich_cvs <- here(input_cvsDir, ct_DARs)
    if (file.exists(DAR_peaks_enrich_cvs)) {
        res_enrich <- read.csv(DAR_peaks_enrich_cvs)
        message("File loaded!")
    } else {
        stop(paste("File not found:", DAR_peaks_enrich_cvs))
    }
    #colnames(res_enrich)

    # get cluster name. ge. "Endo"
    clust_name <- sub("^voomlmFit_DAR_peaks_ALL_in_(.*)\\.csv$", "\\1", ct_DARs)
    
    ## create plot for FDR_thr (s)     
    purrr::map(FDR_thr, ~ {
        create_volcano(
            res_enrich, 
            clust_name,
            FDR_thr = .x,   # current FDR
            lfc_thresh = lfc_thresh,
            plot_Dir
        )
    })
    
    barPlot_list[[clust_name]] <- create_barPlot_significant_DARs(res_enrich, clust_name, FDR_thr, lfc_thresh)
    
    # Identify top up/down peaks based on FDR
    fdr_cols <- grep("^fdr_", colnames(res_enrich), value = TRUE)

    # cluster-specific significant peaks / at 2 FDR thr
    FDR_thr
    # [1] 0.1 0.2
    
    # Loop over each threshold
    purrr::walk(FDR_thr, function(fdr_value) {
        
        # Filter significant peaks for this FDR threshold
        res_sig_ct <- res_enrich |>
            filter(.data[[paste0("fdr_", clust_name)]] < fdr_value) |>
            # abs(.data[[paste0("logFC_", clust_name)]]) > lfc_thresh) |>
            mutate(
                logFC_ct = .data[[paste0("logFC_", clust_name)]],
                direction = case_when(
                    logFC_ct >  0 ~ "Up",    # opening
                    logFC_ct <  0 ~ "Down",  # closing
                    TRUE ~ "NS"
                )
            )
        if (nrow(res_sig_ct) > 0) {
            message(glue::glue(
                "[{clust_name}] {nrow(res_sig_ct)} peaks pass FDR < {fdr_value}"
            ))
            
            # Build output file name
            f_name <- here(processed_Dir, paste0("voomlmFit_DAR_filtered_peaks_", clust_name, "_FDR", fdr_value,".csv"))
            write.csv(res_sig_ct, f_name, row.names = FALSE)
            message("Filtered peaks (cvs) saved [", clust_name, "]")
            
        } else {
            message("No significant peaks for FDR < ", fdr_value, "found")
        }
    })
    
}


# Plot barPlots
if (length(barPlot_list) > 0) {
    
    pdf(here(plot_Dir, "BarPlots_ALL_cellTypes_voomLmFit.pdf"), width = 10, height = 10)
    
    # Remove redundant y-axis labels
    barPlot_list_clean <- lapply(barPlot_list, function(p) {
        p + ylab(NULL) + theme(legend.position = "none")
    })
    
    combined_plot <- wrap_plots(barPlot_list_clean, ncol = 5)

    title_name <- paste0("Number of Up/Down Enriched DARs by CellType (FDR < ", FDR_thr, ")") 
    final_plot <- combined_plot +
        plot_annotation(
            title = title_name,
            theme = theme(
                plot.title = element_text(size = 12, hjust = 0.5),
                axis.title.y = element_text(size = 9)
            )
        )

    print(final_plot)
    
    dev.off()
    
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
