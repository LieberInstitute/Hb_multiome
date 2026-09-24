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
fdr_thr02 = 0.2

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
    # clus = "LHb.7"
    # FDR_thr = 0.2
    
    if (!FDR_thr==0.2) { stop() }
    
    message("Starting Barplot for cluster ", clus)
    
    # Prepare counts of Open and less-open DARs for the specified cluster
    up_down_counts <- res_enrich |>
        mutate(
            logFC = .data[[paste0("logFC_", clus)]], # Access column by dynamic name
            fdr = .data[[paste0("fdr_", clus)]]      # Access column by dynamic name
        ) |>
        filter(
            fdr < FDR_thr                            # Filter for significant DARs
        ) |>
        mutate(
            accessibility = case_when(
                logFC > 0 ~ "Increased Accessibility", # Improved label
                logFC < 0 ~ "Decreased Accessibility", # Improved label
                TRUE ~ "Not"
            )
        ) |>
        filter(accessibility != "Not") |>
        group_by(accessibility) |>
        summarise(n = n(), .groups = "drop")
    
    # Check if there are significant DARs to plot
    if (nrow(up_down_counts) == 0) {
        warning("No significant DARs found for cluster ", clus, " at specified thresholds.")
        return(NULL)
    }
    # Calculate the maximum count for scaling
    max_count <- max(up_down_counts$n, na.rm = TRUE)
    
    g1 <- ggplot(up_down_counts, aes(x = "", y = n, fill = accessibility)) +
        geom_col(position = "dodge") +
        scale_fill_manual(
            values = c("Increased Accessibility" = "red", "Decreased Accessibility" = "blue"),
            breaks = c("Increased Accessibility", "Decreased Accessibility")
        ) +
        labs(
            title = paste0("[", clus, "]"), 
            x = NULL,
            y = NULL,          
            fill = "Differential Accessibility"       
        ) +
        theme_minimal() + 
        geom_text(
            aes(label = n), 
            position = position_dodge(width = 0.9), 
            vjust = -0.75,
            size = 3
        ) + 
        scale_y_continuous(
            expand = c(0, 0), 
            limits = c(0, max_count * 1.15) 
        ) +
        theme(
            plot.title = element_text(size = 10, face = "bold"),
            axis.text.y = element_text(size = 8), 
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        )
    
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
    plot_subtitle <- paste0("Total DARs: ", n_all_dars, " | [FDR < ", FDR_thr, "] = ", n_significant_dars)
    
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


## Parse csv DAR files and plot Volcano and BarPlot for "Up/Down" DARs by cellType

# empty list to store the plots
barPlot_list <- list()

for (ct_DARs in lst_DARs_enrich_cvs) {
    # ct_DARs = lst_DARs_enrich_cvs[2] # "Endo"
    
    message("Processing:\n", ct_DARs)
    
    # load DARs
    DAR_peaks_enrich_cvs <- here(input_cvsDir, ct_DARs)
    if (file.exists(DAR_peaks_enrich_cvs)) {
        res_enrich <- read.csv(DAR_peaks_enrich_cvs)
        message("File loaded!")
    } else {
        stop(paste("File not found:", DAR_peaks_enrich_cvs))
    }
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
    
    # Only for FDR=0.2
    barPlot_list[[clust_name]] <- create_barPlot_significant_DARs(res_enrich, clust_name, fdr_thr02, lfc_thresh)
    
}


# Integrated BarPlots for Open/Close Chromatin
# Only for FDR=0.2

if (length(barPlot_list) > 0) {
    
    combined_plot <- wrap_plots(barPlot_list, ncol = 5)
    
    title_name <- paste0("Number of DARs per Cell Type") 
    
    final_plot <- combined_plot +
        plot_layout(guides = "collect") + 
        plot_annotation(
            title = title_name,
            subtitle =  paste0("FDR < ", fdr_thr02),
            theme = theme(
                # Sets the position of the single merged legend
                legend.position = "bottom",
                element_text(size = rel(10)), 
                plot.title = element_text(size = rel(0.9), face = "bold"), # 90% of base size
                plot.subtitle = element_text(size = rel(0.8))              # 80% of base size
            )
        )
    
    f_name = paste0("BarPlots_ALL_cellTypes_voomLmFit_FDR", fdr_thr02, ".pdf")
    pdf(here(plot_Dir, f_name), width = 8, height = 8)
    print(final_plot)
    dev.off()
    
    message("BarPlots done!")
    
}


message("All done!")


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
