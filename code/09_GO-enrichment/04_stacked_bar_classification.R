library(tidyverse)
library(here)
library(sessioninfo)

## setup variables and file names
FDR = 0.2
inputCSV_Overlaps_Dir <- here(
    "processed-data",
    "06_peak_calling"
)

message("Loading Overlap ...")

## load classified overlaps with 2-shared ct
overlap_df2 <- read.csv(here(inputCSV_Overlaps_Dir, "21_overlaping_FDRscores_TopHeatmap", 
                            "overlaps_linkPeak_DARs_classified_thr_CC0.3_thr_DAR0.1.csv"))

nrow(overlap_df2) #3205
colnames(overlap_df2)
# [1] "source_df"           "peak_id_links"       "CCscore"            
# [4] "gene_name"           "gene_id"             "FDR_CC"             
# [7] "cluster"             "tss"                 "gene_strand"        
# [10] "distance"            "distance_kb"         "signed_distance"    
# [13] "signed_by_strand"    "peak_id"             "cell_type"          
# [16] "FDR_threshold"       "logFC"               "fdr_dars"           
# [19] "overlap_type"        "n_cell_types"        "sig_CC"             
# [22] "sig_DAR"             "category"            "type_classification"

overlap_df2 <- overlap_df2 |>
    mutate(RE_type = case_when(
        abs(distance_kb) <= 2 ~ "Promoter",
        abs(distance_kb) <= 10 ~ "Near-Promoter",
        abs(distance_kb) <= 100 ~ "Enhancer",
        abs(distance_kb) <= 500 ~ "Long-range",
        TRUE ~ "Distal >500kb"
    ))
head(overlap_df2)


plot_df <- overlap_df2 |>
    count(category, RE_type) |>
    group_by(category) |>
    mutate(pct = n / sum(n) * 100)

head(plot_df)
# category                RE_type           n   pct
# <chr>                   <chr>         <int> <dbl>
#     1 Discordant Linked DAR   Enhancer    221 35.8 
# 2 Discordant Linked DAR   Long-range      365 59.2 
# 3 Discordant Linked DAR   Near-Promoter    20  3.24
# 4 Discordant Linked DAR   Promoter         11  1.78
# 5 Linked OCR (+) enriched Enhancer        105 37.8 
# 6 Linked OCR (+) enriched Long-range      155 55.8 


session_info()