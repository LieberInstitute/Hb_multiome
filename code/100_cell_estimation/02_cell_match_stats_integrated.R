#library(slurmjobs)
library("dplyr")
library("readr")
library("here")


########################    Initials ########################  

### working directories
here::here()
main_dir_name <- "100_cell_match_cellranger_gex_atac"
processedDir <- here("processed-data", main_dir_name)
pattern <- "\\.tsv$"
filenames <- list.files(processedDir, pattern=pattern, full.names=TRUE)
basename(filenames)

# all_files <- lapply(filenames, function(x) {
#   read.table(file = x, sep = '\t', 
#              header = FALSE, stringsAsFactors = FALSE)[2,] # only extract 2nd row
#   })
all_files <- lapply(filenames, function(x) {
  read.csv(file = x, 
             header = FALSE, stringsAsFactors = FALSE)[2,] # only extract 2nd row
})


df_all_stats <- bind_rows(all_files)
head(df_all_stats)
df_all_stats$V1 <- NULL

# set column names
# ln_header <- read.table(filenames[1], sep = '\t', header = FALSE)[1,] # only extract 1nd row
ln_header <- read.csv(filenames[1], header = FALSE)[1,] # only extract 1nd row
ln_header$V1 <- NULL
colnames(df_all_stats) <- trimws(c(as.character(ln_header)))

noquote(colnames(df_all_stats))[2:11]
columns_sel <- c(colnames(df_all_stats))[2:11]

# Convert multiple columns to numeric
df_all_stats <- df_all_stats |>
  mutate_at(vars(columns_sel), as.numeric)
str(df_all_stats)

write.csv(df_all_stats, row.names = TRUE, quote = FALSE, 
          here(processedDir, "all_stats.csv"))

message(" Integrated statistics saved on ", here(processedDir, "all_stats.csv"))


# lev <- as.vector(df_all_stats[1])
# 
# df_plt <- df_all_stats |> 
#   gather(Cells, Value, -Sample_id) |>
#   mutate(Sample_id = factor(Sample_id, 
#                         levels = c("10C_Hb_KDM-10A_Hb_KDM", "11C_Hb_KDM-11A_Hb_KDM", "12C_Hb_KDM-12A_Hb_KDM", 
#                                    "1C-Hb-KDM-Hb-1A_Hb_KDM", "2C-Hb-KDM-Hb-2A_Hb_KDM", "3C_Hb_KDM-3A_Hb_KDM", 
#                                    "4C_Hb_KDM-4A_Hb_KDM", "5C_Hb_KDM-5A_Hb_KDM", "6C_Hb_KDM-6A_Hb_KDM", 
#                                    "9C_Hb_KDM-9A_Hb_KDM"))) |>
#   arrange(Sample_id) 
# 
# 
# gex_plt <- df_plt |> filter(Cells == "translated_gex_bc" | Cells == "valid_bc_gex_arc" | Cells == "not_valid_bc_gex_arc")
# 
# bargraph <- ggplot(data = gex_plt) +
#   geom_bar(aes(x = Sample_id,
#                y = Value,
#                fill = Cells,
#                color = Cells),
#            stat = "identity",
#            position = position_dodge()) +
#   theme(legend.title=element_blank())
# bargraph


# job_single(
#   name = "02_cell_match_stats_integrated", memory = "10G", cores = 2, create_shell = TRUE
# )

