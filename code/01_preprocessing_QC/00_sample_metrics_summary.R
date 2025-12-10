## script in progress

## Setup the experiment folder and data info
library(here)
#devtools::install_github("haozhu233/kableExtra")
library(kableExtra)

here::here()

## Directory to save variable features 
csvDir <- here("processed-data", "01_preprocessing_QC", "cellranger_atac")
processedDir <- here("processed-data", "cellrangerATAC_summary_rpts", "csv_web_summary_rpt")

if (!dir.exists(csvDir)) dir.create(csvDir)
if (!dir.exists(processedDir)) dir.create(processedDir)

experiment_name = "Habenula-multiome"
ids <- paste0("S", rep(1:12))


## Read in the cellranger sample metrics csv files
d10x.metrics <- lapply(ids, function(i){
  # remove _Counts is if names don't include them
  #metrics <- read.csv(file.path(processedDir,paste0(i,"_Counts/outs"),"metrics_summary.csv"), colClasses = "character")
  metrics <- read.csv(file.path(processedDir, paste0(i, "_summary.csv")), colClasses = "character")
})
experiment.metrics <- do.call("rbind", d10x.metrics)
rownames(experiment.metrics) <- ids

#sequencing_metrics <- data.frame(t(experiment.metrics[,c(4:17,1,18,2,3,19,20)]))
sequencing_metrics <- data.frame(t(experiment.metrics[,c(2:5,1)]))

row.names(sequencing_metrics) <- gsub("\\."," ", rownames(sequencing_metrics))


## And lets generate a pretty table
sequencing_metrics %>%
  kable(caption = 'Cell Ranger ATAC metrics') %>%
  pack_rows("Sequencing Characteristics", 2, 5, label_row_css = "background-color: #666; color: #fff;") %>%
  kable_styling("striped")



sequencing_metrics %>%
  kable(caption = 'Cell Ranger Results') %>%
  pack_rows("Sequencing Characteristics", 1, 29, label_row_css = "background-color: #666; color: #fff;") %>%
  # pack_rows("Mapping Characteristics", 8, 14, label_row_css = "background-color: #666; color: #fff;") %>%
  # pack_rows("Cell Characteristics", 15, 20, label_row_css = "background-color: #666; color: #fff;") %>%
  kable_styling("striped")


