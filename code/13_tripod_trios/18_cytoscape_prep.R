#   Run locally on a laptop with Cytoscape installed

library(sessioninfo)
library(tidyverse)
library(here)
library(duckplyr)
library(RCy3)
library(RColorBrewer)

trio_path = here(
    "processed-data", "13_tripod_trios", "09_gather_trios",
    "filtered_trios_fine.parquet"
)
cell_map_path = here('raw-data', 'cell_type_map.csv')
cell_type_colors = c(
    MHb_A = "#5e0c01",
    MHb_B = "#943f02",
    MHb_C = "#f67104",
    MHb_D = "#f4d5ab",
    LHb_A = "#ee9630",
    LHb_B = "#306171",
    LHb_C = "#082844",
    GABA_LHb_C.1 = "#9c66c0",
    GABA_LHb_C.2 = "#5e0c56",
    Excit.Thal = "#2e6296",
    Inhib.Thal = "#8DADCA",
    Astrocyte = "#972f2f",
    Endo = "#f65a45",
    Ependymal = "#dbb369",
    Microglia = "#141b02",
    Oligo = "#384a08",
    OPC = "#829454"
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(dirname(out_path), showWarnings = FALSE)

cluster_map = read_csv(cell_map_path, show_col_types = FALSE)
rename_map <- stats::setNames(
    cluster_map$new_cell_type, cluster_map$old_cell_type
)

tf_targets <- read_parquet_duckdb(trio_path, prudence = 'stingy') |>
    filter(is_intersect, stringency_level == 2) |>
    select(gene, TF, cell_type, coef, adj) |>
    collect() |>
    mutate(
        coef = log(coef) - min(log(coef)),
        cell_type = dplyr::coalesce(
            unname(rename_map[cell_type]), cell_type
        )
    )

nodes <- bind_rows(
  tf_targets %>%
    distinct(id = TF) %>%
    mutate(role = "TF"),

  tf_targets %>%
    distinct(id = gene) %>%
    mutate(role = "Gene")
) %>%
  distinct(id, .keep_all = TRUE) %>%
  as.data.frame()

edges <- tf_targets %>%
  transmute(
    source = TF,
    target = gene,
    interaction = "regulates",
    cell_type,
    coef,
    adj
  ) %>%
  as.data.frame()

network_suid <- createNetworkFromDataFrames(
  nodes = nodes,
  edges = edges,
  title = "TF target network"
)

style_name <- "TF target style"

createVisualStyle(
  style.name = style_name,
  defaults = list(
    NODE_SIZE = 35,
    EDGE_TARGET_ARROW_SHAPE = "DELTA"
  )
)

setVisualStyle(style_name)

setNodeLabelMapping("id", style.name = style_name)

setNodeColorMapping(
  table.column = "role",
  table.column.values = c("TF", "Gene"),
  colors = c("#E69F00", "#56B4E9"),
  mapping.type = "d",
  style.name = style_name
)

cell_types <- sort(unique(edges$cell_type))
cell_colors <- unname(cell_type_colors[cell_types])

setEdgeColorMapping(
  table.column = "cell_type",
  table.column.values = cell_types,
  colors = cell_colors,
  mapping.type = "d",
  style.name = style_name
)

coef_range <- range(edges$coef, na.rm = TRUE)

setEdgeLineWidthMapping(
  table.column = "coef",
  table.column.values = coef_range,
  widths = c(1, 8),
  mapping.type = "c",
  style.name = style_name
)

setNodeFontSizeDefault(
  new.size = 12,
  style.name = style_name
)
setNodeSizeDefault(
  new.size = 48,
  style.name = style_name
)
setNodeFontFaceDefault(
  "SansSerif,bold,18",
  style.name = style_name
)

layoutNetwork(
  "force-directed defaultSpringLength=50"
)

fitContent()

session_info()
