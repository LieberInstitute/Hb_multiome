# Run locally with Cytoscape open

library(sessioninfo)
library(tidyverse)
library(here)
library(duckplyr)
library(RCy3)

trio_path = here(
    "processed-data", "13_tripod_trios", "09_gather_trios",
    "filtered_trios_fine.parquet"
)
cell_map_path = here("raw-data", "cell_type_map.csv")
out_dir = here("plots", "13_tripod_trios", "18_cytoscape_prep")
pdf_dir = file.path(out_dir, "pdf_pages")
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

num_cores = Sys.getenv("SLURM_CPUS_PER_TASK") |>
    as.integer() |>
    dplyr::coalesce(1L)
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

cluster_map = read_csv(cell_map_path, show_col_types = FALSE)
rename_map = stats::setNames(
    cluster_map$new_cell_type,
    cluster_map$old_cell_type
)

tf_targets = read_parquet_duckdb(trio_path, prudence = "stingy") |>
    filter(is_intersect, stringency_level == 2) |>
    select(gene, TF, cell_type, coef, adj) |>
    collect() |>
    mutate(
        coef = log(coef) - min(log(coef)),
        neg_log10_adj = -log10(adj),
        cell_type = dplyr::coalesce(
            unname(rename_map[cell_type]),
            cell_type
        )
    )

build_nodes = function(tf_targets) {
    bind_rows(
        tf_targets |>
            distinct(id = TF) |>
            mutate(role = "TF"),
        tf_targets |>
            distinct(id = gene) |>
            mutate(role = "Gene")
    ) |>
        distinct(id, .keep_all = TRUE) |>
        as.data.frame()
}

build_edges = function(tf_targets) {
    tf_targets |>
        transmute(
            source = TF,
            target = gene,
            interaction = "regulates",
            cell_type,
            coef,
            adj,
            neg_log10_adj
        ) |>
        as.data.frame()
}

make_style_name = function(title) {
    paste0("style__", gsub("[^A-Za-z0-9]+", "_", title))
}

make_pdf_path = function(title) {
    file.path(
        pdf_dir,
        paste0(gsub("[^A-Za-z0-9]+", "_", title), ".pdf")
    )
}

reset_style = function(style_name) {
    if (style_name %in% getVisualStyleNames()) {
        deleteVisualStyle(style_name)
    }

    createVisualStyle(
        style.name = style_name,
        defaults = list(
            NODE_SIZE = 35,
            EDGE_TARGET_ARROW_SHAPE = "DELTA"
        )
    )

    setNodeFontSizeDefault(12, style.name = style_name)
    setNodeSizeDefault(48, style.name = style_name)
    setNodeFontFaceDefault("SansSerif,bold,18", style.name = style_name)
}

map_shared_node_style = function(style_name, network_suid) {
    setNodeLabelMapping("id", style.name = style_name, network = network_suid)
    setNodeColorMapping(
        table.column = "role",
        table.column.values = c("TF", "Gene"),
        colors = c("#E69F00", "#56B4E9"),
        mapping.type = "d",
        style.name = style_name,
        network = network_suid
    )
}

map_edge_width = function(style_name, network_suid, coef_range) {
    setEdgeLineWidthMapping(
        table.column = "coef",
        table.column.values = coef_range,
        widths = c(1, 8),
        mapping.type = "c",
        style.name = style_name,
        network = network_suid
    )
}

map_edge_color_by_cell_type = function(style_name, network_suid, edges) {
    cell_types = sort(unique(edges$cell_type))
    cell_colors = unname(cell_type_colors[cell_types])

    setEdgeColorMapping(
        table.column = "cell_type",
        table.column.values = cell_types,
        colors = cell_colors,
        mapping.type = "d",
        style.name = style_name,
        network = network_suid
    )
}

map_edge_color_by_neg_log10_adj = function(
    style_name,
    network_suid,
    neg_log10_adj_range
) {
    setEdgeColorMapping(
        table.column = "neg_log10_adj",
        table.column.values = neg_log10_adj_range,
        colors = c("#d9d9d9", "#b2182b"),
        mapping.type = "c",
        style.name = style_name,
        network = network_suid
    )
}

create_and_export_network = function(
    nodes,
    edges,
    title,
    edge_color_mapper,
    edge_color_args,
    coef_range,
    collection = "TF target networks"
) {
    style_name = make_style_name(title)
    pdf_path = make_pdf_path(title)

    network_suid = createNetworkFromDataFrames(
        nodes = nodes,
        edges = edges,
        title = title,
        collection = collection
    )

    reset_style(style_name)
    map_shared_node_style(style_name, network_suid)
    do.call(
        edge_color_mapper,
        c(
            list(style_name = style_name, network_suid = network_suid),
            edge_color_args
        )
    )
    map_edge_width(style_name, network_suid, coef_range)
    setVisualStyle(style_name, network = network_suid)
    layoutNetwork("force-directed defaultSpringLength=50", network = network_suid)
    fitContent(network = network_suid)
    exportImage(
        filename = pdf_path,
        type = "PDF",
        network = network_suid,
        overwriteFile = TRUE
    )

    tibble(
        network_suid = network_suid,
        title = title,
        style_name = style_name,
        pdf_path = pdf_path
    )
}

overview_edges = build_edges(tf_targets)
overview_nodes = build_nodes(tf_targets)
hb_targets = tf_targets |>
    filter(str_detect(cell_type, "Hb"))
hb_edges = build_edges(hb_targets)
hb_nodes = build_nodes(hb_targets)
coef_range = range(overview_edges$coef, na.rm = TRUE)
neg_log10_adj_range = range(overview_edges$neg_log10_adj, na.rm = TRUE)

bind_rows(
    create_and_export_network(
        nodes = overview_nodes,
        edges = overview_edges,
        title = "TF target network",
        edge_color_mapper = map_edge_color_by_cell_type,
        edge_color_args = list(edges = overview_edges),
        coef_range = coef_range
    ),
    create_and_export_network(
        nodes = hb_nodes,
        edges = hb_edges,
        title = "TF target network - Hb only",
        edge_color_mapper = map_edge_color_by_cell_type,
        edge_color_args = list(edges = hb_edges),
        coef_range = coef_range
    ),
    tf_targets |>
        split(~ cell_type) |>
        imap(function(cell_type_targets, cell_type_name) {
            create_and_export_network(
                nodes = build_nodes(cell_type_targets),
                edges = build_edges(cell_type_targets),
                title = paste("TF target network", cell_type_name, sep = " - "),
                edge_color_mapper = map_edge_color_by_neg_log10_adj,
                edge_color_args = list(
                    neg_log10_adj_range = neg_log10_adj_range
                ),
                coef_range = coef_range
            )
        }) |>
        bind_rows()
)

session_info()
