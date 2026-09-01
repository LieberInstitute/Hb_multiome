# Run locally with Cytoscape open

library(sessioninfo)
library(tidyverse)
library(here)
library(duckplyr)
library(qs2)
library(RCy3)

#   This was manually changed in an interactive R session. Can be "intersect" or
#   "level_2"
mode = "level_2"

trio_path = here(
    "processed-data", "13_tripod_trios", "09_gather_trios",
    "filtered_trios_fine.parquet"
)
risk_path = here(
    "processed-data", "13_tripod_trios", "18_cytoscape_prep", "risk_genes.csv"
)
prep_path = here(
    "processed-data", "13_tripod_trios", "02_tripod_preprocess",
    "preprocessed_objects_%s.qs2"
)
cell_map_path = here("raw-data", "cell_type_map.csv")
pdf_dir = here("plots", "13_tripod_trios", "18_cytoscape_prep", mode)
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

risk_df = read_csv(risk_path, show_col_types = FALSE)
risk_flags = risk_df |>
    distinct(trait, gene) |>
    mutate(
        trait = case_match(
            trait,
            "MDD" ~ "MDD",
            "substance" ~ "substance",
            .default = trait
        ),
        value = TRUE
    ) |>
    pivot_wider(
        id_cols = gene,
        names_from = trait,
        values_from = value,
        values_fill = FALSE,
        names_prefix = "risk_"
    ) |>
    mutate(
        risk_MDD = coalesce(risk_MDD, FALSE),
        risk_substance = coalesce(risk_substance, FALSE),
        is_risk_gene = risk_MDD | risk_substance,
        risk_trait = case_when(
            risk_MDD & risk_substance ~ "MDD + substance",
            risk_MDD ~ "MDD",
            risk_substance ~ "substance",
            TRUE ~ "none"
        )
    )

if (mode == 'intersect') {
    tf_targets = read_parquet_duckdb(trio_path, prudence = "stingy") |>
        filter(is_intersect, stringency_level == 2) |>
        select(gene, TF, cell_type, coef, adj) |>
        collect()
} else {
    tf_targets = read_parquet_duckdb(trio_path, prudence = "stingy") |>
        filter(stringency_level == 2) |>
        select(gene, TF, cell_type, coef, adj) |>
        collect() |>
        group_by(cell_type) |>
        arrange(adj) |>
        slice_head(n = 200) |>
        ungroup()
}

tf_targets = tf_targets |>
    mutate(
        coef = log(coef) - min(log(coef)),
        neg_log10_adj = -log10(adj),
        source_cell_type = cell_type,
        cell_type = dplyr::coalesce(
            unname(rename_map[cell_type]),
            cell_type
        )
    )

get_cell_type_expression = function(cell_type, genes) {
    mats = qs_read(sprintf(prep_path, cell_type))$metacell_seur
    genes_present = intersect(genes, colnames(mats$rna))

    tibble(
        source_cell_type = cell_type,
        id = genes_present,
        mean_expr = colMeans(mats$rna[, genes_present, drop = FALSE])
    )
}

make_node_expression = function(tf_targets) {
    genes_by_ct = tf_targets |>
        transmute(source_cell_type, gene, TF) |>
        pivot_longer(c(gene, TF), values_to = "id") |>
        distinct(source_cell_type, id)

    expr_by_ct = imap_dfr(
        split(genes_by_ct$id, genes_by_ct$source_cell_type),
        \(genes, cell_type) get_cell_type_expression(cell_type, genes)
    )

    expr_by_ct |>
        group_by(id) |>
        summarise(mean_expr = mean(mean_expr), .groups = "drop")
}

make_node_sizes = function(mean_expr, default_size = 48, size_range = c(30, 70)) {
    log_expr = log1p(mean_expr)
    expr_range = range(log_expr, na.rm = TRUE)

    if (!all(is.finite(expr_range)) || diff(expr_range) == 0) {
        return(rep(default_size, length(mean_expr)))
    }

    scales::rescale(log_expr, to = size_range, from = expr_range) |>
        replace_na(default_size)
}

build_nodes = function(tf_targets) {
    node_expression = make_node_expression(tf_targets)

    bind_rows(
        tf_targets |>
            distinct(id = TF) |>
            mutate(role = "TF"),
        tf_targets |>
            distinct(id = gene) |>
            mutate(role = "Gene")
    ) |>
        distinct(id, .keep_all = TRUE) |>
        left_join(risk_flags, by = c("id" = "gene")) |>
        left_join(node_expression, by = "id") |>
        mutate(
            across(
                c(risk_MDD, risk_substance, is_risk_gene),
                \(x) replace_na(x, FALSE)
            ),
            risk_trait = replace_na(risk_trait, "none"),
            node_fill_group = paste(
                role,
                if_else(is_risk_gene, "risk", "nonrisk"),
                sep = "_"
            ),
            node_size = make_node_sizes(mean_expr)
        ) |>
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

map_shared_node_style = function(style_name, network_suid, nodes) {
    setNodeLabelMapping("id", style.name = style_name, network = network_suid)
    setNodeColorMapping(
        table.column = "node_fill_group",
        table.column.values = c("TF_nonrisk", "TF_risk", "Gene_nonrisk", "Gene_risk"),
        colors = c("#E69F00", "#A35F00", "#56B4E9", "#2077A1"),
        mapping.type = "d",
        style.name = style_name,
        network = network_suid
    )
    setNodeShapeMapping(
        table.column = "risk_trait",
        table.column.values = c("none", "MDD", "substance", "MDD + substance"),
        shapes = c("ELLIPSE", "TRIANGLE", "DIAMOND", "HEXAGON"),
        style.name = style_name,
        network = network_suid
    )
    setNodeSizeMapping(
        table.column = "node_size",
        table.column.values = range(nodes$node_size, na.rm = TRUE),
        sizes = range(nodes$node_size, na.rm = TRUE),
        mapping.type = "c",
        default.size = 48,
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
    map_shared_node_style(style_name, network_suid, nodes)
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
    filter(str_detect(cell_type, "Hb")) |>
    arrange(adj) |>
    slice_head(n = 200)
hb_edges = build_edges(hb_targets)
hb_nodes = build_nodes(hb_targets)
focused_hb_targets = read_parquet_duckdb(trio_path, prudence = "stingy") |>
    filter(stringency_level == 2) |>
    select(gene, TF, cell_type, coef, adj) |>
    collect() |>
    mutate(
        coef = log(coef) - min(log(coef)),
        neg_log10_adj = -log10(adj),
        source_cell_type = cell_type,
        cell_type = dplyr::coalesce(
            unname(rename_map[cell_type]),
            cell_type
        )
    ) |>
    filter(str_detect(cell_type, "Hb"))
coef_range = range(overview_edges$coef, na.rm = TRUE)
neg_log10_adj_range = range(overview_edges$neg_log10_adj, na.rm = TRUE)

create_focused_hb_network = function(focus_gene) {
    focus_targets = focused_hb_targets |>
        filter(gene == focus_gene | TF == focus_gene)

    if (nrow(focus_targets) == 0) {
        message(sprintf(
            "No focused Hb level-2 interactions found for %s; skipping network",
            focus_gene
        ))
        return(tibble(
            network_suid = NA_integer_,
            title = sprintf("TF target network - Hb focused - %s", focus_gene),
            style_name = NA_character_,
            pdf_path = NA_character_
        ))
    }

    create_and_export_network(
        nodes = build_nodes(focus_targets),
        edges = build_edges(focus_targets),
        title = sprintf("TF target network - Hb focused - %s", focus_gene),
        edge_color_mapper = map_edge_color_by_cell_type,
        edge_color_args = list(edges = build_edges(focus_targets)),
        coef_range = coef_range
    )
}

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
    create_focused_hb_network("CACNA2D1"),
    create_focused_hb_network("NRXN3"),
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
