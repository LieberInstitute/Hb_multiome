library(shiny)
library(bslib)
library(Seurat)
library(ggplot2)
library(thematic)
library(DT)

validate_color_vars <- function(seur, color_vars, cell_type_var = NULL) {
  if (!is.character(color_vars) || length(color_vars) == 0) {
    stop("`color_vars` must be a non-empty character vector.")
  }

  requested_color_vars <- unique(c(cell_type_var, color_vars))
  meta_vars <- colnames(seur[[]])
  available_color_vars <- intersect(requested_color_vars, meta_vars)

  if (length(available_color_vars) == 0) {
    stop("None of the requested `color_vars` were found in the Seurat metadata.")
  }

  missing_color_vars <- setdiff(requested_color_vars, meta_vars)

  selected_color_by <- cell_type_var
  if (is.null(selected_color_by) || !selected_color_by %in% available_color_vars) {
    selected_color_by <- available_color_vars[[1]]
  }

  list(
    available = available_color_vars,
    missing = missing_color_vars,
    selected = selected_color_by
  )
}

get_matching_colors <- function(values, color_map) {
  if (is.null(color_map) || length(color_map) == 0) {
    return(NULL)
  }

  present_levels <- unique(as.character(values))
  matched_colors <- color_map[names(color_map) %in% present_levels]

  if (length(matched_colors) == 0) {
    return(NULL)
  }

  matched_colors
}

validate_reductions <- function(seur, default_reduction = NULL) {
  available_reductions <- names(seur@reductions)

  if (length(available_reductions) == 0) {
    stop("The Seurat object has no reductions available for plotting.")
  }

  selected_reduction <- default_reduction

  if (is.null(selected_reduction) || !selected_reduction %in% available_reductions) {
    selected_reduction <- if ("wnn_umap" %in% available_reductions) {
      "wnn_umap"
    } else {
      available_reductions[[1]]
    }
  }

  list(
    available = available_reductions,
    selected = selected_reduction
  )
}

validate_metacell_inputs <- function(metacell_seur, trio_df, cell_type_var) {
  if (!cell_type_var %in% colnames(metacell_seur[[]])) {
    stop(paste0("`metacell_seur` must contain the `", cell_type_var, "` metadata column."))
  }

  required_trio_cols <- c("peak", "gene", cell_type_var)
  missing_trio_cols <- setdiff(required_trio_cols, colnames(trio_df))

  if (length(missing_trio_cols) > 0) {
    stop(
      paste(
        "`trio_df` is missing required columns:",
        paste(missing_trio_cols, collapse = ", ")
      )
    )
  }

  valid_rows <- trio_df$gene %in% rownames(metacell_seur[["RNA"]]) &
    trio_df$peak %in% rownames(metacell_seur[["ATAC"]])

  filtered_trio_df <- trio_df[valid_rows, , drop = FALSE]
  filtered_trio_df$trio_id <- seq_len(nrow(filtered_trio_df))

  list(
    trio_df = filtered_trio_df,
    dropped_rows = sum(!valid_rows)
  )
}

build_atlas_panel <- function(reduction_choices, selected_reduction, color_choices, selected_color_by, missing_color_vars) {
  tagList(
    p(
      "Cell-level embeddings colored by cell-type label. By default, the WNN (weighted nearest-neighbor) UMAP, which was built on a consensus of the RNA and ATAC data, is shown for all cells in the dataset."
    ),
    layout_sidebar(
      sidebar = sidebar(
      selectInput(
        inputId = "reduction",
        label = "Reduced dimension",
        choices = reduction_choices,
        selected = selected_reduction
      ),
      selectInput(
        inputId = "color_by",
        label = "Color by",
        choices = color_choices,
        selected = selected_color_by
      ),
      if (length(missing_color_vars) > 0) {
        helpText(
          paste(
            "Configured metadata fields not found in the atlas object:",
            paste(missing_color_vars, collapse = ", ")
          )
        )
      }
    ),
      card(
        full_screen = TRUE,
        card_header(textOutput("plot_title")),
        plotOutput("dim_plot", height = "700px")
      )
    )
  )
}

build_metacell_panel <- function(trio_df, dropped_trio_rows) {
  tagList(
    p(
      "Gene expression, peak accessibility, and optionally transcription factor expression can be explored at the metacell level in this tab. Metacells are aggregates of many individual cells of the same cell type with similar feature expression, used by the TRIPOD method to find gene-peak-TF trios."
    ),
    layout_sidebar(
      sidebar = sidebar(
      radioButtons(
        inputId = "feature_mode",
        label = "Feature selection mode",
        choices = c("Manual feature search" = "manual", "Trio table" = "trio"),
        selected = "trio"
      ),
      conditionalPanel(
        condition = "input.feature_mode === 'manual'",
        selectizeInput(
          inputId = "manual_gene",
          label = "Gene",
          choices = NULL,
          selected = NULL,
          multiple = FALSE
        ),
        selectizeInput(
          inputId = "manual_peak",
          label = "Peak",
          choices = NULL,
          selected = NULL,
          multiple = FALSE
        )
      ),
      conditionalPanel(
        condition = "input.feature_mode === 'trio'",
        helpText("Select one trio row to drive both violin plots.")
      ),
      if (dropped_trio_rows > 0) {
        helpText(sprintf("Dropped %d trio rows with missing gene or peak features.", dropped_trio_rows))
      },
      verbatimTextOutput("feature_selection_text")
    ),
    card(
      full_screen = TRUE,
      card_header("Selected features across cell types"),
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header(textOutput("gene_plot_title")),
          plotOutput("gene_vln_plot", height = "500px")
        ),
        card(
          card_header(textOutput("peak_plot_title")),
          plotOutput("peak_vln_plot", height = "500px")
        )
      )
    ),
    conditionalPanel(
      condition = "input.feature_mode === 'trio'",
      card(
        full_screen = TRUE,
        card_header(textOutput("trio_scatter_title")),
        plotOutput("trio_scatter_plot", height = "550px")
      )
    ),
      card(
        full_screen = TRUE,
        card_header("Trio table"),
        DTOutput("trio_table")
      )
    )
  )
}

build_dar_panel <- function() {
  tagList(
    p(
      "In this tab, differentially accessible regions (DARs), which are ATAC peaks with higher or lower accessibility in one cell type (against all others), can be explored. DARs were computed using the cell-level ATAC data, but for visualization we show data pseudobulked by cell type."
    ),
    layout_sidebar(
      sidebar = sidebar(
        helpText("Select one DAR row to plot its peak accessibility across cell types.")
      ),
      card(
        full_screen = TRUE,
        card_header(textOutput("dar_peak_plot_title")),
        plotOutput("dar_peak_vln_plot", height = "550px")
      ),
      card(
        full_screen = TRUE,
        card_header("DAR table"),
        DTOutput("dar_table")
      )
    )
  )
}

build_app_ui <- function(
  reduction_choices,
  selected_reduction,
  color_choices,
  selected_color_by,
  missing_color_vars,
  trio_df,
  dropped_trio_rows
) {
  page_navbar(
    title = "Habenula Atlas Multiome",
    theme = bs_theme(version = 5),
    nav_panel(
      "Cell-Level Embeddings",
      build_atlas_panel(
        reduction_choices = reduction_choices,
        selected_reduction = selected_reduction,
        color_choices = color_choices,
        selected_color_by = selected_color_by,
        missing_color_vars = missing_color_vars
      )
    ),
    nav_panel(
      "Metacell Features and Trios",
      build_metacell_panel(trio_df = trio_df, dropped_trio_rows = dropped_trio_rows)
    ),
    nav_panel(
      "DAR Exploration",
      build_dar_panel()
    )
  )
}

build_app_server <- function(atlas_seur, metacell_seur, trio_df, dar_df, seur_pb, cell_type_var, cell_type_colors, default_gene = NULL, default_peak = NULL) {
  force(atlas_seur)
  force(metacell_seur)
  force(trio_df)
  force(dar_df)
  force(seur_pb)
  force(cell_type_var)
  force(cell_type_colors)

  rna_assay <- metacell_seur[["RNA"]]
  atac_assay <- metacell_seur[["ATAC"]]
  pb_atac_assay <- seur_pb[["ATAC"]]
  rna_layer <- "data"
  atac_layer <- "data"
  pb_atac_layer <- "data"

  gene_choices <- rownames(rna_assay)
  peak_choices <- rownames(atac_assay)
  pb_peak_choices <- rownames(pb_atac_assay)
  metacell_meta <- metacell_seur[[]]
  metacell_group_var <- cell_type_var
  metacell_cell_type_colors <- get_matching_colors(metacell_meta[[metacell_group_var]], cell_type_colors)
  pb_meta <- seur_pb[[]]
  pb_group_var <- cell_type_var
  pb_cell_type_colors <- get_matching_colors(pb_meta[[pb_group_var]], cell_type_colors)

  function(input, output, session) {
    thematic_shiny()

    updateSelectizeInput(
      session = session,
      inputId = "manual_gene",
      choices = gene_choices,
      selected = if (!is.null(default_gene) && default_gene %in% gene_choices) default_gene else gene_choices[[1]],
      server = TRUE
    )
    updateSelectizeInput(
      session = session,
      inputId = "manual_peak",
      choices = peak_choices,
      selected = if (!is.null(default_peak) && default_peak %in% peak_choices) default_peak else peak_choices[[1]],
      server = TRUE
    )

    output$plot_title <- renderText({
      paste(input$reduction, "colored by", input$color_by)
    })

    output$dim_plot <- renderPlot({
      req(input$reduction, input$color_by)

      validate(
        need(input$reduction %in% names(atlas_seur@reductions), "Selected reduction is not available."),
        need(input$color_by %in% colnames(atlas_seur[[]]), "Selected metadata field is not available.")
      )

      reduction_mat <- Embeddings(atlas_seur[[input$reduction]])

      validate(
        need(ncol(reduction_mat) >= 2, "Selected reduction has fewer than two dimensions.")
      )

      dim_plot <- DimPlot(
        object = atlas_seur,
        reduction = input$reduction,
        group.by = input$color_by,
        raster = TRUE
      )

      if (identical(input$color_by, cell_type_var)) {
        atlas_cell_type_colors <- get_matching_colors(atlas_seur[[]][[cell_type_var]], cell_type_colors)

        if (!is.null(atlas_cell_type_colors)) {
          dim_plot <- dim_plot + scale_color_manual(values = atlas_cell_type_colors)
        }
      }

      dim_plot
    }, res = 110)

    selected_trio_row <- reactive({
      req(input$trio_table_rows_selected)
      trio_df[input$trio_table_rows_selected, , drop = FALSE]
    })

    trio_scatter_df <- reactive({
      req(identical(input$feature_mode, "trio"))
      trio_row <- selected_trio_row()

      validate(
        need("TF" %in% colnames(trio_row), "Selected trio row does not contain a TF column."),
        need(cell_type_var %in% colnames(trio_row), paste0("Selected trio row does not contain the `", cell_type_var, "` column."))
      )

      trio_gene <- trio_row$gene[[1]]
      trio_peak <- trio_row$peak[[1]]
      trio_tf <- trio_row$TF[[1]]
      trio_cell_type <- trio_row[[cell_type_var]][[1]]

      validate(
        need(trio_gene %in% gene_choices, "Selected trio gene is not available in the RNA assay."),
        need(trio_peak %in% peak_choices, "Selected trio peak is not available in the ATAC assay."),
        need(trio_tf %in% gene_choices, "Selected trio TF is not available in the RNA assay.")
      )

      cell_idx <- rownames(metacell_meta)[metacell_meta[[cell_type_var]] == trio_cell_type]

      validate(
        need(length(cell_idx) > 0, "No metacells matched the trio-selected cell type.")
      )

      plot_df <- data.frame(
        cell = cell_idx,
        gene_expr = as.numeric(LayerData(rna_assay, layer = rna_layer, features = trio_gene, cells = cell_idx)),
        peak_expr = as.numeric(LayerData(atac_assay, layer = atac_layer, features = trio_peak, cells = cell_idx)),
        tf_expr = as.numeric(LayerData(rna_assay, layer = rna_layer, features = trio_tf, cells = cell_idx)),
        plot_group = trio_cell_type,
        gene = trio_gene,
        peak = trio_peak,
        TF = trio_tf
      )

      plot_df <- subset(plot_df, tf_expr > 0)

      validate(
        need(nrow(plot_df) > 0, "No metacells had positive TF expression for the selected trio."),
        need(sum(stats::complete.cases(plot_df[, c("gene_expr", "peak_expr")])) >= 2, "Not enough metacells remained to compute the scatter plot.")
      )

      plot_df
    })

    selected_gene <- reactive({
      if (identical(input$feature_mode, "manual")) {
        return(input$manual_gene)
      }

      selected_trio_row()$gene[[1]]
    })

    selected_peak <- reactive({
      if (identical(input$feature_mode, "manual")) {
        return(input$manual_peak)
      }

      selected_trio_row()$peak[[1]]
    })

    output$feature_selection_text <- renderText({
      if (identical(input$feature_mode, "manual")) {
        paste(
          c(
            "Manual selection:",
            paste0("Gene: ", selected_gene() %||% "<gene>"),
            paste0("Peak: ", selected_peak() %||% "<peak>")
          ),
          collapse = "\n"
        )
      } else {
        req(input$trio_table_rows_selected)

        trio_row <- selected_trio_row()
        trio_lines <- c(
          "Trio selection:",
          paste0("Gene: ", trio_row$gene[[1]]),
          paste0("Peak: ", trio_row$peak[[1]])
        )

        if (cell_type_var %in% colnames(trio_row)) {
          trio_lines <- c(trio_lines, paste0(cell_type_var, ": ", trio_row[[cell_type_var]][[1]]))
        }

        if ("TF" %in% colnames(trio_row)) {
          trio_lines <- c(trio_lines, paste0("TF: ", trio_row$TF[[1]]))
        }

        paste(trio_lines, collapse = "\n")
      }
    })

    output$gene_plot_title <- renderText({
      paste("RNA:", selected_gene() %||% "No gene selected")
    })

    output$peak_plot_title <- renderText({
      paste("ATAC:", selected_peak() %||% "No peak selected")
    })

    output$trio_scatter_title <- renderText({
      req(identical(input$feature_mode, "trio"))
      trio_row <- selected_trio_row()
      paste0(
        "RNA vs ATAC in ",
        trio_row[[cell_type_var]][[1]],
        " metacells, colored by TF expression (",
        trio_row$TF[[1]],
        ")"
      )
    })

    output$gene_vln_plot <- renderPlot({
      req(selected_gene())

      validate(
        need(selected_gene() %in% gene_choices, "Selected gene is not available in the RNA assay.")
      )

      gene_plot <- VlnPlot(
        object = metacell_seur,
        features = selected_gene(),
        assay = "RNA",
        group.by = metacell_group_var
      ) +
        guides(fill = "none") +
        labs(x = "Cell Type", y = "Expression", title = NULL) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

      if (!is.null(metacell_cell_type_colors)) {
        gene_plot <- gene_plot + scale_fill_manual(values = metacell_cell_type_colors)
      }

      gene_plot
    }, res = 110)

    output$peak_vln_plot <- renderPlot({
      req(selected_peak())

      validate(
        need(selected_peak() %in% peak_choices, "Selected peak is not available in the ATAC assay.")
      )

      peak_plot <- VlnPlot(
        object = metacell_seur,
        features = selected_peak(),
        assay = "ATAC",
        group.by = metacell_group_var
      ) +
        guides(fill = "none") +
        labs(x = "Cell Type", y = "Accessibility", title = NULL) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

      if (!is.null(metacell_cell_type_colors)) {
        peak_plot <- peak_plot + scale_fill_manual(values = metacell_cell_type_colors)
      }

      peak_plot
    }, res = 110)

    output$trio_scatter_plot <- renderPlot({
      plot_df <- trio_scatter_df()

      global_cor <- stats::cor(
        plot_df$gene_expr,
        plot_df$peak_expr,
        use = "complete.obs"
      )

      x_pos <- min(plot_df$gene_expr, na.rm = TRUE)
      y_pos <- max(plot_df$peak_expr, na.rm = TRUE)

      ggplot(
        plot_df,
        aes(x = gene_expr, y = peak_expr, color = tf_expr)
      ) +
        geom_point(size = 2.5) +
        annotate(
          geom = "text",
          x = x_pos,
          y = y_pos,
          label = sprintf("Global cor: %.3f", global_cor),
          hjust = 0,
          vjust = 1
        ) +
        scale_color_viridis_c() +
        labs(
          x = paste0("RNA: ", unique(plot_df$gene)),
          y = paste0("ATAC: ", unique(plot_df$peak)),
          color = paste0("TF: ", unique(plot_df$TF))
        )
    }, res = 110)

    output$trio_table <- renderDT({
      datatable(
        trio_df,
        rownames = FALSE,
        filter = "top",
        selection = list(mode = "single", selected = 1, target = "row"),
        options = list(
          pageLength = 15,
          lengthMenu = c(15, 30, 50, 100),
          scrollX = TRUE
        )
      )
    })

    selected_dar_row <- reactive({
      req(input$dar_table_rows_selected)
      dar_df[input$dar_table_rows_selected, , drop = FALSE]
    })

    output$dar_peak_plot_title <- renderText({
      req(input$dar_table_rows_selected)
      paste("Peak accessibility:", selected_dar_row()$peak[[1]])
    })

    output$dar_peak_vln_plot <- renderPlot({
      req(input$dar_table_rows_selected)
      selected_peak <- selected_dar_row()$peak[[1]]

      validate(
        need(selected_peak %in% pb_peak_choices, "Selected DAR peak is not available in the pseudobulk ATAC assay."),
        need(pb_group_var %in% colnames(pb_meta), paste0("The pseudobulk object is missing the `", pb_group_var, "` metadata column."))
      )

      dar_plot <- VlnPlot(
        object = seur_pb,
        features = selected_peak,
        assay = "ATAC",
        group.by = pb_group_var
      ) +
        guides(fill = "none") +
        labs(x = "Cell Type", y = "Accessibility", title = NULL) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

      if (!is.null(pb_cell_type_colors)) {
        dar_plot <- dar_plot + scale_fill_manual(values = pb_cell_type_colors)
      }

      dar_plot
    }, res = 110)

    output$dar_table <- renderDT({
      datatable(
        dar_df,
        rownames = FALSE,
        filter = "top",
        selection = list(mode = "single", selected = 1, target = "row"),
        options = list(
          pageLength = 15,
          lengthMenu = c(15, 30, 50, 100),
          scrollX = TRUE
        )
      )
    })
  }
}

run_app <- function(
  atlas_seur,
  color_vars,
  metacell_seur,
  trio_df,
  dar_df,
  seur_pb,
  cell_type_var,
  cell_type_colors = NULL,
  default_reduction = NULL,
  default_gene = NULL,
  default_peak = NULL
) {
  color_info <- validate_color_vars(
    atlas_seur,
    color_vars,
    cell_type_var = cell_type_var
  )
  reduction_info <- validate_reductions(atlas_seur, default_reduction = default_reduction)
  trio_info <- validate_metacell_inputs(
    metacell_seur,
    trio_df,
    cell_type_var = cell_type_var
  )

  app_ui <- build_app_ui(
    reduction_choices = reduction_info$available,
    selected_reduction = reduction_info$selected,
    color_choices = color_info$available,
    selected_color_by = color_info$selected,
    missing_color_vars = color_info$missing,
    trio_df = trio_info$trio_df,
    dropped_trio_rows = trio_info$dropped_rows
  )

  app_server <- build_app_server(
    atlas_seur = atlas_seur,
    metacell_seur = metacell_seur,
    trio_df = trio_info$trio_df,
    dar_df = dar_df,
    seur_pb = seur_pb,
    cell_type_var = cell_type_var,
    cell_type_colors = cell_type_colors,
    default_gene = default_gene,
    default_peak = default_peak
  )

  shinyApp(ui = app_ui, server = app_server)
}
