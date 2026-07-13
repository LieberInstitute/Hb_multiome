library(shiny)
library(bslib)
library(Seurat)
library(ggplot2)
library(thematic)
library(DT)

validate_color_vars <- function(seur, color_vars, default_color_by = NULL) {
  if (!is.character(color_vars) || length(color_vars) == 0) {
    stop("`color_vars` must be a non-empty character vector.")
  }

  meta_vars <- colnames(seur[[]])
  available_color_vars <- intersect(color_vars, meta_vars)

  if (length(available_color_vars) == 0) {
    stop("None of the requested `color_vars` were found in the Seurat metadata.")
  }

  missing_color_vars <- setdiff(color_vars, meta_vars)

  selected_color_by <- default_color_by
  if (is.null(selected_color_by) || !selected_color_by %in% available_color_vars) {
    selected_color_by <- available_color_vars[[1]]
  }

  list(
    available = available_color_vars,
    missing = missing_color_vars,
    selected = selected_color_by
  )
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

get_assay_matrix <- function(seur, assay, preferred_layers = c("data", "counts")) {
  assay_obj <- seur[[assay]]
  assay_layers <- Layers(assay_obj)
  selected_layer <- preferred_layers[preferred_layers %in% assay_layers][[1]]

  if (is.null(selected_layer) || is.na(selected_layer)) {
    stop(
      paste0(
        "Assay `", assay, "` is missing all preferred layers: ",
        paste(preferred_layers, collapse = ", ")
      )
    )
  }

  LayerData(assay_obj, layer = selected_layer)
}

validate_metacell_inputs <- function(metacell_seur, trio_df) {
  if (!"cell_type" %in% colnames(metacell_seur[[]])) {
    stop("`metacell_seur` must contain a `cell_type` metadata column.")
  }

  required_trio_cols <- c("peak", "gene")
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
}

build_metacell_panel <- function(trio_df, dropped_trio_rows) {
  layout_sidebar(
    sidebar = sidebar(
      radioButtons(
        inputId = "feature_mode",
        label = "Feature selection mode",
        choices = c("Manual feature search" = "manual", "Trio table" = "trio"),
        selected = "manual"
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
    title = "Habenula atlas viewer",
    theme = bs_theme(version = 5),
    nav_panel(
      "Atlas embeddings",
      build_atlas_panel(
        reduction_choices = reduction_choices,
        selected_reduction = selected_reduction,
        color_choices = color_choices,
        selected_color_by = selected_color_by,
        missing_color_vars = missing_color_vars
      )
    ),
    nav_panel(
      "Metacell features",
      build_metacell_panel(trio_df = trio_df, dropped_trio_rows = dropped_trio_rows)
    )
  )
}

build_app_server <- function(atlas_seur, metacell_seur, trio_df) {
  force(atlas_seur)
  force(metacell_seur)
  force(trio_df)

  gene_choices <- rownames(metacell_seur[["RNA"]])
  peak_choices <- rownames(metacell_seur[["ATAC"]])
  metacell_meta <- metacell_seur[[]]
  rna_plot_mat <- get_assay_matrix(metacell_seur, assay = "RNA", preferred_layers = c("data", "counts"))
  atac_plot_mat <- get_assay_matrix(metacell_seur, assay = "ATAC", preferred_layers = c("data", "counts"))

  function(input, output, session) {
    thematic_shiny()

    updateSelectizeInput(
      session = session,
      inputId = "manual_gene",
      choices = gene_choices,
      selected = gene_choices[[1]],
      server = TRUE
    )
    updateSelectizeInput(
      session = session,
      inputId = "manual_peak",
      choices = peak_choices,
      selected = peak_choices[[1]],
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

      DimPlot(
        object = atlas_seur,
        reduction = input$reduction,
        group.by = input$color_by,
        raster = TRUE
      )
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
        need("cell_type" %in% colnames(trio_row), "Selected trio row does not contain a cell_type column.")
      )

      trio_gene <- trio_row$gene[[1]]
      trio_peak <- trio_row$peak[[1]]
      trio_tf <- trio_row$TF[[1]]
      trio_cell_type <- trio_row$cell_type[[1]]

      validate(
        need(trio_gene %in% gene_choices, "Selected trio gene is not available in the RNA assay."),
        need(trio_peak %in% peak_choices, "Selected trio peak is not available in the ATAC assay."),
        need(trio_tf %in% gene_choices, "Selected trio TF is not available in the RNA assay.")
      )

      cell_idx <- rownames(metacell_meta)[metacell_meta$cell_type == trio_cell_type]

      validate(
        need(length(cell_idx) > 0, "No metacells matched the trio-selected cell type.")
      )

      plot_df <- data.frame(
        cell = cell_idx,
        gene_expr = as.numeric(rna_plot_mat[trio_gene, cell_idx]),
        peak_expr = as.numeric(atac_plot_mat[trio_peak, cell_idx]),
        tf_expr = as.numeric(rna_plot_mat[trio_tf, cell_idx]),
        cell_type = trio_cell_type,
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

        if ("cell_type" %in% colnames(trio_row)) {
          trio_lines <- c(trio_lines, paste0("Cell type: ", trio_row$cell_type[[1]]))
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
        trio_row$cell_type[[1]],
        " metacells, colored by log TF expression (",
        trio_row$TF[[1]],
        ")"
      )
    })

    output$gene_vln_plot <- renderPlot({
      req(selected_gene())

      validate(
        need(selected_gene() %in% gene_choices, "Selected gene is not available in the RNA assay.")
      )

      VlnPlot(
        object = metacell_seur,
        features = selected_gene(),
        assay = "RNA",
        group.by = "cell_type"
      ) +
        guides(fill = "none") +
        labs(title = NULL) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    }, res = 110)

    output$peak_vln_plot <- renderPlot({
      req(selected_peak())

      validate(
        need(selected_peak() %in% peak_choices, "Selected peak is not available in the ATAC assay.")
      )

      VlnPlot(
        object = metacell_seur,
        features = selected_peak(),
        assay = "ATAC",
        group.by = "cell_type"
      ) +
        guides(fill = "none") +
        labs(title = NULL) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
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
        aes(x = gene_expr, y = peak_expr, color = log1p(tf_expr))
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
          color = paste0("log1p RNA: ", unique(plot_df$TF))
        )
    }, res = 110)

    output$trio_table <- renderDT({
      datatable(
        trio_df,
        rownames = FALSE,
        filter = "top",
        selection = "single",
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
  default_reduction = NULL,
  default_color_by = NULL
) {
  color_info <- validate_color_vars(
    atlas_seur,
    color_vars,
    default_color_by = default_color_by
  )
  reduction_info <- validate_reductions(atlas_seur, default_reduction = default_reduction)
  trio_info <- validate_metacell_inputs(metacell_seur, trio_df)

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
    trio_df = trio_info$trio_df
  )

  shinyApp(ui = app_ui, server = app_server)
}
