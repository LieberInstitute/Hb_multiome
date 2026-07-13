library(shiny)
library(bslib)
library(Seurat)
library(ggplot2)
library(thematic)
library(DT)

validate_color_vars <- function(seur, color_vars) {
  if (!is.character(color_vars) || length(color_vars) == 0) {
    stop("`color_vars` must be a non-empty character vector.")
  }

  meta_vars <- colnames(seur[[]])
  available_color_vars <- intersect(color_vars, meta_vars)

  if (length(available_color_vars) == 0) {
    stop("None of the requested `color_vars` were found in the Seurat metadata.")
  }

  missing_color_vars <- setdiff(color_vars, meta_vars)

  list(
    available = available_color_vars,
    missing = missing_color_vars
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

build_atlas_panel <- function(reduction_choices, selected_reduction, color_choices, missing_color_vars) {
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
        selected = color_choices[[1]]
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
      textOutput("feature_selection_text")
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
          "Manual selection:",
          selected_gene() %||% "<gene>",
          "and",
          selected_peak() %||% "<peak>"
        )
      } else {
        req(input$trio_table_rows_selected)

        trio_row <- selected_trio_row()
        trio_label <- paste(trio_row$gene[[1]], trio_row$peak[[1]], sep = " | ")

        extra_bits <- intersect(c("TF", "cell_type", "trio_cor", "trio_p_adj"), colnames(trio_row))
        extra_text <- paste(
          vapply(extra_bits, function(col) paste0(col, ": ", trio_row[[col]][[1]]), character(1)),
          collapse = " | "
        )

        paste("Trio selection:", trio_label, if (nzchar(extra_text)) paste0(" | ", extra_text) else "")
      }
    })

    output$gene_plot_title <- renderText({
      paste("RNA:", selected_gene() %||% "No gene selected")
    })

    output$peak_plot_title <- renderText({
      paste("ATAC:", selected_peak() %||% "No peak selected")
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
      )
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
  default_reduction = NULL
) {
  color_info <- validate_color_vars(atlas_seur, color_vars)
  reduction_info <- validate_reductions(atlas_seur, default_reduction = default_reduction)
  trio_info <- validate_metacell_inputs(metacell_seur, trio_df)

  app_ui <- build_app_ui(
    reduction_choices = reduction_info$available,
    selected_reduction = reduction_info$selected,
    color_choices = color_info$available,
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
