library(shiny)
library(bslib)
library(Seurat)
library(ggplot2)
library(thematic)

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

build_app_ui <- function(reduction_choices, selected_reduction, color_choices, missing_color_vars) {
  page_sidebar(
    title = "Habenula atlas viewer",
    theme = bs_theme(version = 5),
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
            "Configured metadata fields not found in the Seurat object:",
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

build_app_server <- function(seur) {
  force(seur)

  function(input, output, session) {
    thematic_shiny()

    output$plot_title <- renderText({
      paste(input$reduction, "colored by", input$color_by)
    })

    output$dim_plot <- renderPlot({
      req(input$reduction, input$color_by)

      validate(
        need(input$reduction %in% names(seur@reductions), "Selected reduction is not available."),
        need(input$color_by %in% colnames(seur[[]]), "Selected metadata field is not available.")
      )

      reduction_mat <- Embeddings(seur[[input$reduction]])

      validate(
        need(ncol(reduction_mat) >= 2, "Selected reduction has fewer than two dimensions.")
      )

      DimPlot(
        object = seur,
        reduction = input$reduction,
        group.by = input$color_by,
        raster = TRUE
      )
    }, res = 110)
  }
}

run_app <- function(seur, color_vars, default_reduction = NULL) {
  color_info <- validate_color_vars(seur, color_vars)
  reduction_info <- validate_reductions(seur, default_reduction = default_reduction)

  app_ui <- build_app_ui(
    reduction_choices = reduction_info$available,
    selected_reduction = reduction_info$selected,
    color_choices = color_info$available,
    missing_color_vars = color_info$missing
  )

  app_server <- build_app_server(seur)

  shinyApp(ui = app_ui, server = app_server)
}
