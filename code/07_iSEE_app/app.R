library("SingleCellExperiment")
library("iSEE")
library("shiny")

## Load data
sce <- readRDS("sce_Habenula_iSEE.rds")
#sn_colors <- readRDS("sn_colors.rds")

initial <- list()

################################################################################
# Settings for Row data table 1
################################################################################

initial[["RowDataTable1"]] <- new("RowDataTable", Selected = "POU4F1", Search = "POU4", SearchColumns = "", 
                                  HiddenColumns = character(0), VersionInfo = list(iSEE = structure(list(
                                      c(2L, 18L, 0L)), class = c("package_version", "numeric_version"
                                      ))), PanelId = c(RowDataTable = 1L), PanelHeight = 600L, 
                                  PanelWidth = 4L, SelectionBoxOpen = FALSE, RowSelectionSource = "---", 
                                  ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, 
                                  ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE, 
                                  ColumnSelectionRestrict = FALSE, SelectionHistory = list())

################################################################################
# Settings for Feature assay plot 1
################################################################################

initial[["FeatureAssayPlot1"]] <- new("FeatureAssayPlot", Assay = "logcounts", XAxis = "Column data", 
                                      XAxisColumnData = "cluster_annotation", XAxisFeatureName = "MIR1302-2HG", 
                                      XAxisFeatureSource = "---", XAxisFeatureDynamicSource = FALSE, 
                                      YAxisFeatureName = "POU4F1", YAxisFeatureSource = "RowDataTable1", 
                                      YAxisFeatureDynamicSource = TRUE, FacetRowByColData = "orig.ident", 
                                      FacetColumnByColData = "orig.ident", ColorByColumnData = "seurat_clusters", 
                                      ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000", 
                                      ShapeByColumnData = "orig.ident", SizeByColumnData = "nCount_RNA", 
                                      TooltipColumnData = character(0), FacetRowBy = "None", FacetColumnBy = "None", 
                                      ColorBy = "Column data", ColorByDefaultColor = "#000000", 
                                      ColorByFeatureName = "MIR1302-2HG", ColorByFeatureSource = "---", 
                                      ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "S04_AAACAGCCAGAATGAC-1", 
                                      ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE, 
                                      ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1, 
                                      ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = FALSE, 
                                      VisualChoices = "Color", ContourAdd = FALSE, ContourColor = "#0000FF", 
                                      FixAspectRatio = FALSE, ViolinAdd = TRUE, PointSize = 1, 
                                      PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200, 
                                      CustomLabels = FALSE, CustomLabelsText = "S04_AAACAGCCAGAATGAC-1", 
                                      FontSize = 1, LegendPointSize = 1, LegendPosition = "Bottom", 
                                      HoverInfo = TRUE, LabelCenters = FALSE, LabelCentersBy = "orig.ident", 
                                      LabelCentersColor = "#000000", VersionInfo = list(iSEE = structure(list(
                                          c(2L, 18L, 0L)), class = c("package_version", "numeric_version"
                                          ))), PanelId = c(FeatureAssayPlot = 1L), PanelHeight = 400L, 
                                      PanelWidth = 8L, SelectionBoxOpen = FALSE, RowSelectionSource = "---", 
                                      ColumnSelectionSource = "FeatureAssayPlot3", DataBoxOpen = FALSE, 
                                      RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = TRUE, 
                                      RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
                                      SelectionHistory = list())

################################################################################
# Settings for Feature assay plot 3
################################################################################

initial[["FeatureAssayPlot3"]] <- new("FeatureAssayPlot", Assay = "logcounts", XAxis = "Column data", 
                                      XAxisColumnData = "merged_cluster", XAxisFeatureName = "MIR1302-2HG", 
                                      XAxisFeatureSource = "---", XAxisFeatureDynamicSource = FALSE, 
                                      YAxisFeatureName = "POU4F1", YAxisFeatureSource = "RowDataTable1", 
                                      YAxisFeatureDynamicSource = TRUE, FacetRowByColData = "orig.ident", 
                                      FacetColumnByColData = "orig.ident", ColorByColumnData = "orig.ident", 
                                      ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000", 
                                      ShapeByColumnData = "orig.ident", SizeByColumnData = "nCount_RNA", 
                                      TooltipColumnData = character(0), FacetRowBy = "None", FacetColumnBy = "None", 
                                      ColorBy = "Column data", ColorByDefaultColor = "#000000", 
                                      ColorByFeatureName = "MIR1302-2HG", ColorByFeatureSource = "---", 
                                      ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "S04_AAACAGCCAGAATGAC-1", 
                                      ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE, 
                                      ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1, 
                                      ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = FALSE, 
                                      VisualChoices = "Color", ContourAdd = FALSE, ContourColor = "#0000FF", 
                                      FixAspectRatio = FALSE, ViolinAdd = TRUE, PointSize = 1, 
                                      PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200, 
                                      CustomLabels = FALSE, CustomLabelsText = "S04_AAACAGCCAGAATGAC-1", 
                                      FontSize = 1, LegendPointSize = 1, LegendPosition = "Bottom", 
                                      HoverInfo = TRUE, LabelCenters = FALSE, LabelCentersBy = "orig.ident", 
                                      LabelCentersColor = "#000000", VersionInfo = list(iSEE = structure(list(
                                          c(2L, 18L, 0L)), class = c("package_version", "numeric_version"
                                          ))), PanelId = 3L, PanelHeight = 400L, PanelWidth = 6L, SelectionBoxOpen = FALSE, 
                                      RowSelectionSource = "---", ColumnSelectionSource = "---", 
                                      DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = TRUE, 
                                      RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
                                      SelectionHistory = list())

################################################################################
# Settings for Feature assay plot 4
################################################################################

initial[["FeatureAssayPlot4"]] <- new("FeatureAssayPlot", Assay = "logcounts", XAxis = "Column data", 
                                      XAxisColumnData = "orig.ident", XAxisFeatureName = "MIR1302-2HG", 
                                      XAxisFeatureSource = "---", XAxisFeatureDynamicSource = FALSE, 
                                      YAxisFeatureName = "POU4F1", YAxisFeatureSource = "RowDataTable1", 
                                      YAxisFeatureDynamicSource = TRUE, FacetRowByColData = "orig.ident", 
                                      FacetColumnByColData = "orig.ident", ColorByColumnData = "orig.ident", 
                                      ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000", 
                                      ShapeByColumnData = "orig.ident", SizeByColumnData = "nCount_RNA", 
                                      TooltipColumnData = character(0), FacetRowBy = "None", FacetColumnBy = "None", 
                                      ColorBy = "Column data", ColorByDefaultColor = "#000000", 
                                      ColorByFeatureName = "MIR1302-2HG", ColorByFeatureSource = "---", 
                                      ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "S04_AAACAGCCAGAATGAC-1", 
                                      ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE, 
                                      ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1, 
                                      ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = FALSE, 
                                      VisualChoices = "Color", ContourAdd = FALSE, ContourColor = "#0000FF", 
                                      FixAspectRatio = FALSE, ViolinAdd = TRUE, PointSize = 1, 
                                      PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200, 
                                      CustomLabels = FALSE, CustomLabelsText = "S04_AAACAGCCAGAATGAC-1", 
                                      FontSize = 1, LegendPointSize = 1, LegendPosition = "Bottom", 
                                      HoverInfo = TRUE, LabelCenters = FALSE, LabelCentersBy = "orig.ident", 
                                      LabelCentersColor = "#000000", VersionInfo = list(iSEE = structure(list(
                                          c(2L, 18L, 0L)), class = c("package_version", "numeric_version"
                                          ))), PanelId = 4L, PanelHeight = 400L, PanelWidth = 6L, SelectionBoxOpen = FALSE, 
                                      RowSelectionSource = "---", ColumnSelectionSource = "FeatureAssayPlot1", 
                                      DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = TRUE, 
                                      RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
                                      SelectionHistory = list())

################################################################################
# Settings for Column data plot 1
################################################################################

initial[["ColumnDataPlot1"]] <- new("ColumnDataPlot", XAxis = "Column data", YAxis = "nCount_RNA", 
                                    XAxisColumnData = "orig.ident", FacetRowByColData = "orig.ident", 
                                    FacetColumnByColData = "orig.ident", ColorByColumnData = "orig.ident", 
                                    ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000", 
                                    ShapeByColumnData = "orig.ident", SizeByColumnData = "nCount_RNA", 
                                    TooltipColumnData = character(0), FacetRowBy = "None", FacetColumnBy = "None", 
                                    ColorBy = "Column data", ColorByDefaultColor = "#000000", 
                                    ColorByFeatureName = "MIR1302-2HG", ColorByFeatureSource = "---", 
                                    ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "S04_AAACAGCCAGAATGAC-1", 
                                    ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE, 
                                    ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1, 
                                    ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = FALSE, 
                                    VisualChoices = "Color", ContourAdd = FALSE, ContourColor = "#0000FF", 
                                    FixAspectRatio = FALSE, ViolinAdd = TRUE, PointSize = 1, 
                                    PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200, 
                                    CustomLabels = FALSE, CustomLabelsText = "S04_AAACAGCCAGAATGAC-1", 
                                    FontSize = 1, LegendPointSize = 1, LegendPosition = "Bottom", 
                                    HoverInfo = TRUE, LabelCenters = FALSE, LabelCentersBy = "orig.ident", 
                                    LabelCentersColor = "#000000", VersionInfo = list(iSEE = structure(list(
                                        c(2L, 18L, 0L)), class = c("package_version", "numeric_version"
                                        ))), PanelId = 1L, PanelHeight = 400L, PanelWidth = 4L, SelectionBoxOpen = FALSE, 
                                    RowSelectionSource = "---", ColumnSelectionSource = "---", 
                                    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
                                    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
                                    SelectionHistory = list())

################################################################################
# Settings for Column data plot 6
################################################################################

initial[["ColumnDataPlot6"]] <- new("ColumnDataPlot", XAxis = "Column data", YAxis = "nFeature_RNA", 
                                    XAxisColumnData = "orig.ident", FacetRowByColData = "orig.ident", 
                                    FacetColumnByColData = "orig.ident", ColorByColumnData = "orig.ident", 
                                    ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000", 
                                    ShapeByColumnData = "orig.ident", SizeByColumnData = "nCount_RNA", 
                                    TooltipColumnData = character(0), FacetRowBy = "None", FacetColumnBy = "None", 
                                    ColorBy = "Column data", ColorByDefaultColor = "#000000", 
                                    ColorByFeatureName = "MIR1302-2HG", ColorByFeatureSource = "---", 
                                    ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "S04_AAACAGCCAGAATGAC-1", 
                                    ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE, 
                                    ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1, 
                                    ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = FALSE, 
                                    VisualChoices = "Color", ContourAdd = FALSE, ContourColor = "#0000FF", 
                                    FixAspectRatio = FALSE, ViolinAdd = TRUE, PointSize = 1, 
                                    PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200, 
                                    CustomLabels = FALSE, CustomLabelsText = "S04_AAACAGCCAGAATGAC-1", 
                                    FontSize = 1, LegendPointSize = 1, LegendPosition = "Bottom", 
                                    HoverInfo = TRUE, LabelCenters = FALSE, LabelCentersBy = "orig.ident", 
                                    LabelCentersColor = "#000000", VersionInfo = list(iSEE = structure(list(
                                        c(2L, 18L, 0L)), class = c("package_version", "numeric_version"
                                        ))), PanelId = 6L, PanelHeight = 400L, PanelWidth = 4L, SelectionBoxOpen = FALSE, 
                                    RowSelectionSource = "---", ColumnSelectionSource = "---", 
                                    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
                                    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
                                    SelectionHistory = list())

################################################################################
# Settings for Column data plot 7
################################################################################

initial[["ColumnDataPlot7"]] <- new("ColumnDataPlot", XAxis = "Column data", YAxis = "MTRatio", 
                                    XAxisColumnData = "orig.ident", FacetRowByColData = "orig.ident", 
                                    FacetColumnByColData = "orig.ident", ColorByColumnData = "orig.ident", 
                                    ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000", 
                                    ShapeByColumnData = "orig.ident", SizeByColumnData = "nCount_RNA", 
                                    TooltipColumnData = character(0), FacetRowBy = "None", FacetColumnBy = "None", 
                                    ColorBy = "Column data", ColorByDefaultColor = "#000000", 
                                    ColorByFeatureName = "MIR1302-2HG", ColorByFeatureSource = "---", 
                                    ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "S04_AAACAGCCAGAATGAC-1", 
                                    ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE, 
                                    ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1, 
                                    ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = FALSE, 
                                    VisualChoices = "Color", ContourAdd = FALSE, ContourColor = "#0000FF", 
                                    FixAspectRatio = FALSE, ViolinAdd = TRUE, PointSize = 1, 
                                    PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200, 
                                    CustomLabels = FALSE, CustomLabelsText = "S04_AAACAGCCAGAATGAC-1", 
                                    FontSize = 1, LegendPointSize = 1, LegendPosition = "Bottom", 
                                    HoverInfo = TRUE, LabelCenters = FALSE, LabelCentersBy = "orig.ident", 
                                    LabelCentersColor = "#000000", VersionInfo = list(iSEE = structure(list(
                                        c(2L, 18L, 0L)), class = c("package_version", "numeric_version"
                                        ))), PanelId = 7L, PanelHeight = 400L, PanelWidth = 4L, SelectionBoxOpen = FALSE, 
                                    RowSelectionSource = "---", ColumnSelectionSource = "---", 
                                    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
                                    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
                                    SelectionHistory = list())

################################################################################
# Settings for Column data plot 5
################################################################################

initial[["ColumnDataPlot5"]] <- new("ColumnDataPlot", XAxis = "Column data", YAxis = "nCount_ATAC", 
                                    XAxisColumnData = "orig.ident", FacetRowByColData = "orig.ident", 
                                    FacetColumnByColData = "orig.ident", ColorByColumnData = "orig.ident", 
                                    ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000", 
                                    ShapeByColumnData = "orig.ident", SizeByColumnData = "nCount_RNA", 
                                    TooltipColumnData = character(0), FacetRowBy = "None", FacetColumnBy = "None", 
                                    ColorBy = "Column data", ColorByDefaultColor = "#000000", 
                                    ColorByFeatureName = "MIR1302-2HG", ColorByFeatureSource = "---", 
                                    ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "S04_AAACAGCCAGAATGAC-1", 
                                    ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE, 
                                    ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1, 
                                    ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = FALSE, 
                                    VisualChoices = "Color", ContourAdd = FALSE, ContourColor = "#0000FF", 
                                    FixAspectRatio = FALSE, ViolinAdd = TRUE, PointSize = 1, 
                                    PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200, 
                                    CustomLabels = FALSE, CustomLabelsText = "S04_AAACAGCCAGAATGAC-1", 
                                    FontSize = 1, LegendPointSize = 1, LegendPosition = "Bottom", 
                                    HoverInfo = TRUE, LabelCenters = FALSE, LabelCentersBy = "orig.ident", 
                                    LabelCentersColor = "#000000", VersionInfo = list(iSEE = structure(list(
                                        c(2L, 18L, 0L)), class = c("package_version", "numeric_version"
                                        ))), PanelId = 5L, PanelHeight = 400L, PanelWidth = 3L, SelectionBoxOpen = FALSE, 
                                    RowSelectionSource = "---", ColumnSelectionSource = "---", 
                                    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
                                    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
                                    SelectionHistory = list())

################################################################################
# Settings for Column data plot 2
################################################################################

initial[["ColumnDataPlot2"]] <- new("ColumnDataPlot", XAxis = "Column data", YAxis = "nFeature_ATAC", 
                                    XAxisColumnData = "orig.ident", FacetRowByColData = "orig.ident", 
                                    FacetColumnByColData = "orig.ident", ColorByColumnData = "orig.ident", 
                                    ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000", 
                                    ShapeByColumnData = "orig.ident", SizeByColumnData = "nCount_RNA", 
                                    TooltipColumnData = character(0), FacetRowBy = "None", FacetColumnBy = "None", 
                                    ColorBy = "Column data", ColorByDefaultColor = "#000000", 
                                    ColorByFeatureName = "MIR1302-2HG", ColorByFeatureSource = "---", 
                                    ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "S04_AAACAGCCAGAATGAC-1", 
                                    ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE, 
                                    ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1, 
                                    ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = FALSE, 
                                    VisualChoices = "Color", ContourAdd = FALSE, ContourColor = "#0000FF", 
                                    FixAspectRatio = FALSE, ViolinAdd = TRUE, PointSize = 1, 
                                    PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200, 
                                    CustomLabels = FALSE, CustomLabelsText = "S04_AAACAGCCAGAATGAC-1", 
                                    FontSize = 1, LegendPointSize = 1, LegendPosition = "Bottom", 
                                    HoverInfo = TRUE, LabelCenters = FALSE, LabelCentersBy = "orig.ident", 
                                    LabelCentersColor = "#000000", VersionInfo = list(iSEE = structure(list(
                                        c(2L, 18L, 0L)), class = c("package_version", "numeric_version"
                                        ))), PanelId = 2L, PanelHeight = 400L, PanelWidth = 3L, SelectionBoxOpen = FALSE, 
                                    RowSelectionSource = "---", ColumnSelectionSource = "---", 
                                    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
                                    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
                                    SelectionHistory = list())

################################################################################
# Settings for Column data plot 3
################################################################################

initial[["ColumnDataPlot3"]] <- new("ColumnDataPlot", XAxis = "Column data", YAxis = "nucleosome_group", 
                                    XAxisColumnData = "orig.ident", FacetRowByColData = "orig.ident", 
                                    FacetColumnByColData = "orig.ident", ColorByColumnData = "orig.ident", 
                                    ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000", 
                                    ShapeByColumnData = "orig.ident", SizeByColumnData = "nCount_RNA", 
                                    TooltipColumnData = character(0), FacetRowBy = "None", FacetColumnBy = "None", 
                                    ColorBy = "Column data", ColorByDefaultColor = "#000000", 
                                    ColorByFeatureName = "MIR1302-2HG", ColorByFeatureSource = "---", 
                                    ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "S04_AAACAGCCAGAATGAC-1", 
                                    ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE, 
                                    ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1, 
                                    ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = FALSE, 
                                    VisualChoices = "Color", ContourAdd = FALSE, ContourColor = "#0000FF", 
                                    FixAspectRatio = FALSE, ViolinAdd = TRUE, PointSize = 1, 
                                    PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200, 
                                    CustomLabels = FALSE, CustomLabelsText = "S04_AAACAGCCAGAATGAC-1", 
                                    FontSize = 1, LegendPointSize = 1, LegendPosition = "Bottom", 
                                    HoverInfo = TRUE, LabelCenters = FALSE, LabelCentersBy = "orig.ident", 
                                    LabelCentersColor = "#000000", VersionInfo = list(iSEE = structure(list(
                                        c(2L, 18L, 0L)), class = c("package_version", "numeric_version"
                                        ))), PanelId = 3L, PanelHeight = 400L, PanelWidth = 3L, SelectionBoxOpen = FALSE, 
                                    RowSelectionSource = "---", ColumnSelectionSource = "---", 
                                    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
                                    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
                                    SelectionHistory = list())

################################################################################
# Settings for Column data plot 4
################################################################################

initial[["ColumnDataPlot4"]] <- new("ColumnDataPlot", XAxis = "Column data", YAxis = "high.tss", 
                                    XAxisColumnData = "orig.ident", FacetRowByColData = "orig.ident", 
                                    FacetColumnByColData = "orig.ident", ColorByColumnData = "orig.ident", 
                                    ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000", 
                                    ShapeByColumnData = "orig.ident", SizeByColumnData = "nCount_RNA", 
                                    TooltipColumnData = character(0), FacetRowBy = "None", FacetColumnBy = "None", 
                                    ColorBy = "Column data", ColorByDefaultColor = "#000000", 
                                    ColorByFeatureName = "MIR1302-2HG", ColorByFeatureSource = "---", 
                                    ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "S04_AAACAGCCAGAATGAC-1", 
                                    ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE, 
                                    ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1, 
                                    ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = FALSE, 
                                    VisualChoices = "Color", ContourAdd = FALSE, ContourColor = "#0000FF", 
                                    FixAspectRatio = FALSE, ViolinAdd = TRUE, PointSize = 1, 
                                    PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200, 
                                    CustomLabels = FALSE, CustomLabelsText = "S04_AAACAGCCAGAATGAC-1", 
                                    FontSize = 1, LegendPointSize = 1, LegendPosition = "Bottom", 
                                    HoverInfo = TRUE, LabelCenters = FALSE, LabelCentersBy = "orig.ident", 
                                    LabelCentersColor = "#000000", VersionInfo = list(iSEE = structure(list(
                                        c(2L, 18L, 0L)), class = c("package_version", "numeric_version"
                                        ))), PanelId = 4L, PanelHeight = 400L, PanelWidth = 3L, SelectionBoxOpen = FALSE, 
                                    RowSelectionSource = "---", ColumnSelectionSource = "---", 
                                    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
                                    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
                                    SelectionHistory = list())

## Build the iSEE app
iSEE(sce,
        appTitle = "Habenula_multiome", 
        initial = initial
        # colormap = ExperimentColorMap(
        #     colData = list(
        #          cell_type_fine= function(n) {
        #             return(sn_colors$fine)
        #         }
        #         ,
        #         cell_type_broad = function(n) {
        #             return(bulk_colors$broad)
        #         }
        #     )
        # )
    )
