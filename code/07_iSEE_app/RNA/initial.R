initial <- list()

################################################################################
# Settings for Reduced dimension plot 1
################################################################################

initial[["ReducedDimensionPlot1"]] <- new(
    "ReducedDimensionPlot",
    Type = "wnn_umap",
    XAxis = 1L,
    YAxis = 2L,
    FacetRowByColData = "mid_cluster",
    FacetColumnByColData = "mid_cluster",
    ColorByColumnData = "mid_cluster",
    ColorByFeatureNameAssay = "logcounts",
    ColorBySampleNameColor = "#FF0000",
    ShapeByColumnData = "mid_cluster",
    # SizeByColumnData = "age",
    TooltipColumnData = character(0),
    FacetRowBy = "None",
    FacetColumnBy = "None",
    ColorBy = "Column data",
    ColorByDefaultColor = "#000000",
    ColorByFeatureName = "MBP",
    ColorByFeatureSource = "---",
    ColorByFeatureDynamicSource = FALSE,
    ColorBySampleName = "S03_Hb_r",
    ColorBySampleSource = "---",
    ColorBySampleDynamicSource = FALSE,
    ShapeBy = "None",
    SizeBy = "None",
    SelectionAlpha = 0.1,
    ZoomData = numeric(0),
    BrushData = list(),
    VisualBoxOpen = FALSE,
    VisualChoices = "Color",
    ContourAdd = FALSE,
    ContourColor = "#0000FF",
    PointSize = 1,
    PointAlpha = 1,
    Downsample = FALSE,
    DownsampleResolution = 200,
    CustomLabels = FALSE,
    CustomLabelsText = "S03_Hb_r",
    FontSize = 1,
    LegendPointSize = 1,
    LegendPosition = "Bottom",
    HoverInfo = TRUE,
    LabelCenters = FALSE,
    LabelCentersBy = "mid_cluster",
    LabelCentersColor = "#000000",
    VersionInfo = list(iSEE = structure(
        list(c(2L, 10L, 0L)),
        class = c(
            "package_version",
            "numeric_version"
        )
    )),
    PanelId = c(ReducedDimensionPlot = 1L),
    PanelHeight = 600L,
    PanelWidth = 5L,
    SelectionBoxOpen = FALSE,
    RowSelectionSource = "---",
    ColumnSelectionSource = "---",
    DataBoxOpen = FALSE,
    RowSelectionDynamicSource = FALSE,
    ColumnSelectionDynamicSource = FALSE,
    RowSelectionRestrict = FALSE,
    ColumnSelectionRestrict = FALSE,
    SelectionHistory = list()
)

################################################################################
# Settings for Complex heatmap 1
################################################################################

initial[["ComplexHeatmapPlot1"]] <- new(
    "ComplexHeatmapPlot",
    Assay = "logcounts",
    CustomRows = TRUE,
    CustomRowsText = "SNAP25\nAQP4\nGAD1\nSLC17A6\nSLC32A1\nRELN\nMBP\nGFAP",
    ClusterRows = FALSE,
    ClusterRowsDistance = "spearman",
    ClusterRowsMethod = "ward.D2",
    DataBoxOpen = FALSE,
    VisualChoices = "Annotations",
    ColumnData = c(
        "mid_cluster"#,
        # "subject"
    ),
    RowData = character(0),
    CustomBounds = FALSE,
    LowerBound = NA_real_,
    UpperBound = NA_real_,
    AssayCenterRows = TRUE,
    AssayScaleRows = TRUE,
    DivergentColormap = "blue < white < orange",
    ShowDimNames = "Rows",
    LegendPosition = "Right",
    LegendDirection = "Vertical",
    VisualBoxOpen = FALSE,
    NamesRowFontSize = 10,
    NamesColumnFontSize = 10,
    ShowColumnSelection = TRUE,
    OrderColumnSelection = TRUE,
    VersionInfo = list(iSEE = structure(
        list(c(2L, 10L, 0L)),
        class = c(
            "package_version",
            "numeric_version"
        )
    )),
    PanelId = c(ComplexHeatmapPlot = 1L),
    PanelHeight = 600L,
    PanelWidth = 7L,
    SelectionBoxOpen = FALSE,
    RowSelectionSource = "---",
    ColumnSelectionSource = "---",
    RowSelectionDynamicSource = FALSE,
    ColumnSelectionDynamicSource = FALSE,
    RowSelectionRestrict = FALSE,
    ColumnSelectionRestrict = FALSE,
    SelectionHistory = list()
)

################################################################################
# Settings for Row data table 1
################################################################################

initial[["RowDataTable1"]] <- new(
    "RowDataTable",
    Selected = "SNAP25",
    Search = "",
    SearchColumns = c(
        "",
        "", "", "", "", "", ""
    ),
    HiddenColumns = character(0),
    VersionInfo = list(iSEE = structure(
        list(c(2L, 10L, 0L)),
        class = c(
            "package_version",
            "numeric_version"
        )
    )),
    PanelId = c(RowDataTable = 1L),
    PanelHeight = 600L,
    PanelWidth = 5L,
    SelectionBoxOpen = FALSE,
    RowSelectionSource = "---",
    ColumnSelectionSource = "---",
    DataBoxOpen = FALSE,
    RowSelectionDynamicSource = FALSE,
    ColumnSelectionDynamicSource = FALSE,
    RowSelectionRestrict = FALSE,
    ColumnSelectionRestrict = FALSE,
    SelectionHistory = list()
)

################################################################################
# Settings for Feature assay plot 1
################################################################################

initial[["FeatureAssayPlot1"]] <- new(
    "FeatureAssayPlot",
    Assay = "logcounts",
    XAxis = "Column data",
    XAxisColumnData = "mid_cluster",
    XAxisFeatureName = "MBP",
    XAxisFeatureSource = "---",
    XAxisFeatureDynamicSource = FALSE,
    YAxisFeatureName = "SNAP25",
    YAxisFeatureSource = "RowDataTable1",
    YAxisFeatureDynamicSource = TRUE,
    FacetRowByColData = "mid_cluster",
    FacetColumnByColData = "mid_cluster",
    ColorByColumnData = "mid_cluster",
    ColorByFeatureNameAssay = "logcounts",
    ColorBySampleNameColor = "#FF0000",
    ShapeByColumnData = "mid_cluster",
    # SizeByColumnData = "age",
    TooltipColumnData = character(0),
    FacetRowBy = "None",
    FacetColumnBy = "None",
    ColorBy = "Column data",
    ColorByDefaultColor = "#000000",
    ColorByFeatureName = "MBP",
    ColorByFeatureSource = "---",
    ColorByFeatureDynamicSource = FALSE,
    ColorBySampleName = "S03_Hb_r",
    ColorBySampleSource = "---",
    ColorBySampleDynamicSource = FALSE,
    ShapeBy = "None",
    SizeBy = "None",
    SelectionAlpha = 0.1,
    ZoomData = numeric(0),
    BrushData = list(),
    VisualBoxOpen = FALSE,
    VisualChoices = "Color",
    ContourAdd = FALSE,
    ContourColor = "#0000FF",
    PointSize = 1,
    PointAlpha = 1,
    Downsample = FALSE,
    DownsampleResolution = 200,
    CustomLabels = FALSE,
    CustomLabelsText = "S03_Hb_r",
    FontSize = 1,
    LegendPointSize = 1,
    LegendPosition = "Bottom",
    HoverInfo = TRUE,
    LabelCenters = FALSE,
    LabelCentersBy = "mid_cluster",
    LabelCentersColor = "#000000",
    VersionInfo = list(iSEE = structure(
        list(c(2L, 10L, 0L)),
        class = c("package_version", "numeric_version")
    )),
    PanelId = c(FeatureAssayPlot = 1L),
    PanelHeight = 600L,
    PanelWidth = 7L,
    SelectionBoxOpen = FALSE,
    RowSelectionSource = "---",
    ColumnSelectionSource = "---",
    DataBoxOpen = FALSE,
    RowSelectionDynamicSource = FALSE,
    ColumnSelectionDynamicSource = FALSE,
    RowSelectionRestrict = FALSE,
    ColumnSelectionRestrict = FALSE,
    SelectionHistory = list()
)
