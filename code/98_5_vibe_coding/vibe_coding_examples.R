#Testing out Positron's AI Assistant feature for generating different kinds of plots for single cell data

#############################

#First, orient everyone to the Positron interface
#Show how to check on memory usage

#############################



#Load libraries and data
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(here)

#Mouse habenula data
hashikawa_path = '/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/09_cross_species_analysis/Hashikawa_data'
list.files(hashikawa_path)

#Loads as sce_mouse_sub, list of sce objects
load(paste0(hashikawa_path, '/sce_mouse_habenula.Rdata'))
hashikawa_sce = sce_mouse_sub$all
# Add cpms for later
assay(hashikawa_sce, 'cpm') = MetaMarkers::convert_to_cpm(assay(hashikawa_sce, 'counts'))


hashikawa_sce 

#Swap the rownames to gene symbols
rowData(hashikawa_sce)
rownames(hashikawa_sce) = rowData(hashikawa_sce)$Symbol

#Check out the available metadata

colnames(colData(hashikawa_sce))

table(hashikawa_sce$celltype)

table(hashikawa_sce$stim)


#Make a proportional bar plot of cell types by stimulus condition

hashikawa_sce |>
  colData() |>
  as.data.frame() |>
  group_by(stim, celltype) |>
  summarize(n = n(), .groups = "drop") |>
  group_by(stim) |>
  mutate(proportion = n / sum(n)) |>
  ggplot(aes(x = stim, y = proportion, fill = celltype)) +
  geom_col(position = "fill") +
  labs(
    x = "Condition",
    y = "Proportion",
    fill = "Cell Type",
    title = "Cell Type Proportions by Stimulation Status"
  )


#Try different palettes with subjective descriptions. Try different models
#This is with a cheerful palette
celltype_colors <- c(
  # Astrocytes - bright blues
  "Astrocyte1" = "#0096FF",
  "Astrocyte2" = "#00D9FF",
  
  # Endothelial and Epen - bright greens
  "Endothelial" = "#00D084",
  "Epen" = "#39FF14",
  
  # Microglia - bright reds
  "Microglia" = "#FF1744",
  
  # Mural - bright magentas
  "Mural" = "#FF10F0",
  
  # Neurons - bright yellows and oranges
  "Neuron1" = "#FFB700",
  "Neuron2" = "#FF6B35",
  "Neuron3" = "#FF8C00",
  "Neuron4" = "#FFA500",
  "Neuron5" = "#FFD700",
  "Neuron6" = "#FFDC00",
  "Neuron7" = "#FFED4E",
  "Neuron8" = "#FFF44F",
  
  # Oligodendrocytes - bright cyans
  "Oligo1" = "#00F0FF",
  "Oligo2" = "#0FF0FF",
  "Oligo3" = "#30B0FF",
  
  # OPCs - bright purples and pinks
  "OPC1" = "#D946EF",
  "OPC2" = "#EC4899",
  "OPC3" = "#F472B6"
)

hashikawa_sce |>
  colData() |>
  as.data.frame() |>
  group_by(stim, celltype) |>
  summarize(n = n(), .groups = "drop") |>
  group_by(stim) |>
  mutate(proportion = n / sum(n)) |>
  ggplot(aes(x = stim, y = proportion, fill = celltype)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = celltype_colors) +
  labs(
    x = "Condition",
    y = "Proportion",
    fill = "Cell Type",
    title = "Cell Type Proportions by Stimulation Status"
  )




#Make a bubble plot for a marker gene panel
#script to grab marker genes from: 98_external_Hb_comparisons/03_metaMarkers_wallace_2019

#Try the agent to grab the markers from that script
wallace_markers <- c('Tac2', 'Slc17a7', 'Slc17a6', 'Snap25', 'Gap43', 
                     'Slc6a11', 'Cldn5', 'Abcc9', 'Pdgfrb', 'Cx3cr1', 
                     'Mrc1', 'Col3a1', 'Gpr17', 'Mog', 'Olig1', 'Pdgfra')

#Make the pubble plot
# Filter to Wallace markers and calculate metrics by celltype
bubble_data <- hashikawa_sce[wallace_markers, ] |>
  assay('cpm') |>
  as.matrix() |>
  t() |>
  as.data.frame() |>
  mutate(celltype = hashikawa_sce$celltype) |>
  tidyr::pivot_longer(-celltype, names_to = 'gene', values_to = 'expression') |>
  group_by(celltype, gene) |>
  summarize(
    mean_expression = mean(expression),
    percent_expressing = sum(expression > 0) / n() * 100,
    .groups = 'drop'
  )

ggplot(bubble_data, aes(x = celltype, y = gene, size = percent_expressing, color = mean_expression)) +
  geom_point() +
  scale_size_continuous(name = '% Cells Expressing', range = c(1, 8)) +
  scale_color_gradient(name = 'Mean Expression (CPM)', low = 'lightgray', high = 'darkred') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = 'Cell Type',
    y = 'Gene',
    title = 'Wallace Markers Expression by Cell Type'
  )


# Add a z-scored version
bubble_data_zscore <- hashikawa_sce[wallace_markers, ] |>
  assay('cpm') |>
  as.matrix() |>
  t() |>
  as.data.frame() |>
  mutate(celltype = hashikawa_sce$celltype) |>
  tidyr::pivot_longer(-celltype, names_to = 'gene', values_to = 'expression') |>
  group_by(gene) |>
  mutate(expression_zscore = scale(expression)[,1]) |>
  ungroup() |>
  group_by(celltype, gene) |>
  summarize(
    mean_expression_zscore = mean(expression_zscore),
    percent_expressing = sum(expression > 0) / n() * 100,
    .groups = 'drop'
  )

ggplot(bubble_data_zscore, aes(x = celltype, y = gene, size = percent_expressing, color = mean_expression_zscore)) +
  geom_point() +
  scale_size_continuous(name = '% Cells Expressing', range = c(1, 8)) +
  scale_color_gradient2(name = 'Mean Expression\n(Z-scored)', low = 'blue', mid = 'white', high = 'red', midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = 'Cell Type',
    y = 'Gene',
    title = 'Wallace Markers Expression by Cell Type (Z-scored)'
  )


# Create a matrix for clustering
bubble_matrix <- bubble_data_zscore |>
  dplyr::select(gene, celltype, mean_expression_zscore) |>
  tidyr::pivot_wider(
    names_from = celltype,
    values_from = mean_expression_zscore,
    values_fill = 0
  ) |>
  (\(x) {
    m <- as.matrix(x[, -1, drop = FALSE])
    rownames(m) <- x$gene
    m
  })()

# Cluster genes and cell types
gene_order <- hclust(dist(bubble_matrix))$order
celltype_order <- hclust(dist(t(bubble_matrix)))$order

# Convert to factors with clustered levels
bubble_data_zscore_clustered <- bubble_data_zscore |>
  mutate(
    gene = factor(gene, levels = rownames(bubble_matrix)[gene_order]),
    celltype = factor(celltype, levels = colnames(bubble_matrix)[celltype_order])
  )

ggplot(bubble_data_zscore_clustered, aes(x = celltype, y = gene, size = percent_expressing, color = mean_expression_zscore)) +
  geom_point() +
  scale_size_continuous(name = '% Cells Expressing', range = c(1, 8)) +
  scale_color_gradient2(name = 'Mean Expression\n(Z-scored)', low = 'blue', mid = 'white', high = 'red', midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = 'Cell Type',
    y = 'Gene',
    title = 'Wallace Markers Expression by Cell Type (Z-scored & Clustered)'
  )

#Try a different model

library(ggplot2)

ggplot(
  bubble_data_zscore_clustered,
  aes(x = celltype, y = gene, size = percent_expressing, color = mean_expression_zscore)
) +
  geom_point() +
  scale_size_continuous(name = "% Cells Expressing", range = c(1, 8)) +
  scale_color_gradient2(
    name = "Mean Expression\n(Z-scored)",
    low = "blue", mid = "white", high = "red", midpoint = 0
  ) +
  annotate(
    "segment",
    x = which(levels(bubble_data_zscore_clustered$celltype) == "Neuron4") + 2.5,
    xend = which(levels(bubble_data_zscore_clustered$celltype) == "Neuron4"),
    y = which(levels(bubble_data_zscore_clustered$gene) == "Slc17a6") + 2.5,
    yend = which(levels(bubble_data_zscore_clustered$gene) == "Slc17a6"),
    arrow = arrow(length = unit(0.2, "cm")),
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = which(levels(bubble_data_zscore_clustered$celltype) == "Neuron4") + 3,
    y = which(levels(bubble_data_zscore_clustered$gene) == "Slc17a6") + 2.8,
    label = "Neuron4 / Slc17a6",
    hjust = 0,
    size = 3.5
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Cell Type",
    y = "Gene",
    title = "Wallace Markers Expression by Cell Type (Z-scored & Clustered)"
  )





#Compute cell-type markers

celltype_markers = MetaMarkers::compute_markers(assay(hashikawa_sce, 'cpm'), hashikawa_sce$celltype)

#View the top 10 markers by auroc
top10_per_celltype <- celltype_markers |>
  slice_max(order_by = auroc, n = 10, by = cell_type) |> View()



#Make a heatmap of the top 10 markers per cell type

library(dplyr)
library(tidyr)
library(ggplot2)
library(Matrix)

# 1) Top 10 markers per cell type (by AUROC)
top10_per_celltype <- celltype_markers |>
  slice_max(order_by = auroc, n = 10, by = cell_type) |>
  ungroup()

marker_genes <- top10_per_celltype |>
  pull(gene) |>
  unique()

# 2) Mean expression per gene per cell type from cpm assay
expr_mat <- assay(hashikawa_sce, "cpm")
celltypes <- hashikawa_sce$celltype

heat_df <- lapply(marker_genes, function(g) {
  v <- expr_mat[g, ]
  tibble(
    gene = g,
    cell_type = celltypes,
    expr = as.numeric(v)
  )
}) |>
  bind_rows() |>
  summarize(mean_expr = mean(expr), .by = c(gene, cell_type))

# 3) Z-score per gene across cell types
heat_df_z <- heat_df |>
  group_by(gene) |>
  mutate(z = as.numeric(scale(mean_expr))) |>
  ungroup()

# Optional: order genes by their marker-assigned cell type then AUROC
gene_order <- top10_per_celltype |>
  arrange(cell_type, desc(auroc)) |>
  pull(gene) |>
  unique()

heat_df_z <- heat_df_z |>
  mutate(
    gene = factor(gene, levels = rev(gene_order)),
    cell_type = factor(cell_type, levels = sort(unique(cell_type)))
  )

# 4) Heatmap
ggplot(heat_df_z, aes(x = cell_type, y = gene, fill = z)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#313695",
    mid = "white",
    high = "#A50026",
    midpoint = 0,
    name = "Z-score"
  ) +
  labs(
    x = "Cell type",
    y = "Top markers (10 per cell type)",
    title = "Top marker expression heatmap (z-scored per gene)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )




library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

# 1) Top 10 markers per cell type (by AUROC)
top10_per_celltype <- celltype_markers |>
  slice_max(order_by = auroc, n = 10, by = cell_type) |>
  ungroup()

marker_genes <- top10_per_celltype |>
  pull(gene) |>
  unique()

# 2) Mean expression per gene per cell type from cpm assay
expr_mat <- assay(hashikawa_sce, "cpm")
celltypes <- hashikawa_sce$celltype

heat_df <- lapply(marker_genes, function(g) {
  tibble(
    gene = g,
    cell_type = celltypes,
    expr = as.numeric(expr_mat[g, ])
  )
}) |>
  bind_rows() |>
  summarize(mean_expr = mean(expr), .by = c(gene, cell_type))

# 3) Wide matrix and z-score per gene across cell types
mat <- heat_df |>
  pivot_wider(names_from = cell_type, values_from = mean_expr) |>
  (\(x) {
    m <- as.matrix(x[, -1, drop = FALSE])
    rownames(m) <- x$gene
    m
  })()

mat_z <- t(scale(t(mat)))
mat_z[is.na(mat_z)] <- 0

# 4) Ward.D2 clustering
row_hc <- hclust(dist(mat_z), method = "ward.D2")
col_hc <- hclust(dist(t(mat_z)), method = "ward.D2")

# 5) ComplexHeatmap
Heatmap(
  mat_z,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("#313695", "white", "#A50026")),
  cluster_rows = as.dendrogram(row_hc),
  cluster_columns = as.dendrogram(col_hc),
  row_names_side = "left",
  column_names_rot = 45,
  heatmap_legend_param = list(title = "Z-score")
)


#Highlight specific genes
# Create heatmap with Apoe highlighted




#Maybe try recommendations for dimensionality rededuction plots



#Try visualizing average marker gene expression on a UMAP


##################
# 
# 
# Any other plot suggestions?
#
#  
##################
