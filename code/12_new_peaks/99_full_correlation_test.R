#Testing out time and memory for computing the full peak by gene correlation matrix 


library(Seurat)
library(qs2)
library(dplyr)
library(here)

here::here()
#Peak dataframe
peak_file = here('processed-data','06_peak_calling','01_call_peaks_MACS2','macs_peaks_Mid_resolution.csv')

seurat_file = here('processed-data','12_new_peaks','07_non_pb_seur','non_pb_seur.qs2')

seurat_multiome = qs_read(seurat_file)
seurat_multiome

peaks_df = data.table::fread(peak_file)
peaks_df %>% head() %>% View()

dim(peaks_df)
# 355127 6

#Normalized atac data
dim(seurat_multiome@assays$ATAC@data)
seurat_multiome@assays$ATAC@data[1:10, 1:10]

#Sanity check the peaks match up
seurat_peaks = rownames(seurat_multiome@assays$ATAC@data)
seurat_peak_chrom = sapply(strsplit(seurat_peaks, split = '-', fixed = TRUE ),'[[' ,1 )
seurat_peak_start = sapply(strsplit(seurat_peaks, split = '-', fixed = TRUE ),'[[' ,2 )
seurat_peak_end = sapply(strsplit(seurat_peaks, split = '-', fixed = TRUE ),'[[' ,3 )

table(seurat_peak_chrom == peaks_df$seqnames)
table(seurat_peak_start == peaks_df$start)
table(seurat_peak_end == peaks_df$end)

#  TRUE 
#355127 

#They do match up



#Test out computing pearson correlations with matrix operations
gene_names = rownames(seurat_multiome@assays$RNA@features@.Data)

test_peaks = seurat_multiome@assays$ATAC@data[1:1000,1:100 ]
test_genes = seurat_multiome@assays$RNA@layers$data[1:1000,1:100 ]
rownames(test_genes) = gene_names[1:1000]


dim(test_peaks)
dim(test_genes)

#pearson correlation code, colab with chatGPT

n_cells = ncol(test_peaks)

M <- test_peaks %*% t(test_genes)

mu_peaks = Matrix::rowMeans(test_peaks)
mu_genes = Matrix::rowMeans(test_genes)

#Covariance
cov_peaks_genes = (M / n_cells ) - (mu_peaks %o% mu_genes)

#Standard deviations
# Genes
var_G <- apply(test_genes, 1, var)
std_G <- sqrt(var_G)
std_G[std_G == 0] <- 1e-8

# Peaks
var_P <- apply(test_peaks, 1, var)
std_P <- sqrt(var_P)
std_P[std_P == 0] <- 1e-8

#Compute pearsons
pearson_mat = cov_peaks_genes / (std_P %o% std_G)

as.data.frame(pearson_mat[1:10, 1:10]) %>% View()

#Even just for 1000 genes and peaks with 100 cells, this took like 10 minutes
brute_pearson_mat = matrix(nrow = nrow(test_peaks), ncol = nrow(test_genes))

for(i in 1:nrow(brute_pearson_mat)){
  for(j in 1:ncol(brute_pearson_mat)){
    brute_pearson_mat[i,j] = cor(test_peaks[i,], test_genes[j,], method = 'pearson')
  }
}


brute_pearson_mat[is.na(brute_pearson_mat)] <- 0

pearson_vec <- as.numeric(pearson_mat)
brute_vec <- as.numeric(brute_pearson_mat)

# Plot against each other
plot(pearson_vec, brute_vec, 
  main = 'Comparing Pearson implementations: 1000 genes and peaks, 100 cells', xlab = 'Matrix Pearson', ylab = 'Brute-force Pearson')
abline(a=0, b=1, col='red')

#They are not exactly identical, but we're talking like .00001 off from the identity line



#Try different combinations of genes and peaks with all the cells
#1.8 minutes and ~10gb for 10kpeaks by 10kgenes with all cells
#6.3 minutes and ~15gb for 10k peaks by all genes and all cells
test_peaks = seurat_multiome@assays$ATAC@data[1:20000, ]
test_genes = seurat_multiome@assays$RNA@layers$data
rownames(test_genes) = gene_names


dim(test_peaks)
dim(test_genes)

#pearson correlation code, colab with chatGPT
start_time <- Sys.time()
n_cells = ncol(test_peaks)

M <- test_peaks %*% t(test_genes)

mu_peaks = Matrix::rowMeans(test_peaks)
mu_genes = Matrix::rowMeans(test_genes)

#Covariance
cov_peaks_genes = (M / n_cells ) - (mu_peaks %o% mu_genes)

#Standard deviations
# Peaks
P_sq <- test_peaks
P_sq@x <- P_sq@x^2  # square only the non-zero entries
var_P <- rowMeans(P_sq) - (rowMeans(test_peaks)^2)
std_P <- sqrt(var_P)
std_P[std_P == 0] <- 1e-8  # avoid division by zero
#Genes
G_sq <- test_genes
G_sq@x <- G_sq@x^2  # square only the non-zero entries
var_G <- rowMeans(G_sq) - (rowMeans(test_genes)^2)
std_G <- sqrt(var_G)
std_G[std_G == 0] <- 1e-8  # avoid division by zero

#Compute pearsons
pearson_mat = cov_peaks_genes / (std_P %o% std_G)

Sys.time() - start_time


as.data.frame(pearson_mat[1:10, 1:10]) %>% View()


dim(seurat_multiome@assays$ATAC@data)



#And donor-specific networks
#MHb.2 cell number 2610
rna_meta = as.data.frame(seurat_multiome[[]])
current_donors = rna_meta %>% group_by(orig.ident, mid_cluster) %>% summarise(n_cells = n()) %>% 
  filter(n_cells >= 100 & mid_cluster == 'MHb.2') %>% pull(orig.ident)


gene_index = rowSums(seurat_multiome@assays$RNA@layers$data[, seurat_multiome$mid_cluster == 'MHb.2']) > 0
table(gene_index)

library(Matrix)

rank_normalize_dense <- function(x, ties.method = "average", na.last = "keep") {
  stopifnot(inherits(x, "dgeMatrix") || is.matrix(x))
  
  vals <- as.numeric(x)
  ok <- !is.na(vals)
  
  r <- rank(vals[ok], ties.method = ties.method, na.last = na.last)
  vals[ok] <- r / max(r)
  
  if (inherits(x, "dgeMatrix")) {
    x@x <- vals
    x
  } else {
    matrix(vals, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
  }
}

pearson_mat_rank <- rank_normalize_dense(pearson_mat)

donor_mat_list = vector(mode = 'list', length = length(current_donors))

peak_index = grepl('MHb.2',peaks_df$peak_called_in)
gene_index = rowSums(seurat_multiome@assays$RNA@layers$data[, seurat_multiome$mid_cluster == 'MHb.2']) > 0

for(i in 1:length(current_donors)){
  i = 1
  donor_index = rna_meta$orig.ident == current_donors[i] & rna_meta$mid_cluster == 'MHb.2'

  test_peaks = seurat_multiome@assays$ATAC@data[peak_index, donor_index ]
  test_genes = seurat_multiome@assays$RNA@layers$data[gene_index , donor_index]
  rownames(test_genes) = gene_names[gene_index]


  start_time <- Sys.time()
  n_cells = ncol(test_peaks)

  M <- test_peaks %*% t(test_genes)

  mu_peaks = Matrix::rowMeans(test_peaks)
  mu_genes = Matrix::rowMeans(test_genes)

  #Covariance
  cov_peaks_genes = (M / n_cells ) - (mu_peaks %o% mu_genes)

  #Standard deviations
  # Peaks
  P_sq <- test_peaks
  P_sq@x <- P_sq@x^2  # square only the non-zero entries
  var_P <- rowMeans(P_sq) - (rowMeans(test_peaks)^2)
  std_P <- sqrt(var_P)
  std_P[std_P == 0] <- 1e-8  # avoid division by zero
  #Genes
  G_sq <- test_genes
  G_sq@x <- G_sq@x^2  # square only the non-zero entries
  var_G <- rowMeans(G_sq) - (rowMeans(test_genes)^2)
  std_G <- sqrt(var_G)
  std_G[std_G == 0] <- 1e-8  # avoid division by zero

  #Compute pearsons
  # usage
  pearson_mat = cov_peaks_genes / (std_P %o% std_G)
  pearson_mat_rank <- rank_normalize_dense(pearson_mat)
  donor_mat_list[[i]] = pearson_mat_rank

  rm(pearson_mat, pearson_mat_rank)
  Sys.time() - start_time

}

 pearson_mat[1:10, 1:10]

long_pearson = reshape2::melt(as.matrix(pearson_mat))
dim(long_pearson)

long_pearson %>% filter(value >0.5 & value != 1) %>% head() %>% View()
length(donor_corrs)
hist(sample(donor_corrs, 1000000), breaks = 50)

rownames(seurat_multiome@assays$RNA@layers$data) = gene_names

peak_index = rownames(seurat_multiome@assays$ATAC@data) == 'chr1-17409419-17409777'
gene_index = rownames(seurat_multiome@assays$RNA@layers$data) == 'AL627309.5'
atac_signal = seurat_multiome@assays$ATAC@data[peak_index, donor_index ]
gene_signal = seurat_multiome@assays$RNA@layers$data[gene_index, donor_index]

plot(atac_signal, gene_signal)

