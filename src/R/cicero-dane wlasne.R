library(cicero)
install.packages("hdf5r")
h5ls("/home/julia/Desktop/uni/enhancer-promoter-interactions/muszka/filtered_feature_bc_matrix.h5")
install.packages("remotes")
remotes::install_github("scfurl/m3addon", build = FALSE)
set.seed(736)
library(rhdf5)
library(m3addon)
library(monocle3)


path='/home/julia/Desktop/uni/enhancer-promoter-interactions/muszka/filtered_feature_bc_matrix.h5'
result <- load_cellranger_data_h5(
  files = path)
class(result)
setClass(
  "input_cds",
  contains = "input_cds",
  slots = c(lowerDetectionLimit = "numeric")
)
input_cds$lowerDetectionLimit <- 0
result <- detectGenes(result)

traceh5data = load_cellranger_data_h5(folders='/home/julia/Desktop/uni/enhancer-promoter-interactions/muszka')
mydata <- h5read("/home/julia/Desktop/uni/enhancer-promoter-interactions/muszka/filtered_feature_bc_matrix.h5", "/matrix")
library(readr)
start_data <- read.table(gzfile("/home/julia/Desktop/strony/muszka/GSM6614577_fragments.tsv.gz"), nrows=10000)   

colnames(start_data) <- c('chrom','start','end','cell', 'count') 
start_data$name <- paste(start_data$chrom, start_data$start, start_data$end, sep="_")
start_data = start_data[,c('name', 'cell', 'count')]

input_cds <- make_atac_cds(start_data, binarize = TRUE)
class(input_cds)
rm(start_data)

input_cds <- detectGenes(input_cds)
input_cds <- estimateSizeFactors(input_cds)

input_cds <- reduceDimension(input_cds, max_components = 2, num_dim=6,
                             reduction_method = 'tSNE', norm_method = "none")

tsne_coords <- t(reducedDimA(input_cds))
row.names(tsne_coords) <- row.names(pData(input_cds))
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = tsne_coords)
rm(input_cds)
gc()

chrom_sizes <- read.table(file("/home/julia/Desktop/strony/muszka/dmel-all-chromosome-r6.36.chrom.sizes")) 
sample_genome <- subset(chrom_sizes, V1=='2L')
dist_p = estimate_distance_parameter(cicero_cds, genomic_coords = sample_genome, s=0.85, distance_constraint = 50000, window = 100000)
#conns <- run_cicero(cicero_cds, genomic_coords=sample_genome, sample_num=5, max_elements=1000) # Takes a few minutes to run
head(conns)

library(cicero)
# library(igraph)
# install.packages("remotes")
library(rhdf5)
path='/home/julia/Desktop/uni/enhancer-promoter-interactions/muszka/filtered_feature_bc_matrix.h5'
# Read the datasets from the HDF5 file
data <- h5read(path, "/matrix/data")
indices <- h5read(path, "/matrix/indices")
indptr <- h5read(path, "/matrix/indptr")
barcodes <- h5read(path, "/matrix/barcodes")
feature_ids <- h5read(path, "/matrix/features/id")
feature_names <- h5read(path, "/matrix/features/name")

# Combine the datasets into a data frame
df <- data.frame(
  site_name = feature_ids[indices + 1],
  cell_name = rep(barcodes, diff(indptr)),
  read_count = data
)

# Use the make_atac_cds function to create a CDS object
cds <- make_atac_cds(df)

# Use the detectGenes function
cds <- detectGenes(cds)
# result <- load_cellranger_data_h5(
#   files = path)
# input_cds=result
# class(input_cds)
# setClass(
#   "input_cds",
#   contains = "input_cds",
#   slots = c(lowerDetectionLimit = "numeric")
# )
# input_cds$lowerDetectionLimit <- 0
# input_cds <- detectGenes(input_cds)
cds <- estimateSizeFactors(cds)

cds <- reduceDimension(cds, max_components = 2, num_dim=6,
                       reduction_method = 'tSNE', norm_method = "none")

tsne_coords <- t(reducedDimA(cds))
row.names(tsne_coords) <- row.names(pData(cds))
cicero_cds <- make_cicero_cds(cds, reduced_coordinates = tsne_coords)
