library(cicero)
library(igraph)
set.seed(736)

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
dist_p = estimate_distance_parameter(cicero_cds, genomic_coords = sample_genome, max_elements = 1000)
#conns <- run_cicero(cicero_cds, genomic_coords=sample_genome, sample_num=5, max_elements=1000) # Takes a few minutes to run
head(conns)