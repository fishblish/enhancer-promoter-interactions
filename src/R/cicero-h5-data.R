library(cicero)
library(rhdf5)
library(GenomicFeatures)
library(rtracklayer)
library(dplyr)

# próba użycia funkcji load_cellranger_data_h5 zakończyła się niepomyślnie - w uzyskanym obiekcie cds brak niektórych pól, np lowerDetectionLimit

path='/home/julia/Desktop/uni/enhancer-promoter-interactions/muszka/filtered_feature_bc_matrix.h5'

fileh5 <- h5read(path, 'matrix')
data <- h5read(path, "/matrix/data")
indices <- h5read(path, "/matrix/indices")
indptr <- h5read(path, "/matrix/indptr")
barcodes <- h5read(path, "/matrix/barcodes")
feature_ids <- h5read(path, "/matrix/features/id")
feature_names <- h5read(path, "/matrix/features/name")


df <- data.frame(
  site_name = feature_ids[indices + 1],
  cell_name = rep(barcodes, diff(indptr)),
  read_count = data
)
df = df[1:500000,]
rm(data, indices, indptr, barcodes, feature_ids, feature_names)

cds <- make_atac_cds(df)
rm(df)
gc()

cds <- detectGenes(cds)
cds <- estimateSizeFactors(cds)
gc()

cds <- reduceDimension(cds, max_components = 2, num_dim = 6,
                       reduction_method = 'tSNE', norm_method = "none",
                       perplexity = 3, theta = 0.2)

tsne_coords <- t(reducedDimA(cds))
row.names(tsne_coords) <- row.names(pData(cds))
cicero_cds <- make_cicero_cds(cds, reduced_coordinates = tsne_coords)

chrom_sizes <- read.table(file("/home/julia/Desktop/uni/enhancer-promoter-interactions/muszka/dmel-all-chromosome-r6.36.chrom.sizes")) 
sample_genome <- subset(chrom_sizes, V1=='2L')
conns <- run_cicero(cicero_cds, genomic_coords=sample_genome) # Takes a few minutes to run

# plik z annotacjami genów
gtf_file <- "/home/julia/Desktop/uni/enhancer-promoter-interactions/muszka/Drosophila_melanogaster.BDGP6.32.57.gtf.gz"
gr <- readGFF(gtf_file)

colnames(gr)[which(names(gr) %in% c('seqid', 'transcript_id', 'gene_name'))] <- c('chromosome', 'transcript', 'symbol')
gr <- gr %>%
  rename(chromosome = seqid, transcript=transcript_id, symbol=gene_name)

plot_connections(conns, "2L", 2288, 1000000, 
                 gene_model = gr, 
                 coaccess_cutoff = .6, 
                 connection_width = .5, 
                 collapseTranscripts = "longest" )




              