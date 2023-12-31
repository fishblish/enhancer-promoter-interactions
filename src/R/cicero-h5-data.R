library(cicero)
library(rhdf5)
library(GenomicFeatures)
library(rtracklayer)
library(dplyr)
library(m3addon)
path = '/home/julia/Desktop/uni/enhancer-promoter-interactions'
#path = sub('enhancer-promoter-interactions.*', 'enhancer-promoter-interactions', getwd())
path_muszka=paste0(path, '/data/muszka')
path_h5=paste0(path_muszka, '/filtered_feature_bc_matrix.h5')

# próba użycia funkcji load_cellranger_data_h5 zakończyła się niepomyślnie - w uzyskanym obiekcie cds brak niektórych pól, np lowerDetectionLimit


fileh5 <- h5read(path_h5, 'matrix')
data <- h5read(path_h5, "/matrix/data")
indices <- h5read(path_h5, "/matrix/indices")
indptr <- h5read(path_h5, "/matrix/indptr")
barcodes <- h5read(path_h5, "/matrix/barcodes")
feature_ids <- h5read(path_h5, "/matrix/features/id")
feature_names <- h5read(path_h5, "/matrix/features/name")

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

                       reduction_method = 'tSNE', norm_method = "none")


tsne_coords <- t(reducedDimA(cds))
row.names(tsne_coords) <- row.names(pData(cds))
cicero_cds <- make_cicero_cds(cds, reduced_coordinates = tsne_coords)

chrom_sizes <- read.table(file(paste0(path_muszka,'/dmel-all-chromosome-r6.36.chrom.sizes'))) 
sample_genome <- subset(chrom_sizes, V1=='2L')
#conns <- run_cicero(cicero_cds, genomic_coords=sample_genome, s=) # Takes a few minutes to run

distance_parameters = estimate_distance_parameter(cicero_cds, s=0.85, distance_constraint = 5e4, window = 1e5, genomic_coords = sample_genome )
mean_distance_parameter <- mean(unlist(distance_parameters))
cicero_out <- generate_cicero_models(cicero_cds, distance_parameter = mean_distance_parameter, s=0.85, window = 1e5, genomic_coords = sample_genome)
conns <- assemble_connections(cicero_out)

saveRDS(conns, file = paste0(path, '/results/saved_conns_params.rds'))
saveRDS(cicero_cds, file = paste0(path, '/results/saved_cds_params.rds'))
conns = readRDS(file = paste0(path, '/results/saved_conns.rds'))
cds = readRDS(file = paste0(path, '/results/saved_cds_params.rds'))

# plik z annotacjami genów
gtf_file <- paste0(path_muszka, '/dmel-all-r6.36.gtf.gz')
gr <- readGFF(gtf_file)

colnames(gr)[which(names(gr) %in% c('seqid', 'transcript_id', 'gene_id'))] <- c('chromosome', 'transcript', 'symbol')


plot_connections(conns, "2L", 22e6, 22.1e6, 
                 gene_model = gr, 
                 coaccess_cutoff = 0.4, 
                 connection_width = .5, 
                 collapseTranscripts = "longest" )

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



library(rtracklayer)
library(Signac)

conns = na.omit(conns)

color = ifelse(conns$coaccess < 0, '255,0,0', '0,0,255')
conns.pairs <- Pairs(first  = StringToGRanges(conns$Peak1, sep=c(":", "-")),
                     second = StringToGRanges(conns$Peak2, sep=c(":", "-")),
                     score = conns$coaccess,
                     color = color)

rtracklayer::export(object = conns.pairs, format = 'bedpe', con =  paste0(path,'/results/out_params.bedpe'))



              
