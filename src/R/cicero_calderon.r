library(cicero)
path = sub('enhancer-promoter-interactions.*', 'enhancer-promoter-interactions', getwd())
#path= '/home/julia/Desktop/uni/enhancer-promoter-interactions'
data_path = paste0(path, '/data/muszka/calderon_data')
print(path)
print(paste0(Sys.time(), ' Ładownie macierzy'))

gz= gzfile(paste0(data_path, "/GSE190130_exp1_hrs06-10_b1.peak_matrix.mtx.gz"))
indata <- Matrix::readMM(gz)
indata@x[indata@x > 0] <- 1
print(paste0(Sys.time(), ' Przygotowywanie input_cds'))

cellinfo <- read.table(paste0(data_path,"/GSE190130_exp1_hrs06-10_b1.peak_matrix.columns.txt"))
row.names(cellinfo) <- cellinfo$V1
names(cellinfo) <- "cells"

peakinfo <- read.table(paste0(data_path,"/GSE190130_exp1_hrs06-10_b1.peak_matrix.rows.txt"))
names(peakinfo) <- 'site_name'
row.names(peakinfo) <- peakinfo$site_name

row.names(indata) <- row.names(peakinfo)
colnames(indata) <- row.names(cellinfo)

#load rds file with cell types
annot = readRDS(file  = paste0(data_path, '/atac_meta.rds'))
cell_types = c('Neural', 'Ventral nerve cord', 'Brain', 'Ventral nerve cord prim.', 'Brain prim.')
cells_subset = annot[annot$refined_annotation %in% cell_types,]$cell
indata = indata[, colnames(indata) %in% cells_subset]
cellinfo = cellinfo[cellinfo$cells %in% cells_subset,]
cellinfo = data.frame(cells = cellinfo)

# Set the row names
row.names(cellinfo) = cellinfo$cells

# # Find indices where row names contain '2L'
# indices <- grepl("2L", rownames(indata))
# indata <- indata[indices, ]

# make CDS
input_cds <-  suppressWarnings(new_cell_data_set(indata,
                                            cell_metadata = cellinfo,
                                            gene_metadata = peakinfo  ))
rm(indata)
input_cds <- monocle3::detect_genes(input_cds)


#Ensure there are no peaks included with zero reads
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]

set.seed(2017)
input_cds <- detect_genes(input_cds)
print(paste0(Sys.time(), ' Input_cds gotowe. EstimateSizeFactors'))
input_cds <- estimate_size_factors(input_cds)
print(paste0(Sys.time(), ' preprocess_cds i reduce_dimention'))
input_cds <- preprocess_cds(input_cds, method = 'LSI')

input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP',
                                preprocess_method = 'LSI')

umap_coords <- reducedDims(input_cds)$UMAP

print(paste0(Sys.time(), ' make_cicero_cds'))

cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)
saveRDS(cicero_cds, file = paste0(path, '/results/saved_cicero_cds_calderon.rds'))
cicero_cds = readRDS(file = paste0(path, '/results/saved_cicero_cds_calderon.rds'))

chrom_sizes <- read.table(file(paste0(path,'/data/muszka/dmel-all-chromosome-r6.36.chrom.sizes')))
print(paste0(Sys.time(), ' wyznaczanie distance parameter'))

distance_parameters = estimate_distance_parameter(cicero_cds, s=0.85, distance_constraint = 5e4, window = 1e5, genomic_coords = chrom_sizes )
mean_distance_parameter <- mean(unlist(distance_parameters))
print(paste0(Sys.time(), ' wyznaczanie generate_cicero_models'))

cicero_out <- generate_cicero_models(cicero_cds, distance_parameter = mean_distance_parameter, s=0.85, window = 1e5, genomic_coords = chrom_sizes)
print(paste0(Sys.time(), ' wyznaczanie assemble_connections'))

conns <- assemble_connections(cicero_out)
head(conns)
saveRDS(conns, file = paste0(path, '/results/saved_conns_calderon.rds'))
conns = readRDS(file = paste0(path, '/results/saved_conns_calderon.rds'))

library(rtracklayer)
library(Signac)
print(paste0(Sys.time(), ' Zapisywanie wyników'))

conns = na.omit(conns)

color = ifelse(conns$coaccess < 0, '255,0,0', '0,0,255')
conns.pairs <- Pairs(first  = StringToGRanges(conns$Peak1, sep=c("_", "_")),
                     second = StringToGRanges(conns$Peak2, sep=c("_", "_")),
                     score = conns$coaccess,
                     color = color)

rtracklayer::export(object = conns.pairs, format = 'bedpe', con =  paste0(path,'/results/calderon/links_calderon.bedpe'))

print(paste0(Sys.time(), ' koniec'))



