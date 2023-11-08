library(cicero)
library(igraph)
set.seed(736)

data(cicero_data)

input_cds <- make_atac_cds(cicero_data, binarize = TRUE)
class(input_cds)


input_cds <- detectGenes(input_cds)
input_cds <- estimateSizeFactors(input_cds)

input_cds <- reduceDimension(input_cds, max_components = 2, num_dim=6,
                             reduction_method = 'tSNE', norm_method = "none")

tsne_coords <- t(reducedDimA(input_cds))
row.names(tsne_coords) <- row.names(pData(input_cds))
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = tsne_coords)

data("human.hg19.genome")
sample_genome <- subset(human.hg19.genome, V1 == "chr18")
conns <- run_cicero(cicero_cds, sample_genome) # Takes a few minutes to run

data(gene_annotation_sample)
plot_connections(conns, "chr18", 8575097, 8839855, 
                 gene_model = gene_annotation_sample, 
                 coaccess_cutoff = .25, 
                 connection_width = .5, 
                 collapseTranscripts = "longest" )
