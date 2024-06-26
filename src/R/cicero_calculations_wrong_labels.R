library(cicero)
library(rtracklayer)

path = sub('enhancer-promoter-interactions.*', 'enhancer-promoter-interactions', getwd())
path_muszka=paste0(path, '/data/muszka')
conns = readRDS(file = paste0(path, '/results/saved_conns.rds'))

gtf_file <- paste0(path_muszka, '/dmel-all-r6.36.gtf.gz')
gr <- readGFF(gtf_file)

colnames(gr)[which(names(gr) %in% c('seqid', 'transcript_id', 'gene_id'))] <- c('chromosome', 'transcript', 'symbol')


plot_connections(conns, "2L", 22e6, 22.1e6, 
                 gene_model = gr, 
                 coaccess_cutoff = 0.4, 
                 connection_width = .5, 
                 collapseTranscripts = "longest" )
