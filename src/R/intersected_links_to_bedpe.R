#save results to bedpe to visualize them in IGV browser
library(rtracklayer)

#path = '/home/julia/Desktop/uni/enhancer-promoter-interactions'
path = sub('enhancer-promoter-interactions.*', 'enhancer-promoter-interactions', getwd())
print(path)
#read csv
finds = read.csv(file = paste0(path, '/results/whole_data/intersection_results_bigger_set.csv'), sep = ',')
cicero = read.csv(file = paste0(path, '/results/whole_data/cicero_conns.csv'), sep = ',')
cicero = cicero[unlist(finds['cicero_index']), ]
finds$score = cicero[unlist(finds['cicero_index']), ]$score

peak1 = GRanges(seqnames=cicero$chr, ranges=IRanges(start=cicero$start1, end=cicero$end1))
peak2 = GRanges(seqnames=cicero$chr, ranges=IRanges(start=cicero$start2, end=cicero$end2))

pairs <- Pairs(first  = peak1,
                     second = peak2,
                     score = cicero$score)

rtracklayer::export(object = pairs, format = 'bedpe', con =  paste0(path,'/results/whole_data/intersected_links.bedpe'))

loops_df = read.csv(file = paste0(path, '/data/muszka/long_and_short_range_loops_D_mel_formatted.csv'), sep = ',')
anchor1 = GRanges(seqnames=loops_df$chr, ranges=IRanges(start=loops_df$start1, end=loops_df$end1))
anchor2 = GRanges(seqnames=loops_df$chr, ranges=IRanges(start=loops_df$start2, end=loops_df$end2))

loop_pairs <- Pairs(first  = anchor1,
                     second = anchor2)
                     

rtracklayer::export(object = loop_pairs, format = 'bedpe', con =  paste0(path,'/results/whole_data/long_and_short_range_loops_D_mel.bedpe'))
