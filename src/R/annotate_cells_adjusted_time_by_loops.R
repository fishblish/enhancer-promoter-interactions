# The script to annotate the cells in the context of loops from long_and_short_range_loops_D_mel.tsv file.
# Each cell gets annotation xy, where x and y are indicators of non-zero ATAC-seq signal in anchor 1 and 2 of the loop
# The script uses time windows from the NNv1_time.new column in the metadata file named 'atac.meta.rds'.


library(GenomicRanges)
library(Matrix)

# Read chromatin loops
loops <- read.table('data/long_and_short_range_loops_D_mel.tsv', header = TRUE, sep = '\t')
loops_A1_gr <- with(loops, GRanges(chr1, IRanges(x1 + 1L, x2)))
loops_A2_gr <- with(loops, GRanges(chr2, IRanges(y1 + 1L, y2)))

# Read scATAC-seq data
data_path <- 'data/calderon_data'
print(data_path)

annot <- readRDS(file.path(data_path, 'atac_meta.rds'))

time_windows <- unique(annot$NNv1_time.new)
time_windows <- sort(time_windows)

for (window in time_windows) {
  print(window)
  cells <- annot$cell[annot$NNv1_time.new == window]
  files <- unique(annot$sample[annot$NNv1_time.new == window])
  new_matrix_list <- list()
  
  for (file in files) {
    print(file)

    # Load .mtx.gz file
    mtx_path = paste0(data_path, '/GSE190130_', file, '.peak_matrix.mtx.gz')
    mtx = readMM(mtx_path)
    print("Read mtx")

    columns_path = paste0(data_path, '/GSE190130_', file, '.peak_matrix.columns.txt.gz')
    colnames(mtx) <- readLines(columns_path)
    mtx <- mtx[, colnames(mtx) %in% cells, drop = FALSE]

    rows_path = paste0(data_path, '/GSE190130_', file, '.peak_matrix.rows.txt.gz')
    rownames(mtx) <- readLines(rows_path)

    # Note on coordinates: chr2L_5534_5763 corresponds to line "chr2L   5534    5763" in BED file,
    # i.e. coordinates chr2L:5535-5763 inclusive.
    regions_split <- strsplit(rownames(mtx), '_')
    regions_chrom <- sub("^chr", "", sapply(regions_split, `[`, 1))
    regions_start <- as.integer(sapply(regions_split, `[`, 2))
    regions_end <- as.integer(sapply(regions_split, `[`, 3))
    regions_gr <- GRanges(regions_chrom, IRanges(regions_start + 1L, regions_end))

    loop_cell_list <- list()

    for (i in seq_len(nrow(loops)))
    {
      loop_id <- loops$loop_id[i]
      # print(loop_id)

      mtx_A1 <- mtx[overlapsAny(regions_gr, loops_A1_gr[i]), , drop = FALSE]
      sum_A1 <- colSums(mtx_A1)

      mtx_A2 <- mtx[overlapsAny(regions_gr, loops_A2_gr[i]), , drop = FALSE]
      sum_A2 <- colSums(mtx_A2)

      # cells will have numeric annotations 00, 01, 10, 11
      loop_cell_list[[loop_id]] <- 10 * pmin(sum_A1, 1) + pmin(sum_A2, 1)
    }

    loop_cell_matrix <- do.call(rbind, loop_cell_list)
    colnames(loop_cell_matrix) <- colnames(mtx)
    new_matrix_list[[file]] <- loop_cell_matrix
  }
  
  new_matrix <- do.call(cbind, new_matrix_list)
  stopifnot(setequal(cells, colnames(new_matrix)))
  saveRDS(new_matrix, file = file.path('results/calderon/new_time', paste0('hrs', window, '_NNv1_time_matrix_loops.rds')))
  write.table(new_matrix, file = gzfile(file.path('results/calderon/new_time', paste0('hrs', window, '_NNv1_time_matrix_loops.tsv.gz'))),
    row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)


}
