# The script to prepare the Calderon motif data with adjusted time windows. 
# It uses time windows from the NNv1_time.new column in the metadata file named 'atac.meta.rds'.
# The script reads the RDS files with motif activity for each cell, filters the cells that belong to the time window, 
# and saves the new motif files in the 'results/calderon/new_time' folder in RDS format.


library(Matrix)
library(R.utils)
library(plyr)

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
    file_name <- paste0(file, '_motif_activity')
    mat <- readRDS(file.path(data_path, 'motifs', paste0(file_name, '.rds')))
    mat <- mat[, colnames(mat) %in% cells, drop = FALSE]
    new_matrix_list[[file]] <- mat
  }
  
  new_matrix <- do.call(cbind, new_matrix_list)
  stopifnot(setequal(cells, colnames(new_matrix)))
  saveRDS(new_matrix, file = file.path('results/calderon/new_time', paste0('hrs', window, '_NNv1_time_matrix_motifs.rds')))
  write.table(new_matrix, file = gzfile(file.path('results/calderon/new_time', paste0('hrs', window, '_NNv1_time_matrix_motifs.tsv.gz'))),
    row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)


}
