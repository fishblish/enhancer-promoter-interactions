# The script to prepare the Calderon motif data with adjusted time windows. 
# It uses time windows from the NNv1_time.new column in the metadata file named 'atac.meta.rds'.
# The script reads the RDS files with motif activity for each cell, filters the cells that belong to the time window, 
# and saves the new motif files in the 'data/muszka/calderon_data/motifs/new_time/NN' folder in RDS format.


library(Matrix)
library(R.utils)
library(plyr)

data_path <- '/home/jbartczak/enhancer-promoter-interactions/data/muszka/calderon_data'
print(data_path)

annot <- readRDS(file.path(data_path, 'atac_meta.rds'))

time_windows <- unique(annot$NNv1_time.new)
time_windows <- sort(time_windows)

for (window in c('18-20')) {
  print(window)
  cells <- annot$cell[annot$NNv1_time.new == window]
  files <- unique(annot$sample[annot$NNv1_time.new == window])
  new_matrix_list <- list()
  cells_to_save <- list()
  
  for (file in files) {
    print(file)
    file_name <- paste0(file, '_motif_activity')
    mat <- readRDS(file.path(data_path, 'motifs', paste0(file_name, '.rds')))
    mat <- mat[, colnames(mat) %in% cells]
    cells_to_save <- c(cells_to_save, colnames(mat))
    new_matrix_list[[file]] <- as(mat, "CsparseMatrix")
  }
  
  new_matrix <- do.call(cbind, new_matrix_list)
  saveRDS(new_matrix, file = file.path(data_path, paste0('motifs/new_time/NN/GSE190130_', window, '_new_timeNN_matrix_motifs.rds')))


}
