# A script to prepare the Calderon ATAC-seq data with adjusted time windows. 
# It uses time windows from the NNv1_time.new column in the metadata file named 'atac.meta.rds'.
# The script reads the peak matrix files, filters the cells that belong to the time window, 
# and saves the new peak matrix files in the 'data/muszka/calderon_data/new_time/NN' folder.

import os
import scipy.io
import scipy.sparse
from scipy.io import mmwrite
import pyreadr
import gzip
import pandas as pd
import numpy as np
data_path = os.getcwd().split('src')[0] + 'data/muszka/calderon_data/'
print(data_path)

annot = pyreadr.read_r(data_path + 'atac_meta.rds')[None]

time_windows = set(annot['NNv1_time.new'])
time_windows = sorted(time_windows)
for window in time_windows:
    print(window)
    cells = annot[annot['NNv1_time.new']==window]['cell']
    files = set(annot[annot['NNv1_time.new']==window]['sample'])
    new_matrix = []
    cells_to_save=[]
    for file in files:
        print(file)
        file_name = 'GSE190130_' + file
        peakinfo = pd.read_table(gzip.open(data_path+file_name+'.peak_matrix.rows.txt.gz', 'rt'), header=None)
        cellinfo = pd.read_table(gzip.open(data_path+file_name+'.peak_matrix.columns.txt.gz', 'rt'), header=None)
        
        mat = scipy.io.mmread(data_path+file_name+'.peak_matrix.mtx.gz')
        mat = mat.tocsr()
        mat = pd.DataFrame.sparse.from_spmatrix(mat)
        mat.columns = cellinfo.iloc[:, 0].values
        mat.index = peakinfo.iloc[:, 0].values
        mat = mat[np.intersect1d(cells, mat.columns)]
        cells_to_save.extend(mat.columns.tolist())
        new_matrix.append(mat)
    new_matrix = pd.concat(new_matrix, axis=1)
    cells_to_save = pd.concat([pd.DataFrame(cells_to_save)], axis=1)
    sparse_matrix = scipy.sparse.csr_matrix(new_matrix)

    # Save to gzipped mtx
    with gzip.open(data_path+'new_time/NN/GSE190130_'+window+'_new_timeNN.peak_matrix.mtx.gz', 'wb') as f:
        mmwrite(f, sparse_matrix)
        
    # save peakinfo and cells_to_save to txt.gz
    peakinfo.to_csv(data_path+'new_time/NN/GSE190130_'+window+'_new_timeNN.peak_matrix.rows.txt.gz', sep='\t', header=False, index=False, compression='gzip')
    pd.DataFrame(cells_to_save).to_csv(data_path+'new_time/NN/GSE190130_'+window+'_new_timeNN.peak_matrix.columns.txt.gz', sep='\t', header=False, index=False, compression='gzip')
    

