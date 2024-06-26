# The script to select cells that are considered as open for each loop from long_and_short_range_loops_D_mel.tsv file.
# Cell is open if it has non-zero signal in both anchors of the loop.
# More precisely, for both anchors of the loop there is a region in the peak matrix that overlaps with the anchor and has non-zero signal in the cell.
# It works on calderon data with adjusted time windows.

library(R.utils)
library(Matrix)
library(parallel)  
library(pbmcapply)

hrs = "18-20"
print(hrs)
path = '/home/jbartczak/enhancer-promoter-interactions/data/muszka/calderon_data/new_time/NN/'

# Load .mtx.gz file
file_name = paste0('GSE190130_', hrs, '_new_timeNN')
print(file_name)
mtx_path = paste0(path, file_name, '.peak_matrix.mtx.gz')
mtx = readMM(mtx_path)
print("Read mtx")

columns = read.table(gzfile(paste0(path, file_name, '.peak_matrix.columns.txt.gz')), sep='\t', header=FALSE)
columns_vector = na.omit(unlist(columns, use.names=FALSE))
columns_vector = columns_vector[columns_vector != ""]
colnames(mtx) = columns_vector

rownames(mtx) = readLines(paste0(path, file_name, '.peak_matrix.rows.txt.gz'))

loops = read.table('/home/jbartczak/enhancer-promoter-interactions/data/muszka/long_and_short_range_loops_D_mel.tsv', header = TRUE, sep = '\t')

regions = data.frame(region = rownames(mtx))
regions_split = strsplit(regions$region, split = '_')
regions$chromosome = sapply(regions_split, `[`, 1)
regions$start = as.numeric(sapply(regions_split, `[`, 2))
regions$end = as.numeric(sapply(regions_split, `[`, 3))

if_intersect = function(region1, region2) {
    region1$chromosome == paste0('chr', region2['chromosome']) && 
    region1$start <= as.numeric(region2['end']) && 
    region1$end >= as.numeric(region2['start'])
}

# For each loop create list of cells which have non-zero signal in both anchors
loops_ids = loops$loop_id
len_loops = length(loops_ids)

process_loop = function(loop_id, loops, regions, mtx) {
    loop = loops[loops$loop_id == loop_id, ]
    anchor1 = c(loop$x1, loop$x2, loop$chr1)
    anchor2 = c(loop$y1, loop$y2, loop$chr2)
    names(anchor1) = c('start', 'end', 'chromosome')
    names(anchor2) = c('start', 'end', 'chromosome')

    left_region = regions[sapply(1:nrow(regions), function(i) if_intersect(regions[i, ], anchor1)), ]
    right_region = regions[sapply(1:nrow(regions), function(i) if_intersect(regions[i, ], anchor2)), ]

    if (nrow(left_region) > 0 && nrow(right_region) > 0) {
        left_cells = colnames(mtx)[apply(mtx[left_region$region, , drop = FALSE], 2, function(x) any(x > 0))]
        right_cells = colnames(mtx)[apply(mtx[right_region$region, , drop = FALSE], 2, function(x) any(x > 0))]
        common_cells = intersect(left_cells, right_cells)
    } else {
        common_cells = character(0)
    }

    return(list(loop_id = loop_id, common_cells = common_cells))
}

# Wrapper function to include print statements
wrapper_process_loop = function(loop_id, loops, regions, mtx, index, len_loops, hrs) {
    result = process_loop(loop_id, loops, regions, mtx)
    print(paste0("Step ", index, " out of ", len_loops, ": ", loop_id, ' ', hrs))
    return(result)
}

# Use parallel processing
num_cores = 8
print(paste0('There will be used ', num_cores, ' cores.'))
results = pbmclapply(seq_along(loops_ids), function(i) {
    wrapper_process_loop(loops_ids[i], loops, regions, mtx, i, len_loops, hrs)
}, mc.cores = num_cores)

loops_cells = setNames(lapply(results, `[[`, "common_cells"), loops_ids)
saveRDS(loops_cells, file = paste0("/home/jbartczak/enhancer-promoter-interactions/results/calderon/new_time/NN/open_cells_all_for_loops_", hrs, "_new_timeNN_faster.rds"))

