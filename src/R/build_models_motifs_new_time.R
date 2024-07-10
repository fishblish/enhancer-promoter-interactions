library(R.utils)
library(Matrix)
library(heatmaply)

loops = read.table(paste0(getwd(),'/data/muszka/long_and_short_range_loops_D_mel.tsv'), header = TRUE, sep = '\t')
loops_ids = loops$loop_id

## add ground true loops
loops = loops[c('loop_id', "Dmel_6.8h_Neuroblasts", "Dmel_6.8h_Neurons",   
                "Dmel_6.8h_Glia", "Dmel_10.12h_Neuroblasts",
                "Dmel_10.12h_Neurons", "Dmel_10.12h_Glia",    
                "Dmel_14.16h_Neuroblasts", "Dmel_14.16h_Neurons", 
                "Dmel_14.16h_Glia", "Dmel_16.18h_whole_embryo",
                "Dmel_larval_brain", "Dmel_adult_brain")]
rownames(loops) <- loops$loop_id
loops <- loops[,-1]
loops_data = t(loops)

## prepare foldchanges for motifs
fold_change_dict = {}
motifs_path = paste0(getwd(),'/data/muszka/calderon_data/motifs/new_time/')
calderon_path = paste0(getwd(),'/data/muszka/calderon_data/new_time/')
hrs_list = c('00-02', '02-04', '04-06', '06-08', '08-10', '10-12', '12-14', '14-16', '16-18', '18-20')
for (hrs in hrs_list){
    print(paste0('Computing fold change for ', hrs))
    loops_cells = readRDS(paste0("results/calderon/new_time/open_cells_all_for_loops_",hrs,".rds"))
    file_name = paste0('GSE190130_', hrs, '_new_time')
    motifs <- readRDS(paste0(motifs_path, file_name,'_matrix_motifs.rds'))
    motifs <- motifs[, colSums(is.na(motifs)) == 0]
    file_name = paste0('GSE190130_', hrs, '_new_time.peak_matrix')
    #cell_names = readLines(paste0(calderon_path,file_name, '.columns.txt.gz'))
    columns = read.table(gzfile(paste0(calderon_path,file_name, '.columns.txt.gz')), sep='\t', header=FALSE)
    columns_vector = na.omit(unlist(columns, use.names=FALSE))
    empty_indices = grep("^$", columns_vector)
    cell_names = columns_vector[-empty_indices]

    fold_change = data.frame()
    for (loop_id in loops_ids){
        open_cells = loops_cells[[loop_id]]
        if (length(open_cells) == 0){
            next
        }
        motifs_for_open_cells = motifs[,intersect(open_cells, colnames(motifs))]
        motifs_for_closed_cells = motifs[,intersect(setdiff(cell_names, open_cells), colnames(motifs))]
        mean_open = rowMeans(motifs_for_open_cells)
        mean_closed = rowMeans(motifs_for_closed_cells)
        diff = mean_open-mean_closed
        fold_change <- rbind(fold_change, data.frame(t(diff)))
        rownames(fold_change)[nrow(fold_change)] <- loop_id
    }
    fold_change_dict[[hrs]] = t(fold_change)
}

backup = fold_change_dict
saveRDS(fold_change_dict, paste0(getwd(),'/results/calderon/motifs/new_time/fold_change_dict_10.rds'))
#saveRDS(fold_change_dict, paste0(getwd(),'/results/calderon/motifs/new_time/fold_change_dict.rds'))
## prepare distance
loops = read.table(paste0(getwd(),'/data/muszka/long_and_short_range_loops_D_mel.tsv'), header = TRUE, sep = '\t')
dist = ((loops$x2 - loops$x1) - (loops$y1 - loops$y1))/2
dist = log10((dist))
dist = data.frame(loop_id = loops$loop_id, dist = dist, dist1 = dist)

rownames(dist) = dist$loop_id
dist = dist[,-1]
dist = t(dist)
#dist <- matrix(as.numeric(unlist(dist)), nrow = nrow(dist))

## prepare difference in fold change
fold_change_diff_dict = {}
for (idx in 1:(length(hrs_list)-1)){
    hrs1 = hrs_list[[idx]]
    hrs2 = hrs_list[[idx+1]]
    fold_change_diff_dict[[hrs2]] = fold_change_dict[[hrs2]] - fold_change_dict[[hrs1]]
}

## prepare features
## limit dist to loops that are in the fold change
dist = dist[,match(colnames(fold_change_dict[[hrs_list[[1]]]]), colnames(dist))]

## change features names and order
for (idx in 1:length(hrs_list)){
    hrs = hrs_list[[idx]]
    hrs_name = gsub('-','.',hrs)
    fold_change_dict[[hrs]] = fold_change_dict[[hrs]][,match(colnames(dist), colnames(fold_change_dict[[hrs]]))]
    rownames(fold_change_dict[[hrs]]) = gsub(".1.02$", paste0('_',hrs_name), rownames(fold_change_dict[[hrs]]))
    if (idx>1){
        fold_change_diff_dict[[hrs]] = fold_change_diff_dict[[hrs]][,match(colnames(dist), colnames(fold_change_diff_dict[[hrs]]) )]
        rownames(fold_change_diff_dict[[hrs]]) = gsub(".1.02$", paste0('_',hrs_name,'_diff'), rownames(fold_change_diff_dict[[hrs]]))
    }
}

fold_change_matrix = do.call(rbind, fold_change_dict)
fold_change_diff_matrix = do.call(rbind, fold_change_diff_dict)
str(fold_change_matrix)
str(fold_change_diff_matrix)

features_to_use = rbind(fold_change_matrix, fold_change_diff_matrix, dist)
features_to_use = t(features_to_use)

# build models
library(mlr3)
library(mlr3learners)
library(xgboost)
library(mlr3measures)
library(mlr3viz)
library(ggplot2)
library(precrec)
library(gridExtra)

plots = list()


for (i in 1:dim(loops_data)[1]){
    loops_data_row = loops_data[i,]
    loops_data_row = loops_data_row[!is.na(loops_data_row)]
    loops_data_row = loops_data_row[intersect(names(loops_data_row), rownames(features_to_use))]
    features_to_use = features_to_use[match(names(loops_data_row), rownames(features_to_use)),]
    features_to_use = features_to_use[match(names(loops_data_row), rownames(features_to_use)),]
    features_to_use = as.data.frame(features_to_use)
    loops_data_row = as.factor(loops_data_row)
    features_to_use$target <- loops_data_row

    task <- TaskClassif$new(id = "task", backend = features_to_use, target = "target")

    learners <- list(
    lrn("classif.xgboost", predict_type = "prob"),
    lrn("classif.ranger", predict_type = "prob"), #random forest
    lrn("classif.svm", predict_type = "prob")
    )
    resampling <- rsmp("cv", folds = 5)

    bmr <- benchmark(benchmark_grid(task, learners, resampling))


    auc_scores = paste(round(bmr$aggregate(msr("classif.auc"))[['classif.auc']],3), collapse = ", ")

    p1 <- autoplot(bmr, type = "roc", measure = msr("classif.auc")) +
    ggtitle(paste(rownames(loops_data)[i],"\nAUC ROC score ", auc_scores)) +
    theme(plot.title = element_text(hjust = 0.5, size = 14))

    plots[[length(plots)+1]] = p1
}
p=grid.arrange(grobs = plots, ncol = 2)

# Save the combined plot to a file
ggsave(filename = "results/calderon/motifs/new_time/mlr3_models_benchmark_new_time_10_windows.pdf", plot = p, width = 10, height = 20, dpi = 300)


