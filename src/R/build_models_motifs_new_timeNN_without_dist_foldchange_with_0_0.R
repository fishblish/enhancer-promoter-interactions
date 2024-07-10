# The script with models for predicting enhancer-promoter interactions based on motif activities and distance
# It is adapted to the new time windows from the Calderon data.
# It consist of preparing features, building models with mlr3, plotting results
# and checking model performance by distance and tissue type.
# It has also section which generates boruta importance heatmaps.

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
motifs_path = paste0(getwd(),'/data/muszka/calderon_data/motifs/new_time/NN/')
calderon_path = paste0(getwd(),'/data/muszka/calderon_data/new_time/NN/')
#hrs_list = c('06-08', '10-12', '14-16')


# this loop takes some time so result is saved in rds file
hrs_list = c('00-02', '02-04', '04-06', '06-08', '08-10', '10-12', '12-14', '14-16', '16-18')
for (hrs in hrs_list){
    print(paste0('Computing fold change for ', hrs))
    loops_cells = readRDS(paste0("results/calderon/new_time/NN/open_cells_all_for_loops_",hrs,"_new_timeNN_faster.rds"))
    loops_cells_0_0 = readRDS(paste0("results/calderon/new_time/NN/0_0_cells_all_for_loops_",hrs,"_new_timeNN.rds"))
    file_name = paste0('GSE190130_', hrs, '_new_timeNN')
    motifs <- readRDS(paste0(motifs_path, file_name,'_matrix_motifs.rds'))
    motifs <- motifs[, colSums(is.na(motifs)) == 0]
    file_name = paste0('GSE190130_', hrs, '_new_timeNN.peak_matrix')
    #cell_names = readLines(paste0(calderon_path,file_name, '.columns.txt.gz'))
    columns = read.table(gzfile(paste0(calderon_path,file_name, '.columns.txt.gz')), sep='\t', header=FALSE)
    cell_names = na.omit(unlist(columns, use.names=FALSE))

    fold_change = data.frame()
    for (loop_id in loops_ids){
        open_cells = loops_cells[[loop_id]]
        closed_cells = loops_cells_0_0[[loop_id]]
        if (length(open_cells) == 0){
            next
        }
        motifs_for_open_cells = motifs[,intersect(open_cells, colnames(motifs))]
        if (is.vector(motifs_for_open_cells)) {
            motifs_for_open_cells = data.frame(motifs_for_open_cells)
        }
        motifs_for_closed_cells = motifs[,intersect(closed_cells, colnames(motifs))]
        mean_open = rowMeans(motifs_for_open_cells)
        mean_closed = rowMeans(motifs_for_closed_cells)
        diff = mean_open-mean_closed
        fold_change <- rbind(fold_change, data.frame(t(diff)))
        rownames(fold_change)[nrow(fold_change)] <- loop_id
    }
    fold_change_dict[[hrs]] = t(fold_change)
}

backup = fold_change_dict
fold_change_dict = fold_change_dict[hrs_list]
saveRDS(fold_change_dict, paste0(getwd(),'/results/calderon/motifs/new_time/fold_change_with_0_0_dict_NN_9windows.rds'))

fold_change_dict = readRDS(paste0(getwd(),'/results/calderon/motifs/new_time/fold_change_with_0_0_dict_NN_9windows.rds'))
loops = read.table(paste0(getwd(),'/data/muszka/long_and_short_range_loops_D_mel.tsv'), header = TRUE, sep = '\t')

## prepare difference in fold change
fold_change_diff_dict = {}
for (idx in 1:(length(hrs_list)-1)){
    hrs1 = hrs_list[[idx]]
    hrs2 = hrs_list[[idx+1]]
    print(paste0('Computing fold change diff for ', hrs1, ' and ', hrs2))
    avail_loops = intersect(colnames(fold_change_dict[[hrs1]]), colnames(fold_change_dict[[hrs2]]))
    fold_change_diff_dict[[hrs2]] = fold_change_dict[[hrs2]][,avail_loops] - fold_change_dict[[hrs1]][,avail_loops]
}

## prepare features

## change features names and order
for (idx in 1:length(hrs_list)){
    hrs = hrs_list[[idx]]
    hrs_name = gsub('-','.',hrs)
    fold_change_dict[[hrs]] = fold_change_dict[[hrs]][,match(colnames(fold_change_dict[['00-02']]), colnames(fold_change_dict[[hrs]]))]
    rownames(fold_change_dict[[hrs]]) = gsub(".1.02$", paste0('_',hrs_name), rownames(fold_change_dict[[hrs]]))
    if (idx>1){
        fold_change_diff_dict[[hrs]] = fold_change_diff_dict[[hrs]][,match(colnames(fold_change_dict[['00-02']]), colnames(fold_change_diff_dict[[hrs]]) )]
        rownames(fold_change_diff_dict[[hrs]]) = gsub(".1.02$", paste0('_',hrs_name,'_diff'), rownames(fold_change_diff_dict[[hrs]]))
    }
}

fold_change_matrix = do.call(rbind, fold_change_dict)
fold_change_diff_matrix = do.call(rbind, fold_change_diff_dict)
str(fold_change_matrix)
str(fold_change_diff_matrix)

features_to_use = rbind(fold_change_matrix, fold_change_diff_matrix)
features_to_use = t(features_to_use)



# change motif names
dict_motifs = read.csv('/home/jbartczak/enhancer-promoter-interactions/data/muszka/calderon_data/motifs/dict_motifs.csv')
temp = dict_motifs$names
dict_motifs = dict_motifs$values
names(dict_motifs) = temp
features = features_to_use[,which(colnames(features_to_use) != 'target')]
colnames(features)
split_names = strsplit(colnames(features), split = '_')

translated_names = sapply(split_names, function(x) {
  if (x[1] %in% names(dict_motifs)) {
    x[1] = dict_motifs[x[1]]
  }
  return(x)
})

new_colnames = sapply(translated_names, function(x) paste(x, collapse = '_'))

colnames(features) = new_colnames
str(features)

# omit other species motifs
features = features[,!grepl('other_species', colnames(features))]


# build models
library(mlr3)
library(mlr3learners)
library(xgboost)
library(mlr3measures)
library(mlr3viz)
library(ggplot2)
library(precrec)
library(gridExtra)
library(stringr)


colnames(features) = make.names(colnames(features))
features_to_use_dict = list()

results_to_compare = data.frame()

feature_combinations = list()
feature_combinations[['2prev']] = "features[,grepl(paste0(time_window,'$'), colnames(features)) | grepl(paste0(prev_time_window,'$'), colnames(features))]"
feature_combinations[['2diff']] = "features[,grepl(paste0(time_window,'$'), colnames(features)) | grepl(paste0(prev_time_window,'_diff$'), colnames(features))]"
feature_combinations[['3prevprev']] = "features[,grepl(paste0(time_window,'$'), colnames(features)) | grepl(paste0(prev_time_window,'$'), colnames(features)) | grepl(paste0(prev_prev_time_window,'$'), colnames(features)) ]"
feature_combinations[['3']] = "features[,grepl(paste0(time_window,'$'), colnames(features)) | grepl(paste0(time_window,'_diff$'), colnames(features)) | grepl(paste0(prev_time_window,'$'), colnames(features)) ]"
feature_combinations[['4']] = "features[,grepl(paste0(time_window,'$'), colnames(features)) | grepl(paste0(time_window,'_diff$'), colnames(features)) | grepl(paste0(prev_time_window,'$'), colnames(features)) | grepl(paste0(prev_time_window,'_diff$'), colnames(features)) ]"
feature_combinations[['5']] = "features[,grepl(paste0(time_window,'$'), colnames(features)) | grepl(paste0(time_window,'_diff$'), colnames(features)) | grepl(paste0(prev_time_window,'$'), colnames(features)) | grepl(paste0(prev_time_window,'_diff$'), colnames(features)) | grepl(paste0(prev_prev_time_window,'$'), colnames(features)) ]"
feature_combinations[['6']] = "features[,grepl(paste0(time_window,'$'), colnames(features)) | grepl(paste0(time_window,'_diff$'), colnames(features)) | grepl(paste0(prev_time_window,'$'), colnames(features)) | grepl(paste0(prev_time_window,'_diff$'), colnames(features)) | grepl(paste0(prev_prev_time_window,'$'), colnames(features)) | grepl(paste0(prev_prev_time_window,'_diff$'), colnames(features)) ]"
feature_combinations[['all']] = "features"

for (comb in names(feature_combinations)){
    plots = list()
    print(comb)
    for (i in 1:nrow(loops_data)){
        if (grepl(".+_.+\\.+.+h_.+", rownames(loops_data)[i])){
            time_window = str_extract(rownames(loops_data)[i], "\\d+\\.\\d+h")
            print(time_window)
            if (time_window=='6.8h'){
                time_window = '06.08'
                prev_time_window = '04.06'
                prev_prev_time_window = '02.04'
            }
            if (time_window=='10.12h'){
                time_window = '10.12'
                prev_time_window = '08.10'
                prev_prev_time_window = '06.08'
            }
            if (time_window=='14.16h'){
                time_window = '14.16'
                prev_time_window = '12.14'
                prev_prev_time_window = '10.12'
            }
            if (time_window=='16.18h'){
                time_window = '16.18'
                prev_time_window = '14.16'
                prev_prev_time_window = '12.14'
            }
            print(time_window)
            features_to_use = eval(parse(text = feature_combinations[[comb]]))
            #features_to_use_dict[[time_window]] = features_to_use
        }
        else{
            next
        }
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
        results_to_compare[comb, rownames(loops_data)[i]]=round(bmr$aggregate(msr("classif.auc"))[['classif.auc']],3)[2]
    }
    
    p=grid.arrange(grobs = plots, ncol = 2)

    ggsave(filename = paste0("results/calderon/motifs/without_dist_0_0/mlr3_models_benchmark_omitted_other_species_NN_9windows_different_time_",comb,"_without_dist_0_0.pdf"), plot = p, width = 10, height = 20, dpi = 300)
    }
# module with concatenated predicted vectors among time windows
# one model per tissue
saveRDS(results_to_compare, paste0(getwd(),'/results/calderon/motifs/without_dist_0_0/results_to_compare_omitted_other_species_NN_9windows_without_dist_0_0.rds'))
results_to_compare = readRDS(paste0(getwd(),'/results/calderon/motifs/without_dist_0_0/results_to_compare_omitted_other_species_NN_9windows_without_dist_0_0.rds'))

# plot results comparison as bar plot
transposed_data <- t(results_to_compare)
transposed_data <- data.frame(Category = rownames(transposed_data), transposed_data)
library(reshape2)
melted_data <- melt(transposed_data, id.vars = "Category")

# Create the barplot
library(ggplot2)
colors <- c("X2prev" = rgb(1, 0, 0),  # red
            "X2diff" = rgb(0, 0.5, 0),  # green
            "X3prevprev" = rgb(0, 0, 1),  # blue
            "X3" = rgb(0.75, 0.75, 0),  # yellow
            "X4" = rgb(0.75, 0, 0.75),  # magenta
            "X5" = rgb(0, 0.75, 0.75),  # cyan
            "X6" = rgb(0.75, 0.75, 0.75),  # grey
            "all" = rgb(0, 0, 0))  # black

# Create the barplot
p <- ggplot(melted_data, aes(x = Category, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Comparison of AUC ROC scores for different feature combinations, models without distance") +
  labs(x = "Category", y = "Value", fill = "Group") +
  coord_cartesian(ylim = c(0.65, 0.85)) +
  scale_fill_manual(values = colors)

ggsave(filename = paste0("results/calderon/motifs/without_dist_0_0/plotted_results_to_compare_omitted_other_species_NN_9windows_without_dist_0_0.pdf"), plot = p, width = 15, height = 10, dpi = 300)


# module with concatenated predicted vectors among time windows
# one model per tissue

neuroblasts = c('Dmel_6.8h_Neuroblasts','Dmel_10.12h_Neuroblasts','Dmel_14.16h_Neuroblasts')
neurons = c('Dmel_6.8h_Neurons','Dmel_10.12h_Neurons','Dmel_14.16h_Neurons')
glia = c('Dmel_6.8h_Glia','Dmel_10.12h_Glia','Dmel_14.16h_Glia')
tissue_names = c('neuroblast', 'neuron', 'glia')

neuroblast_data = loops_data[neuroblasts,]
neuron_data = loops_data[neurons,]
glia_data = loops_data[glia,]

neuroblast_vector = c(apply(neuroblast_data, 1, as.vector))
neuron_vector = c(apply(neuron_data, 1, as.vector))
glia_vector = c(apply(glia_data, 1, as.vector))

col_names = colnames(neuroblast_data)

col_names_1 = paste('06.08', col_names, sep="_")
col_names_2 = paste('10.12', col_names, sep="_")
col_names_3 = paste('14.16', col_names, sep="_")

names(neuroblast_vector) = c(col_names_1, col_names_2, col_names_3)
names(neuron_vector) = c(col_names_1, col_names_2, col_names_3)
names(glia_vector) = c(col_names_1, col_names_2, col_names_3)

for (name in names(features_to_use_dict)){
    rownames(features_to_use_dict[[name]]) = paste(name, rownames(features_to_use_dict[[name]]), sep="_")
}

features_to_use_dict = list()
for (comb in names(feature_combinations)){
    plots = list()
    print(comb)
    for (i in 1:nrow(loops_data)){
        if (grepl(".+_.+\\.+.+h_.+", rownames(loops_data)[i])){
            time_window = str_extract(rownames(loops_data)[i], "\\d+\\.\\d+h")
            print(time_window)
            if (time_window=='6.8h'){
                time_window = '06.08'
                prev_time_window = '04.06'
                prev_prev_time_window = '02.04'
            }
            if (time_window=='10.12h'){
                time_window = '10.12'
                prev_time_window = '08.10'
                prev_prev_time_window = '06.08'
            }
            if (time_window=='14.16h'){
                time_window = '14.16'
                prev_time_window = '12.14'
                prev_prev_time_window = '10.12'
            }
            if (time_window=='16.18h'){
                time_window = '16.18'
                prev_time_window = '14.16'
                prev_prev_time_window = '12.14'
            }
            print(time_window)
            features_to_use = eval(parse(text = feature_combinations[[comb]]))
            features_to_use_dict[[comb]][[time_window]] = features_to_use
            rownames(features_to_use_dict[[comb]][[time_window]]) = paste(time_window, rownames(features_to_use_dict[[comb]][[time_window]]), sep="_")
        }
        else{
            next
        }
    }
}

results_to_compare_per_tissue = data.frame()
for (comb in names(feature_combinations)){
    plots = list()
    print(comb)
    for (i in 1:3){
        print(tissue_names[i])
        loops_data_row = list(neuroblast_vector, neuron_vector, glia_vector)[[i]]
        features_to_use = rbind(features_to_use_dict[[comb]][['06.08']], features_to_use_dict[[comb]][['10.12']], features_to_use_dict[[comb]][['14.16']])
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
        ggtitle(paste(tissue_names[i],"\nAUC ROC score ", auc_scores)) +
        theme(plot.title = element_text(hjust = 0.5, size = 14))

        plots[[length(plots)+1]] = p1
        results_to_compare_per_tissue[comb, tissue_names[i]]=round(bmr$aggregate(msr("classif.auc"))[['classif.auc']],3)[2]
    }
    p=grid.arrange(grobs = plots, ncol = 1)

    ggsave(filename = paste0("results/calderon/motifs/without_dist_0_0/mlr3_models_benchmark_omitted_other_species_NN_9windows_different_time_",comb,"_without_dist_0_0_per_tissue.pdf"), plot = p, width = 15, height = 10, dpi = 300)
}
library(ggplot2)
saveRDS(results_to_compare_per_tissue, paste0(getwd(),'/results/calderon/motifs/without_dist_0_0/results_to_compare_per_tissue_omitted_other_species_NN_9windows_without_dist_0_0.rds'))
results_to_compare_per_tissue_x = readRDS(paste0(getwd(),'/results/calderon/motifs/without_dist_0_0/results_to_compare_per_tissue_omitted_other_species_NN_9windows_without_dist_0_0.rds'))
# plot results comparison as bar plot
transposed_data <- t(results_to_compare_per_tissue)
transposed_data <- data.frame(Category = rownames(transposed_data), transposed_data)
melted_data <- melt(transposed_data, id.vars = "Category")
melted_data$Category <- factor(melted_data$Category, levels = colnames(results_to_compare_per_tissue))


p <- ggplot(melted_data, aes(x = Category, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Comparison of AUC ROC scores for different feature combinations, models without distance, per tissue") +
  labs(x = "Category", y = "Value", fill = "Group") +
  coord_cartesian(ylim = c(0.65, 0.85)) +
  scale_fill_manual(values = colors)

ggsave(filename = paste0("results/calderon/motifs/without_dist_0_0/plotted_results_to_compare_per_tissue_omitted_other_species_NN_9windows_without_dist_0_0.pdf"), plot = p, width = 10, height = 10, dpi = 300)


## boruta importance heatmaps
# boruta importance - for set of 4 features
selected_comb = '4'
library(Boruta)
set.seed(111)
boruta_results = list()
for (i in 1:3){
    print(i)
    loops_data_row = list(neuroblast_vector, neuron_vector, glia_vector)[[i]]
    features_to_use = rbind(features_to_use_dict[[selected_comb]][['06.08']], features_to_use_dict[[selected_comb]][['10.12']], features_to_use_dict[[selected_comb]][['14.16']])
    loops_data_row = loops_data_row[!is.na(loops_data_row)]
    loops_data_row = loops_data_row[intersect(names(loops_data_row), rownames(features_to_use))]
    features_to_use = features_to_use[match(names(loops_data_row), rownames(features_to_use)),]
    features_to_use = features_to_use[match(names(loops_data_row), rownames(features_to_use)),]
    features_to_use = as.data.frame(features_to_use)
    loops_data_row = as.factor(loops_data_row)
    res = Boruta(features_to_use, loops_data_row, doTrace = 1, maxRuns = 200) 
    importances = res$ImpHistory
    mean = importances[,'shadowMean']
    importances=importances-mean
    importances = colMeans(importances)
    boruta_results[[i]] = importances[!grepl('shadow', names(importances))]

}

# and for one model for all tissues 
loops_data_row = c(neuroblast_vector, neuron_vector, glia_vector)
features_to_use = rbind(features_to_use_dict[[selected_comb]][['06.08']], features_to_use_dict[[selected_comb]][['10.12']], features_to_use_dict[[selected_comb]][['14.16']],
                        features_to_use_dict[[selected_comb]][['06.08']], features_to_use_dict[[selected_comb]][['10.12']], features_to_use_dict[[selected_comb]][['14.16']],
                        features_to_use_dict[[selected_comb]][['06.08']], features_to_use_dict[[selected_comb]][['10.12']], features_to_use_dict[[selected_comb]][['14.16']])
loops_data_row = loops_data_row[!is.na(loops_data_row)]
loops_data_row = loops_data_row[intersect(names(loops_data_row), rownames(features_to_use))]
features_to_use = features_to_use[match(names(loops_data_row), rownames(features_to_use)),]
features_to_use = features_to_use[match(names(loops_data_row), rownames(features_to_use)),]
features_to_use = as.data.frame(features_to_use)
loops_data_row = as.factor(loops_data_row)
res = Boruta(features_to_use, loops_data_row, doTrace = 1, maxRuns = 200) 
importances = res$ImpHistory
mean = importances[,'shadowMean']
importances=importances-mean
importances = colMeans(importances)
boruta_results[[4]] = importances[!grepl('shadow', names(importances))]

# prepare heatmap
data_matrix <- do.call(rbind, lapply(boruta_results, function(x) as.numeric(x)))

data_df <- as.data.frame(data_matrix)
rownames(data_df) <- c('neuroblast', 'neuron', 'glia', 'all')
colnames(data_df) <- colnames(features_to_use)
saveRDS(data_df, paste0(getwd(),'/results/calderon/motifs/without_dist_0_0/boruta_results_omitted_other_species_NN_features4.rds'))
data_df = readRDS(paste0(getwd(),'/results/calderon/motifs/without_dist_0_0/boruta_results_omitted_other_species_NN_features4.rds'))

a=colnames(data_df)[grepl('_04.06$', colnames(data_df))]
b=colnames(data_df)[grepl('_06.08$', colnames(data_df))]
c=colnames(data_df)[grepl('_04.06_diff$', colnames(data_df))]
d=colnames(data_df)[grepl('_06.08_diff$', colnames(data_df))]

#change colnames, they comes from first row of features_to_use

colnames(data_df) = gsub("_04.06", "_prev", colnames(data_df))
colnames(data_df) = gsub("_06.08", "_current", colnames(data_df))

boruta_results_prev = data_df[,grepl("_prev$", colnames(data_df)) | grepl("dist", colnames(data_df))]
boruta_results_curr = data_df[,grepl("_current$", colnames(data_df)) | grepl("dist", colnames(data_df))]
boruta_results_prev_diff = data_df[,grepl("_prev_diff$", colnames(data_df)) | grepl("dist", colnames(data_df))]
boruta_results_curr_diff = data_df[,grepl("_current_diff$", colnames(data_df)) | grepl("dist", colnames(data_df))]

colnames(boruta_results_prev) = gsub("_prev", "", colnames(boruta_results_prev))
colnames(boruta_results_curr) = gsub("_current", "", colnames(boruta_results_curr))
colnames(boruta_results_prev_diff) = gsub("_prev_diff", "", colnames(boruta_results_prev_diff))
colnames(boruta_results_curr_diff) = gsub("_current_diff", "", colnames(boruta_results_curr_diff))

rownames(boruta_results_prev) = c('Neuroblasts_prev', 'Neurons_prev', 'Glia_prev', 'All_prev')
rownames(boruta_results_curr) = c('Neuroblasts_curr', 'Neurons_curr', 'Glia_curr', 'All_curr')
rownames(boruta_results_prev_diff) = c('Neuroblasts_prev_diff', 'Neurons_prev_diff', 'Glia_prev_diff', 'All_prev_diff')
rownames(boruta_results_curr_diff) = c('Neuroblasts_curr_diff', 'Neurons_curr_diff', 'Glia_curr_diff', 'All_curr_diff')

heatmap_data = rbind(boruta_results_prev, boruta_results_curr, boruta_results_prev_diff, boruta_results_curr_diff)

my_palette <- colorRampPalette(c("#f9f8f7","#ADD8E6", "#1ebb1e"))(n = 100) 

library(gplots)
df=heatmap_data

# Remove columns where all values are -Inf
df = df[, !apply(heatmap_data == -Inf, 2, all)]
df = df[rownames(df)[order(rownames(df))],]
df[df == -Inf] <- 0
# Sort columns by their maximum values
df <- df[, order(sapply(df, max), decreasing = TRUE)]
pdf("results/calderon/motifs/without_dist_0_0/boruta_heatmap_omitted_other_species_NN_features4.pdf", width = 35, height = 10) 



# Create the heatmap
heatmap.2(as.matrix(df), 
        Rowv = NA, 
        Colv = NA, 
        scale = "none", 
        trace = "none", 
        dendrogram = "none", 
        col = my_palette, 
        keysize = 1.5, 
        symkey = FALSE, 
        density.info = "none", 
        tracecol = NA,
        key = FALSE,
        labRow = rownames(df),
        labCol = colnames(df),
        cexRow = 0.8,
        cexCol = 0.8,
        margins = c(15, 15))

legend("left", legend = color_key$labels, fill = colorRampPalette(c("#f9f8f7","#ADD8E6", "#1ebb1e"))(length(breaks)))

dev.off()
