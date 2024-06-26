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


hrs_list = c('00-02', '02-04', '04-06', '06-08', '08-10', '10-12', '12-14', '14-16', '16-18')
for (hrs in hrs_list){
    print(paste0('Computing fold change for ', hrs))
    loops_cells = readRDS(paste0("results/calderon/new_time/NN/open_cells_all_for_loops_",hrs,"_new_timeNN_faster.rds"))
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
        if (length(open_cells) == 0){
            next
        }
        motifs_for_open_cells = motifs[,intersect(open_cells, colnames(motifs))]
        if (is.vector(motifs_for_open_cells)) {
            motifs_for_open_cells = data.frame(motifs_for_open_cells)
        }
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
fold_change_dict = fold_change_dict[hrs_list]
#saveRDS(fold_change_dict, paste0(getwd(),'/results/calderon/motifs/new_time/fold_change_dict_NN_9windows.rds'))
#saveRDS(fold_change_dict, paste0(getwd(),'/results/calderon/motifs/new_time/fold_change_dict.rds'))
## prepare distance
fold_change_dict = readRDS(paste0(getwd(),'/results/calderon/motifs/new_time/fold_change_dict_NN_9windows.rds'))
loops = read.table(paste0(getwd(),'/data/muszka/long_and_short_range_loops_D_mel.tsv'), header = TRUE, sep = '\t')
dist = ((loops$y1 + loops$y2) - (loops$x1 + loops$x2))/2
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
    print(paste0('Computing fold change diff for ', hrs1, ' and ', hrs2))
    avail_loops = intersect(colnames(fold_change_dict[[hrs1]]), colnames(fold_change_dict[[hrs2]]))
    fold_change_diff_dict[[hrs2]] = fold_change_dict[[hrs2]][,avail_loops] - fold_change_dict[[hrs1]][,avail_loops]
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

plots = list()
colnames(features) = make.names(colnames(features))
features_to_use_dict = list()

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
        features_to_use = features[,grepl(paste0(time_window,'$'), colnames(features)) | grepl("dist", colnames(features)) | grepl(paste0(time_window,"_diff$"), colnames(features)) | grepl(paste0(prev_time_window,'$'), colnames(features)) | grepl(paste0(prev_time_window,"_diff$"), colnames(features)) ]
        features_to_use_dict[[time_window]] = features_to_use
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
}
p=grid.arrange(grobs = plots, ncol = 2)

ggsave(filename = "results/calderon/motifs/mlr3_models_benchmark_omitted_other_species_NN_9windows_different_time_features4.pdf", plot = p, width = 10, height = 20, dpi = 300)

# module with concatenated predicted vectors among time windows
# one model per tissue

neuroblasts = c('Dmel_6.8h_Neuroblasts','Dmel_10.12h_Neuroblasts','Dmel_14.16h_Neuroblasts')
neurons = c('Dmel_6.8h_Neurons','Dmel_10.12h_Neurons','Dmel_14.16h_Neurons')
glia = c('Dmel_6.8h_Glia','Dmel_10.12h_Glia','Dmel_14.16h_Glia')

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


plots = list()
tissue_names = c('neuroblast', 'neuron', 'glia')
for (i in 1:3){
    loops_data_row = list(neuroblast_vector, neuron_vector, glia_vector)[[i]]
    features_to_use = rbind(features_to_use_dict[['06.08']], features_to_use_dict[['10.12']], features_to_use_dict[['14.16']])
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
}
p=grid.arrange(grobs = plots, ncol = 1)

ggsave(filename = "results/calderon/motifs/mlr3_models_benchmark_omitted_other_species_NN_9windows_different_time_features4_per_tissue.pdf", plot = p, width = 10, height = 20, dpi = 300)



## boruta importance heatmaps
# boruta importance - for set of 4 features
library(Boruta)
set.seed(111)
boruta_results = list()
for (i in 1:3){
    print(i)
    loops_data_row = list(neuroblast_vector, neuron_vector, glia_vector)[[i]]
    features_to_use = rbind(features_to_use_dict[['06.08']], features_to_use_dict[['10.12']], features_to_use_dict[['14.16']])
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
features_to_use = rbind(features_to_use_dict[['06.08']], features_to_use_dict[['10.12']], features_to_use_dict[['14.16']],
                        features_to_use_dict[['06.08']], features_to_use_dict[['10.12']], features_to_use_dict[['14.16']],
                        features_to_use_dict[['06.08']], features_to_use_dict[['10.12']], features_to_use_dict[['14.16']])
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


boruta_results_f1 = data_df[,grepl("_04.06$", colnames(data_df)) | grepl("dist", colnames(data_df))]
boruta_results_f2 = data_df[,grepl("_06.08$", colnames(data_df)) | grepl("dist", colnames(data_df))]
boruta_results_f3 = data_df[,grepl("_04.06_diff$", colnames(data_df)) | grepl("dist", colnames(data_df))]
boruta_results_f4 = data_df[,grepl("_06.08_diff$", colnames(data_df)) | grepl("dist", colnames(data_df))]

#change colnames by replacement
for (i in 1:4){
    df = get(paste0("boruta_results_f", i))
    colnames(df) = gsub("04.06", "prev", colnames(df))
    colnames(df) = gsub("06.08", "current", colnames(df))
    assign(paste0("boruta_results_f", i), df)
}
my_palette <- colorRampPalette(c("#f9f8f7","#ADD8E6", "#1ebb1e"))(n = 100) 

pdf("results/calderon/motifs/boruta_heatmap_omitted_other_species_NN_features4.pdf", width = 50, height = 8) 
library(gplots)
for(i in 1:4) {
  df <- get(paste0("boruta_results_f", i))
  df = df[,colnames(df)[order(colnames(df))]]
  df[df == -Inf] <- 0

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
  
  breaks <- seq(min(df), max(df), length.out = 5)
  
  labels <- round(breaks, 2)  
  
  color_key <- list(at = breaks, labels = labels, space = "bottom")
  
  legend("left", legend = color_key$labels, fill = colorRampPalette(c("#f9f8f7","#ADD8E6", "#1ebb1e"))(length(breaks)))
}

dev.off()



# check model performance by distance
# one model per tissue
# plot dist
loops = read.table(paste0(getwd(),'/data/muszka/long_and_short_range_loops_D_mel.tsv'), header = TRUE, sep = '\t')
dist = ((loops$y2 + loops$y1) - (loops$x2 + loops$x1))/2
dist = log10((dist))
dist_distribution = hist(dist, breaks=40,main = 'Distance distribution', xlab = 'Distance')
#save it
pdf("results/calderon/motifs/distance_distribution.pdf", width = 10, height = 10)
plot(dist_distribution)
dev.off()
short_loops = loops[dist<5.5,'loop_id']
long_loops = loops[dist>=5.5,'loop_id']
# 347 short loops, 70 long loops


neuroblasts = c('Dmel_6.8h_Neuroblasts','Dmel_10.12h_Neuroblasts','Dmel_14.16h_Neuroblasts')
neurons = c('Dmel_6.8h_Neurons','Dmel_10.12h_Neurons','Dmel_14.16h_Neurons')
glia = c('Dmel_6.8h_Glia','Dmel_10.12h_Glia','Dmel_14.16h_Glia')

neuroblast_data_short = loops_data[neuroblasts,short_loops]
neuroblast_data_long = loops_data[neuroblasts,long_loops]
neuron_data_short = loops_data[neurons,short_loops]
neuron_data_long = loops_data[neurons,long_loops]
glia_data_short = loops_data[glia,short_loops]
glia_data_long = loops_data[glia,long_loops]

neuroblast_short_vector = c(apply(neuroblast_data_short, 1, as.vector))
neuroblast_long_vector = c(apply(neuroblast_data_long, 1, as.vector))
neuron_short_vector = c(apply(neuron_data_short, 1, as.vector))
neuron_long_vector = c(apply(neuron_data_long, 1, as.vector))
glia_short_vector = c(apply(glia_data_short, 1, as.vector))
glia_long_vector = c(apply(glia_data_long, 1, as.vector))


col_names_short = colnames(neuroblast_data_short)
col_names_long = colnames(neuroblast_data_long)

col_names_short_1 = paste('06.08', col_names_short, sep="_")
col_names_short_2 = paste('10.12', col_names_short, sep="_")
col_names_short_3 = paste('14.16', col_names_short, sep="_")

col_names_long_1 = paste('06.08', col_names_long, sep="_")
col_names_long_2 = paste('10.12', col_names_long, sep="_")
col_names_long_3 = paste('14.16', col_names_long, sep="_")

names(neuroblast_short_vector) = c(col_names_short_1, col_names_short_2, col_names_short_3)
names(neuroblast_long_vector) = c(col_names_long_1, col_names_long_2, col_names_long_3)
names(neuron_short_vector) = c(col_names_short_1, col_names_short_2, col_names_short_3)
names(neuron_long_vector) = c(col_names_long_1, col_names_long_2, col_names_long_3)
names(glia_short_vector) = c(col_names_short_1, col_names_short_2, col_names_short_3)
names(glia_long_vector) = c(col_names_long_1, col_names_long_2, col_names_long_3)


plots = list()
tissue_names = c('neuroblast short loops', 'neuroblast long loops', 'neuron short loops', 'neuron long loops', 'glia short loops', 'glia long loops')
for (i in 1:6){
    loops_data_row = list(neuroblast_short_vector, neuroblast_long_vector, neuron_short_vector, neuron_long_vector, glia_short_vector, glia_long_vector)[[i]]
    features_to_use = rbind(features_to_use_dict[['06.08']], features_to_use_dict[['10.12']], features_to_use_dict[['14.16']])
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
}
p=grid.arrange(grobs = plots, ncol = 2)

ggsave(filename = "results/calderon/motifs/mlr3_models_benchmark_omitted_other_species_NN_9windows_different_time_features4_per_tissue_distance.pdf", plot = p, width = 10, height = 20, dpi = 300)



#one model
#now let's create one common model for all tissues
loops_data_row = c(neuroblast_vector, neuron_vector, glia_vector)
features_to_use = rbind(features_to_use_dict[['06.08']], features_to_use_dict[['10.12']], features_to_use_dict[['14.16']],
                        features_to_use_dict[['06.08']], features_to_use_dict[['10.12']], features_to_use_dict[['14.16']],
                        features_to_use_dict[['06.08']], features_to_use_dict[['10.12']], features_to_use_dict[['14.16']])
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
ggtitle(paste("All tissues", "\nAUC ROC score ", auc_scores)) +
theme(plot.title = element_text(hjust = 0.5, size = 14))

# Save the plot to a file
ggsave(filename = "results/calderon/motifs/mlr3_models_benchmark_omitted_other_species_NN_9windows_different_time_features4_all_tissues.pdf", plot = p1, width = 10, height = 20, dpi = 300)




#boruta discrete
# run boruta analysis old - all features
library(Boruta)
set.seed(111)
boruta_results = list()
features = t(features)
for (i in 1:dim(loops_data)[1]){
    print(i)
    loops_data_row = loops_data[i,]
    loops_data_row = loops_data_row[!is.na(loops_data_row)]
    features_to_use = features[,intersect(colnames(features),names(loops_data_row))]
    loops_data_row = loops_data_row[intersect(names(loops_data_row), colnames(features_to_use))]
    features_to_use = features_to_use[,match(names(loops_data_row), colnames(features_to_use))]
    features_to_use = t(features_to_use)
    res = Boruta(features_to_use, loops_data_row, doTrace = 1, maxRuns = 200) 
    boruta_results[[i]] = res$finalDecision
}

data_matrix <- do.call(rbind, lapply(boruta_results, function(x) as.numeric(x)))

data_df <- as.data.frame(data_matrix)
rownames(data_df) <- rownames(loops_data)
colnames(data_df) <- names(boruta_results[[1]])


boruta_results_f1 = data_df[,grepl("_00.02$", colnames(data_df)) | grepl("dist", colnames(data_df))]
boruta_results_f2 = data_df[,grepl("_02.04$", colnames(data_df)) | grepl("dist", colnames(data_df))]
boruta_results_f3 = data_df[,grepl("_04.06$", colnames(data_df)) | grepl("dist", colnames(data_df))]
boruta_results_f4 = data_df[,grepl("_06.08$", colnames(data_df)) | grepl("dist", colnames(data_df))]
boruta_results_f5 = data_df[,grepl("_08.10$", colnames(data_df)) | grepl("dist", colnames(data_df))]
boruta_results_f6 = data_df[,grepl("_10.12$", colnames(data_df)) | grepl("dist", colnames(data_df))]
boruta_results_f7 = data_df[,grepl("_12.14$", colnames(data_df)) | grepl("dist", colnames(data_df))]
boruta_results_f8 = data_df[,grepl("_14.16$", colnames(data_df)) | grepl("dist", colnames(data_df))]
boruta_results_f9 = data_df[,grepl("_16.18$", colnames(data_df)) | grepl("dist", colnames(data_df))]
boruta_results_f10 = data_df[,grepl("_02.04_diff$", colnames(data_df)) | grepl("dist", colnames(data_df))]
boruta_results_f11 = data_df[,grepl("_04.06_diff$", colnames(data_df)) | grepl("dist", colnames(data_df))]
boruta_results_f12 = data_df[,grepl("_06.08_diff$", colnames(data_df)) | grepl("dist", colnames(data_df))]
boruta_results_f13 = data_df[,grepl("_08.10_diff$", colnames(data_df)) | grepl("dist", colnames(data_df))]
boruta_results_f14 = data_df[,grepl("_10.12_diff$", colnames(data_df)) | grepl("dist", colnames(data_df))]
boruta_results_f15 = data_df[,grepl("_12.14_diff$", colnames(data_df)) | grepl("dist", colnames(data_df))]
boruta_results_f16 = data_df[,grepl("_14.16_diff$", colnames(data_df)) | grepl("dist", colnames(data_df))]
boruta_results_f17 = data_df[,grepl("_16.18_diff$", colnames(data_df)) | grepl("dist", colnames(data_df))]



my_palette <- colorRampPalette(c("#6969ff", "#1ebb1e", "#e4dfdf"))(n = 3)
# Open a PDF device
pdf("results/calderon/motifs/boruta_heatmap_omitted_other_species_NN_9windows.pdf", width = 50, height = 8)  # Increase the height to fit all heatmaps
library(gplots)
# Loop to generate 7 heatmaps
for(i in 1:17) {
  df <- get(paste0("boruta_results_f", i))
  df = df[,colnames(df)[order(colnames(df))]]

  # Generate the heatmap
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
            labRow = rownames(loops_data),
            labCol = colnames(df),
            cexRow = 0.8,
            cexCol = 0.8,
            margins = c(15, 15))
  
  breaks <- seq(min(df), max(df), length.out = 5)
  
  labels <- c("Tentative", "Confirmed", "Rejected")
  
  color_key <- list(at = breaks, labels = labels, space = "bottom")
  
  legend("left", legend = color_key$labels, fill = my_palette)
}

# Close the PDF device
dev.off()



# tries with shapley values
library(iml)
library(shapr)
plots = list()
tissue_names = c('neuroblast', 'neuron', 'glia')

for (i in 1:3) {
    loops_data_row = list(neuroblast_vector, neuron_vector, glia_vector)[[i]]
    features_to_use = rbind(features_to_use_dict[['06.08']], features_to_use_dict[['10.12']], features_to_use_dict[['14.16']])
    loops_data_row = loops_data_row[!is.na(loops_data_row)]
    loops_data_row = loops_data_row[intersect(names(loops_data_row), rownames(features_to_use))]
    features_to_use = features_to_use[match(names(loops_data_row), rownames(features_to_use)),]
    features_to_use = features_to_use[match(names(loops_data_row), rownames(features_to_use)),]
    features_to_use = as.data.frame(features_to_use)
    loops_data_row = as.factor(loops_data_row)
    features_to_use$target <- loops_data_row

    task <- TaskClassif$new(id = "task", backend = features_to_use, target = "target")

    learner <- lrn("classif.ranger", predict_type = "prob")
    learner$train(task)
    model=learner$model

    #explainer <- shapr(features_to_use, model, n_combinations=10000)  #nie działa, za dużo cech



    # Compute Shapley values for all samples
    predictor_ranger <- Predictor$new(learners[[1]], data = features_to_use, y = features_to_use$target)
    shapley_values_list <- lapply(1:2, function(row) {
        shapley_ranger <- Shapley$new(predictor_ranger, x.interest = features_to_use[row, ], sample.size=1)
        shapley_ranger$results
    })

    # Combine Shapley values for all samples
    shapley_values <- do.call(rbind, shapley_values_list)
    shapley_summary <- aggregate(shapley_values$value, by = list(feature = shapley_values$feature), FUN = mean)
    colnames(shapley_summary) <- c("feature", "mean_shapley_value")

    # Plot the mean Shapley values
    ggplot(shapley_summary, aes(x = reorder(feature, mean_shapley_value), y = mean_shapley_value)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        ggtitle(paste("Mean Shapley values for", tissue_names[i])) +
        xlab("Feature") +
        ylab("Mean Shapley value") +
        theme_minimal() -> p2

    plots[[length(plots) + 1]] = p2
}
