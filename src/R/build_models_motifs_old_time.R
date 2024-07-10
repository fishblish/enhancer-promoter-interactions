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
motifs_path = paste0(getwd(),'/data/muszka/calderon_data/motifs/')
calderon_path = paste0(getwd(),'/data/muszka/calderon_data/')
hrs_list = c('exp1_hrs06-10_b1','exp1_hrs10-14_b1','exp1_hrs14-18_b1')
for (hrs in hrs_list){
    print(paste0('Computing fold change for ', hrs))
    loops_cells = readRDS(paste0("results/calderon/motifs/open_cells_all_for_loops_",hrs,".rds"))
    file_name = paste0('GSE190130_', hrs, '')
    motifs <- readRDS(paste0(motifs_path, hrs,'_motif_activity.rds'))
    motifs <- motifs[, colSums(is.na(motifs)) == 0]
    file_name = paste0('GSE190130_', hrs, '.peak_matrix')
    #cell_names = readLines(paste0(calderon_path,file_name, '.columns.txt.gz'))
    cell_names = read.table(gzfile(paste0(calderon_path,file_name, '.columns.txt.gz')), sep='\t', header=FALSE)
    cell_names = cell_names$V1

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
saveRDS(fold_change_dict, paste0(getwd(),'/results/calderon/motifs/fold_change_dict_3.rds'))
#saveRDS(fold_change_dict, paste0(getwd(),'/results/calderon/motifs/new_time/fold_change_dict.rds'))
fold_change_dict = readRDS(paste0(getwd(),'/results/calderon/motifs/fold_change_dict_3.rds'))
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

#motif names
names = readRDS(paste0(getwd(),'/data/muszka/calderon_data/motif_names.rds'))
good_motifs <- !grepl('.*\\(.*\\).*.*\\(.*\\)', names)
names[!good_motifs] <- 'other_species'

# a bunch of missing genes - E2f, CG31782 are ambiguous
motifs_swaps_this <- paste0('^', c('vfl', # 'Her', # this is already in motifs
    'FBgn0054031', 'FBpp0292044', 'CG42234',
    'FBpp0073832', 'CG7056', 'Lag1', 'FBgn0261930', 'STAT92E', 'FBgn0262975',
    'Bteb2', 'CG10267', 'CG11071', 'FBgn0263239', 'CG13424', 'CG13897', 'CG14962',
    'FBgn0263240', 'FBgn0262477', 'CG17181', 'CG31670',
    # 'CG32830', # this is ab which already in motifs
    'CG5669',
    'CG9437', 'FBgn0261705', 'FBgn0262582', 'FBgn0263118', 'FBgn0262656', 'HLHm3',
    'HLHm5', 'HLHm7', 'HLHmbeta', 'HLHmdelta', 'HLHmgamma', 'FBgn0263112',
	'Side', 'CG32830', 'CG3838', 'CG7928', 'CG6272', 'Dip3', 'dys', 'Hr46',
	'fd64A', 'Mio', 'pfk'),'$')
motifs_swaps_that <- c('zld', # 'her',
    'CG34031', 'pdm3', 'Dbx', 'acj6', 'HHEX',
    'schlank', 'vnd', 'Stat92E', 'cnc', 'Klf15', 'Zif', 'mamo', 'dar1', 'lms', 'hng3',
    'Asciz', 'Coop', 'FoxP', 'Kah', 'erm', # 'ab',
    'Spps', 'hng1', 'CG42741', 'cic',
    'tx', 'Myc', 'E(spl)m3-HLH', 'E(spl)m5-HLH', 'E(spl)m7-HLH', 'E(spl)mbeta-HLH',
    'E(spl)mdelta-HLH', 'E(spl)mgamma-HLH', 'Mitf',
	'side', 'Lmx1a', 'brwl', 'ZIPIC', 'Irbp18', 'Dlip3', 'Dys', 'Hr3', 'FoxL1',
	'mio', 'Pfk')

for (i in 1:length(motifs_swaps_this)) {
    names <- gsub(motifs_swaps_this[i], motifs_swaps_that[i], names)
}

old_names = rownames(fold_change_dict[[1]])
old_names = sapply(old_names, strsplit, split = '_')
old_names = sapply(old_names, function(x) x[1])
dict_motifs = names
names(dict_motifs) = old_names
write.csv(dict_motifs, file = paste0(calderon_path,'/motifs/dict_motifs.csv'))

# translate names from old to new in features_to_use

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
new_colnames = sub('_exp1_hrs', '_', colnames(features))

colnames(features) = new_colnames
str(features)

# omit other species motifs
features = features[,!grepl('other_species', colnames(features))]

#build models

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
colnames(features) = make.names(colnames(features))
for (i in 1:dim(loops_data)[1]){
    loops_data_row = loops_data[i,]
    loops_data_row = loops_data_row[!is.na(loops_data_row)]
    loops_data_row = loops_data_row[intersect(names(loops_data_row), rownames(features_to_use))]
    features_to_use = features[match(names(loops_data_row), rownames(features)),]
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
ggsave(filename = "results/calderon/motifs/mlr3_models_benchmark_3_omitted_other_species.pdf", plot = p, width = 10, height = 20, dpi = 300)



# run boruta analysis
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

boruta_results_f1 = data_df[,grepl("_06.10_b1$", colnames(data_df)) | grepl("dist", colnames(data_df))]
boruta_results_f2 = data_df[,grepl("_10.14_b1$", colnames(data_df)) | grepl("dist", colnames(data_df))]
boruta_results_f3 = data_df[,grepl("_14.18_b1$", colnames(data_df)) | grepl("dist", colnames(data_df))]
boruta_results_f4 = data_df[,grepl("_10.14_b1_diff$", colnames(data_df)) | grepl("dist", colnames(data_df))]
boruta_results_f5 = data_df[,grepl("_14.18_b1_diff$", colnames(data_df)) | grepl("dist", colnames(data_df))]



my_palette <- colorRampPalette(c("#6969ff", "#1ebb1e", "#e4dfdf"))(n = 3)
# Open a PDF device
pdf("results/calderon/motifs/boruta_heatmap_omitted_other_species.pdf", width = 50, height = 8)  # Increase the height to fit all heatmaps
library(gplots)
# Loop to generate 7 heatmaps
for(i in 1:5) {
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

