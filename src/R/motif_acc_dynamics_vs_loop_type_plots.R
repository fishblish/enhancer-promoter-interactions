library(ggplot2)
library(cowplot)
library(Matrix)

loops = read.table(paste0(getwd(),'/data/muszka/long_and_short_range_loops_D_mel.tsv'), header = TRUE, sep = '\t')
loops=na.omit(loops)
annot=readRDS('/home/jbartczak/enhancer-promoter-interactions/data/muszka/calderon_data/atac_meta.rds')
motifs_acc_time=annot[,c('NNv1_age','cell')]
str(motifs_acc_time)

dict_motifs=read.csv('/home/jbartczak/enhancer-promoter-interactions/data/muszka/calderon_data/motifs/dict_motifs.csv')
temp=dict_motifs$names
dict_motifs=dict_motifs$values
names(dict_motifs)=temp

heatmap_data=readRDS(paste0(getwd(),'/results/calderon/motifs/without_dist/heatmap_data_omitted_other_species_NN_features4.rds'))
important_motifs=colnames(heatmap_data)

plot_combined_list=list()
for (selected_motif in important_motifs){
  print(selected_motif)
  plot_list=list()
  for (tissue in c('Neurons','Neuroblasts','Glia')){
    print(tissue)
    boruta_values=heatmap_data[grepl(tissue,rownames(heatmap_data)),selected_motif]
    if (any(boruta_values>0)){
      importance_flag='***'
    } else {
      importance_flag=''
    }
    loops_types = list()

    cols = colnames(loops)[grepl(tissue, colnames(loops))]
    loops_types[['000']] = loops[loops[cols[1]] == 0 & loops[cols[2]] == 0 & loops[cols[3]] == 0, 'loop_id']
    loops_types[['001']] = loops[loops[cols[1]] == 0 & loops[cols[2]] == 0 & loops[cols[3]] == 1, 'loop_id']
    loops_types[['010']] = loops[loops[cols[1]] == 0 & loops[cols[2]] == 1 & loops[cols[3]] == 0, 'loop_id']
    loops_types[['100']] = loops[loops[cols[1]] == 1 & loops[cols[2]] == 0 & loops[cols[3]] == 0, 'loop_id']
    loops_types[['011']] = loops[loops[cols[1]] == 0 & loops[cols[2]] == 1 & loops[cols[3]] == 1, 'loop_id']
    loops_types[['101']] = loops[loops[cols[1]] == 1 & loops[cols[2]] == 0 & loops[cols[3]] == 1, 'loop_id']
    loops_types[['110']] = loops[loops[cols[1]] == 1 & loops[cols[2]] == 1 & loops[cols[3]] == 0, 'loop_id']
    loops_types[['111']] = loops[loops[cols[1]] == 1 & loops[cols[2]] == 1 & loops[cols[3]] == 1, 'loop_id']

    time_window_size=0.2
    time_window_start = 0 #adjusted start and end to avoid windows with small number of cells
    time_window_end = time_window_size+time_window_start
    time_windows=data.frame()
    while (time_window_end <= 18) { 
      time_windows=rbind(time_windows, data.frame(start=time_window_start, end=time_window_end))
      time_window_start <- time_window_end
      time_window_end <- round(time_window_end + time_window_size,1)
    }

    i=0
    for (hrs in c('00-02', '02-04', '04-06', '06-08', '08-10', '10-12', '12-14', '14-16', '16-18')) {
      print(hrs)
      open_cells_loops = readRDS(paste0('/home/jbartczak/enhancer-promoter-interactions/results/calderon/new_time/NN/open_cells_all_for_loops_', hrs, '_new_timeNN_faster.rds'))
      motifs=readRDS(paste0(getwd(),'/data/muszka/calderon_data/motifs/new_time/NN/GSE190130_', hrs, '_new_timeNN_matrix_motifs.rds'))
      
      split_names=strsplit(rownames(motifs),split='-')
      split_names=sapply(split_names,function(x) x[-length(x)])
      translated_names=sapply(split_names,function(x) {
        if (x[1] %in% names(dict_motifs)) {
          x[1]=dict_motifs[x[1]]
        }
        return(x[1])
      })
      rownames(motifs)=translated_names
      motifs=motifs[,!grepl('other_species',colnames(motifs))]
      rownames(motifs) = gsub('[()\\-]', '.', rownames(motifs))

      one_motif=motifs[selected_motif,]
      one_motif_df=data.frame(cell=names(one_motif),score=one_motif,stringsAsFactors=FALSE)
      motif_score=merge(one_motif_df,motifs_acc_time,by='cell')


      for (j in 1:(2/time_window_size)){
        time_window=time_windows[i+j,c('start','end')]
        motif_score_window=motif_score[motif_score$NNv1_age>=time_window$start & motif_score$NNv1_age<time_window$end,]
        for (loop in loops$loop_id){
          open_cells = open_cells_loops[[loop]]
          open_mean_score=mean(motif_score_window[motif_score_window$cell %in% open_cells,]$score, na.rm=TRUE)
          rest_score=mean(motif_score_window[!motif_score_window$cell %in% open_cells,]$score, na.rm=TRUE)
          mean_diff=open_mean_score-rest_score
          time_windows[i+j,loop]=mean_diff
          

        }
      }
      i=i+(2/time_window_size)
    }
    time_windows=time_windows[time_windows$start>=0.4,]
    data_to_plot=data.frame(time_window = paste0(time_windows$start,'-',time_windows$end))
  
    time_stamp_means <- rowMeans(time_windows[, !(names(time_windows) %in% c("start", "end"))], na.rm = TRUE)

    for (group in names(loops_types)){
      data_to_plot[group]=rowMeans(time_windows[,loops_types[[group]]], na.rm=TRUE)-time_stamp_means
    } 

    #Create line plot
    data_to_plot_long <- tidyr::pivot_longer(data_to_plot, -time_window, names_to = "group", values_to = "value")
    data_to_plot_long$time_window <- factor(data_to_plot_long$time_window, levels = data_to_plot$time_window)

    custom_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
                   "#FF7F00", "#00CED1", "#A65628", "#F781BF")

    
    p <- ggplot(data_to_plot_long, aes(x = time_window, y = value, color = group, group = group)) +
      geom_line() +  
      theme_bw() +
      labs(title = paste0("Motif ",selected_motif,", accessibility score differences by loop type for ", tissue,'\n',importance_flag),
          x = "Time Window",
          y = "Mean Difference",
          color = "Group") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))+
      scale_x_discrete(breaks = data_to_plot$time_window[seq(1, nrow(data_to_plot), by = 5)])+
      scale_color_manual(values = custom_colors)+
      geom_rect(data = data.frame(xmin = c(6/time_window_size+1,10/time_window_size+1,14/time_window_size+1), xmax = c(8/time_window_size+1,12/time_window_size+1,16/time_window_size+1), ymin = -Inf, ymax = Inf), 
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
                alpha = 0.1, inherit.aes = FALSE)

    plot_list=c(plot_list,list(p))
  

    # generate transparent line for each loop
    # data_to_plot_transparent=time_windows
    # data_to_plot_transparent$time_window=paste0(time_windows$start,'-',time_windows$end)
    # data_to_plot_transparent=data_to_plot_transparent[,!(names(data_to_plot_transparent) %in% c('start', 'end'))]
    # long_transparent_data=tidyr::pivot_longer(data_to_plot_transparent,-c(time_window),names_to='loop', values_to='value')
    # long_transparent_data$time_window <- factor(long_transparent_data$time_window, levels = data_to_plot$time_window)

    # for (group in names(loops_types)){
    #   long_transparent_data[long_transparent_data$loop %in% loops_types[[group]],'group']=group
    # }
    # p <- ggplot(long_transparent_data, aes(x = time_window, y = value, color = group, group = loop)) +
    #   geom_line(alpha=0.1) +  
    #   theme_bw() +
    #   labs(title = paste0("Motif ",selected_motif,", accessibility score fold changes by loop type for ", tissue),
    #       x = "Time Window",
    #       y = "Mean Difference",
    #       color = "Group") +
    #   theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    #   scale_x_discrete(breaks = data_to_plot$time_window[seq(1, nrow(data_to_plot), by = 5)])+
    #   scale_color_manual(values = custom_colors)+
    #   guides(color = guide_legend(override.aes = list(alpha = 1)))

    # plot_list=c(plot_list,list(p))

  }
  
  plot_combined = plot_grid(plotlist = plot_list, ncol = 3)
  plot_combined_list=c(plot_combined_list,list(plot_combined))
}
pdf('/home/jbartczak/enhancer-promoter-interactions/results/calderon/motifs/motif_accessibility_dynamics/all_important_motif_acc_dynamics_vs_loop_type.pdf', width = 24, height = 7)
for (i in 1:length(plot_combined_list)){
  print(plot_combined_list[[i]])
}
dev.off()

#saveRDS(plot_list, '/home/jbartczak/enhancer-promoter-interactions/results/calderon/motifs/motif_accessibility_dynamics/motif_acc_dynamics_vs_loop_type.rds')
#plot_combined = plot_grid(plotlist = plot_list, ncol = 3)
#ggsave(paste0('/home/jbartczak/enhancer-promoter-interactions/results/calderon/motifs/motif_accessibility_dynamics/motif_acc_dynamics_vs_loop_type_time_window_',time_window_size,'_standarized.png'), plot_combined, width = 24, height = 20)
#ggsave(paste0('/home/jbartczak/enhancer-promoter-interactions/results/calderon/motifs/motif_accessibility_dynamics/',selected_motif,'_motif_acc_dynamics_vs_loop_type_',tissue,'_time_window_',time_window_size,'.png'), p, width = 20, height = 10)

#plot histogram: number of cells in time window

# annot=readRDS('/home/jbartczak/enhancer-promoter-interactions/data/muszka/calderon_data/atac_meta.rds')
# cell_age=annot[,'NNv1_age']

# png(filename = '/home/jbartczak/enhancer-promoter-interactions/results/calderon/motifs/motif_accessibility_dynamics/cell_age_histogram.png', width = 800, height = 600)
# hist(cell_age, breaks = 100, main = 'Number of cells in time window', xlab = 'Time window', xlim = c(0, 20), xaxt = 'n')
# x_ticks <- seq(0, 20, by = 1) 
# axis(1, at = x_ticks, labels = x_ticks)
# dev.off()

# hist_data = hist(cell_age, breaks = 100, main = 'Number of cells in time window', xlab = 'Time window', xlim = c(0, 20), xaxt = 'n', plot = TRUE)
# print(paste0('Start points of time windows with number of cells lower then 500: ', list(hist_data$breaks[which(hist_data$counts < 500)])))
