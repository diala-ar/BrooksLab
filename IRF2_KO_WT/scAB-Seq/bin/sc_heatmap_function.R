sc_heatmap = function(htp_data, meta_data, cluster_by, use_yellow_purple=F, 
                      min_cutoff=NULL, max_cutoff=NULL, fontsize=10, scale_data=T, clust_rows=T) {
  source('~/projects/bin/functions/custom_pheatmap_colors.R')
  
  meta_data$cluster_by = meta_data[, cluster_by]
  
  # order meta_data by grp to order then the htp_data columns
  meta_data  = meta_data[order(meta_data$cluster_by), ] 
  htp_data   = htp_data[, row.names(meta_data)]
  if (scale_data==T) htp_data = t(scale(t(htp_data)))
  htp_data   = na.omit(htp_data)
  if (!is.null(min_cutoff))
    htp_data[htp_data < min_cutoff] = min_cutoff
  if (!is.null(max_cutoff)) 
    htp_data[htp_data > max_cutoff] = max_cutoff
  
  color_n_breaks = custom_pheatmap_colors(min(htp_data), max(htp_data), use_yellow_purple)
  gaps_nb = as.array(table(meta_data[cluster_by[1]]))
  gaps_cols_index = sapply(1:length(gaps_nb), function(i) sum(gaps_nb[1:i]))
  p = pheatmap::pheatmap(as.matrix(htp_data), 
                         #color = PurpleAndYellow(20), #RColorBrewer::brewer.pal(9, "Blues"), #colorRampPalette(c("black", "yellow", "blue"))(50),
                         color  = color_n_breaks$my_colors,
                         breaks = unique(color_n_breaks$breaks),
                         border_color = NA,
                         cluster_rows = clust_rows, 
                         cluster_cols = F,
                         show_rownames = T,
                         show_colnames = F,
                         scale='none',
                         annotation_col = select(meta_data, all_of(rev(cluster_by))),
                         gaps_col = gaps_cols_index,
                         fontsize = fontsize)
  p
}
