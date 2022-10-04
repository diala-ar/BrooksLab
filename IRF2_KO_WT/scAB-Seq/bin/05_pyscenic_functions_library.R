# Identify differentially active regulons
da_regulons = function(data, min_clust_size, clust_1, clust_2) {
  clust_regulons = NULL
  row.names(data) = data$Cell
  cells_clust = data[, 'cluster', drop=F]
  data        = t(data[, -1:-3]) 
  cells_1     = row.names(cells_clust)[cells_clust[, 1]==clust_1]
  cells_2     = row.names(cells_clust)[cells_clust[, 1]==clust_2]
  if (length(cells_1)<min_clust_size | length(cells_2)<min_clust_size) next
  mn1        = rowMeans(data[, cells_1])
  mn2        = rowMeans(data[, cells_2])
  pval       = sapply(1:nrow(data), function(j) wilcox.test(data[j, cells_1], data[j, cells_2])$p.value)
  e          = min(setdiff(c(mn1, mn2), 0))/10#  quantile(c(mn1, mn2), .05)
  log2_avgFC = log2(mn1 + e) - log2(mn2+e)
  da_regulons = data.frame(regulon=row.names(data), cluster_1=clust_1, cluster_2=clust_2, 
                          pval=pval, mn1=mn1, mn2=mn2, FC=mn1/mn2, log2_avgFC=log2_avgFC)
  row.names(da_regulons) = NULL
  da_regulons$adj_pval = p.adjust(da_regulons$pval, method='fdr')

  da_regulons
}  

# generate regulons activation score dotplots
auc_dotplot = function(data, clust_1, clust_2, cond, FC_cutoff=0.25) {
  library(scales)
  data = data[data$adj_pval<0.05,]
  data$regulon = paste0(data$regulon, ' (', data$size, 'g)')
  if (FC_cutoff>0)   data = data[abs(data$log2_avgFC) > .25,]
  data$log_adjp = -log10(data$adj_pval)
  
  data = data[, c('regulon', 'size', 'log2_avgFC', 'log_adjp')]
  data = data[order(data$log2_avgFC), ]
  data$regulon = factor(data$regulon, levels=data$regulon)
  
  if (all(data$log2_avgFC>0)) {
    values  = c(0, mean(data$log2_avgFC), max(data$log2_avgFC))
    colours = c('#FFF0EC', '#FFA594', 'firebrick2')
  } else if (all(data$log2_avgFC<0)) {
    values  = c(0, mean(data$log2_avgFC), max(data$log2_avgFC))
    colours = c('#E5EEF7', '#BCB3D9', 'steelblue3')
  } else {
    values  = rescale(c(min(data$log2_avgFC), 0, max(data$log2_avgFC)))
    colours = c('steelblue3', 'white', 'firebrick2')
  }

  x_lab = paste0('log2(', clust_1, '_RAS/', clust_2, '_RAS)')
  ggplot(data, aes(log2_avgFC, regulon)) +
    geom_point(aes(size=log_adjp, color=log2_avgFC)) +
    scale_color_gradientn(colours=colours, values=rescale(values)) + 
    xlab(x_lab) + ylab('') + labs(color='log2(FC)', size='-log10(FDR))') +
    theme_bw() + ggtitle(cond)  
}


