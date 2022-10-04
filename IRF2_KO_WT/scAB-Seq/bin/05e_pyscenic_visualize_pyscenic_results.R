setwd('~/projects/sabelo/AB_seq_v2/')
source('bin/05_pyscenic_functions_library.R')

library(ggplot2)
library(Seurat)
library(SCENIC)
library(dplyr)
library(grid)


conds = c('Tum_CD8', 'Tum_CD8_KO.C3_WT.C0')

# # print AUC scores with their threshold and update the thresholds in the excel file on need
# for (cond in conds) {
#   reg_thresh = read.csv(paste0('results/Seurat_AB_seq/Tum_CD8_v2/res_0.5/pyscenic_100_runs/aucs_thresholds_',cond,'.csv'))
#   names(reg_thresh) = c('reg', 'thresh')
#   auc_score = read.csv(paste0('results/Seurat_AB_seq/Tum_CD8_v2/res_0.5/pyscenic_100_runs/aucs_', cond, '.csv'), row.names=NULL)
#   
#   pdf(paste0('results/Seurat_AB_seq/Tum_CD8_v2/res_0.5/pyscenic_100_runs/aucs_histogram_',cond,'.pdf'))
#   for (reg in reg_thresh$reg) {
#     hist_data = data.frame(auc=auc_score[, reg])
#     p = ggplot(hist_data, aes(auc)) +
#       geom_histogram(fill='red') + 
#       geom_vline(xintercept= reg_thresh[reg_thresh$reg==reg, 'thresh']) +
#       theme_bw() + ggtitle(reg)
#     print(p)
#   }
#   dev.off()
# }

for (i in 1:length(conds)) {
  cond = conds[i]
  seur = readRDS(paste0('results/Seurat_AB_seq/Tum_CD8_v3/Seurat/Seur_', cond, '_clustered.rds'))
  if (cond=='Tum_CD8') seur$seurat_clusters = seur$sample
  
  min_sample_size = min(table(seur$sample))
  Idents(seur) = 'seurat_clusters'
  sub_sampled = subset(x=seur, downsample=min_sample_size)
  sub_cells   = Cells(sub_sampled)
  sub_sampled
  
  meta_data = data.frame(Cell=row.names(seur@meta.data),
                         sample=seur@meta.data$sample, cluster=seur@meta.data$seurat_clusters)
  
  # get regulon targets and size
  reg      = readLines(paste0('results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/pyscenic_100_runs/high_confidence_reg_', cond,'.gmt'))
  reg_list = lapply(1:length(reg), function(i) setdiff(strsplit(reg[[i]], '\t')[[1]], 'na')[-1])
  TFs      = sapply(1:length(reg), function(i) strsplit(reg[[i]], '\t')[[1]][1])
  names(reg_list) = TFs
  reg_size = sapply(reg_list, length)
  
  clust_1 = ifelse(cond=='Tum_CD8', 'Tum_KO', 'Tum_KO.C3')
  clust_2 = ifelse(cond=='Tum_CD8', 'Tum_WT', 'Tum_WT.C0')
  cells_1 = meta_data[meta_data$cluster==clust_1, 'Cell'] 
  cells_2 = meta_data[meta_data$cluster==clust_2, 'Cell']
  
  if (cond=='Tum_CD8_KO.C3_WT.C0') {
    degs   = read.csv('results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/Tum_CD8_DEGs_C3.KO_vs_C0_WT_without_cutoff.csv')
  } else {
    degs = read.csv('results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/Tum_CD8_DEGs_without_cutoff.csv')
    degs = degs[degs$cluster_1=='Tum_KO.All', c('symbol', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj')]
    names(degs)[1] = 'gene'
  }
  degs$rank = sign(degs$avg_log2FC)*(-log10(degs$p_val))
  degs = degs[order(degs$rank, decreasing=T), ]
  leading_edge = c(degs[1:round(nrow(degs)*.05), 'gene'], degs[nrow(degs):(nrow(degs)-round(nrow(degs)*.05)), 'gene'])

  # read AUC scores
  auc_score = read.csv(paste0('results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/pyscenic_100_runs/aucs_', cond, '.csv'), row.names=NULL)
  row.names(auc_score) = auc_score$Cell
  auc_score = t(auc_score[, -1])
  colnames(auc_score) = sub('.', '-', colnames(auc_score), fixed=T) 
  mn1        = round(rowMeans(auc_score[, cells_1]),3)
  mn2        = round(rowMeans(auc_score[, cells_2]),3)
  pval       = sapply(1:nrow(auc_score), function(j) wilcox.test(auc_score[j, cells_1], auc_score[j, cells_2])$p.value)
  e          = min(setdiff(c(mn1, mn2), 0))/10
  log2_avgFC = round(log2(mn1 + e) - log2(mn2+e),2)
  
  # read binary AUC scores
  bin_auc = read.csv(paste0('results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/pyscenic_100_runs/binary_aucs_', cond,'.csv'))
  row.names(bin_auc) = bin_auc$Cell
  bin_auc = t(bin_auc[, -1])
  colnames(bin_auc) = sub('.', '-', colnames(bin_auc), fixed=T)
  active_prct_1 = round(rowSums(bin_auc[, cells_1])/length(cells_1) * 100, 2)
  active_prct_2 = round(rowSums(bin_auc[, cells_2])/length(cells_2) * 100, 2)
  

  regulons   = data.frame(regulon=names(reg_list), size=reg_size, cond=cond, clust_1=clust_1, clust_2=clust_2,
                          active_prct_1=active_prct_1, active_prct_2=active_prct_2, 
                          mn1=mn1, mn2=mn2, mn_FC=round(mn1/mn2,2),
                          pval=pval, log2_avgFC=log2_avgFC)
  regulons$adj_pval = p.adjust(regulons$pval, method='fdr')
  selected_reg = regulons[regulons$size>=5 & abs(regulons$log2_avgFC)>.25 & regulons$adj_pval<0.05 &
                          (regulons$active_prct_1>10 | regulons$active_prct_2>10), ]
  selected_reg = selected_reg[order(-sign(selected_reg$log2_avgFC)*(-log10(selected_reg$pval))), ]

  write.csv(selected_reg, 
            paste0('results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/pyscenic_100_runs/differentially_activated_reg_', cond,'.csv'),
            row.names=F)
  
  # print the dotplot
  ras_dtp = auc_dotplot(selected_reg, clust_1, clust_2, cond) +
    theme(axis.text.y = element_text(size=11))
  ggsave(paste0('results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/pyscenic_100_runs/differentially_activated_reg_', cond,'_RAS_FC_dotplot.pdf'), 
         plot=ras_dtp, width=4, height=6)
  
  meta_data$cluster = factor(meta_data$cluster, levels=c(clust_2, clust_1))
  row.names(meta_data) = meta_data$Cell
  meta_data = meta_data[sub_cells,]
  htp_data = as.matrix(auc_score[selected_reg$regulon, sub_cells])
  ras_htmp = sc_heatmap(htp_data, meta_data, cluster_by='cluster', title=paste0('Regulon RAS score in ', cond), 
                        clust_rows=F, use_yellow_purple=T, min_cutoff=-4, max_cutoff=4)
  ggsave(paste0('results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/pyscenic_100_runs/differentially_activated_reg_',cond,'_RAS_heatmaps.pdf'), ras_htmp)
  dev.off()
}

