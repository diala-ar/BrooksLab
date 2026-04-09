library(Seurat)
library(ggplot2)
library(openxlsx)
library(pheatmap)
library(ggrastr)


source('bin/00_functions.R')

res_dir  = file.path('results')
dir_seur = file.path(res_dir, 'seurat_object') 
dir_degs = file.path(res_dir, 'DEGs')
dir_gsea = file.path(res_dir, 'GSEA')
dir_gsea_clones = file.path(res_dir, 'Top_clones_GSEA')


seurs_ls = readRDS(file.path(dir_seur, '04_integrated_seurs_ls_n_clustered_CD8_T_cells.rds'))
groups   = names(seurs_ls)

genes_ls = list('D8' = c('Gzma', 'Cx3cr1', 'S1pr5', 'Zeb2', 'Klrg1', 'Tbx21', 'Klf2', 
                         'Klf3', 'Sidt1', 'Gapdh', 'Pdcd1', 'Igkc', 'Cxcr3', 'Cxcr6', 
                         'Il7r', 'Tcf7', 'Id2', 'Id3', 'Sell', 'Ccr7'),
                'D50' = c('Gzma', 'Cx3cr1', 'S1pr5', 'Zeb2', 'Klrg1', 'Tbx21', 'Klf2', 
                          'Klf3', 'Cd7', 'Irf7', 'mt-Nd5', 'mt-Co3', 'Il2rb', 'Ifng', 
                          'Ifngr1', 'Cxcr3', 'Cxcr6', 'Il7r', 'Tcf7', 'Id2', 'Id3', 
                          'Sell', 'Ccr7', 'Uba52', 'Birc5', 'Cenpa'))

## visualize specific DEGs (log2FC and FDR on the level of each cluster and bulk) 
pdf_heights = c('D8'=2.5, 'D50'=4)
for (group_i in groups) {
  seur     = seurs_ls[[group_i]]
  
  degs_file = file.path(dir_degs, paste0(group_i, '_DEGs_ON_vs_OFF.xlsx'))
  genes     = genes_ls[[group_i]]
  ident_1   = paste0('ON_', group_i)
  ident_2   = paste0('OFF_', group_i)
  
  htp_data = get_logFC_n_FDR(seur, degs_file, genes, significant_rows=T,
                             ident_1, ident_2)
  
  p = logFC_n_FDR_htmp(htp_data, fontsize=7, use_yellow_purple=F, 
                       abs_logFC_cutoff=0, cluster_rows=T, title=group_i) 
  ggsave(file.path(dir_degs, paste0(group_i, '_Heatmap_log2FC_selected_DEGs_per_cluster.pdf')),
         plot=p, width=3, height=pdf_heights[[group_i]])
}  


## visualize all significant DEGs (log2FC and FDR on the level of each cluster and bulk) 
pdf_heights = c('D8'=5, 'D50'=40)
for (group_i in groups) {
  seur     = seurs_ls[[group_i]]
  
  degs_file = file.path(dir_degs, paste0(group_i, '_DEGs_ON_vs_OFF.xlsx'))
  ident_1   = paste0('ON_', group_i)
  ident_2   = paste0('OFF_', group_i)
  genes = NULL
  
  htp_data = get_logFC_n_FDR(seur, degs_file=degs_file, significant_rows=F,
                             ident_1=ident_1, ident_2=ident_2)
  
  p = logFC_n_FDR_htmp(htp_data, fontsize=7, use_yellow_purple=F, 
                       abs_logFC_cutoff=0, cluster_rows=T, title=group_i) 
  ggsave(file.path(dir_degs, paste0(group_i, '_Heatmap_log2FC_all_DEGs_per_cluster.pdf')),
         plot=p, width=3.5, height=pdf_heights[[group_i]])
}  



# visualise mean of gene expression per cluster across groups
pdf_heights = c('D8'=3, 'D50'=4)
for(group_i in groups) {
  seur = seurs_ls[[group_i]]
  
  p = mean_htp(seur, cluster_by='seurat_clusters', features=genes_ls[[group_i]],
               fontsize=7, title=group_i)
  
  ggsave(file.path(dir_degs, paste0(group_i, '_Heatmap_mean_selected_DEGs_per_cluster.pdf')),
         plot=p, width=3.5, height=pdf_heights[[group_i]])
}



# visualise single cells heatmap of all degs in bulk
pdf_heights = c('D8'=2.5, 'D50'=12)
for(group_i in groups) {
  seur = seurs_ls[[group_i]]
  degs_file = file.path(dir_degs, paste0(group_i, '_DEGs_ON_vs_OFF.xlsx'))
  degs = read.xlsx(degs_file, sheet='Bulk')
  
  font_size = ifelse(group_i=='D8', 7, 5)
  p = sc_heatmap(seur, cluster_by='sample', use_yellow_purple=T, title=group_i,
                 min_cutoff=-2, max_cutoff=4, fontsize=font_size, clust_rows=F, 
                 features=degs$gene)
  ggsave(file.path(dir_degs, paste0(group_i, '_Heatmap_SC_Zscore_all_DEGs_per_cluster.pdf')),
         plot=p, width=4, height=pdf_heights[[group_i]])
}





#### Visualize feature plots for TCR gene ####
out_dir_res = file.path(res_dir, 'TCRs_FeaturePlots')
for (group_i in groups) {
  seurs = seurs_ls[[group_i]]
  TCRs  = sort(c(grep('Trav', row.names(seurs), value=T),
                 grep('Trbv', row.names(seurs), value=T)))
  p_ls = featureplots_faceted(seurs, split_by='sample', features=TCRs, 
                              ident_1=paste0('ON_', group_i), ident_2=paste0('OFF_', group_i), 
                              sep='_', alpha=0.5)
  pdf(file.path(out_dir_res, paste0(group_i, '_Featureplots_TCRs_with_alpha.pdf')), 
      width=6, height=2.5) 
  for (i in 1:length(TCRs)) {
    print(p_ls[[i]])
  }
  dev.off()
}



#### Visualize feature plots for other genes ####
genes = sort(c("Mki67", "Cenpa", "Pdia6", "Cx3cr1", "Klrg1",  "Il7r", "Tcf7", 
               "Sell", "Cd27", "Id3", "Cxcr3", "Cxcr6"))
out_dir_res = file.path(res_dir, 'umaps_n_markers')
pdf(file.path(out_dir_res, paste0('D8_n_D50_All_cells_featureplots_without_alpha.pdf')), 
    width=10, height=8) 
for (group_i in groups) {
  seurs = seurs_ls[[group_i]]
  p_ls = FeaturePlot(seurs, features=genes, max.cutoff='q99', combine=F, pt.size=.3)
  p_ls = lapply(p_ls, function(p) {
    p = p +
      theme(legend.key.size = unit(.3, 'cm'),
            legend.text = element_text(size=5),
            axis.text = element_text(size=7),
            title = element_text(size=8)) 
    # p$layers[[1]]$aes_params$alpha = 0.5
    p
  })
  p = ggpubr::ggarrange(plotlist=p_ls, nrow=3, ncol=4)
  p = ggpubr::annotate_figure(p, top=group_i)
  print(p)
}
dev.off()


## Visualize GSEA results as barplots
pdf_heights = c(D8=1.65, D50=5)
for (group_i in groups) {
  data = read.xlsx(file.path(dir_gsea, paste0(group_i, '_GSEA_ON_vs_OFF_Hallmark.xlsx')), sheet='Bulk')
  data = data[data$padj < 0.05, ]
  p    = gsea_barplot(data)
  ggsave(file.path(dir_gsea, paste0(group_i, '_Barplot_activated_n_inhibited_pathways_ON_vs_OFF.pdf')), 
         plot=p, width=5, height=pdf_heights[[group_i]])  
}


## Visualize GSEA results as dotplots
pdf_heights = c(D8=4, D50=6)
for (group_i in groups) {
  file_name = file.path(dir_gsea, paste0(group_i, '_GSEA_ON_vs_OFF_Hallmark.xlsx'))
  p = gsea_dotplot(file_name, title=group_i)
  
  ggsave(file.path(dir_gsea, paste0(group_i, '_Dotplot_activated_n_inhibited_pathways_ON_vs_OFF.pdf')), 
         plot=p, width=8, height=pdf_heights[[group_i]])  
}




## Visualize GSEA of top3 clones results as barplots
pdf_heights = c(D8=4.7, D50=6)
for (group_i in groups) {
  p_ls = lapply(1:3, function(sheet_i) {
    data = read.xlsx(file.path(dir_gsea_clones, paste0(group_i, '_GSEA_ON_vs_OFF_Hallmark.xlsx')), sheet=sheet_i)
    data = data[data$padj < 0.05, ]
    gsea_barplot(data) + ggtitle(paste0(group_i, ', Top ', sheet_i,' clones'))
  })
  if (group_i=='D8') {
    p = ggarrange(plotlist = p_ls[c(1,3)], ncol=1, nrow=2, heights=c(1, .7))
  } else if (group_i=='D50') {
    p = ggarrange(plotlist = p_ls[c(2,3)], ncol=1, nrow=2, heights=c(1, .75))
  }
  ggsave(file.path(dir_gsea_clones, paste0(group_i, '_Barplot_activated_n_inhibited_pathways_ON_vs_OFF.pdf')), 
         plot=p, width=6, height=pdf_heights[[group_i]])  
}
