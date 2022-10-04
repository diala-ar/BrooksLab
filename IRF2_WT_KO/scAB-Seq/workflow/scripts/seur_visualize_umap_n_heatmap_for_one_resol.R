library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)


seur  = readRDS(file=snakemake@input[['seur']]) # seur = readRDS("results/Seurat_AB_seq/Tum_CD8_v3/Seur_Tum_CD8_clustered.rds")
resol = snakemake@params[['resolution']] # resol = 0.4


# generate UMAP split by sample (condition) for this resolution
Idents(seur) = 'sample'
min_sample_size = min(table(seur$sample))
sub_sampled = subset(x=seur, downsample=min_sample_size) 
sub_sampled

this_idents    = ifelse('ADT' %in% Assays(seur), paste0('wsnn_res.', resol), paste0('RNA_snn_res.', resol))
this_reduction = ifelse('ADT' %in% Assays(seur), 'wnn.umap', 'umap')
p1 = DimPlot(sub_sampled, reduction=this_reduction, group.by=this_idents, 
             label=T, repel=T, label.size=3.5, pt.size=2.5, shuffle=T) +
  theme(legend.text=element_text(size=7)) +
  guides(color=guide_legend(override.aes=list(size=2), keywidth=.2, keyheight=.6, ncol=1) )
p1$layers[[1]]$aes_params$alpha = 0.6
p2 = DimPlot(sub_sampled, reduction=this_reduction, split.by='sample', group.by=this_idents, 
             label=T, repel=T, label.size=3.5, pt.size=2.5, shuffle=T) + NoLegend() + ggtitle(NULL)  
p2$layers[[1]]$aes_params$alpha = 0.6
print(ggarrange(p1, p2, ncol=2, nrow=1, widths=c(1.7, 2.5)))
p = ggarrange(p1, p2, ncol=2, nrow=1, widths=c(1.7, 2.5))

p1 = DimPlot(sub_sampled, reduction='wnn.umap', 
             group.by=paste0('wsnn_res.', resol), label=T, repel=T, label.size=5, pt.size=2.5) + NoLegend()
p1$layers[[1]]$aes_params$alpha = 0.6
p2 = DimPlot(sub_sampled, reduction='wnn.umap', split.by='sample',
             group.by=paste0('wsnn_res.', resol), pt.size=2.5, label=T, repel=T, label.size=5) + NoLegend() + ggtitle(NULL)  
p2$layers[[1]]$aes_params$alpha = 0.6
p = ggarrange(p1, p2, ncol=2, nrow=1, widths=c(1.5, 2.5))

ggsave(snakemake@output[['umap_split_by_sample']], plot=p, height=3, width=7)


# write csv files of adt and rna markers of this resolution
adt_markers = read.csv(snakemake@input[['all_resol_adt_markers']])
adt_markers = adt_markers[adt_markers$resolution==resol, -1]
rna_markers = read.csv(snakemake@input[['all_resol_rna_markers']])
rna_markers = rna_markers[rna_markers$resolution==resol, -1]

write.csv(adt_markers, file=snakemake@output[['adt_markers']], row.names=F)
write.csv(rna_markers, file=snakemake@output[['rna_markers']], row.names=F)


# generate adt markers heatmap
Idents(seur) = paste0('wsnn_res.', resol)
p = DoHeatmap(seur, features=unique(adt_markers$gene), assay="ADT", angle=90, size=4) + NoLegend() + 
  ggtitle(paste0('Res: ', resol, ', ADT markers'))

ggsave(snakemake@output[['adt_heatmap']], plot=p, height=2)


# generate rna markers heatmap
top_nbs     = c(10, 20)
pdf_heights = c(6, 12)
for (j in 1:length(top_nbs)) {
  top_nb = top_nbs[j]
  Idents(seur) = ifelse('ADT' %in% Assays(seur), paste0('wsnn_res.', resol), paste0('RNA_snn_res.', resol))
  rna_markers %>%
    group_by(cluster) %>%
    top_n(n=top_nb, wt=avg_log2FC) -> top_markers
  p = DoHeatmap(seur, features=unique(top_markers$gene), assay="RNA", angle=90) +
    theme(text = element_text(size = 10)) + #NoLegend() + 
    ggtitle(paste0('Res: ', resol, ', Top ', top_nb, ' up-regulated RNA markers'))
  ggsave(filename=snakemake@output[[paste0('rna_top_', top_nb, '_heatmap')]], plot=p, height=pdf_heights[j])  
}



# visualization percentage of cells in each cluster for a specific resolution
clusters_prct = data.frame(table(seur$sample, seur$wsnn_res.0.4))
names(clusters_prct) = c('sample', 'cluster', 'cells_nb')
total = clusters_prct %>%
  group_by(sample) %>%
  summarize(total_cells=sum(cells_nb))
clusters_prct = merge(clusters_prct, total)
clusters_prct$percentage = round(clusters_prct$cells_nb / clusters_prct$total_cells * 100)

p = ggplot(clusters_prct, aes(x=sample, y=percentage, fill=cluster)) +
  geom_bar(stat='identity', width=0.7) +
  ggtitle('Res 0.4') + theme_bw()
ggsave(snakemake@output[['clust_prct_barplot']], p, height=3, width=2.5)
write.csv(clusters_prct, snakemake@output[['clust_prct_csv']], row.names=F)

# cells number per cluster
cells_nb = aggregate(cells_nb, list(as.numeric(cells_nb$cluster_id)), length)
cells_nb = cells_nb[, -3]
names(cells_nb) = c('cluster_id', 'cells_nb')
write.csv(cells_nb, "results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/cells_nb_per_cluster.csv", row.names=F)

cells_nb = data.frame(sample=sub_sampled$sample, cluster_id=sub_sampled$wsnn_res.0.4)
cells_nb = cells_nb %>%
  group_by(cluster_id, sample) %>%
  summarise(n = n())
names(cells_nb)[3] = 'cells_nb'
write.csv(cells_nb, "results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/cells_nb_per_cluster_n_sample.csv", row.names=F)
