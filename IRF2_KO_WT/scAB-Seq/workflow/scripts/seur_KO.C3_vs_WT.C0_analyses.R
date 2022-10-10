library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(pheatmap)
library(msigdbr)
library(fgsea)

source(snakemake@input[['sc_heatmap_function']])


seur       = readRDS(snakemake@input[['seur']])
cluster_by = snakemake@params[['cluster_by']]
clust_1    = snakemake@params[['clust_1']]
clust_2    = snakemake@params[['clust_2']]

meta_data  = seur@meta.data
meta_data$sample.cluster = paste(meta_data$sample, meta_data[, cluster_by], sep='.')
seur@meta.data = meta_data
rm(meta_data)
seur = subset(seur, sample.cluster %in% c(clust_1, clust_2))
seur@meta.data$sample.cluster = factor(seur@meta.data$sample.cluster, levels=c('Tum_WT.0', 'Tum_KO.3'))


Idents(seur) = 'sample.cluster'
degs = FindMarkers(seur, ident.1=clust_1, ident.2=clust_2, logfc.threshold=0, min.pct=0, 
                   min.diff.pct=0)
degs = data.frame(gene=row.names(degs), degs)
degs$p_val_adj = p.adjust(degs$p_val, method='fdr')
write.csv(degs, file=snakemake@output[['degs_without_cutoff']], row.names=F)
filt_degs = degs[degs$p_val_adj <= 0.05 & abs(degs$avg_log2FC) >= 0.25 & (degs$pct.1 >= 0.1 | degs$pct.2 >= 0.1), ]
dim(filt_degs)
write.csv(filt_degs, file=snakemake@output[['degs_with_cutoff']], row.names=F)



# generate sc_heatmap of DEGs for KO.C3_vs_WT.C0
# sub sample seurat to smallest size
Idents(seur) = 'sample.cluster'
min_sample_size = min(table(seur$sample.cluster))
sub_sampled = subset(x=seur, downsample=min_sample_size) 
sub_cells   = Cells(sub_sampled)
sub_sampled
filt_degs = filt_degs[order(filt_degs$avg_log2FC), ]
norm_expr = as.matrix(sub_sampled[['RNA']]@data)
norm_expr = norm_expr[filt_degs$gene, ]
p = sc_heatmap(norm_expr, sub_sampled@meta.data, cluster_by='sample', fontsize=6.5, use_yellow_purple=T,
               clust_rows=F, max_cutoff=3, min_cutoff=-3)
ggsave(snakemake@output[['htmp_DEGs']], p, height=8)


my_genes_up = c('Tox', 'Crbn', 'Klre1', 'Eif4b', 'Eif3e', 'Nt5e', 'Irf8', 'Nr4a2', 'Eef1a1', 'Cd3e', 'Il2rb', 
                'Cd3g', 'Rpl31', 'Rps12')
my_genes_dn = c('Batf', 'Irf4', 'Prdm1', 'Cd274', 'Srgn', 'Stat5a', 'Nfkbiz', 'Nfkbid', 'Il12rb2', 'Stx11', 
                'Gzmk', 'Gbp7', 'Gbp5', 'Ttc39c', 'Traf1', 'Ccl5', 'Tnfrsf4', 'Il2ra', 'Ifitm1', 'Ifit1',
                'Gzma', 'Gzmb', 'Ccl3')
my_genes  = c(my_genes_up, rev(my_genes_dn))
norm_expr = as.matrix(sub_sampled[['RNA']]@data)
my_degs   = norm_expr[my_genes, ]
p = sc_heatmap(my_degs, sub_sampled@meta.data, cluster_by='sample', fontsize=8, use_yellow_purple=T, 
               clust_rows=F, max_cutoff=3, min_cutoff=-3)
ggsave(snakemake@output[['htmp_DEGs_of_interest']], p, height=4.2, width=5)