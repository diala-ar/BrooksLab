library(Seurat)
library(annotables)


resol = snakemake@params[['resolution']]
seur  = readRDS(file=snakemake@input[['seur']])
conditions = read.csv(snakemake@params[['conds']])$conditions
org=snakemake@params[['org']]

metadata   = seur@meta.data
this_idents = ifelse('ADT' %in% Assays(seur), paste0('wsnn_res.', resol), paste0('RNA_snn_res.', resol))
metadata$sample.cluster = paste(metadata$sample, metadata[, this_idents], sep='.')
seur@meta.data = metadata


# DEGs between the 2 conditions within each cluster
degs = NULL
Idents(seur) = 'sample.cluster'
clusters = sort(as.character(unique(seur@meta.data[, this_idents])))
for (i in 1:length(clusters)) {
  print(i)
  temp      = NULL
  cluster   = as.numeric(clusters[i])
  cluster_1 = paste0(conditions[2], '.', cluster)
  cluster_2 = paste0(conditions[1], '.', cluster)
  if (any(cluster_1 %in% seur$sample.cluster) && any(cluster_2 %in% seur$sample.cluster) && 
      length(colnames(subset(seur, sample.cluster==cluster_1)))>3 && 
      length(colnames(subset(seur, sample.cluster==cluster_2)))>3) {
    
    temp = FindMarkers(seur, ident.1=cluster_1, ident.2=cluster_2, logfc.threshold=0, min.pct=0, min.diff.pct=0)
    temp$avg_FC = (2^temp$avg_log2FC)
    temp = data.frame(gene=row.names(temp), cluster_1=cluster_1, cluster_2=cluster_2, avg_FC=temp$avg_FC,
                      avg_log2FC=temp$avg_log2FC, p_val=temp$p_val, 
                      pct.1=temp$pct.1, pct.2=temp$pct.2)
    degs[[i]] = temp
  }
}

# DEGs between the 2 conditions (total cells)
Idents(seur) = 'sample'
ind = length(degs)+1
degs[[ind]] = FindMarkers(seur, ident.1=conditions[2], ident.2=conditions[1], logfc.threshold=0, min.pct=0, 
                          min.diff.pct=0)
degs[[ind]] = data.frame(gene=row.names(degs[[ind]]),
                         cluster_1=paste0(conditions[2], '.All'), cluster_2=paste0(conditions[1], '.All'),
                         avg_FC=(2^degs[[ind]]$avg_log2FC), avg_log2FC=degs[[ind]]$avg_log2FC, 
                         p_val=degs[[ind]]$p_val, 
                         pct.1=degs[[ind]]$pct.1, pct.2=degs[[ind]]$pct.2)

degs = do.call(rbind, degs)
degs$p_val_adj = p.adjust(degs$p_val, method='fdr')
if (org=='Mouse') {
  annot_genes = as.data.frame(grcm38)
  annot_genes = annot_genes[annot_genes$chr %in% c(1:19, 'X', 'Y'), ]
  annot_genes = annot_genes[, c('ensgene', 'symbol', 'description')]
  annot_genes = unique(annot_genes[, -1])
  if (any(duplicated(annot_genes$symbol))) {
    annot_genes = annot_genes[!duplicated(annot_genes$symbol, fromLast=TRUE), ]
  }
} else if (org=='Human'){
  annot_genes = as.data.frame(grch38[, c('ensgene', 'symbol', 'description')])
  annot_genes = annot_genes[, -1]
}
degs = merge(annot_genes, degs, by.x='symbol', by.y='gene', all.y=T)
degs = degs[order(degs$cluster_1, degs$p_val, -degs$avg_log2FC), ]
write.csv(degs, file=snakemake@output[['degs_without_cutoff']], row.names=F)

filt_degs = degs[degs$p_val_adj <= 0.05 & abs(degs$avg_log2FC) >= 0.25 & (degs$pct.1 >= 0.1 | degs$pct.2 >= 0.1), ]
write.csv(filt_degs, file=snakemake@output[['degs_with_cutoff']], row.names=F)


