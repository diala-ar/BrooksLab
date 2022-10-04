#library(dplyr)
library(Seurat)
library(ggplot2)

conds = read.csv(snakemake@params[['conds']], row.names=NULL, header=T)$conditions

merged_seur = readRDS(snakemake@input[['seur']])
metadata    = merged_seur@meta.data

nFeature_plot = ggplot(metadata, aes(x=nFeature_RNA, color=sample, fill=sample)) + 
  geom_histogram(bins=100) + 
  facet_wrap(~sample, scales='free') +
  theme_bw()
nCount_plot = ggplot(metadata, aes(x=nCount_RNA, color=sample, fill=sample)) + 
  geom_density() + 
  facet_wrap(~sample, scales='free') +
  theme_bw()
mito_plot = ggplot(metadata, aes(x=pc_mito, color=sample, fill=sample)) + 
  geom_density() + 
  facet_wrap(~sample, scales='free') +
  theme_bw()
ribo_plot = ggplot(metadata, aes(x=pc_ribo, color=sample, fill=sample)) + 
  geom_density() + 
  facet_wrap(~sample, scales='free') +
  theme_bw()
nGenesPerUMI_plot = ggplot(metadata, aes(x=log10GenesPerUMI, color=sample, fill=sample)) + 
  geom_histogram(bins=100) + 
  facet_wrap(~sample, scales='free') +
  theme_bw()

pdf(file=snakemake@output[['hist']], height=4.5, width=8)
print(nFeature_plot)
print(nCount_plot)
print(mito_plot)
print(ribo_plot)
print(nGenesPerUMI_plot)
dev.off()

pdf(file=snakemake@output[['violon']], height=5)
print(VlnPlot(merged_seur, features=c("nFeature_RNA"), group.by="sample")) 
print(VlnPlot(merged_seur, features=c("pc_mito"), group.by="sample")) 
print(VlnPlot(merged_seur, features=c("pc_ribo"), group.by="sample")) 
print(VlnPlot(merged_seur, features=c("nCount_RNA"), group.by="sample")) 
print(VlnPlot(merged_seur, features=c("log10GenesPerUMI"), group.by="sample")) 
Idents(merged_seur) = 'sample'
for (i in 1:length(conds)) {
  plot1 = FeatureScatter(subset(merged_seur, idents=conds[i]), feature1="nCount_RNA", feature2="pc_mito") + ggtitle(conds[i]) + theme(legend.position='none')
  plot2 = FeatureScatter(subset(merged_seur, idents=conds[i]), feature1="nCount_RNA", feature2="nFeature_RNA") + ggtitle(conds[i]) + theme(legend.position='none')
  print(plot1 + plot2)
}
dev.off()
