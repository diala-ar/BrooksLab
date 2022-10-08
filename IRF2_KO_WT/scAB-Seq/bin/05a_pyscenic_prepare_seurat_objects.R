# Preprocess data: 
# 1. Rename some meta_data columns to prepare data for pyscenic analyses
# 2. Subset data to include only clusters of interest (cluster 3 from KO cells and cluster 0 from WT cells)

setwd('~/projects/sabelo/AB_seq_v2/')
library(Seurat)

seur  = readRDS("results/Seurat_AB_seq/Tum_CD8_v3/Seur_Tum_CD8_clustered.rds")
seur$seurat_clusters = seur$wsnn_res.0.4
saveRDS(seur, file="results/Seurat_AB_seq/Tum_CD8_v3/Seurat/Seur_Tum_CD8_clustered.rds")


seur$seurat_clusters = paste(seur$sample, seur$wsnn_res.0.4, sep='.C')
seur  = subset(seur, seurat_clusters %in% c('Tum_WT.C0', 'Tum_KO.C3'))
saveRDS(seur, file=paste0('results/Seurat_AB_seq/Tum_CD8_v3/Seurat/Seur_Tum_CD8_KO.C3_WT.C0_clustered.rds'))