# convert RDS files containing seurat objects to loom file which will be used by arboreto_with_multiprocessing
# generate tsv files that will be used by aucell

library(Seurat)
library(SCopeLoomR)
library(SCENIC)

setwd('~/projects/sabelo/AB_seq_v2/')
system('mkdir -vp results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/pyscenic_100_runs')

convert_seur_to_loom = function(seur_obj, file_name) {
  exp_mat  = as.data.frame(GetAssayData(seur_obj, assay='RNA', slot='counts'))
  dim(exp_mat)
  
  # create loom object
  loom <- build_loom(file_name, dgem=exp_mat)
  cellInfo <- data.frame(seuratCluster=Idents(seur_obj))
  row.names(cellInfo) = colnames(exp_mat)
  loom <- add_cell_annotation(loom, cellInfo)
  close_loom(loom)
}


conds = c('Tum_CD8', 'Tum_CD8_KO.C3_WT.C0')
for (cond in conds) {
  seur = readRDS(paste0('results/Seurat_AB_seq/Tum_CD8_v3/Seurat/Seur_', cond, '_clustered.rds'))
  convert_seur_to_loom(seur, paste0('results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/pyscenic_100_runs/Seur_', cond, '.loom'))
  write.table(seur[['RNA']]@counts, file=paste0('results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/pyscenic_100_runs/Seur_', cond, '_counts.tsv'), sep='\t')
}


