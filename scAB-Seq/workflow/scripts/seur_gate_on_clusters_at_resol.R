library(Seurat)

seur       = readRDS(snakemake@input[['seur']])
seur
cluster_id = snakemake@params[['cluster_id']]
resol      = snakemake@params[['resolution']]
cluster_id
class(cluster_id)
resol

meta_data   = seur@meta.data
this_idents = ifelse('ADT' %in% Assays(seur), paste0('wsnn_res.', resol), paste0('RNA_snn_res.', resol))
these_cells = row.names(meta_data)[meta_data[, this_idents] %in% cluster_id]
seur        = subset(seur, cells %in% these_cells)

saveRDS(seur, file=snakemake@output[['seur']])
