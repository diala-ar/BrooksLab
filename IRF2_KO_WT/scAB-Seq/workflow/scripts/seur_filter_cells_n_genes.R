library(Seurat)

conds = read.csv(snakemake@params[['conds']], row.names=NULL, header=T)$conditions
criteria_cutoffs = read.csv(snakemake@params[['criteria_cutoffs_file']], row.names=1)


# Filter out low quality cells using selected thresholds
merged_seur = readRDS(snakemake@input[['seur']])
seur_list   = SplitObject(merged_seur, split.by="sample")


print("Filter out cells with low quality reads using selected thresholds")
filt_seur = NULL
for (i in 1:length(seur_list)) {
  filt_seur[[i]] <- subset(x=seur_list[[conds[i]]],
                           subset=(nCount_RNA >= criteria_cutoffs[conds[i], 'nCount_RNA']) &
                             (nFeature_RNA >= criteria_cutoffs[conds[i], 'nFeature_low']) &
                             (nFeature_RNA < criteria_cutoffs[conds[i], 'nFeature_up']) &
                             (pc_mito < criteria_cutoffs[conds[i], 'pc_mito']))

  print(paste0("Filtered ", conds[i], ': ', ncol(filt_seur[[i]]), ' cells (', round(ncol(filt_seur[[i]])/ncol(seur_list[[i]]),2)*100, '%)'))
}  
names(filt_seur) = conds

filt_seur = merge(x=filt_seur[[1]], y=c(filt_seur[[2]]))
metadata  = filt_seur@meta.data
metadata$orig.ident = factor(metadata$orig.ident, levels=conds)  # keep to render plot in WT, KO order 
metadata$sample     = metadata$orig.ident                        # keep to render plot in WT, KO order 
filt_seur@meta.data = metadata
print(paste0("Merging '", paste(conds, collapse=', '), "' done! Dimensions of filtered seurat object: ", nrow(filt_seur), ' genes, ', ncol(filt_seur), ' cells'))


# filter out genes expressed in less than 10 cells
rna_counts  = GetAssayData(object=filt_seur, assay='RNA', slot="counts")
adt_counts  = GetAssayData(object=filt_seur, assay='ADT', slot="counts")
nonzero     = rna_counts > 0
keep_genes  = Matrix::rowSums(nonzero) >= 10
filt_counts = rna_counts[keep_genes, ]
filt_seur   = CreateSeuratObject(filt_counts, meta.data=metadata)
filt_seur[['ADT']] = CreateAssayObject(counts=adt_counts[, row.names(metadata)])


date()
saveRDS(filt_seur, file=snakemake@output[['seur']])
date()
##************************************End STEP4:Quality control of filtered cells******************************

