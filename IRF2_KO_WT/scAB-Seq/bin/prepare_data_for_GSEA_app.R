library(msigdbr)
library(Seurat)
library(fgsea)
library(msigdbi)

# write gmt files
Hallmark_df = msigdbr(species='Mus musculus', category='H')
gene_sets = Hallmark_df %>% split(x=.$gene_symbol, f=.$gs_name)
gene_sets_output = lapply(1:length(gene_sets), function(i) paste0(names(gene_sets)[i], '\t', 'na\t', paste(gene_sets[[i]], collapse = '\t'),'\t'))
sink("results/Seurat_AB_seq/Tum_CD8_v2/res_0.5/GSEA_app/Hallmark.gmt")
writeLines(unlist(lapply(gene_sets_output, paste, collapse=" ")))
sink()

GO_BP_df = msigdbr(species='Mus musculus', category='C5', subcategory='GO:BP')
gene_sets = GO_BP_df %>% split(x=.$gene_symbol, f=.$gs_name)
gene_sets_output = lapply(1:length(gene_sets), function(i) paste0(names(gene_sets)[i], '\t', 'na\t', paste(gene_sets[[i]], collapse = '\t'),'\t'))
sink("results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/GSEA_app/GO_BP.gmt")
writeLines(unlist(lapply(gene_sets_output, paste, collapse=" ")))
sink()

sel_go_bp = unique(grep('INFLAMMATORY_RESPONSE|INTERFERON_ALPHA|INTERFERON_GAMMA|APOPTOTIC_PROCESS|NECROTIC_CELL_DEATH|CYTOKINE_MEDIATED_SIGNALING|TORC1_SIGNALING|HYPOXIA|COMPLEMENT', GO_BP_df$gs_name, value = T))
gene_sets = gene_sets[sel_go_bp]
length(gene_sets)
gene_sets_output = lapply(1:length(gene_sets), function(i) paste0(names(gene_sets)[i], '\t', 'na\t', paste(gene_sets[[i]], collapse = '\t'),'\t'))
sink("results/Seurat_AB_seq/Tum_CD8_v2/res_0.5/GSEA_app/GO_BP_subset.gmt")
writeLines(unlist(lapply(gene_sets_output, paste, collapse=" ")))
sink()

# write expression data file
seur = readRDS("results/Seurat_AB_seq/Tum_CD8_v3/Seur_Tum_CD8_clustered.rds")
cluster_by = 'wsnn_res.0.4'
clust_1    = 'Tum_KO.3'
clust_2    = 'Tum_WT.0'
meta_data = seur@meta.data
meta_data$sample_cluster = paste(meta_data$sample, meta_data[, cluster_by], sep='.')
seur@meta.data = meta_data
seur = subset(seur, sample_cluster %in% c(clust_1, clust_2))
expr_data = GetAssayData(seur, assay='RNA', slot='data')
expr_data = data.frame(Name=row.names(expr_data), expr_data, check.names = F)
expr_data_output = expr_data
names(expr_data_output) = gsub('-', '_', names(expr_data_output))
write.table(expr_data_output, file='results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/GSEA_app/expr_data.txt', 
            sep='\t', row.names=F, quote=F)


# write pheno data file
cells_clust = sapply(Cells(seur), function(cell)
  ifelse(cell %in% Cells(subset(seur, sample_cluster=='Tum_KO.3')), 'KO_C3', 'WT_CO'))
pheno_data = list(paste(c(nrow(seur@meta.data), 2, 1), collapse='\t'),
                  c('#\tWT_C0\tKO_C3'),
                  #c(as.numeric(names(expr_data)[-1] %in% Cells(subset(seur, sample_cluster=='Tum_KO.3')))))
                  paste(cells_clust, collapse='\t'))
sink("results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/GSEA_app/pheno_data.CLS")
writeLines(unlist(lapply(pheno_data, paste, collapse=" ")))
sink()


# write ranked gene list file
degs  = read.csv('results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/Tum_CD8_DEGs_C3.KO_vs_C0_WT_without_cutoff.csv')
degs$rank = sign(degs$avg_log2FC)*(-log10(degs$p_val))
degs = degs[order(degs$rank, decreasing=T), ]
ranks = degs[, c('gene', 'rank')]
names(ranks)[1] = 'GeneName'
write.table(ranks, file="results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/GSEA_app/ranked_gene_list.RNK", 
            sep='\t', row.names=F, quote=F)

