# Analyse D8 dataset apart and D50 dataset apart and integrate each dataset using FastMNN.
# remove unwanted genes from variable features in RNA assay
# unwanted genes are TCRs in D8


library(Seurat)
library(ggplot2)
library(scater)
library(batchelor)
library(SeuratWrappers)
library(ggpubr)
library(SingleR)
library(dplyr)
library(openxlsx)
library(msigdbr)
library(fgsea)
library(magrittr)
library(tidyr)




in_dir = 'resources'
out_dir = file.path('results')
out_dir_seur = file.path(out_dir, 'seurat_object') 
out_dir_QC   = file.path(out_dir, 'QC') 

dir.create(out_dir_seur, recursive=T)
dir.create(out_dir_QC, recursive=T)

samples = c('ON_D8', 'OFF_D8', 'ON_D50', 'OFF_D50')
groups  = c('D8', 'D50')

source('bin/00_functions.R')

##### 01. Create seurat object ####
seur_ls = sapply(samples, function(sample_i) {
  dir_name = file.path(in_dir, sample_i, 'filtered_feature_bc_matrix')
  seur     = create_RNA_seurat(dir_name=dir_name, specie='mouse', 
                               project_name='scRNA_Rbpj_ON_OFF', sample_id=sample_i)
})
saveRDS(seur_ls, file.path(out_dir_seur, '01_seur_ls_before_filtering_cells_n_genes.rds'))


##### 02. Quality control ####
if (!exists('seur_ls'))
  seur_ls = readRDS(file.path(out_dir_seur, '01_seur_ls_before_filtering_cells_n_genes.rds'))
seurs = merge(seur_ls[[1]], seur_ls[-1])
seurs$sample = factor(seurs$sample, levels=samples)

## plot RNA features before filtering poor quality cells
feats = c("nFeature_RNA", "nCount_RNA", "mito_prct", "ribo_prct", "hmgb_prct", "log10GenesPerUMI")
print("Visualize QC plots and metrics before filtering")
p = VlnPlot(object=seurs, features=feats, pt.size=0, ncol=length(feats), group.by='sample') + xlab('')
ggsave(file.path(out_dir_QC, 'ViolonPlot_before_dead_cells_filtering.pdf'), plot=p, height=4.5, width=18)

## perform quality control 
filt_seur_ls = lapply(seur_ls, function(seur) {
  kept_cells = keep_cells(seur, nmads_count=3, nmads_feature=4, nmads_mito=5)
  selected_f = rownames(seur)[Matrix::rowSums(seur) > 10]
  
  subset(seur, cells %in% kept_cells, features=selected_f)
})

for (sample_i in samples) {
  print(prct_filtered_cells_n_genes(seur_ls[[sample_i]], filt_seur_ls[[sample_i]], sample_i))
}

## plot rna features after filtering
filt_seurs = merge(filt_seur_ls[[1]], filt_seur_ls[-1])
filt_seurs$sample = factor(filt_seurs$sample, levels=samples)
feats = c("nFeature_RNA", "nCount_RNA", "mito_prct", "ribo_prct", "hmgb_prct", "log10GenesPerUMI")
print("Visualize QC plots and metrics after filtering")
p = VlnPlot(object=filt_seurs, features=feats, pt.size=0, ncol=length(feats), group.by='sample') + xlab('')
ggsave(file.path(out_dir_QC, 'ViolonPlot_after_dead_cells_filtering.pdf'), plot=p, height=4.5, width=18)

saveRDS(filt_seur_ls, file.path(out_dir_seur, '02_seur_ls_after_filtering_cells_n_genes.rds'))



##### 03. Integrate samples using FastMNN ####
## integrate samples, remove TCR genes from D8 VariableFeatures
seur_ls  = readRDS(file.path(out_dir_seur, '02_seur_ls_after_filtering_cells_n_genes.rds'))
seurs_ls = list(merge(seur_ls[[1]], seur_ls[[2]]),
                merge(seur_ls[[3]], seur_ls[[4]]))
names(seurs_ls) = groups
seurs_ls[['D8']]$sample  = factor(seurs_ls[['D8']]$sample, c('ON_D8', 'OFF_D8'))
seurs_ls[['D50']]$sample = factor(seurs_ls[['D50']]$sample, c('ON_D50', 'OFF_D50'))

## Preprocess the RNA assay to remove TCR genes from Variable features
## and visualize UMAP before integration
umap_p_ls = NULL
for (group_i in groups) {
  seurs_ls[[group_i]] = NormalizeData(seurs_ls[[group_i]]) %>%
    FindVariableFeatures(nfeatures=2000) 
  if (group_i == 'D8') {
    TCR_genes  = grep('^Tr[abd]v', row.names(seurs_ls[[group_i]]), value=T)
    VariableFeatures(seurs_ls[[group_i]]) = setdiff(VariableFeatures(seurs_ls[[group_i]]), TCR_genes)
  }
  seurs_ls[[group_i]] = ScaleData(seurs_ls[[group_i]], features=row.names(seurs_ls[[group_i]])) %>%  
    RunPCA() %>%
    RunUMAP(dims=1:50)
  umap_p_ls[[group_i]]  = DimPlot(seurs_ls[[group_i]], group.by="sample", reduction='umap')
}

## integrate samples using FastMNN
integ_mnn_ls = NULL
for (group_i in groups) {
  seur_ls = SplitObject(seurs_ls[[group_i]], split.by='sample') # don't omit this line otherwise will have problem getting less integrated features than it should be 
  ## integrate samples using fastmnn
  integ_mnn_ls[[group_i]] = RunFastMNN(lapply(seur_ls, DietSeurat),
                                       reduction.name = 'integrated_mnn',
                                       reduction.key = 'integrated_mnn_')
  
  integ_mnn_ls[[group_i]] = RunUMAP(integ_mnn_ls[[group_i]], reduction="integrated_mnn", dims=1:50)
  # feed in variable features and scaled data for RNA assay from RNA assay of seurs_ls
  VariableFeatures(integ_mnn_ls[[group_i]]) = VariableFeatures(seurs_ls[[group_i]])
  integ_mnn_ls[[group_i]]@assays$RNA@scale.data = seurs_ls[[group_i]]@assays$RNA@scale.data
}
integ_mnn_ls[['D8']]@meta.data$sample = factor(integ_mnn_ls[['D8']]$sample, levels=c('ON_D8', 'OFF_D8'))
integ_mnn_ls[['D50']]@meta.data$sample = factor(integ_mnn_ls[['D50']]$sample, levels=c('ON_D50', 'OFF_D50'))

saveRDS(integ_mnn_ls, file.path(out_dir_seur, '03_integrated_seurs_ls_FastMNN.rds'))


## Visualize before and after FastMNN integration UMAPs
integ_umap_p_ls = NULL
pdf(file.path(out_dir_QC, 'UMAP_before_n_after_samples_integration_FastMNN.pdf'), width=11, height=5)
for (group_i in groups) {
  integ_umap_p_ls[[group_i]] = DimPlot(integ_mnn_ls[[group_i]], group.by="sample")
  p = ggarrange(umap_p_ls[[group_i]] + ggtitle('Before integration'),
                integ_umap_p_ls[[group_i]] + ggtitle('After integration'),
                nrow=1)
  print(p)
}
dev.off()



##### 04. cluster all cells, identify CD8 T cells and recluster CD8 T cells ######
## 04.a clustering all cells
resol  = 0.1
bed.se = readRDS("/Users/dabdrabb/OneDrive - UHN/projects/ref/scrna/scRNAseq_datasets/immgen.rds") 
if (!exists('integ_mnn_ls')) 
  integ_mnn_ls = readRDS(file.path(out_dir_seur, '03_integrated_seurs_ls_FastMNN.rds'))

out_dir_res = file.path(out_dir, 'umaps_n_markers')
dir.create(out_dir_res, recursive=T)

graph_name  = 'RNA_snn'
for (group_i in groups) {
  integ_seurs = cluster_cells(seur=integ_mnn_ls[[group_i]], dims=1:50, 
                              reduction="integrated_mnn", graph_name=graph_name, resolution=resol, 
                              reduction_name="umap", reduction_key="UMAP_",
                              assays='RNA', clust_algorithm=1)
  integ_seurs = annotate_cells(seur=integ_seurs, graph_name=graph_name, 
                               bed.se=bed.se, resolution=resol, assay='RNA')
  top_markers = get_plots(seur=integ_seurs, sample_col='sample', resolution=resol, out_dir=out_dir_res,
                          graph_name=graph_name, reduction='umap', suffix=paste0(group_i, '_All_cells'))
  integ_mnn_ls[[group_i]] = integ_seurs
}


#### 04.b filtering out Monocytes and reclustering
for (group_i in groups) {
  integ_mnn_ls[[group_i]] = subset(integ_mnn_ls[[group_i]], immgen_res_0.1_label.main != 'Monocytes') #%in% c('T cells', 'NKT'))
  integ_mnn_ls[[group_i]]@meta.data = integ_mnn_ls[[group_i]]@meta.data[, -10:-13]
}
integ_mnn_ls_o = integ_mnn_ls


## 04.c re-clustering CD8 T cells
integ_mnn_ls2 = NULL
graph_name  = 'RNA_snn'
resols = c(0.1, 0.15, 0.2)
for (group_i in groups) {
  for (resol in resols) {
    integ_seurs = cluster_cells(seur=integ_mnn_ls[[group_i]], dims=1:50, 
                                reduction="integrated_mnn", graph_name=graph_name, resolution=resol, 
                                reduction_name="umap", reduction_key="UMAP_",
                                assays='RNA', clust_algorithm=1)
    integ_seurs = annotate_cells(seur=integ_seurs, graph_name=graph_name, 
                                 bed.se=bed.se, resolution=resol, assay='RNA')
    top_markers = get_plots(seur=integ_seurs, sample_col='sample', resolution=resol, out_dir=out_dir_res,
                            graph_name=graph_name, reduction='umap', suffix=paste0(group_i, '_CD8_T_cells'))
    ind = paste0(group_i, '_resol_', resol)
    integ_mnn_ls2[[ind]] = integ_seurs
  }
}

integ_seurs_ls = integ_mnn_ls2[c('D8_resol_0.15', 'D50_resol_0.2')]
names(integ_seurs_ls) = c('D8', 'D50')
saveRDS(integ_seurs_ls, file.path(out_dir_seur, '04_integrated_seurs_ls_n_clustered_CD8_T_cells.rds'))


## 04.d visualize clusters percentage
out_dir_res = file.path(out_dir, 'clusters_percentage')
dir.create(out_dir_res)
p_ls = NULL
for (group_i in groups) {
  res_vis = viz_clust_prop(integ_seurs_ls[[group_i]], 'sample', 'seurat_clusters') 
  write.csv(res_vis$df, row.names = F,
            file=file.path(out_dir_res, paste0(group_i, '_clusters_percentage.csv')))
  p_ls[[group_i]] = res_vis$plot
}
p = ggarrange(plotlist = p_ls, nrow=1)
ggsave(file.path(out_dir_res, 'clusters_precentage.pdf'), plot=p, height=6)




##### 05. Identify differentially expressed genes (DEGs) ######
if (!exists('integ_seurs_ls')) 
  integ_seurs_ls = readRDS(file.path(out_dir_seur, '04_integrated_seurs_ls_n_clustered_CD8_T_cells.rds'))
out_dir_res = file.path(out_dir, 'DEGs')
dir.create(out_dir_res)
for (group_i in groups) {
  seur     = integ_seurs_ls[[group_i]]
  clusters = sort(as.numeric(as.character(unique(seur$seurat_clusters))))
    
  signif_DEGs = lapply(clusters, function(clust) {
      degs = get_DEGs(seur, ident_1=paste0('ON_', group_i), group_by='sample', subset_ident=clust,
                      significant=T, degs_by='cluster', idents_col='seurat_clusters')
    })
  names(signif_DEGs) = paste0('C', clusters)
  signif_DEGs[['Bulk']] = get_DEGs(seur, ident_1=paste0('ON_', group_i), 
                                   ident_2=paste0('OFF_', group_i), significant=T,
                                   degs_by='sample', idents_col='sample')
  write.xlsx(signif_DEGs,
             file=file.path(out_dir_res, paste0(group_i, '_DEGs_ON_vs_OFF.xlsx')),
             rowNames=F)
}



##### 06. GSEA analyses ########
## 06.a Prepare data for GSEA
if (!exists('integ_seurs_ls')) 
  integ_seurs_ls = readRDS(file.path(out_dir_seur, '04_integrated_seurs_ls_n_clustered_CD8_T_cells.rds'))
out_dir_res = file.path(out_dir, 'GSEA')
dir.create(out_dir_res)

for (group_i in groups) {
  seur     = integ_seurs_ls[[group_i]]
  clusters = sort(as.numeric(as.character(unique(seur$seurat_clusters))))
  
  all_DEGs_ls = lapply(clusters, function(clust) {
    degs = get_DEGs(seur, ident_1=paste0('ON_', group_i), group_by='sample', subset_ident=clust,
                    significant=F, degs_by='cluster', idents_col='seurat_clusters')
  })
  names(all_DEGs_ls) = paste0('C', clusters)
  all_DEGs_ls[['Bulk']] = get_DEGs(seur, ident_1=paste0('ON_', group_i), 
                                   ident_2=paste0('OFF_', group_i), significant=F,
                                   degs_by='sample', idents_col='sample')
  write.xlsx(all_DEGs_ls,
             file=file.path(out_dir_res, paste0(group_i, '_all_DEGs_ON_vs_OFF.xlsx')),
             rowNames=F)
}


## 06.b GSEA analyses
# HALLMARK gene sets
hallmark_df = msigdbr(species = "Mus musculus", category = "H")
hallmark_df = hallmark_df[, c('gs_name', 'gene_symbol')]
hallmark_ls = hallmark_df %>% split(x=.$gene_symbol, f=.$gs_name)

# C7 immunologic signature gene sets
C7_df = msigdbr(species='Mus musculus', category='C7', subcategory='IMMUNESIGDB')
C7_df = data.frame(gs_name=C7_df$gs_name, gene_symbol=C7_df$gene_symbol)
immune_ls = C7_df %>% split(x=.$gene_symbol, f=.$gs_name)

gene_sets_ls = list(Hallmark=hallmark_ls, Immune=immune_ls)

for (group_i in groups) {
  for (j in 1:length(gene_sets_ls)) {
    gene_set_j = gene_sets_ls[[j]]
    # for (gene_set_j in gene_sets_ls) {
    clusters = getSheetNames(file.path(out_dir_res, paste0(group_i, '_all_DEGs_ON_vs_OFF.xlsx')))
    # standard way
    gsea_res = lapply(clusters, function(clust) {
      degs = read.xlsx(file.path(out_dir_res, paste0(group_i, '_all_DEGs_ON_vs_OFF.xlsx')), sheet=clust) 
      degs$rank = sign(degs$avg_log2FC)*(-log10(degs$p_val))
      ranks = setNames(degs$rank, degs$gene) 
      ranks = sort(ranks, decreasing=T)
      ranks[ranks==Inf] = max(setdiff(ranks, Inf)) + 10
      ranks[ranks==-Inf] = min(setdiff(ranks, -Inf)) - 10
      set.seed(27912)
      res = fgsea(pathways=gene_set_j, stats=ranks, minSize=15, maxSize=500, nproc=1)#, gseaParam=1.5)
      res$leadingEdge = sapply(res$leadingEdge, function(x) paste0(x, collapse=', '))
      res = res[order(res$pval), ]
      res = res[res$padj<0.05, ]
    })
    names(gsea_res) = clusters
    write.xlsx(gsea_res, file.path(out_dir_res, paste0(group_i, '_GSEA_ON_vs_OFF_', names(gene_sets_ls)[j], '.xlsx')))
  }
}



