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
library(SCPA)
library(magrittr)
library(tidyr)
library(scRepertoire)


source('bin/00_functions.R')


out_dir = 'results_C0_in_D8'
out_dir_QC = file.path(out_dir, 'QC')
out_dir_seur = file.path(out_dir, 'seurat_object')
dir.create(out_dir_QC, recursive=T)
dir.create(out_dir_seur, recursive = T)

comb_seurs_ls = readRDS(file=file.path('results', 'seurat_object', '05_seurs_ls_with_scTCR_seq_data.rds'))

seurs = comb_seurs_ls$D8

# subset cluster 0 and re-cluster it to identify transcriptome associated with knocking out RBPJ 
seurs = subset(seurs, RNA_snn_res.0.15==0)
keep_genes = row.names(seurs@assays$RNA@counts)[rowSums(seurs@assays$RNA@counts>0) > 3]
seurs = subset(seurs, features=keep_genes)

## Preprocess the RNA assay to remove TCR genes from Variable features
## and visualize UMAP before integration
umap_p_ls = NULL
seurs = NormalizeData(seurs) %>%
  FindVariableFeatures(nfeatures=2000) 

TCR_genes  = grep('^Tr[abd]v', row.names(seurs), value=T)
seurs <- AddModuleScore(object=seurs, features=list(TCR_genes), name="TCR_Score")
seurs <- ScaleData(seurs, vars.to.regress = "TCR_Score1", features=row.names(seurs)) %>%  
  RunPCA() %>%
  RunUMAP(dims=1:30)
umap_p  = DimPlot(seurs, group.by="sample", reduction='umap')
pc_scores <- seurs@reductions$pca@cell.embeddings[, 1:2]
cor(pc_scores, seurs@meta.data[, c("TCR_Score1")])

## integrate samples using FastMNN
seur_ls = SplitObject(seurs, split.by='sample') # don't omit this line otherwise will have problem getting less integrated features than it should be 
## integrate samples using fastmnn
integ_mnn = RunFastMNN(lapply(seur_ls, DietSeurat),
                       reduction.name = 'integrated_mnn',
                       reduction.key = 'integrated_mnn_')
integ_mnn = RunUMAP(integ_mnn, reduction="integrated_mnn", dims=1:30)
# feed in variable features and scaled data for RNA assay from RNA assay of seurs_ls
VariableFeatures(integ_mnn) = VariableFeatures(seurs)
integ_mnn@assays$RNA@scale.data = seurs@assays$RNA@scale.data
integ_mnn@meta.data$sample = factor(integ_mnn$sample, levels=c('ON_D8', 'OFF_D8'))

## Visualize before and after FastMNN integration UMAPs
pdf(file.path(out_dir_QC, 'UMAP_before_n_after_D8_C0_integration_FastMNN.pdf'), width=11, height=5)
integ_umap_p = DimPlot(integ_mnn, group.by="sample")
p = ggarrange(umap_p + ggtitle('Before integration'),
              integ_umap_p + ggtitle('After integration'),
              nrow=1)
print(p)
dev.off()



##### 04. cluster all cells, identify CD8 T cells and recluster CD8 T cells ######
## 04.a clustering all cells
resols = c(0.1, 0.3)
bed.se = readRDS("/Users/dabdrabb/projects/ref/scrna/scRNAseq_datasets/immgen.rds") 

out_dir_res = file.path(out_dir, 'umaps_n_markers')
dir.create(out_dir_res, recursive=T)

graph_name  = 'RNA_snn'
for (resol in resols) {
  integ_seurs = cluster_cells(seur=integ_mnn, dims=1:30, 
                              reduction="integrated_mnn", graph_name=graph_name, resolution=resol, 
                              reduction_name="umap", reduction_key="UMAP_",
                              assays='RNA', clust_algorithm=1)
  integ_seurs = annotate_cells(seur=integ_seurs, graph_name=graph_name,
                               bed.se=bed.se, resolution=resol, assay='RNA')
  top_markers = get_plots(seur=integ_seurs, sample_col='sample', resolution=resol, out_dir=out_dir_res,
                          graph_name=graph_name, reduction='umap', suffix='C0_in_D8', top_nb=20)
  integ_mnn = integ_seurs
}


integ_seurs = integ_mnn
integ_seurs@meta.data = integ_seurs@meta.data[, -24:-26]
saveRDS(integ_seurs, file.path(out_dir_seur, '05_integrated_seurs_n_clustered_C0_in_D8.rds'))
# load here integ_seurs = readRDS(file.path(out_dir_seur, '05_integrated_seurs_n_clustered_C0_in_D8.rds'))

## 04.b visualize clusters percentage
out_dir_res = file.path(out_dir, 'clusters_percentage')
dir.create(out_dir_res)
res_vis = viz_clust_prop(integ_seurs, 'sample', 'seurat_clusters') 
write.csv(res_vis$df, row.names = F,
          file=file.path(out_dir_res, 'clusters_percentage.csv'))
p = res_vis$plot
ggsave(file.path(out_dir_res, 'clusters_precentage.pdf'), plot=p, height=6, width=3.5)




##### 05. Identify differentially expressed genes (DEGs) ######
## 05.a Identify DEGs
if (!exists('integ_seurs')) 
  integ_seurs = readRDS(file.path(out_dir_seur, '05_integrated_seurs_n_clustered_C0_in_D8.rds'))
out_dir_res = file.path(out_dir, 'DEGs')
dir.create(out_dir_res)
seurs     = integ_seurs
clusters = sort(as.numeric(as.character(unique(seurs$seurat_clusters))))

signif_DEGs = lapply(clusters, function(clust) {
  degs = get_DEGs(seurs, ident_1='ON_D8', group_by='sample', subset_ident=clust,
                  significant=T, degs_by='cluster', idents_col='seurat_clusters')
})
names(signif_DEGs) = paste0('C', clusters)
signif_DEGs[['Bulk']] = get_DEGs(seurs, ident_1='ON_D8', 
                                 ident_2='OFF_D8', significant=T,
                                 degs_by='sample', idents_col='sample')
write.xlsx(signif_DEGs,
           file=file.path(out_dir_res, 'C0_in_D8_DEGs_ON_vs_OFF.xlsx'), rowNames=F)


## 05.b. visualize all significant DEGs (log2FC and FDR on the level of each cluster and bulk) 
seurs      = integ_seurs
degs_file  = file.path(out_dir_res, 'C0_in_D8_DEGs_ON_vs_OFF.xlsx')
ident_1    = 'ON_D8'
ident_2    = 'OFF_D8'
genes = NULL

htp_data = get_logFC_n_FDR(seurs, degs_file=degs_file, significant_rows=F,
                           ident_1=ident_1, ident_2=ident_2)

p = logFC_n_FDR_htmp(htp_data, fontsize=7, use_yellow_purple=F, 
                     abs_logFC_cutoff=0, cluster_rows=T, title='Re-clustered cells of C0 in D8') 
ggsave(file.path(out_dir_res, 'C0_in_D8_Heatmap_log2FC_all_DEGs_per_cluster.pdf'),
       plot=p, width=2.5, height=5.5)





##### 06. GSEA analyses ########
## 06.a Prepare data for GSEA
if (!exists('integ_seurs')) 
  integ_seurs = readRDS(file.path(out_dir_seur, '05_integrated_seurs_n_clustered_C0_in_D8.rds'))
out_dir_res = file.path(out_dir, 'GSEA')
dir.create(out_dir_res)

seurs    = integ_seurs
clusters = sort(as.numeric(as.character(unique(seurs$seurat_clusters))))

all_DEGs_ls = lapply(clusters, function(clust) {
  degs = get_DEGs(seurs, ident_1='ON_D8', group_by='sample', subset_ident=clust,
                  significant=F, degs_by='cluster', idents_col='RNA_snn_res.0.3')
})
names(all_DEGs_ls) = paste0('C', clusters)
all_DEGs_ls[['Bulk']] = get_DEGs(seurs, ident_1='ON_D8', ident_2='OFF_D8', significant=F,
                                 degs_by='sample', idents_col='sample')
write.xlsx(all_DEGs_ls,
           file=file.path(out_dir_res, 'C0_in_D8_all_DEGs_ON_vs_OFF.xlsx'),
           rowNames=F)


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

for (j in 1:length(gene_sets_ls)) {
  gene_set_j = gene_sets_ls[[j]]
  # for (gene_set_j in gene_sets_ls) {
  clusters = getSheetNames(file.path(out_dir_res, 'C0_in_D8_all_DEGs_ON_vs_OFF.xlsx'))
  # standard way
  gsea_res = lapply(clusters, function(clust) {
    degs = read.xlsx(file.path(out_dir_res, 'C0_in_D8_all_DEGs_ON_vs_OFF.xlsx'), sheet=clust) 
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
  write.xlsx(gsea_res, file.path(out_dir_res, paste0('C0_in_D8_GSEA_ON_vs_OFF_', names(gene_sets_ls)[j], '.xlsx')))
}


## 06.c. Visualize GSEA results as dotplots
pdf_heights = c(D8=4, D50=6)
file_name = file.path(out_dir_res, 'C0_in_D8_GSEA_ON_vs_OFF_Hallmark.xlsx')
p = gsea_dotplot(file_name, title='Re-clustered cells of C0 in D8')

ggsave(file.path(out_dir_res, 'C0_in_D8_Dotplot_activated_n_inhibited_pathways_ON_vs_OFF.pdf'), 
       plot=p, width=8, height=2.5)  





## Abundance of top 10 clones with each sub-cluster
comb_TCR_ls = readRDS(file.path('results', 'rds_files', 'comb_TCR_ls.rds'))
comb_TCR_ls = list(ON_D8 = comb_TCR_ls$ON_D8[comb_TCR_ls$ON_D8$barcode %in% integ_seurs$cells, ],
                   OFF_D8 = comb_TCR_ls$OFF_D8[comb_TCR_ls$OFF_D8$barcode %in% integ_seurs$cells, ])

my_cols = rainbow(20)
out_dir_res = file.path('results_C0_in_D8', 'TCRs_percentage')
dir.create(out_dir_res, recursive = T)
pdf(file.path(out_dir_res, 'C0_in_D8_Top_10_clones_in_abundance.pdf'), width=5.5)
p = clonalCompare(comb_TCR_ls, 
                  top.clones = 10, 
                  samples = c("ON_D8", "OFF_D8"),
                  order.by = c("ON_D8", "OFF_D8"),
                  cloneCall="gene", 
                  graph = "alluvial") +
  scale_fill_manual(values=my_cols)
print(p)
dev.off()
