# perform the analyses on 2 combinated contigs (1. default (all False), 2. all True) 
# show top 10 clones

library(scRepertoire)
library(Seurat)
library(openxlsx)
library(ggpubr)
library(fgsea)
library(SCPA)
library(msigdbr)

source('bin/00_functions.R')


in_dir = file.path('resources', 'scTCR_data')
out_dir = file.path('results')
out_dir_seur = file.path(out_dir, 'seurat_object') 

samples = c('ON_D8', 'OFF_D8', 'ON_D50', 'OFF_D50')


seurs_ls = readRDS(file.path(out_dir_seur, '04_integrated_seurs_ls_n_clustered_CD8_T_cells.rds'))
seurs    = merge(seurs_ls$D8, seurs_ls$D50)
seur_ls  = SplitObject(seurs, 'sample') 


##### 01. read in contig annotations for each sample ####
contig_ls = lapply(samples, function(sample_i) {
  # day = sub('(.+)_(.+)', '\\2', sample_i)
  dir_name = file.path(in_dir, sample_i)
  contig_annot = read.csv(file.path(dir_name, "filtered_contig_annotations.csv"))
  contig_annot$barcode = gsub('-1', '', contig_annot$barcode)
  contig_annot$cells   = paste0(sample_i, '_', contig_annot$barcode)
  contig_annot = contig_annot[contig_annot$cells %in% Cells(seur_ls[[sample_i]]), ]
})
names(contig_ls) = samples


#### 02. Combine the list of T cell receptor contigs into clones ####
comb_TCR_ls = combineTCR(contig_ls, 
                         samples = samples,
                         removeNA = T, 
                         removeMulti = F, 
                         filterMulti = T)
saveRDS(comb_TCR_ls, file.path('results', 'rds_files', 'comb_TCR_ls.rds'))

#### 03. combine seurat and TCR data
comb_seur_ls <- lapply(samples, function(sample_i) {
  seur = combineExpression(comb_TCR_ls[[sample_i]], 
                           seur_ls[[sample_i]], 
                           group.by = 'sample',
                           cloneCall = "gene", 
                           proportion = TRUE,
                           filterNA = TRUE)
})
names(comb_seur_ls) = samples

### merge combined seurs for same day, add their umap coordinates from seurs_ls to generate UMAP of the top3 clones 
comb_seurs_ls = list(D8 = merge(comb_seur_ls$ON_D8, comb_seur_ls$OFF_D8), 
                     D50 = merge(comb_seur_ls$ON_D50, comb_seur_ls$OFF_D50))
for (day in c('D8', 'D50')) {
  seurs_i = seurs_ls[[day]]
  comb_seurs_i = comb_seurs_ls[[day]]
  cells_i = Cells(comb_seurs_i) 
  comb_seurs_i[['umap']] = CreateDimReducObject(seurs_i[['umap']]@cell.embeddings[cells_i, ],
                                                key="UMAP_", global=T, assay="RNA")
  comb_seurs_ls[[day]] = comb_seurs_i
}
## clean seurat object before saving it
comb_seurs_ls$D8@meta.data = comb_seurs_ls$D8@meta.data[, -14:-16]
comb_seurs_ls$D50@meta.data = comb_seurs_ls$D50@meta.data[, -10:-13]
saveRDS(comb_seurs_ls, file='results/seurat_object/05_seurs_ls_with_scTCR_seq_data.rds')
## load here comb_seurs_ls = readRDS(file='results/seurat_object/05_seurs_ls_with_scTCR_seq_data.rds')

my_cols = rainbow(20)
out_dir_perct = file.path('results', 'TCRs_percentage')
pdf(file.path(out_dir_perct, 'D8_n_D50_Top_10_clones_in_abundance.pdf'), width=5.5)
for (tp in c('D8', 'D50')) {
  p = clonalCompare(comb_TCR_ls, 
                    top.clones = 10, 
                    samples = paste0(c("ON_", "OFF_"), tp), 
                    cloneCall="gene", 
                    graph = "alluvial") +
    scale_fill_manual(values=my_cols)
  print(p)
}
dev.off()

## generate clones proportions and percentages per sample
clonal_percent_ls = lapply(comb_seur_ls, function(seur) {
  clonal_percent = unique(data.frame(clone = seur$CTgene, 
                                     clone_proportion = seur$clonalProportion,
                                     clone_percentage = round(seur$clonalProportion*100, 2),
                                     clone_frequency = seur$clonalFrequency,
                                     cluster_id = seur$seurat_clusters,
                                     sample = seur$sample))
  clonal_percent = clonal_percent[order(clonal_percent$clone_proportion, decreasing = T), ]
})
names(clonal_percent_ls) = samples

out_dir_perct = 'results/TCRs_percentage'
dir.create(out_dir_perct, recursive = T)
write.xlsx(clonal_percent_ls, file.path(out_dir_perct, 'TCRs_percentage.xlsx'))


#### generate umap of top 3 clones
top_TCRs_ls = lapply(samples, function(sample_i) {
  top_TCRs_i = read.xlsx(file.path(out_dir_perct, 'TCRs_percentage.xlsx'), sheet = sample_i)
  unlist(unique(top_TCRs_i$clone)[1:3])
})
names(top_TCRs_ls) = samples
p_ls = NULL
for (day in c('D8', 'D50')) {
  comb_seurs_i = comb_seurs_ls[[day]] 
  comb_seurs_i = highlightClones(comb_seurs_i, 
                                 cloneCall = 'gene', 
                                 sequence = unlist(top_TCRs_ls[grep(day, names(top_TCRs_ls))]))
  
  p = DimPlot(comb_seurs_i, 
          group.by = "highlight", 
          pt.size = 3, 
          split.by = 'sample',
          shuffle = T,
          raster = T) + 
    theme(legend.text = element_text(size=9)) + # unit(10, 'cm'))
    ggtitle('Top 3 clones per sample')
  p_ls[[day]] = p + 
    scale_color_manual(values = ggplot_build(p)$data[[1]]$colour %>% unique(), na.value = "grey") # change NA color from dark to light grey
}
p = ggarrange(plotlist = p_ls, ncol=1, nrow=2)
ggsave('results/Top_TCRs/UMAP_of_Top_3_clones_per_sample.pdf', plot=p, width=10, dpi=600)



### Identify top 3 clones in frequency in ON and top 3 clones in OFF and identify their DEGs comparing:
### Top1 clone in ON_D8 vs Top1 clone in OFF_D8, do the same for Top2 and Top3 clones in percentage
out_dir_degs = 'results/Top_TCRs'
dir.create(out_dir_degs, recursive = T)
top_TCRs_ls = lapply(samples, function(sample_i) {
  top_TCRs_i = read.xlsx(file.path(out_dir_perct, 'TCRs_percentage.xlsx'), sheet = sample_i)
  unlist(unique(top_TCRs_i$clone)[1:3])
})
names(top_TCRs_ls) = samples

# perform DEGs analyses
for (day in c('D8', 'D50')) {
  wb <- createWorkbook()
  for (i in 1:3) {
    sheet_i = paste('Top', i, 'clones')
    addWorksheet(wb, sheet_i)

    seur     = comb_seurs_ls[[day]]
    sample_1 = paste0('ON_', day)
    sample_2 = paste0('OFF_', day)
    seur$sample_clone = paste0(seur$sample, '_', seur$CTgene)
    Idents(seur) = 'sample_clone'
    ident_1 = paste0(sample_1, '_', top_TCRs_ls[[sample_1]][i]) 
    ident_2 = paste0(sample_2, '_', top_TCRs_ls[[sample_2]][i]) 
 
    degs        = FindMarkers(seur, ident.1=ident_1, ident.2=ident_2)
    degs$gene   = row.names(degs)
    degs$avg_FC = 2^(degs$avg_log2FC)
    degs        = degs[degs$p_val_adj < 0.05, c('gene', 'avg_FC', 'avg_log2FC', 'p_val_adj', 'pct.1', 'pct.2', 'p_val')]
    degs        = degs[order(degs$avg_log2FC, decreasing=T), ]
    comparison  = paste0('Cells in ', sample_1,' expressing ', top_TCRs_ls[[sample_1]][i], 
                         ' clone     VERSUS     cells in ', 
                         sample_2,' expressing ', top_TCRs_ls[[sample_2]][i], ' clone')
    
    writeData(wb, sheet_i, comparison, startCol = 1, startRow = 1)
    writeData(wb, sheet_i, degs, startCol = 1, startRow = 3, rowNames = F)
  }
  saveWorkbook(wb, file.path(out_dir_degs, paste0('DEGs_of_Top3_clones_in_', day, "_ON_vs_OFF.xlsx")), overwrite = TRUE)
}





#### assess diversity ####
res_ls = lapply(samples,  function(sample_i) 
  clonalDiversity(comb_seur_ls[[sample_i]], 
                  cloneCall = "gene", 
                  n.boots = 1000, 
                  metrics = "inv.simpson", 
                  group.by = 'seurat_clusters',
                  x.axis = 'sample',
                  exportTable = T)
)
diversity_df = dplyr::bind_rows(res_ls)
diversity_df$sample = factor(diversity_df$sample, levels=samples)
diversity_df$seurat_clusters = factor(diversity_df$seurat_clusters, levels=0:4)

p = ggplot(diversity_df, aes(seurat_clusters, inv.simpson, group=sample, colour = sample)) +
  geom_line() + 
  geom_point() +
  theme_bw() +
  labs(title='Clonal diversity', x='clusters', y='Inverse simpson score')
ggsave('results/TCRs_clonal_diversity.pdf', plo=p, height=5)



##### 06. GSEA analyses ########
## 06.a Prepare data for GSEA
if (!exists('comb_seurs_ls')) 
  comb_seurs_ls = readRDS(file.path(out_dir_seur, '05_seurs_ls_with_scTCR_seq_data.rds'))
out_dir_res = file.path(out_dir, 'Top_clones_GSEA')
dir.create(out_dir_res)


top_TCRs_ls = lapply(samples, function(sample_i) {
  top_TCRs_i = read.xlsx(file.path(out_dir_perct, 'TCRs_percentage.xlsx'), sheet = sample_i)
  unlist(unique(top_TCRs_i$clone)[1:3])
})
names(top_TCRs_ls) = samples

# identify DEGs including non-significant genes to perform GSEA analyses
for (day in c('D8', 'D50')) {
  wb <- createWorkbook()
  for (i in 1:3) { # 
    sheet_i = paste('Top', i, 'clones')
    addWorksheet(wb, sheet_i)
    
    seur     = comb_seurs_ls[[day]]
    sample_1 = paste0('ON_', day)
    sample_2 = paste0('OFF_', day)
    seur$sample_clone = paste0(seur$sample, '_', seur$CTgene)
    Idents(seur) = 'sample_clone'
    ident_1 = paste0(sample_1, '_', top_TCRs_ls[[sample_1]][i]) 
    ident_2 = paste0(sample_2, '_', top_TCRs_ls[[sample_2]][i]) 
    degs        = FindMarkers(seur, ident.1=ident_1, ident.2=ident_2, min.pct=0,
                              logfc.threshold=0)
    degs$gene   = row.names(degs)
    degs$avg_FC = 2^(degs$avg_log2FC)
    degs        = degs[, c('gene', 'avg_FC', 'avg_log2FC', 'p_val_adj', 'pct.1', 'pct.2', 'p_val')]
    degs        = degs[order(degs$avg_log2FC, decreasing=T), ]
    comparison  = paste0('Cells in ', sample_1,' expressing ', top_TCRs_ls[[sample_1]][i], 
                         ' clone     VERSUS     cells in ', 
                         sample_2,' expressing ', top_TCRs_ls[[sample_2]][i], ' clone')
    
    writeData(wb, sheet_i, comparison, startCol = 1, startRow = 1)
    writeData(wb, sheet_i, degs, startCol = 1, startRow = 3, rowNames = F)
  }
  saveWorkbook(wb, file.path(out_dir_res, paste0('DEGs_of_Top3_clones_in_', day, "_ON_vs_OFF_including_non_significant.xlsx")), overwrite = TRUE)
}


## 06.b GSEA analyses
# HALLMARK gene sets
hallmark_df = msigdbr(species = "Mus musculus", category = "H")
hallmark_df = hallmark_df[, c('gs_name', 'gene_symbol')]
hallmark_ls = hallmark_df %>% split(x=.$gene_symbol, f=.$gs_name)

kegg_df = msigdbr(species='Mus musculus', category = 'C2', subcategory = 'CP:KEGG')
kegg_df = kegg_df[, c('gs_name', 'gene_symbol')]
kegg_ls = kegg_df %>% split(x=.$gene_symbol, f=.$gs_name)

# C7 immunologic signature gene sets
C7_df = msigdbr(species='Mus musculus', category='C7', subcategory='IMMUNESIGDB')
C7_df = data.frame(gs_name=C7_df$gs_name, gene_symbol=C7_df$gene_symbol)
immune_ls = C7_df %>% split(x=.$gene_symbol, f=.$gs_name)

gene_sets_ls = list(Hallmark=hallmark_ls, KEGG=kegg_ls, Immune=immune_ls)

for (day in c('D8', 'D50')) {
  for (j in 1:length(gene_sets_ls)) {
    gene_set_j = gene_sets_ls[[j]]
    clones = getSheetNames(file.path(out_dir_res, paste0('DEGs_of_Top3_clones_in_', day, "_ON_vs_OFF_including_non_significant.xlsx")))
    gsea_res = lapply(clones, function(clone) {
      degs = read.xlsx(file.path(out_dir_res, 
                                 paste0('DEGs_of_Top3_clones_in_', day, "_ON_vs_OFF_including_non_significant.xlsx")), 
                       sheet=clone, startRow = 3) 
      degs$rank = sign(degs$avg_log2FC)*(-log10(degs$p_val))
      ranks = setNames(degs$rank, degs$gene) 
      ranks = sort(ranks, decreasing=T)
      ranks[ranks==Inf] = max(setdiff(ranks, Inf)) + 10
      ranks[ranks==-Inf] = min(setdiff(ranks, -Inf)) - 10
      set.seed(1234567)
      res = fgsea(pathways=gene_set_j, stats=ranks, minSize=15, maxSize=500, nproc=1)#, gseaParam=1.5)
      res$leadingEdge = sapply(res$leadingEdge, function(x) paste0(x, collapse=', '))
      res = res[order(res$pval), ]
      res = res[res$padj<0.1, ]
    })
    names(gsea_res) = clones
    write.xlsx(gsea_res, file.path(out_dir_res, paste0(day, '_GSEA_ON_vs_OFF_', names(gene_sets_ls)[j], '.xlsx')))
  }
}