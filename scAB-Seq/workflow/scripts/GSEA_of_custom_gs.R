setwd('~/projects/sabelo/AB_seq_v2/')
library(msigdbr)
library(Seurat)
library(ggpubr)
library(fgsea)

# read in and preprocess the gene sets of interest
C7_df = msigdbr(species='Mus musculus', category='C7', subcategory='IMMUNESIGDB')
C7_df = C7_df[grep('EFFECTOR_VS_EXHAUSTED_CD8_TCELL', C7_df$gs_name), ]
msigbdr_gs = data.frame(gs_name=C7_df$gs_name, gene_symbol=C7_df$gene_symbol)
msigbdr_gs[msigbdr_gs$gs_name=='GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN', 'gs_name'] = 'EFFector VS Exhausted_DN'
msigbdr_gs[msigbdr_gs$gs_name=='GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP', 'gs_name'] = 'EFFector VS Exhausted_UP'

# Read in data taken from Khan et al., 2019, Nature, PMID: 31207603, tox_degs taken from Table S2
tox_degs                    = read.csv(snakemake@input[['tox_degs']], skip=2)
tox_up_with_adj_pval_cutoff = tox_degs[tox_degs$Adjusted.p.value<0.00001 & tox_degs$Log2.Fold.Change>0, ]
tox_up_with_adj_pval_cutoff = data.frame(gs_name='Tox UP DEGs', gene_symbol=tox_up_with_adj_pval_cutoff$Gene.ID)
tox_dwn_with_adj_pval_cutoff = tox_degs[tox_degs$Adjusted.p.value<0.00001 & tox_degs$Log2.Fold.Change<0, ]
tox_dwn_with_adj_pval_cutoff = data.frame(gs_name='Tox DN DEGs', gene_symbol=tox_dwn_with_adj_pval_cutoff$Gene.ID)

gene_sets = rbind(msigbdr_gs, tox_up_with_adj_pval_cutoff, 
                  tox_dwn_with_adj_pval_cutoff)
gene_sets_list = gene_sets %>% split(x=.$gene_symbol, f=.$gs_name)

degs       = read.csv(snakemake@input[['degs']])


degs$rank = sign(degs$avg_log2FC)*(-log10(degs$p_val))
degs = degs[order(degs$rank, decreasing=T), ]
ranks = setNames(degs$rank, degs$gene) 
set.seed = 5779
res   = fgsea(pathways=gene_sets_list, stats=ranks, minSize=5, maxSize=3500, 
              nPermSimple=1000, eps=0, gseaParam=1, nproc=1)
res = res[c(2,4), ]
gene_sets_list = gene_sets_list[c(2, 4)]

p_list = lapply(1:nrow(res), function(i) {
  title = paste0(res$pathway[i], ' (NES=', round(res$NES[i], digits=2), ', FDR=',format(res$padj[i], scientific=T, digits=2), ')')
  plotEnrichment(gene_sets_list[[i]], ranks) + labs(title=title) + theme(title=element_text(size=7))
})
p = ggarrange(plotlist=p_list, ncol=1, nrow=2)
ggsave(snakemake@output[["gsea_plots"]], plot=p, width=5, height=6) 


