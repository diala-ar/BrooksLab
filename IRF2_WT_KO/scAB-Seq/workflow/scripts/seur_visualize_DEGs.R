library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(pheatmap)
library(ggiraph)
library(gridExtra)

source('~/projects/bin/functions/heatmaps.R')


conditions = read.csv(snakemake@params[['conds']])$conditions
degs     = read.csv(snakemake@input[['degs_without_cutoff']])
resol    = snakemake@params[['resolution']]
seur     = readRDS(file=snakemake@input[['seur']])
isgs     = read.csv(snakemake@input[['DE_ISGs']], skip=19, sep='\t', header=T)
genes_df = read.csv(snakemake@input[['genes_of_interest']])


# sub sample seurat to smallest size
Idents(seur) = 'sample'
min_sample_size = min(table(seur$sample))
sub_sampled = subset(x=seur, downsample=min_sample_size) 
sub_cells   = Cells(sub_sampled)
sub_sampled


# generate volcano plot: static and interactive plots
degs = degs[degs$cluster_1==paste0(conditions[2],'.All'), ]
degs$is.significant = degs$p_val_adj < 0.05
degs$label = ifelse(degs$p_val_adj<0.05, degs$symbol, '')
head(degs)
# static plot
pdf(snakemake@output[['volcano_plot']])
ggplot(degs, aes(x=avg_log2FC, y=-log10(p_val_adj), label=label)) + 
  geom_text_repel(size=5) +
  geom_point(aes(color=is.significant)) +
  scale_color_manual(values = c("black", "red")) +
  geom_hline(yintercept=1.3, color='red') +
  geom_vline(xintercept=-1, color='blue3') +
  geom_vline(xintercept=1, color='blue3') +
  xlab('log2 fold change (KO vs WT)') +
  ylab('-log10 q-value') +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))
dev.off()
# interactive plot
degs$tooltip_text = paste0(degs$symbol, ', FC=', round(degs$avg_FC,2), ', fdr=', format(degs$p_val_adj,scientific=T, digits=2))
p = ggplot(degs, aes(x=avg_log2FC, y=-log10(p_val_adj), tooltip=tooltip_text)) + 
  geom_point_interactive(aes(color=is.significant), size=.3) +
  scale_color_manual(values = c("black", "red")) +
  geom_hline(yintercept=1.3, color='red') +
  geom_vline(xintercept=-1, color='blue3') +
  geom_vline(xintercept=1, color='blue3') +
  xlab('log2 fold change (KO vs WT)') +
  ylab('-log10 q-value') +
  theme_bw() +
  theme(legend.position = "none", #aspect.ratio=1/1.5, 
        text=element_text(size=6))
my_widget = girafe(ggobj=p, width_svg=3, height_svg=2,
                   options = list(opts_sizing(width=.7), opts_zoom(max=5)))
htmlwidgets::saveWidget(my_widget, snakemake@output[['interactive_volcano_plot']], selfcontained=TRUE)
print('interactive volcano plot done')


# generate heatmap of DE ISGs
my_genes  = isgs$Gene.Name
norm_expr = as.matrix(sub_sampled[['RNA']]@data)
my_genes  = intersect(my_genes, row.names(norm_expr))
norm_expr = norm_expr[my_genes, ]
p = sc_heatmap(norm_expr, sub_sampled@meta.data, cluster_by='sample', fontsize=9, use_yellow_purple = T)
ggsave(snakemake@output[['heatmap_DE_ISGs']], p, height=4, width=5)



# generate tiles for genes of interest
row.names(genes_df) = genes_df$gene
genes_df   = genes_df[-1, ] # remove IRF2
categories = unique(genes_df$category)
abbr_categ = c('TF', 'IR', 'EM', 'ISGs', 'Misc')

degs = read.csv(snakemake@input[['degs_without_cutoff']])
degs = degs[degs$cluster_1==paste0(conditions[2],'.All'), ]

norm_expr_WT_mean = log2(rowMeans(subset(seur, sample=='Tum_WT')[['RNA']]@counts)+1)
norm_expr_KO_mean = log2(rowMeans(subset(seur, sample=='Tum_KO')[['RNA']]@counts)+1)
norm_expr_mean = cbind(Tum_WT=norm_expr_WT_mean, Tum_KO=norm_expr_KO_mean)
norm_expr_mean = norm_expr_mean[genes_df$gene, ]
fdr            = sapply(row.names(norm_expr_mean), function(g) min(degs[degs$symbol==g, 'p_val_adj']))
labels = sapply(fdr, function(x) ifelse(x>0.05, '', ifelse(x<0.001, '***', ifelse(x<0.01, '**', '*')))) 
colfunc <- colorRampPalette(c("white", "firebrick2"))
p_list = NULL
for (i in 1:length(categories)) {
  these_genes = genes_df[genes_df$category %in% categories[i], ]
  this_norm_expr = norm_expr_mean[these_genes$gene, ]
  row.names(this_norm_expr) = paste(row.names(this_norm_expr), labels[these_genes$gene])
  annotation_row = genes_df[genes_df$category==categories[i], 'category', drop=F]
  p = pheatmap(this_norm_expr, 
           color = colfunc(100),
           cluster_rows = F, 
           cluster_cols = F,
           show_rownames = T,
           show_colnames=T,
           scale='none',
           main=abbr_categ[i],
           angle_col=45)
  p_list[[i]] = p[[4]]
}

g = grid.arrange(arrangeGrob(grobs=p_list, nrow=1))
ggsave(snakemake@output[['selected_genes_tiles']], g, width=8, height=5)
