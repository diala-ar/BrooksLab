library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)


seur      = readRDS(snakemake@input[['seur']]) #to browse data of CD8 clustered seur=readRDS("results/Seurat_AB_seq/Tum_CD8_v2/Seur_Tum_CD8_clustered.rds")
pca_dims  = snakemake@params[['pca_dims']]     #pca_dims = list(rna=20, adt=12)
resols    = snakemake@params[['resolutions']]  #resols = seq(0.1, 1.5, .2)
recluster = snakemake@params[['recluster']]    

#############################################
# renormalize, rerun pca and umap
#############################################
if (recluster==T) {
  seur = NormalizeData(seur) %>%
    FindVariableFeatures() %>%
    ScaleData(features=row.names(seur)) %>%
    RunPCA()
  
  # Normalize AB-seq data
  if ('ADT' %in% Assays(seur)) {
    seur <-  NormalizeData(seur,
                           assay="ADT",
                           normalization.method = "CLR") %>%  
      ScaleData(assay="ADT") %>% 
      RunPCA(reduction.name='apca')
  }
}

#############################################
# (re)cluster
#############################################
if ('ADT' %in% Assays(seur)) {
  seur = FindMultiModalNeighbors(seur, 
                                 reduction.list = list("pca", "apca"), 
                                 dims.list = list(1:pca_dims$rna, 1:pca_dims$adt), 
                                 modality.weight.name = c("RNA.weight"))
  seur = FindClusters(seur, graph.name="wsnn", algorithm=3, resolution=resols)
  seur = RunUMAP(seur, nn.name="weighted.nn", reduction.name="wnn.umap", reduction.key="wnnUMAP_")
} else {
  seur = FindNeighbors(seur, dims=1:pca_dims$rna)
  seur = FindClusters(seur, resolution=resols, graph.name="RNA_snn")
  seur = RunUMAP(seur, dims=1:pca_dims$rna, reduction.name="umap", reduction.key="UMAP_")
}
saveRDS(seur, file=snakemake@output[['seur']])


#############################################
# generate all umaps
#############################################
# UMAP all resolutions
p = NULL
this_reduction = ifelse('ADT' %in% Assays(seur), 'wnn.umap', 'umap')
for (i in 1:length(resols)) {
  resol = resols[i]
  this_idents = ifelse('ADT' %in% Assays(seur), paste0('wsnn_res.', resol), paste0('RNA_snn_res.', resol))
  this_title  = ifelse('ADT' %in% Assays(seur), paste0('WNN, Res: ', resol), paste0('SNN, Res: ', resol))
  Idents(seur) = this_idents
  p[[i]] = DimPlot(seur, reduction=this_reduction, label=T, repel=T, label.size=5.5, pt.size=2.5, shuffle=T) + 
    ggtitle(this_title) +
    theme(legend.text=element_text(size=11)) +
    guides(color=guide_legend(override.aes=list(size=3), keywidth=.6, keyheight=.9, ncol=1) )
  p[[i]]$layers[[1]]$aes_params$alpha = 0.6
}

x = length(p)
cols = round(sqrt(x),0)
rows = ceiling(x/cols)
pdf(paste0(snakemake@output[['umap_all_resol']]), width=20, height=17)
ggarrange(plotlist=p, ncol=cols, nrow=rows)
dev.off()
print('***umap_all_resol done***')



# UMAP all resolutions split by sample
# downsample cells to have the same nb of cells in WT and KO
Idents(seur) = 'sample'
min_sample_size = min(table(seur$sample))
sub_sampled = subset(x=seur, downsample=min_sample_size) 
sub_sampled

pdf(snakemake@output[['umap_all_resol_split_by_sample']], width=8.5, height=4)
for (i in 1:length(resols)) {
  resol = resols[i]
  this_idents = ifelse('ADT' %in% Assays(seur), paste0('wsnn_res.', resol), paste0('RNA_snn_res.', resol))
  p1 = DimPlot(sub_sampled, reduction=this_reduction, group.by=this_idents, 
               label=T, repel=T, label.size=3.5, pt.size=2.5, shuffle=T) +
               theme(legend.text=element_text(size=7)) +
               guides(color=guide_legend(override.aes=list(size=2), keywidth=.2, keyheight=.6, ncol=1) )
  p1$layers[[1]]$aes_params$alpha = 0.6
  p2 = DimPlot(sub_sampled, reduction=this_reduction, split.by='sample', group.by=this_idents, 
               label=T, repel=T, label.size=3.5, pt.size=2.5, shuffle=T) + NoLegend() + ggtitle(NULL)  
  p2$layers[[1]]$aes_params$alpha = 0.6
  print(ggarrange(p1, p2, ncol=2, nrow=1, widths=c(1.7, 2.5)))
}
dev.off()
print('***umap_all_resot_split_by_sample done***')


markers     = read.csv(snakemake@input[['markers']])
rna_markers = markers$RNA[markers$RNA!='']
if ('ADT' %in% Assays(seur)) {
  adt_markers = paste0(markers$ADT[markers$ADT!=''], '-adt')
  markers     = sort(c(adt_markers, rna_markers))
} else {
  markers = sort(rna_markers)
}

# plot adt and rna UMAP markers split by sample
pdf(snakemake@output[['umap_markers']], width=8.5, height=4)
all_cells = Cells(sub_sampled)
wt_cells  = grep('Tum_WT', all_cells, value=T)
ko_cells  = grep('Tum_KO', all_cells, value=T)
expr_markers = rbind(GetAssayData(seur, assay='ADT', slot='data'),
                     GetAssayData(seur, assay='RNA', slot='data'))
umap_coords  = sub_sampled@reductions$wnn.umap@cell.embeddings
umap_coords  = data.frame(cell=row.names(umap_coords), umap_coords)
for (i in 1:length(markers)) {
  p = NULL
  expr_marker = expr_markers[markers[i], ]
  expr_marker = data.frame(cell=names(expr_marker), expr_mark=expr_marker)
  temp_1 = temp_2 = merge(umap_coords, expr_marker, x.all=T)
  temp_1$cell  = markers[i]
  temp_2$cell  = gsub('(.+)-(.+)', '\\1', temp_2$cell)
  plot_data    = rbind(temp_1, temp_2)
  plot_data$cell = factor(plot_data$cell, levels=c(markers[i], 'Tum_WT', 'Tum_KO'))
  q1_expr_mark = quantile(temp_1$expr_mark, 0.05)
  q2_expr_mark = quantile(temp_1$expr_mark, probs=0.95)
  plot_data[plot_data$expr_mark<=q1_expr_mark, 'expr_mark'] = q1_expr_mark
  plot_data[plot_data$expr_mark>=q2_expr_mark, 'expr_mark'] = q2_expr_mark
  
  if (length(unique(plot_data$expr_mark))==1 && unique(plot_data$expr_mark)==q1_expr_mark){
    p = ggplot(plot_data, aes(x=wnnUMAP_1, y=wnnUMAP_2)) + theme_bw() +
      geom_point(aes(color=expr_mark), alpha=.5, size=2.5) +
      facet_wrap(~cell) +
      scale_color_gradient(low='grey85', high='grey85')
  } else {
    p = ggplot(plot_data, aes(x=wnnUMAP_1, y=wnnUMAP_2)) + theme_bw() +
      geom_point(aes(color=expr_mark), alpha=.5, size=2.5) +
      facet_wrap(~cell) +
      scale_color_gradient(low='grey85', high='blue')
  }
  print(p)
}
dev.off()

print('***umap_markers done***')



#############################################
# Identify markers for each cluster
#############################################
# Find protein markers for all clusters, and draw a heatmap
if ('ADT' %in% Assays(seur)) {
  all_adt_markers = NULL
  pdf(snakemake@output[['adt_heatmap']], height=2.5)
  for (i in 1:length(resols)) {
    resol = resols[i]
    Idents(seur) = paste0('wsnn_res.', resol)
    adt_markers  = FindAllMarkers(seur, assay="ADT")
    adt_markers$avg_FC = (2^adt_markers$avg_log2FC)
    adt_markers = adt_markers[, c('gene', 'cluster', 'avg_FC', 'avg_log2FC', 'p_val_adj', 'pct.1', 'pct.2')]
    p = DoHeatmap(seur, features=unique(adt_markers$gene), assay="ADT", angle=90, size=4) + NoLegend() + 
      ggtitle(paste0('Res: ', resol, ', ADT markers'))
    print(p)
    
    adt_markers = data.frame(resolution=resol, adt_markers)
    all_adt_markers = rbind(all_adt_markers, adt_markers)
  }
  dev.off()
  write.csv(all_adt_markers, file=snakemake@output[['adt_markers']], row.names=F)
}


# Find RNA markers for all clusters, and draw a heatmap of top 20/10 and botttom 20/10 RNA markers per cluster
rna_markers = NULL
for (i in 1:length(resols)) {
  resol = resols[i]
  Idents(seur) = ifelse('ADT' %in% Assays(seur), paste0('wsnn_res.', resol), paste0('RNA_snn_res.', resol))
  rna_markers[[i]]  = FindAllMarkers(seur, assay="RNA")
  rna_markers[[i]]$avg_FC = (2^rna_markers[[i]]$avg_log2FC)
  rna_markers[[i]] = rna_markers[[i]][, c('gene', 'cluster', 'avg_FC', 'avg_log2FC', 'p_val_adj', 'pct.1', 'pct.2')]
  rna_markers[[i]] = data.frame(resolution=resol, rna_markers[[i]])
}

top_nbs     = c(10, 20)
pdf_heights = c(25, 45)
for (j in 1:length(top_nbs)) {
  top_nb = top_nbs[j]
  pdf(snakemake@output[[paste0('rna_top_', top_nb, '_heatmap')]], height=pdf_heights[j])
  for (i in 1:length(resols)) {
    resol = resols[i]
    Idents(seur) = ifelse('ADT' %in% Assays(seur), paste0('wsnn_res.', resol), paste0('RNA_snn_res.', resol))
    rna_markers[[i]] %>%
      group_by(cluster) %>%
      top_n(n=top_nb, wt=avg_log2FC) -> top_markers
    p1 = DoHeatmap(seur, features=unique(top_markers$gene), assay="RNA", angle=90) +
      theme(text = element_text(size = 13)) + NoLegend() + 
      ggtitle(paste0('Res: ', resol, ', Top ', top_nb, ' up-regulated RNA markers'))
    print(p1)
    
    rna_markers[[i]] %>%
      group_by(cluster) %>%
      top_n(n=-top_nb, wt=avg_log2FC) -> bottom_markers
    p2 = DoHeatmap(seur, features=unique(bottom_markers$gene), assay="RNA", angle=90) +
      theme(text = element_text(size = 13)) + NoLegend() +
      ggtitle(paste0('Res: ', resol, ', Top ', top_nb, ' down-regulated RNA markers'))
    print(p2)
  }
  dev.off()
}
all_rna_markers = do.call(rbind, rna_markers)
write.csv(all_rna_markers, file=snakemake@output[['rna_markers']], row.names=F)