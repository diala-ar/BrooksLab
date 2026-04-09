#### Create seurat object ####
create_RNA_seurat = function(dir_name, specie='mouse', project_name, sample_id) {
  if (specie=='mouse') {
    mito_pattern = '^mt'
    ribo_pattern = '^Rp[sl]'
    hmgb_patter  = '^Hb[^(eps)]'
  } else if (specie=='human') {
    mito_pattern = '^MT'
    ribo_pattern = '^RP[SL]'
    hmgb_patter  = '^HB[^(EPS)]'
  }
  ncounts  = Read10X(dir_name)
  colnames(ncounts) = sub('-1', '', colnames(ncounts))
  seur = CreateSeuratObject(counts=ncounts, project=project_name, min.cells=10, min.features=200)
  seur$mito_prct = PercentageFeatureSet(object=seur, pattern=mito_pattern)
  seur$ribo_prct = PercentageFeatureSet(object=seur, pattern=ribo_pattern)
  seur$hmgb_prct = PercentageFeatureSet(object=seur, pattern=hmgb_patter) # hymoglobin percentage indicates blood contamination
  seur$log10GenesPerUMI = log10(seur$nFeature_RNA) / log10(seur$nCount_RNA)
  seur = RenameCells(seur, add.cell.id=sample_id)
  
  seur$cells <- Cells(seur)
  seur$sample = sample_id
  
  seur
}



keep_cells = function(seur, nmads_count=3, nmads_feature=3, nmads_mito=5) {
  seur$nCount_is_outlier = isOutlier(seur$nCount_RNA, nmads=nmads_count, log=T, type='both')
  seur$nFeature_is_outlier = isOutlier(seur$nFeature_RNA, nmads=nmads_feature, log=T, type='both')
  seur$mito_is_outlier = isOutlier(seur$mito_prct, nmads=nmads_mito, log=F, type='higher')
  Cells(subset(seur, nCount_is_outlier==F & nFeature_is_outlier==F & mito_is_outlier==F & 
               ribo_prct>5 & hmgb_prct<0.00025 & log10GenesPerUMI>0.8))
}



#### Get % of remaining cells and genes after filtering them out ####
prct_filtered_cells_n_genes = function(seur, filt_seur, sample_id) {
  paste0("Filtered ", sample_id, ': ', ncol(filt_seur), ' out of ', ncol(seur), ' cells (', 
         round(ncol(filt_seur)/ncol(seur),2)*100, '%); ',
         nrow(filt_seur), ' out of ', nrow(seur), ' genes (', 
         round(nrow(filt_seur)/nrow(seur),2)*100, '%)')
}




#### cluster cells in a unimodal or multimodal seurat object ####
cluster_cells = function(seur, dims=1:50, reduction, graph_name, resolution, reduction_name, 
                         reduction_key, assays, clust_algorithm=1, 
                         modality_weight_name=NULL, nn_name=NULL) {
  if (length(assays) == 1) { # unimodal clustering
    seur = FindNeighbors(seur, dims=dims, reduction=reduction)
    seur = FindClusters(seur, graph.name=graph_name, resolution=resolution, algorithm=clust_algorithm)
    seur = RunUMAP(seur, dims=dims, reduction=reduction,
                   reduction.name=reduction_name, reduction.key=reduction_key)
    
    
  } else if (length(assays > 1)) { # multimodal clustering using the nearest neighbor analysis
    seur = FindMultiModalNeighbors(seur, reduction.list=reduction, 
                                   dims.list=dims, modality.weight.name=modality_weight_name)
    seur = FindClusters(seur, graph.name=graph_name, algorithm=clust_algorithm, resolution=resol)
    seur = RunUMAP(seur, nn.name=nn_name, reduction.name=reduction_name, 
                   reduction.key=reduction_key)
  }
}



#### Annotate cells using SingleR and reference (bed.se) ####
annotate_cells = function(seur, graph_name, bed.se, resolution, assay) {
  clust_col = paste0(graph_name, '_res.', resolution)
  ## annotate cell types using singleR
  for (label_type in c('label.main', 'label.fine')) {
    immgen_anno = SingleR(test=GetAssayData(seur, assay=assay, slot='data'), 
                          ref=bed.se, labels=bed.se[[label_type]],
                          clusters=seur@meta.data[, clust_col])   
    cluster_ids = setNames(immgen_anno$labels, as.character(rownames(immgen_anno)))
    seur@meta.data[, paste0('immgen_res_', resol, '_', label_type)] = cluster_ids[as.character(seur@meta.data[, clust_col])]
  }
  
  seur
}




#### Generate UMAPs (by sample and cell type annotated), cluster markers as heatmap and excel file ####
get_plots = function(seur, sample_col, resolution=resol, out_dir, suffix, top_nb=10,
                     graph_name, reduction, widths=c(1.7, 2.5), width=8.5) {
  # UMAP, label by cluster number, split by sample
  DefaultAssay(seur) = 'RNA'
  seur@meta.data[, paste0(graph_name, '_res.', resol)] = factor(seur@meta.data[, paste0(graph_name, '_res.', resol)],
                                                                levels = sort(unique(as.numeric(as.character(seur@meta.data[, paste0(graph_name, '_res.', resol)])))))
  
  seur = SetIdent(seur, value=sample_col)
  min_sample_size = min(table(seur@meta.data[, sample_col]))
  sub_sampled     = subset(x=seur, downsample=min_sample_size) 
  sub_sampled
  
  this_idents  = paste0(graph_name, '_res.', resol)
  p1 = DimPlot(sub_sampled, reduction=reduction, group.by=this_idents, 
               label=T, repel=T, label.size=4, pt.size=1, shuffle=T) +
    theme(legend.text=element_text(size=7)) +
    guides(color=guide_legend(override.aes=list(size=2), keywidth=.2, keyheight=.6, ncol=1)) +
    ggtitle(paste0('Resol: ', resol, ', cluster #'))
  p1$layers[[1]]$aes_params$alpha = 0.6
  p2 = DimPlot(sub_sampled, reduction=reduction, split.by=sample_col, group.by=this_idents, 
               label=T, repel=T, label.size=4, pt.size=1, shuffle=T) + NoLegend() + ggtitle(NULL)  
  p2$layers[[1]]$aes_params$alpha = 0.6
  umap_split_by_sample = ggarrange(p1, p2, ncol=2, nrow=1, widths=widths)
  ggsave(file.path(out_dir, paste0(suffix, '_resol_', resolution,'_UMAP_split_by_sample.pdf')),
         plot=umap_split_by_sample, width=width, height=3.5)
  print('***umap_all_resot_split_by_sample done***')
  
  # UMAP, label by cluster number and cell type annotations
  p1 = DimPlot(seur, reduction=reduction, group.by=this_idents, 
               label=T, repel=T, label.size=4, pt.size=1, shuffle=T) +
    theme(legend.position='bottom', legend.direction='horizontal', legend.text=element_text(size=7)) +
    guides(color=guide_legend(override.aes=list(size=2), keywidth=.2, keyheight=.6) ) + 
    ggtitle(paste0('Resol: ', resol, ', cluster #')) 
  p1$layers[[1]]$aes_params$alpha = 0.6
  p2 = DimPlot(seur, reduction=reduction, label=T, repel=T, 
               group.by=paste0('immgen_res_', resol, '_label.main'), 
               label.size=3, pt.size=1, shuffle=T) +
    theme(legend.position='bottom', legend.direction='horizontal', legend.text=element_text(size=7)) +
    guides(color=guide_legend(override.aes=list(size=2), keywidth=.2, keyheight=.6) ) + 
    ggtitle('Main label')   
  p2$layers[[1]]$aes_params$alpha = 0.6
  width_scale = 6
  p3 = DimPlot(seur, reduction=reduction, label=T, repel=T, 
               group.by=paste0('immgen_res_', resol, '_label.fine'), 
               label.size=2, pt.size=1, shuffle=T) +
    theme(legend.position='bottom', legend.direction='horizontal', legend.text=element_text(size=5)) +
    guides(color=guide_legend(override.aes=list(size=2), keywidth=.2, keyheight=.6, ncol=3) ) + 
    theme(legend.margin = margin(.1, 1, .1, .1, "cm")) +
    ggtitle('Fine label') 
  p3$layers[[1]]$aes_params$alpha = 0.6
  umap_annotated = ggarrange(p1, p2, p3, ncol=3, nrow=1, widths=c(1, 1.2, 1.5))
  ggsave(file.path(out_dir, paste0(suffix, '_resol_', resolution,'_UMAP_with_cell_type_annotations.pdf')),
         plot=umap_annotated, width=14, height=5)
  
  # get heatmap of top 20 up-regulated RNA markers
  DefaultAssay(seur) = 'RNA'
  Idents(seur) = paste0(graph_name, '_res.', resol)
  markers = FindAllMarkers(seur)
  markers$avg_FC = exp(markers$avg_log2FC)
  markers = markers[, c('gene', 'cluster', 'avg_FC', 'avg_log2FC', 'p_val_adj', 'pct.1', 'pct.2')]
  markers = data.frame(resolution=resol, markers)
  markers = markers[markers$p_val_adj<0.05, ]
  markers %>%
    group_by(cluster) %>%
    top_n(n=top_nb, wt=avg_log2FC) -> top_markers
  markers_heatmap = DoHeatmap(seur, features=unique(top_markers$gene), assay="RNA", angle=90, size=3) +
    theme(text = element_text(size = 9)) + #NoLegend() + 
    ggtitle(paste0('Res: ', resol, ', Top 10 up-regulated RNA markers'))
  ggsave(file.path(out_dir, paste0(suffix, '_resol_', resolution,'_Heatmap_top_', top_nb, '_up_regulated_markers.pdf')),
         plot=markers_heatmap, width=9, height=11)
  
  markers = markers[order(markers$cluster, markers$avg_log2FC, decreasing=c(F, T)), ]
  markers_ls = split(markers, f=as.factor(markers$cluster))
  names(markers_ls) = paste0('C', names(markers_ls))
  write.xlsx(markers_ls, 
             file=file.path(out_dir, paste0(suffix, '_resol_', resolution,'_RNA_markers.xlsx')), 
             rowNames=F)
  
  top_markers
}



get_DEGs = function(seur, ident_1, ident_2, group_by=NULL, subset_ident=NULL, 
                    significant, degs_by='sample', idents_col) {
  # This function generates DEGs by cluster (degs_by='cluster') or by sample (degs_by='sample') see example 2
  # for by cluster, set ident_1 to the sample of interest, group_by to the column having the sample of interes,
  #                 subset_ident to the cluster of interest, degs_by to 'cluster', 
  #                 idents_col to the column having the clusters
  # for sample, set ident_1 to the first sample of interest, ident_2 to the second sample of interes,
  #                 group_by to NULL, subset_ident to NULL, degs_by to 'sample', idents_col to NULL
  # significant = F returns all DEGs, significant = T returns significant DEGs only (those with adj_pval<0.05)
  DefaultAssay(seur) = 'RNA'
  seur = SetIdent(seur, value=idents_col)
  if (significant==F) {
    logfc_threshold=0; min_pct=0
  } else {
    logfc_threshold=0.25; min_pct=0.1
  }
  if (degs_by=='cluster') {
    DEGs = FindMarkers(seur, ident.1=ident_1, group.by=group_by, subset.ident=subset_ident,
                       logfc.threshold=logfc_threshold, min.pct=min_pct)
    DEGs$cluster = subset_ident
  } else if (degs_by=='sample') {
    DEGs = FindMarkers(seur, ident.1=ident_1, ident.2=ident_2,
                       logfc.threshold=logfc_threshold, min.pct=min_pct)
    DEGs$cluster = 'Bulk'
  }
  DEGs$gene    = row.names(DEGs)
  DEGs$avg_FC  = 2^(DEGs$avg_log2FC)
  DEGs = DEGs[, c('gene', 'avg_FC', 'avg_log2FC', 'p_val_adj', 'pct.1', 'pct.2', 'cluster', 'p_val')]
  if (significant==T) {
    DEGs = DEGs[DEGs$p_val_adj<0.05, ]
    print('filter')
  }
  DEGs = DEGs[order(DEGs$avg_log2FC, decreasing=T), ]
}




get_cells_prct = function(seur, genes, rna_count_cutoff) {
  UMIs_count  = GetAssayData(seur, assay='RNA', slot='counts')[genes, ]
  cells_count = apply(UMIs_count, 1, function(x) sum(x >= rna_count_cutoff))
  cells_prct  = round(cells_count / ncol(seur) * 100, 2)
}


chi_square_prct_test = function(seurs_i, tcr, rna_count_cutoff, clust) {
  if (clust != 'Bulk') {
    seurs_i = subset(seurs_i, seurat_clusters==clust)
  }
  UMIs_count = GetAssayData(seurs_i, assay='RNA', slot='counts')[tcr, ]
  seurs_i$express_tcr = (UMIs_count >= rna_count_cutoff)
  
  conting_tbl = table(seurs_i$sample, seurs_i$express_tcr)
  correct     = any(conting_tbl < 5)
  pval = chisq.test(seurs_i$sample, seurs_i$express_tcr, correct=correct)$p.value
}




#############################################
#### visualize cluster proportions
## ex.: group_col = 'sample', clust_col='integrated_snn_res.0.2'
viz_clust_prop = function(seur, group_col, clust_col) {
  # plot cluster percentage per samples
  seurs_ls = SplitObject(seur, split.by=group_col)
  df_ls    = lapply(seurs_ls, function(seu) seu@meta.data)
  cells_percent = do.call(rbind, lapply(df_ls, function(df) {
    temp = data.frame(table(cluster=df[, clust_col]), group=unique(df[, group_col]))
    names(temp)[2] = 'cells_nb'
    total = temp %>%
      group_by(group) %>%
      summarize(total_cells=sum(cells_nb))
    temp = merge(temp, total)
    temp$percentage = round(temp$cells_nb / temp$total_cells * 100, 1)
    temp$cluster = factor(temp$cluster, levels=sort(temp$cluster))
    
    temp
  }))
  p2 = ggplot(cells_percent, aes(x=group, y=percentage, fill=cluster, group=group, label=round(percentage))) +
    geom_bar(stat="identity", color='black', linewidth=.5, width=.3) + 
    geom_text(size=4, position=position_stack(vjust = 0.5)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=11)) +
    #scale_fill_brewer(palette = "Set3") +
    xlab('') + ylab('Cells percentage per cluster')
  
  return(list(df=cells_percent, plot=p2))
}





custom_pheatmap_colors = function(min_data, max_data, use_yellow_purple=F) {
  paletteLength <- 200
  if (use_yellow_purple==F) {
    myColors = colorRampPalette(c("#4675B3", "#74A9CC", "#ABD0E4", "#E5F5ED", "white", '#FEF5B0', '#FFDA8B', '#FF935F', '#DB362D'))(paletteLength)
  } else {
    myColors = PurpleAndYellow(paletteLength)
  }
  # use floor and ceiling to deal with even/odd length pallettelengths
  myBreaks <- c(seq(min_data, 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max_data/paletteLength, max_data, length.out=floor(paletteLength/2)))
  return(list(breaks=myBreaks, my_colors=myColors))
}




get_logFC_n_FDR = function(seur, degs_file, genes=NULL, significant_rows=T,
                                ident_1, ident_2) {
  # seur: seurat object
  # degs_file: (string) the file name that contains the DEGs including the full path
  # genes: (vector) if NULL all DEGs will be visualized, if genes is a vector of genes symbol, then those genes will be visaulized
  # significant_rows: (logical) if T only genes with at least 1 significant p-value will be visualized, if F all genes will be visualized including none-significant ones.
  # ident_1: (string) the condition that is in the numerator of the fold change (this value should be in a column called sample in meta.data of seur)
  # ident_2: (string) the condition that is in the denominator of the fold change (this value should be in a column called sample in meta.data of seur)
  
  clusters = getSheetNames(degs_file)
  fdr_ls = lapply(clusters, function(clust) {
    fdr = read.xlsx(degs_file, sheet=clust)
    fdr = fdr[, c('gene', 'p_val_adj')]
    names(fdr)[2] = clust
    if (!is.null(genes)) { # heatmap is for specific genes
      missing_genes = setdiff(genes, fdr$gene)
      fdr_missing_genes = data.frame(gene=missing_genes, pval=1)
      names(fdr_missing_genes)[2] = clust
      fdr = rbind(fdr[fdr$gene %in% genes, ], fdr_missing_genes)
    }
    return(fdr)
  })
  fdr_df = Reduce(function(x, y) merge(x, y, all=TRUE), fdr_ls) 
  row.names(fdr_df) = fdr_df$gene
  fdr_df = fdr_df[, -1]
  if (significant_rows) { # remove all genes that are not significant in any cluster
    fdr_df = fdr_df[rowSums(fdr_df) < ncol(fdr_df), ] 
  }
  fdr_df = fdr_df[order(row.names(fdr_df)), ]
  fdr_df[is.na(fdr_df)] = 1
  
  
  clusters2 = sort(as.numeric(as.character(unique(seur$seurat_clusters))))
  fc_df = do.call(cbind, lapply(clusters2, function(clust) {
    cells_1 = Cells(subset(seur, seurat_clusters==clust & sample==ident_1))
    cells_2 = Cells(subset(seur, seurat_clusters==clust & sample==ident_2))
    fc = FoldChange(seur, ident.1=cells_1, ident.2=cells_2, 
                    features=row.names(fdr_df))
    fc = data.frame(avg_log2FC=fc$avg_log2FC)
  }))
  # add bulk
  cells_1 = Cells(subset(seur, sample==ident_1))
  cells_2 = Cells(subset(seur, sample==ident_2))
  fc = FoldChange(seur, ident.1=cells_1, ident.2=cells_2, 
                  features=row.names(fdr_df))
  fc = data.frame(Bulk=fc$avg_log2FC)
  fc_df = cbind(fc_df, fc)
  row.names(fc_df) = row.names(fdr_df)
  names(fc_df) = clusters

  return(list(fdr_df = fdr_df, fc_df = fc_df))
}



logFC_n_FDR_htmp = function(htp_data, fontsize=10, use_yellow_purple=F, abs_logFC_cutoff=0, 
                            cluster_rows=T, title, angle_col=0) {
  library(tidyr)
  # plot heatmap on cluster level
  
  # prepare data for heatmap, transform from long to wide format
  fc_df = htp_data$fc_df
  fdr_df = htp_data$fdr_df
  stars_df = matrix('', nrow=nrow(fdr_df), ncol=ncol(fdr_df))
  for (i in 1:nrow(fdr_df)) {
    for (j in 1:ncol(fdr_df)) {
      stars_df[i, j] = ifelse(fdr_df[i, j] < 0.001, '***',
                                ifelse(fdr_df[i, j] < 0.01 & fdr_df[i, j] > 0.001, '**',
                                       ifelse(fdr_df[i, j] < 0.05 & fdr_df[i, j] > 0.01, '*', '')))
    }
  }

  # winsorize for color contrast
  if (abs_logFC_cutoff>0) {
    fc_df[fc_df > abs_logFC_cutoff] = abs_logFC_cutoff
    fc_df[fc_df < -abs_logFC_cutoff] = -abs_logFC_cutoff
  }
  
  color_n_breaks = custom_pheatmap_colors(min(fc_df), max(fc_df), use_yellow_purple=F)
  p = pheatmap::pheatmap(as.matrix(fc_df), 
                         color  = color_n_breaks$my_colors,
                         breaks = unique(color_n_breaks$breaks),
                         cluster_rows = cluster_rows,
                         cluster_cols = F,
                         show_rownames = T,
                         show_colnames = T,
                         scale='none',
                         display_numbers=as.matrix(stars_df),
                         angle_col=angle_col,
                         fontsize=fontsize,
                         main=title)
  
  p
}


sc_heatmap = function(seur, cluster_by, use_yellow_purple=F, title='', features,
                      min_cutoff=NULL, max_cutoff=NULL, fontsize=10, clust_rows=T) {
  ## htp_data (data.frame or matrix): the expression data from the data assay, meta_data
  ## meta_data (data.frame): seur@meta.data containing the column to cluster cells by
  ## cluster_by (character): name of the column in meta_data to cluster cells with
  ## min_cutoff and max_cutoff (numeric): the min and max cutoffs to winsorize expression with.
  meta_data = seur@meta.data
  htp_data  = GetAssayData(seur, assay='RNA', slot='data')
  htp_data  = t(scale(t(as.matrix(htp_data))))
  htp_data  = htp_data[features, ]
  library(dplyr)
  meta_data$cluster_by = meta_data[, cluster_by]
  
  # order meta_data by grp to order then the htp_data columns
  meta_data  = meta_data[order(meta_data$cluster_by), ] 
  htp_data   = as.matrix(htp_data[, row.names(meta_data)])
  if (!is.null(min_cutoff))
    htp_data[htp_data < min_cutoff] = min_cutoff
  if (!is.null(max_cutoff)) 
    htp_data[htp_data > max_cutoff] = max_cutoff
  
  
  color_n_breaks = custom_pheatmap_colors(min(htp_data), max(htp_data), use_yellow_purple)
  gaps_nb = as.array(table(meta_data[cluster_by[1]]))
  gaps_cols_index = sapply(1:length(gaps_nb), function(i) sum(gaps_nb[1:i]))
  p = pheatmap::pheatmap(as.matrix(htp_data), 
                         #color = PurpleAndYellow(20), #RColorBrewer::brewer.pal(9, "Blues"), #colorRampPalette(c("black", "yellow", "blue"))(50),
                         color  = color_n_breaks$my_colors,
                         breaks = unique(color_n_breaks$breaks),
                         border_color = NA,
                         cluster_rows = clust_rows, 
                         cluster_cols = F,
                         show_rownames = T,
                         show_colnames = F,
                         scale='none',
                         annotation_col = select(meta_data, all_of(rev(cluster_by))),
                         gaps_col = gaps_cols_index,
                         main = title,
                         fontsize = fontsize)
  p
}



mean_htp = function(seur, cluster_by, features, fontsize, title) {
  library(viridis)
  meta_data = seur@meta.data
  expr_mat  = seur[['RNA']]@data
  clusters = sort(as.numeric(as.character(unique(meta_data[, cluster_by]))))
  mean_df = do.call(cbind, lapply(clusters, function(clust) {
    cells_1 = row.names(meta_data[meta_data$sample==paste0('ON_', group_i) & meta_data[, cluster_by]==clust, ])
    cells_2 = row.names(meta_data[meta_data$sample==paste0('OFF_', group_i) & meta_data[, cluster_by]==clust, ])
    mn_df = data.frame(log2(rowMeans(expm1(expr_mat[, cells_1])) + 1),
                       log2(rowMeans(expm1(expr_mat[, cells_2])) + 1))
    names(mn_df) = paste0(clust, c('_ON', '_OFF'))
    mn_df
  }))
  mean_df = t(scale(t(as.matrix(mean_df))))
  htp_data = mean_df[features, ]
  
  gaps_nb = length(clusters) - 1
  gaps_cols_index = seq(2, (gaps_nb)*2, 2) 
  p = pheatmap::pheatmap(as.matrix(htp_data), 
                         color  = inferno(100),
                         border_color = 'black',
                         cluster_rows = T, 
                         cluster_cols = F,
                         show_rownames = T,
                         show_colnames = T,
                         scale='none',
                         gaps_col = gaps_cols_index,
                         main = title,
                         fontsize = fontsize)
  p
  
}




featureplots_faceted = function(seurs, split_by='sample', features, 
                                ident_1, ident_2, sep, alpha) {
  library(ggrastr)
  library(ggplot2)
  Idents(seurs) = split_by
  min_sample_size = min(table(seurs$sample))
  sub_sampled = subset(x=seurs, downsample=min_sample_size) 
  sub_sampled
  
  expr_mat  = GetAssayData(sub_sampled, assay='RNA', slot='data')
  umap_coords  = sub_sampled@reductions$umap@cell.embeddings
  umap_coords  = data.frame(cell=row.names(umap_coords), umap_coords)
  p_ls = lapply(1:length(features), function(i) {
    expr_feat = expr_mat[features[i], ]
    expr_feat = data.frame(cell=names(expr_feat), expr=expr_feat)
    temp_1 = temp_2 = merge(umap_coords, expr_feat, x.all=T)
    temp_1$cell  = features[i]
    temp_2$cell  = gsub(paste0('(.+)', sep, '(.+)'), '\\1', temp_2$cell)
    plot_data    = rbind(temp_1, temp_2)
    plot_data$cell = factor(plot_data$cell, levels=c(features[i], ident_1, ident_2))
    q_expr = quantile(temp_1$expr, probs=0.99)
    plot_data[plot_data$expr >= q_expr, 'expr'] = q_expr
    plot_data = plot_data[order(plot_data$expr), ]
    if (q_expr > 0){
      p = ggplot(plot_data, aes(x=UMAP_1, y=UMAP_2)) + theme_bw() +
        ggrastr::rasterise(geom_point(aes(color=expr), size=10^-12, alpha=alpha), scale=.4) +
        facet_wrap(~cell) +
        scale_color_gradient(low='grey85', high='blue') +
        theme(legend.key.size = unit(.3, 'cm'),
              legend.text = element_text(size=5)) 
    }
  })
}




gsea_barplot = function(data) {
  library(scales)
  data$score = -log10(data$pval) * sign(data$NES)
  if ( min(data$score) < 0 & max(data$score) > 0 ) {
    values  = c(min(data$score), 0, max(data$score))
    colours = c('steelblue3', 'white', 'firebrick2')
    data = data[order(data$NES), ]
  } else if (min(data$score) > 0) {
    values  = c( 0, max(data$score))
    colours = c('white', 'firebrick2')
    data = data[order(data$NES), ]
  } else if (max(data$score) < 0 ) {
    values  = c(min(data$score), 0)
    colours = c('steelblue3', 'white')
    data = data[order(data$NES, decreasing=T), ]
  }
  data$pathway = sub('HALLMARK_', '', data$pathway)
  data$pathway = factor(data$pathway, levels=data$pathway)
  
  p = ggplot(data, aes(pathway, NES, fill=score)) +
    geom_bar(stat='identity', width=.6, colour='black', size=.2) +
    coord_flip() +
    scale_fill_gradientn(colours=colours, values=rescale(values), 
                         limits=values[c(1, length(values))], name='-log10(pval) * sign(NES)') +
    theme_bw() + xlab('') +
    theme(legend.direction = 'horizontal', 
          legend.position = 'bottom',
          legend.key.height = unit(.3, 'cm'), 
          legend.key.width = unit(.5, 'cm'), 
          legend.text=element_text(size=6),
          legend.title = element_text(size=7)) +
    ggtitle(group_i)
}



gsea_dotplot = function(file_name, title) {
  library(dplyr)
  library(scales)
  
  # read in fgsea results from each excel sheet
  clusters  = getSheetNames(file_name)
  all_data = do.call(rbind, lapply(clusters, function(clust)  {
    df = read.xlsx(file_name, sheet=clust)
    if (nrow(df) > 0) df$cluster = clust
    df
  }))
  all_data$pathway = sub('HALLMARK_', '', all_data$pathway)
  
  # calculate gene ratio and score
  all_data$gene_count = sapply(all_data$leadingEdge, function(x)
    length(unlist(strsplit(x, ', '))))
  all_data$gene_ratio = all_data$gene_count / all_data$size
  all_data$score = -log10(all_data$pval) * sign(all_data$NES)
  
  # prepare for color mapping
  if ( min(all_data$score) < 0 & max(all_data$score) > 0 ) {
    values  = c(min(all_data$score), 0, max(all_data$score))
    colours = c('steelblue3', 'white', 'tomato1')
  } else if (min(all_data$score) > 0) {
    values  = c( 0, max(all_data$score))
    colours = c('white', 'tomato1')
  } else if (max(data$score) < 0 ) {
    values  = c(min(all_data$score), 0)
    colours = c('steelblue3', 'white')
  }
  
  # order columns by clusters and rows by pathways
  pathways_order = all_data %>% 
    group_by(pathway) %>% 
    summarise(mean_score=mean(score)) %>% 
    arrange(mean_score) %>% 
    select(pathway) %>% 
    unlist()
  all_data$pathway = factor(all_data$pathway, levels=pathways_order)
  all_data$cluster = factor(all_data$cluster, levels=unique(all_data$cluster))
  
  ggplot(all_data, aes(cluster, pathway)) + 
    geom_point(aes(size = gene_ratio, color = score)) +
    theme_bw(base_size = 14) +
    scale_color_gradientn(colours=colours, values=rescale(values), 
                          limits=range(values), 
                          name='-log10(pval) * sign(NES)') +
    ylab(NULL) +
    ggtitle(title)
}