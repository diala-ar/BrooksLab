library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(patchwork)
library(GenomicRanges)
library(future)
library(SingleR)
library(dplyr)
library(ggpubr)

set.seed(1234)



wrk_dir    = '/cluster/projects/mcgahalab/data/brookslab/sabelo/2022_03_scATAC/' # working directory
in_dir     = 'resources/scATAC_data'                                             # directory containing input files
out_dir    = 'results'
macs2_path = '/cluster/home/dabdrabb/miniconda3/envs/r413/bin/macs2'             # directory containing MACS2 path
setwd(wrk_dir)

groups = c('CD8_cre_pos', 'CD8_cre_neg', 'IRF2_KO', 'WT')
genome = seqinfo(BSgenome.Mmusculus.UCSC.mm10) # to avoid R searching for genome info on UCSC utl (no internet on server)
blacklist = blacklist_mm10
# annotations was created using the following code
# annotations = GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
# annotations = renameSeqlevels(annotations, c(paste0("chr", seqlevels(annotations))))
annotations = readRDS(file.path('resources', 'ens_Mmusculus_v79_ucsc_annotations.rds'))


###########################
##### read in grouped peaks and create a merged seurat object for all groups
peaks_ls = lapply(groups, function(group) {
   peaks = read.table(file = file.path(in_dir, group, "peaks.bed"), check.names=F, sep="\t",
                      col.names = c("chr", "start", "end"))
   gr    = makeGRangesFromDataFrame(peaks)
})
# Create a unified set of peaks to quantify in each dataset
combined_peaks = reduce(x=c(peaks_ls[[1]], peaks_ls[[2]], peaks_ls[[3]], peaks_ls[[4]]))
# Filter out bad peaks based on length
peak_widths    = width(combined_peaks)
combined_peaks = combined_peaks[peak_widths<10000 & peak_widths>20]
combined_peaks

# create seurat object
seur_ls = lapply(groups, function(group) {
   metadata = read.csv(
     file = file.path(in_dir, group, "singlecell.csv"),
     header = TRUE,
     row.names = 1
   )[-1, ]
   metadata$group = group
   print(paste(date(), 'read metadata'))

   # assess blacklist regions
   frag_filename = file.path(in_dir, group, "fragments.tsv.gz")
   frag          = import.bed(frag_filename)
   overlap       = findOverlaps(frag, blacklist)
   overlap_df    = as.data.frame(frag[overlap@from, ])
   overlap_df    = overlap_df[,c("name","score")]
   overlap_df[]  = lapply(overlap_df, function(x) type.convert(as.character(x)))
   aggr_overlap_df = aggregate(.~name, overlap_df, FUN="sum")
   row.names(aggr_overlap_df) = aggr_overlap_df$name
   blacklist_region_fragments = list()
   blacklist_frags = sapply(row.names(metadata), function(x) 
      ifelse(x %in% aggr_overlap_df$name, unlist(append(blacklist_region_fragments, aggr_overlap_df[x,'score'])), unlist(append(blacklist_region_fragments, 0) )))
   metadata$blacklist_region_fragments = unlist(blacklist_frags)
   print(paste(date(), 'identified blacklist fragments'))

   # create fragment objects
   frags = CreateFragmentObject(
      path = frag_filename,
      cells = rownames(metadata), 
      max.lines = NULL
   )
   print(paste(date(), 'created fragment object'))


   # quantify peaks in fragments
   frag_counts = FeatureMatrix(
      fragments = frags,
      features = combined_peaks,
      cells = rownames(metadata)
   )
   print(paste(date(), 'fragment counts done'))

   chrom_assay = CreateChromatinAssay(
     counts = frag_counts,
     sep = c(":", "-"),
     genome = genome,
     fragments = frags,
     min.cells = 10,
     min.features = 200
   )
   print(paste(date(), 'chromatinassay done'))

   seur = CreateSeuratObject(
     counts = chrom_assay,
     assay = "scATAC",
     meta.data = metadata
   )
   print(paste(date(), 'create seurat object done'))

   # call peaks using MACS2
   peaks_cells = CallPeaks(seur, macs2.path=macs2_path,
                           outdir=file.path(out_dir, "macs2"), 
                           name=group, cleanup=F, format='BEDPE', 
                           effective.genome.size='mm')
   
   # remove peaks on nonstandard chromosomes and in genomic blacklist regions
   peaks = keepStandardChromosomes(peaks_cells, pruning.mode="coarse")
   peaks = subsetByOverlaps(x=peaks, ranges=blacklist, invert=TRUE) #ENCODE blacklist
   print(paste(date(), 'MACS2 done'))

   # quantify counts in each peak
   macs2_counts = FeatureMatrix(fragments = Fragments(seur),
                                features = peaks,
                                cells = colnames(seur))

   print(paste(date(), 'MACS2 counts done'))
   # create a new assay using the MACS2 peak set and add it to the Seurat object
   seur[['peaks']] = CreateChromatinAssay(counts = macs2_counts,
                                          fragments = Fragments(seur),
                                          annotation = Annotation(seur))
   DefaultAssay(seur) = 'peaks'
   print(paste(date(), 'Added peaks assay done'))


   # annotate peaks: extract gene annotations from EnsDb
   Annotation(seur) = annotations       # add the gene information to the object
   #### compute QC metrics
   seur$pct_reads_in_peaks = seur$nCount_peaks / seur$passed_filters * 100
   seur$blacklist_ratio = seur$blacklist_region_fragments / seur$peak_region_fragments
   seur = NucleosomeSignal(object=seur)
   seur$nucleosome_group = ifelse(seur$nucleosome_signal>4, 'NS > 4', 'NS < 4')
   seur = TSSEnrichment(seur, fast=FALSE)
   seur$high.tss = ifelse(seur$TSS.enrichment>2, 'High', 'Low')
   
   
   seur = RenameCells(seur, 
                      new.names=gsub("-1$", "", Cells(seur)))
   saveRDS(seur, file=file.path(out_dir, "seurat_object", paste0("01_seur_", group, ".rds")))

   print(paste(group, 'Seurat object!'))
   seur
})

names(seur_ls) = groups
saveRDS(seur_ls, file.path(out_dir, 'seurat_object', '01_seurs_list_before_filtering.rds'))



################################################
#### Quality control
# print violon plots of features before filtering
pdf(file.path(out_dir, 'QC', 'QC_before_filtering.pdf'), width=13, height=4.5)
my_features = c('nCount_peaks', 'pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal', 'TSS.enrichment', 'nFeature_peaks')
for (this_group in groups) {
   print(seur_ls[[this_group]])
   p = VlnPlot(object=seur_ls[[this_group]], features=my_features, pt.size=0, ncol=6, combine=F) 
   p = lapply(1:length(p), function(i) 
      p[[i]] + labs(x=my_features[i], title=ifelse(i==1, this_group, '')) + NoLegend() +
               theme(axis.text.x=element_blank()))
   print(p[[1]] | p[[2]] | p[[3]] | p[[4]] | p[[5]] | p[[6]])
}
dev.off()

# determine criteria cutoffs based on the features violon plots
criteria_cutoffs = rbind(
   CD8_cre_pos=c(nCount_peaks_low=2000,
                 nCount_peaks_up=quantile(seur_ls[['CD8_cre_pos']]$nCount_peaks, .99), 
                 pct_reads_in_peaks=35, blacklist_ratio=0.1, nucleosome_signal=4, TSS.enrichment=2),
   CD8_cre_neg=c(nCount_peaks_low=2000, 
                 nCount_peaks_up=quantile(seur_ls[['CD8_cre_neg']]$nCount_peaks, .99), 
                 pct_reads_in_peaks=37, blacklist_ratio=0.1, nucleosome_signal=4, TSS.enrichment=2),
   IRF2_KO=c(nCount_peaks_low=800, 
             nCount_peaks_up=quantile(seur_ls[['IRF2_KO']]$nCount_peaks, .99), 
             pct_reads_in_peaks=30, blacklist_ratio=0.12, nucleosome_signal=4, TSS.enrichment=2),
   WT=c(nCount_peaks_low=800, 
        nCount_peaks_up=quantile(seur_ls[['WT']]$nCount_peaks, .99), 
        pct_reads_in_peaks=30, blacklist_ratio=0.12, nucleosome_signal=4, TSS.enrichment=2))
colnames(criteria_cutoffs)[c(2,4)] = c('nCount_peaks_up', 'blacklist_ratio')
criteria_cutoffs


# filter out low quality cells
print("Filter out cells with low quality reads using selected thresholds")
filt_seur_ls = NULL
for (i in 1:length(groups)) {
   group = groups[i]
   filt_seur_ls[[i]] = subset(x=seur_ls[[group]],
                            subset=(nCount_peaks > criteria_cutoffs[group, 'nCount_peaks_low']) &
                               (nCount_peaks < criteria_cutoffs[group, 'nCount_peaks_up']) &
                               (pct_reads_in_peaks > criteria_cutoffs[group, 'pct_reads_in_peaks']) &
                               (blacklist_ratio < criteria_cutoffs[group, 'blacklist_ratio']) &
                               (nucleosome_signal < criteria_cutoffs[group, 'nucleosome_signal']) &
                               (TSS.enrichment > criteria_cutoffs[group, 'TSS.enrichment']))
   
   print(paste0("Filtered ", group, ': ', ncol(filt_seur_ls[[i]]), ' cells (', round(ncol(filt_seur_ls[[i]])/ncol(seur_ls[[i]]),2)*100, '%)'))
}  
names(filt_seur_ls) = groups


# print QC metrics after filtering
pdf(file.path(out_dir, 'QC', 'QC_after_filtering.pdf'), width=13, height=4.5)
for (this_group in groups) {
   print(filt_seur_ls[[this_group]])
   p = VlnPlot(object=filt_seur_ls[[this_group]], features=my_features, pt.size=0, ncol=6, combine=F) 
   p = lapply(1:length(p), function(i) 
      p[[i]] + labs(x=my_features[i], title=ifelse(i==1, this_group, '')) + NoLegend() +
         theme(axis.text.x=element_blank()))
   print(p[[1]] | p[[2]] | p[[3]] | p[[4]] | p[[5]] | p[[6]])
}
dev.off()

saveRDS(filt_seur_ls, file.path(out_dir, 'seurat_object', '02_seurs_list_after_filtering.rds'))



# ##############################
# #### normalization, idenitfication of variable features and dimensional reduction 
# create a common peak set
combined_peaks_ls   = lapply(filt_seur_ls, granges)
combined_peaks_cre  = reduce(x=c(combined_peaks_ls[['CD8_cre_pos']], combined_peaks_ls[['CD8_cre_neg']]))
combined_peaks_KOWT = reduce(x=c(combined_peaks_ls[['IRF2_KO']], combined_peaks_ls[['WT']]))
combined_peaks_ls   = list(combined_peaks_all, combined_peaks_cre, combined_peaks_KOWT)
# Filter out bad peaks based on length
peakwidths_ls     = lapply(combined_peaks_ls, width)
combined_peaks_ls = lapply(1:length(combined_peaks_ls), function(i) combined_peaks_ls[[i]][peakwidths_ls[[i]]  < 10000 & peakwidths_ls[[i]] > 20])
names(combined_peaks_ls) = c('cre', 'KOWT')


# create seurat object with merged peaks from CD8_cre_pos/neg samples
pre_merge_seur_ls  = lapply(c('CD8_cre_pos', 'CD8_cre_neg'), function(group) {
   meta_data = filt_seur_ls[[group]]@meta.data 
   frag_obj  = Fragments(filt_seur_ls[[group]])
   seur_counts = FeatureMatrix(
      fragments = frag_obj,
      features = combined_peaks_ls[['cre']],
      cells = rownames(meta_data)
   )
   seur_assay = CreateChromatinAssay(seur_counts, fragments=frag_obj)
   CreateSeuratObject(seur_assay, assay="peaks", meta.data=meta_data)
})
merged_seurs_cre = merge(pre_merge_seur_ls[[1]], pre_merge_seur_ls[[2]], add.cell.ids=c('CD8_cre_pos', 'CD8_cre_neg'))

# create seurat object with merged peaks from IRF2_KO and WT samples
pre_merge_seur_ls  = lapply(c('IRF2_KO', 'WT'), function(group) {
   meta_data = filt_seur_ls[[group]]@meta.data 
   frag_obj  = Fragments(filt_seur_ls[[group]])
   seur_counts = FeatureMatrix(
      fragments = frag_obj,
      features = combined_peaks_ls[['KOWT']],
      cells = rownames(meta_data)
   )
   seur_assay = CreateChromatinAssay(seur_counts, fragments=frag_obj)
   CreateSeuratObject(seur_assay, assay="peaks", meta.data=meta_data)
})
merged_seurs_KOWT = merge(pre_merge_seur_ls[[1]], pre_merge_seur_ls[[2]], add.cell.ids=c('IRF2_KO', 'WT'))

merged_seurs_ls = list(cre=merged_seurs_cre, KOWT=merged_seurs_KOWT)
saveRDS(merged_seurs_ls, file.path(out_dir, 'seurat_object', '03_merged_seurs_list.rds'))


# ##############################
# #### normalization, idenitfication of variable features and dimensional reduction 
norm_seurs_ls_for_plot = lapply(merged_seurs_ls, function(seur)
   seur %>%
   RunTFIDF() %>%
   FindTopFeatures(min.cutoff=5) %>%
   RunSVD() %>%
   RunUMAP(reduction='lsi', dims=2:30))
norm_seurs_ls_for_integ = lapply(merged_seurs_ls, function(merged_seurs) {
   seur_ls = SplitObject(merged_seurs, split.by="group")
   for (i in 1:length(seur_ls)) {
      seur_ls[[i]] = seur_ls[[i]] %>%
         RunTFIDF() %>%
         FindTopFeatures(min.cutoff=5) %>%
         RunSVD() %>%
         RunUMAP(reduction='lsi', dims=2:30)  
   }
   seur_ls
})

saveRDS(norm_seurs_ls_for_integ, file.path(out_dir, 'seurat_object', '04_normalized_seurs_list.rds'))
saveRDS(norm_seurs_ls_for_plot, file.path(out_dir, 'seurat_object', '04_normalized_seurs_list_for_plot.rds'))
p1_ls = lapply(norm_seurs_ls_for_plot, function(seur) DimPlot(object=seur, label=TRUE, group.by='group'))

pdf(file.path(out_dir, 'QC', 'cor_of_depth_n_lsi.pdf'))
for (x in c('all', 'cre', 'KOWT')) {
   print(DepthCor(norm_seurs_ls_for_plot[[x]]) + ggtitle(x))
}
dev.off()
print('normalization and other processing steps done!')



##########################
#### integrate cre samples (CD8_cre_pos/neg) and IRF2_KO with WT (KOWT)
#### cell type annotation
#### add the gene activity matrix to the Seurat object as a new assay and normalize it
#### annotate cells using 2 types of labels (main and fine labels)
# find integration anchors
integ_seurs_ls = p2_ls = clust_seurs_ls = NULL
resols = seq(.3, 1.3, .2)
pdf(file.path(out_dir, 'UMAP_groups.pdf'), width=12, height=5)
bed.se = readRDS("/cluster/projects/mcgahalab/ref/scrna/scRNAseq_datasets/immgen.rds") 
lbls = c('label.main', 'label.fine')
for (x in c('cre', 'KOWT')) {
   # find anchors among samples
   integration_anchors = FindIntegrationAnchors(
      object.list = norm_seurs_ls_for_integ[[x]],
      anchor.features = rownames(norm_seurs_ls_for_integ[[x]][[1]]),
      reduction = "rlsi",
      dims = 2:30
   )
   
   # integrate LSI embeddings
   integ_seurs_ls[[x]] = IntegrateEmbeddings(
      anchorset  = integration_anchors,
      reductions = norm_seurs_ls_for_plot[[x]][["lsi"]],
      new.reduction.name = "integrated_lsi",
      dims.to.integrate  = 1:30
   )
   
   # create a new UMAP using the integrated embeddings
   integ_seurs_ls[[x]] = RunUMAP(integ_seurs_ls[[x]], reduction="integrated_lsi", dims=2:30)
   print(paste('*******', date(), x, 'samples integrated done!'))
   integ_seurs_ls[[x]]
   
   p2_ls[[x]] = DimPlot(integ_seurs_ls[[x]], group.by="group")
   p = (p1_ls[[x]] + ggtitle("Unintegrated")) | (p2_ls[[x]] + ggtitle("Integrated"))
   print(p)
   
   # find neighbors and cluster
   clust_seurs_ls[[x]] = FindNeighbors(object=integ_seurs_ls[[x]], reduction='integrated_lsi', dims=2:30)
   clust_seurs_ls[[x]] = FindClusters(object=clust_seurs_ls[[x]], algorithm=3, resolution=resols)
   print(paste('*******', date(), x, 'samples clustered done!'))
   
   # add gene information to the seurat object and calculate gene activity
   Annotation(clust_seurs_ls[[x]]) = annotations       
   gene_activities = GeneActivity(clust_seurs_ls[[x]])
   clust_seurs_ls[[x]][['gene_activity']] = CreateAssayObject(counts=gene_activities)
   
   DefaultAssay(clust_seurs_ls[[x]]) = 'gene_activity'
   clust_seurs_ls[[x]] = clust_seurs_ls[[x]] %>%
      NormalizeData(scale.factor=median(clust_seurs_ls[[x]]$nCount_gene_activity)) %>%
      FindVariableFeatures() %>%
      ScaleData(features=row.names(clust_seurs_ls[[x]]))
   print(paste('*******', date(), x, 'samples gene activity assay added and preprocessed done!'))
   
   # annotate cells using SingleR and gene activity data
   for (resol in resols) {
      for (lbl in lbls) {
         immgen_anno = SingleR(test=GetAssayData(clust_seurs_ls[[x]]), ref=bed.se, 
                               assay.type.test=1, labels=bed.se[[lbl]],
                               clusters=clust_seurs_ls[[x]]@meta.data[, paste0('peaks_snn_res.', resol)])   
         
         cluster_ids = setNames(immgen_anno$labels, 
                                as.character(rownames(immgen_anno)))
         clust_seurs_ls[[x]]@meta.data[, paste0('immgen_res_', resol, '_', lbl)] = cluster_ids[as.character(clust_seurs_ls[[x]]@meta.data[, paste0('peaks_snn_res.', resol)])]
         print(paste("*********", date(), 'Resolution', resol, lbl, x, 'done!'))
      }
   }
   
}
dev.off()
saveRDS(integ_seurs_ls, file.path(out_dir, 'seurat_object', '05_seurs_list_integrated_peaks.rds'))
saveRDS(clust_seurs_ls, file.path(out_dir, 'seurat_object', '06_seurs_list_clustered_n_annotated.rds'))




# UMAP all resolutions split by sample
# downsample cells to have the same nb of cells in WT and KO
pdf(file.path(out_dir, 'UMAP_all_resol.pdf'), width=13, height=4)
for (x in c('cre', 'KOWT')) {
   min_sample_size = min(table(clust_seurs_ls[[x]]$group))
   sub_sampled = subset(x=clust_seurs_ls[[x]], downsample=min_sample_size) 
   sub_sampled
   
   for (resol in resols) {
      # print with number cell labels
      p1 = DimPlot(sub_sampled, reduction='umap', group.by=paste0('peaks_snn_res.', resol), 
                   label=T, repel=T, label.size=3, pt.size=2, shuffle=T) +
         ggtitle(paste(x, 'samples - Cluster #, Res:', resol)) +
         theme(plot.title=element_text(size=11), legend.text=element_text(size=6)) +
         guides(color=guide_legend(override.aes=list(size=2), keywidth=.2, keyheight=.6, ncol=1) )
      p1$layers[[1]]$aes_params$alpha = 0.6
      p2 = DimPlot(sub_sampled, reduction='umap', split.by='group', group.by=paste0('peaks_snn_res.', resol), 
                   label=T, repel=T, label.size=3, pt.size=2, shuffle=T) + NoLegend() + ggtitle(NULL)  
      p2$layers[[1]]$aes_params$alpha = 0.6
      print(ggarrange(p1, p2, ncol=2, nrow=1, widths=c(1.5, 4))) #c(1.7, 2.5)))
      
      # print with main cell type
      p3 = DimPlot(sub_sampled, reduction='umap', group.by=paste0('immgen_res_', resol, '_label.main'), 
                   label=T, repel=T, label.size=2, pt.size=2, shuffle=T) +
         ggtitle(paste(x, 'samples - Main label, Res:', resol)) +
         theme(plot.title=element_text(size=11), legend.text=element_text(size=4)) +
         guides(color=guide_legend(override.aes=list(size=2), keywidth=.2, keyheight=.6, ncol=1) )
      p3$layers[[1]]$aes_params$alpha = 0.6
      p4 = DimPlot(sub_sampled, reduction='umap', split.by='group', group.by=paste0('immgen_res_', resol, '_label.main'), 
                   label=T, repel=T, label.size=2, pt.size=2, shuffle=T) + NoLegend() + ggtitle(NULL)  
      p4$layers[[1]]$aes_params$alpha = 0.6
      print(ggarrange(p3, p4, ncol=2, nrow=1, widths=c(1.5, 4))) #c(1.7, 2.5)))
      
      # print with fine cell type
      p5 = DimPlot(sub_sampled, reduction='umap', group.by=paste0('immgen_res_', resol, '_label.fine'), 
                   label=T, repel=T, label.size=2, pt.size=2, shuffle=T) +
         ggtitle(paste(x, 'samples - Fine label, Res:', resol)) +
         theme(plot.title=element_text(size=11), legend.text=element_text(size=4)) +
         guides(color=guide_legend(override.aes=list(size=2), keywidth=.2, keyheight=.6, ncol=1) )
      p5$layers[[1]]$aes_params$alpha = 0.6
      p6 = DimPlot(sub_sampled, reduction='umap', split.by='group', group.by=paste0('immgen_res_', resol, '_label.fine'), 
                   label=T, repel=T, label.size=2, pt.size=2, shuffle=T) + NoLegend() + ggtitle(NULL)  
      p6$layers[[1]]$aes_params$alpha = 0.6
      print(ggarrange(p5, p6, ncol=2, nrow=1, widths=c(1.7, 4))) #c(1.7, 2.5)))
   }
}
dev.off()
print('***umap_all_resot_split_by_sample done***')



# Find RNA markers for clusters of Cre samples and KOWT samples, and draw a heatmap of top 20/10 and botttom 20/10 RNA markers per cluster
lbls  = c('peaks_snn', 'label.main', 'label.fine')
resol = 1.3
for (x in c('cre', 'KOWT')) {
   for (lbl in lbls) {
     rna_markers_ls = NULL
     file_suffix = ifelse(lbl=='peaks_snn', 'cluster_nb', ifelse(lbl=='label.main', 'main_label', 'fine_label'))
     i = 1
     this_ident  = ifelse(lbl=='peaks_snn', paste0(lbl, '_res.', resol), 
                                            paste0('immgen_res_', resol, '_', lbl))
     Idents(clust_seurs_ls[[x]]) = this_ident
     rna_markers_ls[[i]]         = FindAllMarkers(clust_seurs_ls[[x]], assay="gene_activity")
     rna_markers_ls[[i]]$avg_FC  = (2^rna_markers_ls[[i]]$avg_log2FC)
     rna_markers_ls[[i]] = rna_markers_ls[[i]][, c('gene', 'cluster', 'avg_FC', 'avg_log2FC', 'p_val', 'p_val_adj', 'pct.1', 'pct.2')]
     rna_markers_ls[[i]] = data.frame(resolution=resol, rna_markers_ls[[i]])
     rna_markers = do.call(rbind, rna_markers_ls)
     write.csv(rna_markers, file=file.path(out_dir, paste0('gene_activity_markers_', file_suffix, '_', x,'.csv')), row.names=F)
         
   
      top_nbs     = c(10, 20)
      pdf_heights = c(25, 45)
      for (j in 1:length(top_nbs)) {
         top_nb = top_nbs[j]
         pdf(file.path(out_dir, paste0('Heatmap_gene_activity_top_', top_nb, '_', file_suffix, '_', x, '.pdf')), height=pdf_heights[j])
         Idents(clust_seurs_ls[[x]]) = ifelse(lbl=='peaks_snn', paste0(lbl, '_res.', resol), 
                                              paste0('immgen_res_', resol, '_', lbl))
         rna_markers_ls[[i]] %>%
            group_by(cluster) %>%
            top_n(n=top_nb, wt=avg_log2FC) -> top_markers
         p1 = DoHeatmap(clust_seurs_ls[[x]], features=unique(top_markers$gene), assay="gene_activity", angle=90) +
            theme(text = element_text(size = 13), title=element_text(size=11)) + NoLegend() + 
            ggtitle(paste0(x, ' samples, Res: ', resol, ', Top ', top_nb, ' up-regulated gene activity markers'))
         print(p1)
        
         rna_markers_ls[[i]] %>%
            group_by(cluster) %>%
            top_n(n=-top_nb, wt=avg_log2FC) -> bottom_markers
         p2 = DoHeatmap(clust_seurs_ls[[x]], features=unique(bottom_markers$gene), assay="gene_activity", angle=90) +
            theme(text = element_text(size = 13), title=element_text(size=11)) + NoLegend() +
            ggtitle(paste0(x, ' samples, Res: ', resol, ', Top ', top_nb, ' down-regulated gene activity markers'))
         print(p2)
         dev.off()
      }
   }
}




#######################################
#### remove cell types that do not represent CD8 T cells for resolution 0.3
resol = 0.3
clust_seurs_ls = readRDS(file.path(out_dir, 'seurat_object', '06_seurs_list_clustered_n_annotated.rds'))

keep_cells = list(all=c('T cells (T.4MEM44H62L)', 'T cells (T.8EFF.OT1.D5.VSVOVA)', 'Tgd (Tgd.vg2-)'),
                  cre=c('T cells (T.4MEM44H62L)', 'T cells (T.8EFF.OT1.D5.VSVOVA)', 'T cells (T.CD8.1H)'),
                  KOWT=c('ILC (LPL.NCR+CNK)', 'T cells (T.4MEM44H62L)', 'T cells (T.8EFF.OT1.D5.VSVOVA)', 'Tgd (Tgd.vg2-.TCRbko)'))
clean_seurs_ls = lapply(c('all', 'cre', 'KOWT'), function(x) {
   subset(clust_seurs_ls[[x]], immgen_res_0.3_label.fine %in% keep_cells[[x]])  
})
names(clean_seurs_ls) = c('all', 'cre', 'KOWT')



#######################################################
pdf(file.path(out_dir, 'violon_plot_removed_cells_features.pdf'), width=25, height=8)
for (x in c('all', 'cre', 'KOWT')) {
   clean_seurs_ls[[x]]$cells = row.names(clean_seurs_ls[[x]]@meta.data)
   umap_df  = clean_seurs_ls[[x]]@reductions$umap@cell.embeddings
   if (x=='all') rm_cells = c(row.names(umap_df[umap_df[,1] < -9, , drop=F]), row.names(umap_df[umap_df[,2] > 4, , drop=F]))
   if (x=='cre') rm_cells = row.names(umap_df[umap_df[,1] < -10, , drop=F])
   if (x=='KOWT') rm_cells = row.names(umap_df[umap_df[,2] > 2, , drop=F])
   clean_seurs_ls[[x]]$rm_cell = ifelse(clean_seurs_ls[[x]]$cells %in% rm_cells, 'remove', 'keep')
   
   # print with fine cell type
   p = VlnPlot(
      object = clean_seurs_ls[[x]],
      features = c('nCount_peaks', 'nFeature_peaks', 'duplicate', 'chimeric', 'unmapped', 'lowmapq',
                   'TSS_fragments', 'on_target_fragments', 'peak_region_cutsites',
                   'mitochondrial', 'nonprimary', 'blacklist_region_fragments', 'peak_region_fragments',
                   'pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_percentile', 'TSS.enrichment', 'TSS.percentile', 'nCount_gene_activity'),
      pt.size = 0.1, group.by='rm_cell', 
      ncol = 10
   )
   
   print(p)
}
dev.off()


#######################################################
# Need more cleaning to remove some cells that are appearing far away from most of cells in their cluster
for (x in c('all', 'cre', 'KOWT')) {
   clean_seurs_ls[[x]]$cells = row.names(clean_seurs_ls[[x]]@meta.data)
   umap_df  = clean_seurs_ls[[x]]@reductions$umap@cell.embeddings
   if (x=='all') rm_cells = c(row.names(umap_df[umap_df[,1] < -9, , drop=F]), row.names(umap_df[umap_df[,2] > 4, , drop=F]))
   if (x=='cre') rm_cells = row.names(umap_df[umap_df[,1] < -10, , drop=F])
   if (x=='KOWT') rm_cells = row.names(umap_df[umap_df[,2] > 2, , drop=F])
   clean_seurs_ls[[x]] = subset(clean_seurs_ls[[x]], cells %in% rm_cells, invert=T)
}
saveRDS(clean_seurs_ls, file=file.path(out_dir, 'seurat_object', '07_clean_seurs_list.rds'))
##### end #####







