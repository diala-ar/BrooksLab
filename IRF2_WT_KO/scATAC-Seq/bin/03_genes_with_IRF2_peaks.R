library(dplyr)
library(GenomicFeatures)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(Signac)
library(Seurat)
library(BSgenome.Mmusculus.UCSC.mm10)
library(org.Mm.eg.db)
library(openxlsx)



wrk_dir    = '/cluster/projects/mcgahalab/data/brookslab/sabelo/2022_03_scATAC/'
in_dir     = 'resources/scATAC_data'
out_dir    = 'results_v2'
setwd(wrk_dir)

mouse_txdb = '/cluster/projects/mcgahalab/ref/genomes/mouse/GRCm38/GTF/genome.gtf'
motifs     = list("IRF2"=c('Mus musculus'='MA0051.1', 'Homo sapiens'='MA0051.1'))


########################
#### Reference Data ####
# Create a reference map of ENSEMBL to SYMBOL
txdb     = makeTxDbFromGFF(file = mouse_txdb, format = "gtf")
  
bsgenome = BSgenome.Mmusculus.UCSC.mm10
chrs     = GRanges(seqinfo(bsgenome))
chrs     = keepStandardChromosomes(chrs, pruning.mode='coarse')

genome    = org.Mm.eg.db
txby      = keys(genome, 'ENSEMBL')
gene_ids  = mapIds(genome, keys=txby, column='SYMBOL',
                   keytype='ENSEMBL', multiVals="first")



#### Get all predicted IRF2 sites ####
idx = 1 # Set Motif of interest
jaspar_species = 'Mus musculus'
motif_interest = setNames(motifs[[idx]][jaspar_species],
                          names(motifs)[idx]) 
## Get Motif data from JASPAR2020
pfm = getMatrixSet(x = JASPAR2020,
                   opts = list(collection = "CORE", ID=motif_interest))
# motif inference from JASPAR2020 to genome
motif_ix <- matchMotifs(pwms = pfm,
                        subject = chrs,
                        out = 'positions',
                        bg='genome',
                        p.cutoff = 0.00005, #0.0001, # for lenient_cutoff files
                        genome = bsgenome)
irf2_motif <- sort(motif_ix[[1]])


# get promoters from genome
dtss = 3000
promoter = promoters(x=genes(txdb), upstream=dtss, downstream=dtss) %>%  
   sort %>%
   keepStandardChromosomes(., pruning.mode='coarse')
seqlevelsStyle(promoter) = 'UCSC'
promoter$gene_name = gene_ids[promoter$gene_id]
promoter$gene_name[is.na(promoter$gene_name)] = promoter$gene_id[is.na(promoter$gene_name)]


# Intersect motif with promoter (Not Saved)
ov_idx        = findOverlaps(irf2_motif, promoter, type='within')
ov_irf2_motif = data.frame(irf2_motif[queryHits(ov_idx)])[, -5]
ov_promoter   = data.frame(promoter[subjectHits(ov_idx)])
names(ov_promoter) = paste0('promoter_', names(ov_promoter))
irf2_motif_in_promoter = cbind(ov_irf2_motif, ov_promoter)
irf2_motif_in_promoter = makeGRangesFromDataFrame(irf2_motif_in_promoter, keep.extra.columns=T)


annot_seurs_ls = readRDS(file.path(out_dir, 'seurat_object', '09_seurs_list_clustered_n_annotated.rds'))
# my_genes is a list of interest defined by Sabelo
my_genes       = read.csv(file.path('resources', 'gene_list_investigate_irf2_motifs_in_promoter.csv'))[, 1] 


for (x in c('cre', 'KOWT')) {
   DefaultAssay(annot_seurs_ls[[x]]) = 'peaks'
   Idents(annot_seurs_ls[[x]]) = 'group'
   samples  = unique(annot_seurs_ls[[x]]$group) 
   
   annot_seurs_ls[[x]] = AddMotifs(object=annot_seurs_ls[[x]], genome=bsgenome, pfm=pfm)
   
   peaks_df = data.frame(chr=paste0('chr', gsub('chr(.+)-([0-9]+)-([0-9]+)', '\\1', row.names(annot_seurs_ls[[x]]))), 
                         start=gsub('chr(.+)-([0-9]+)-([0-9]+)', '\\2', row.names(annot_seurs_ls[[x]])), 
                         end=gsub('chr(.+)-([0-9]+)-([0-9]+)', '\\3', row.names(annot_seurs_ls[[x]])),
                         peak=row.names(annot_seurs_ls[[x]]))
   peaks_gr = makeGRangesFromDataFrame(peaks_df, keep.extra.columns=T)
   
   # write excel file with coordinates of irf2 motifs, promoter and peaks
   ov_idx   = findOverlaps(peaks_gr, irf2_motif_in_promoter)
   ov_peaks = data.frame(peaks_gr[queryHits(ov_idx)])
   names(ov_peaks) = paste0('peak_', names(ov_peaks))
   ov_irf2_motif_in_promoter = data.frame(irf2_motif_in_promoter[subjectHits(ov_idx)])
   names(ov_irf2_motif_in_promoter)[1:6] = paste0('irf2_mot_', names(ov_irf2_motif_in_promoter)[1:6])
   genes_with_irf2_motif_in_promoter = cbind(samples=x, ov_peaks, ov_irf2_motif_in_promoter)

   temp = NULL
   for (spl in samples) {
      open_peaks  = AccessiblePeaks(annot_seurs_ls[[x]], idents=spl)
      temp[[spl]] = genes_with_irf2_motif_in_promoter[genes_with_irf2_motif_in_promoter$peak_peak %in% open_peaks, ]
      temp[[spl]]$sample = spl
   }
   write.xlsx(temp, file = file.path(out_dir, paste0('genes_with_irf2_motif_in_their_promoter_in_', x, '_samples.xlsx')))#_lenient_cutoff.xlsx')))
   
   
   # print coverage plot for peaks of certain genes of interest
   pdf(file.path(out_dir, paste0('coverage_plot_irf2_motifs_in_promoter_region_in_', x, '_samples.pdf')), width=5.5)#_lenient_cutoff.pdf')), width=5.5)
   for (my_gene in my_genes) {
      regions = unique(genes_with_irf2_motif_in_promoter[genes_with_irf2_motif_in_promoter$promoter_gene_name==my_gene, 'peak_peak'])
      strand  = unique(genes_with_irf2_motif_in_promoter[genes_with_irf2_motif_in_promoter$promoter_gene_name==my_gene, 'promoter_strand'])
      
      for (reg in regions) {
         if (length(reg)==0) next
         chr   = gsub('(.+)-(.+)-(.+)', '\\1', reg)
         start = as.numeric(gsub('(.+)-(.+)-(.+)', '\\2', reg))
         end   = as.numeric(gsub('(.+)-(.+)-(.+)', '\\3', reg))
         
         reg_gr    = StringToGRanges(reg)
         reg_width = width(reg_gr)
         extend_by = 3000 - reg_width
         
         new_start = ifelse(strand=='+', start, start-extend_by)
         new_end   = ifelse(strand=='+', end+extend_by, end)
         reg_gr    = StringToGRanges(paste(chr, new_start, new_end, sep = '-'))
         ov_idx    = findOverlaps(reg_gr, irf2_motif)
         irf2_mot  = irf2_motif[subjectHits(ov_idx)]
         irf2_mot$color = 'black'
         
         cov_plot = CoveragePlot(
            object = annot_seurs_ls[[x]], 
            idents = samples,
            region = paste(chr, new_start, new_end, sep = '-'),
            extend.downstream = 0,
            extend.upstream = 0, 
            annotation = F,
            peaks = F,
            region.highlight = irf2_mot
         )
         peak_plot = PeakPlot(
            object = annot_seurs_ls[[x]], 
            region = paste(chr, new_start, new_end, sep = '-')
         )
         gene_plot = AnnotationPlot(
            object = annot_seurs_ls[[x]], 
            region = paste(chr, new_start, new_end, sep = '-')
         )
         p = CombineTracks(
            plotlist = list(cov_plot, peak_plot, gene_plot),
            heights = c(.8, .05, .08)
         )
         print(p)
      }
   }
   dev.off()
}




#########################################
#### does my_genes contain peaks in their promoters?
for (x in c('cre', 'KOWT')) {
  DefaultAssay(annot_seurs_ls[[x]]) = 'peaks'
  Idents(annot_seurs_ls[[x]]) = 'group'
  samples  = unique(annot_seurs_ls[[x]]$group) 
  
  annot_seurs_ls[[x]] = AddMotifs(object=annot_seurs_ls[[x]], genome=bsgenome, pfm=pfm)
  
  peaks_df = data.frame(chr=paste0('chr', gsub('chr(.+)-([0-9]+)-([0-9]+)', '\\1', row.names(annot_seurs_ls[[x]]))), 
                        start=gsub('chr(.+)-([0-9]+)-([0-9]+)', '\\2', row.names(annot_seurs_ls[[x]])), 
                        end=gsub('chr(.+)-([0-9]+)-([0-9]+)', '\\3', row.names(annot_seurs_ls[[x]])),
                        peak=row.names(annot_seurs_ls[[x]]))
  peaks_gr = makeGRangesFromDataFrame(peaks_df, keep.extra.columns=T)
  
  # write excel file with coordinates of irf2 motifs, promoter and peaks
  ov_idx   = findOverlaps(peaks_gr, promoter)
  ov_peaks = data.frame(peaks_gr[queryHits(ov_idx)])
  names(ov_peaks) = paste0('peak_', names(ov_peaks))
  ov_promoter = data.frame(promoter[subjectHits(ov_idx)])
  names(ov_promoter)[1:6] = paste0('promoter_', names(ov_promoter)[1:6])
  genes_with_peaks_in_promoter = cbind(samples=x, ov_peaks, ov_promoter)
  
  genes_with_peaks_in_promoter = genes_with_peaks_in_promoter[genes_with_peaks_in_promoter$gene_name %in% my_genes, ]

  print(paste(x, 'sample: those genes has peaks in their promoter region', 
               paste(sort(unique(genes_with_peaks_in_promoter$gene_name)), collapse=', ')))

}


#########################################
#### does my_genes contain peaks in their irf2_motif not necessarily in promoter?
for (x in c('cre', 'KOWT')) {
  DefaultAssay(annot_seurs_ls[[x]]) = 'peaks'
  Idents(annot_seurs_ls[[x]]) = 'group'
  samples  = unique(annot_seurs_ls[[x]]$group) 
  
  annot_seurs_ls[[x]] = AddMotifs(object=annot_seurs_ls[[x]], genome=bsgenome, pfm=pfm)
  
  peaks_df = data.frame(chr=paste0('chr', gsub('chr(.+)-([0-9]+)-([0-9]+)', '\\1', row.names(annot_seurs_ls[[x]]))), 
                        start=gsub('chr(.+)-([0-9]+)-([0-9]+)', '\\2', row.names(annot_seurs_ls[[x]])), 
                        end=gsub('chr(.+)-([0-9]+)-([0-9]+)', '\\3', row.names(annot_seurs_ls[[x]])),
                        peak=row.names(annot_seurs_ls[[x]]))
  peaks_gr = makeGRangesFromDataFrame(peaks_df, keep.extra.columns=T)
  
  # write excel file with coordinates of irf2 motifs, promoter and peaks
  ov_idx   = findOverlaps(peaks_gr, irf2_motif)
  ov_peaks = data.frame(peaks_gr[queryHits(ov_idx)])
  names(ov_peaks) = paste0('peak_', names(ov_peaks))
  ov_irf2_motif = data.frame(irf2_motif[subjectHits(ov_idx)])
  names(ov_irf2_motif)[1:6] = paste0('irf2_mot_', names(ov_irf2_motif)[1:6])
  genes_with_irf2_mot_in_peaks = cbind(samples=x, ov_peaks, ov_irf2_motif)
  
  genes_with_irf2_mot_in_peaks = genes_with_irf2_mot_in_peaks[genes_with_irf2_mot_in_peaks$gene_name %in% my_genes, ]
  
  print(paste(x, 'sample - the following genes have IRF2 motif in their peaks', 
              paste(sort(unique(genes_with_peaks_in_promoter$gene_name)), collapse=', ')))
  
}