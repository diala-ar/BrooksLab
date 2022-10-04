log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("ChIPseeker")
library("clusterProfiler")
library("enrichplot")
library("ggplot2")
library("dplyr")
library("reshape2")
library("GenomicFeatures")

library("KEGG.db")
library("GO.db")
library("Category")
library("GOstats")

setwd('/cluster/projects/mcgahalab/data/mcgahalab/sabelo_irf2/cutandrun/results/peaks/seacr/IgG')
gbuild <- 'GRCm38'
txdb <- makeTxDbFromGFF(file = '/cluster/projects/mcgahalab/ref/genomes/mouse/GRCm38/GTF/genome.gtf', format = "gtf")

# Create a reference map of ENSEMBL to SYMBOL
gbuild <- snakemake@params[["build"]]
txdb <- makeTxDbFromGFF(file = snakemake@params[['gtf']], format = "gtf")
if(gbuild %in% c('GRCh38', 'GRCh37', 'hg19', 'hg38')){
  library("org.Hs.eg.db")
  species <- 'org.Hs.eg.db'
} else if(gbuild %in% c('mm10', 'mm9', 'GRCm38')){
  library("org.Mm.eg.db")
  species <- 'org.Mm.eg.db'
}


## Read a file containing all the files and groupings
# files <- strsplit(paste(list.files(pattern='bed$'), collapse=","), ",")[[1]]
files <- strsplit(snakemake@params[['files']], ',')[[1]]
file_list <- data.frame('file'=files,
                        'group'=gsub("^.*_([a-zA-Z0-9]+).stringent.*$", "\\1", files))
print(file_list)
# colnames(file_list) <- c('file', 'group')
# file_list <- data.frame("file"=c("HD125_Mono_K27.stringent.bed",
#                                  "HD33_Mono_K27.stringent.bed",
#                                  "HD125_Mono_Me3.stringent.bed"),
#                         "group"=c("K27", "K27", "Me3"))


## Read in all the peak files based on their group ID
file_list_split <- split(file_list, file_list$group)
promoter <- getPromoters(TxDb=txdb, upstream=as.integer(snakemake@params[['promoter_range']]), 
                         downstream=as.integer(snakemake@params[['promoter_range']]))
# promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
ggps <- lapply(file_list_split, function(files_df){
  gg_group <- list()
  
  file_sizes <- unlist(sapply(files_df$file, file.info)['size',])
  if(any(file_sizes == 0)) files_df <- files_df[-which(file_sizes<1000),]
  files <- as.list(files_df$file)
  names(files) <- gsub("\\.stringent.bed", "", files_df$file) %>%
    gsub("^.*?\\/", "", .)
  tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)
  
  
  ## Coverage around TSS
  gg_group[['tss_cov']] <- plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
  # pdf("~/xfer/tss_cov.pdf"); print(gg_group[['tss_cov']]); dev.off()
  
  ## Average coverage around TSS per sample
  tryCatch({
    gg_group[['tss_cov_avg']] <- plotAvgProf(tagMatrixList, xlim=c(-3000, 3000),
                                          conf=0.95,resample=500, facet="row")
  }, error=function(e){
    print("Parallelization failed, too much memory required")
    gg_group[['tss_cov_avg']] <- ggplot(data.frame(1))
  })
  
  ## Plot Deeptools style barplots around the TSS
  pdf("~/xfer/tagmat_deep.pdf")
  gg_group[['tagmat_deep']] <- tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)
  dev.off()
  # pdf("~/xfer/tagmat_deep.pdf"); print(gg_group[['tagmat_deep']]); dev.off()
  
  ## Categorize the peaks based on their genomic annotation (CEAS style plots)
  peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                         tssRegion=c(-3000, 3000), verbose=FALSE)
  gg_group[['peak_anno']][['a']] <- plotAnnoBar(peakAnnoList)
  gg_group[['peak_anno']][['b']] <- plotDistToTSS(peakAnnoList)
  # pdf("~/xfer/peak_anno_bar.pdf"); print(gg_group[['peak_anno']][['a']]); dev.off()
  # pdf("~/xfer/peak_anno_dist.pdf"); print(gg_group[['peak_anno']][['b']]); dev.off()

  # Write out the peak annotation to a flat file
  lapply(names(peakAnnoList), function(sampleid){
    sampleid <- gsub("^.*\\/", "", sampleid)
    # outid <- gsub("\\.[a-zA-Z]*$", paste0(".", sampleid, ".tsv"), snakemake@output[['peak_anno']])
    outid <- gsub("\\.[a-zA-Z]*$", paste0(".", sampleid, ".tsv"), file.path(getwd(), "anno.txt"))

    write.table(as.data.frame(peakAnnoList[[sampleid]]),
                file = outid, quote = FALSE, row.names = FALSE, col.names = TRUE)
  })
  
  ## Functional analysis on the genes associated with the peaks
  genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
  if(all(grepl("^ENS", head(genes[[1]])))){
    keytype <- 'ENSEMBL'
  } else {
    keytype <- 'ENTREZID'
  }
  if(gbuild %in% c('GRCh38', 'GRCh37', 'hg19', 'hg38')){
    compKEGG <- compareCluster(geneCluster   = genes,
                               fun           = "enrichKEGG",
                               pvalueCutoff  = 0.05,
                               use_internal_data = TRUE,
                               pAdjustMethod = "fdr")
    gg_group[['functional_anno']][['a']] <- dotplot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")
    # pdf("~/xfer/functional_kegg.pdf"); print(gg_group[['functional_anno']][['a']]); dev.off()
  }
  
  compGO <- compareCluster(geneCluster   = genes,
                           fun           = "enrichGO",
                           OrgDb         = species,
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "fdr",
                           keyType       = keytype,
                           ont           = 'ALL')
  gg_group[['functional_anno']][['b']] <- dotplot(compGO, showCategory = 15, title = "GO Ontology Enrichment Analysis")
  # pdf("~/xfer/functional_go.pdf", width=14); print(gg_group[['functional_anno']][['b']]); dev.off()
  
  ## Enrichment of peak overlap between files
  # Permutation Test for overlapping peaks betwen the query peak and all target peaks
  # shuffle bootstraps
  if(length(files)>1){
    enr_overlap <- lapply(files, function(query_f){
      target_f <- unlist(files)
      enrichPeakOverlap(queryPeak     = query_f,
                        targetPeak    = target_f,
                        TxDb          = txdb,
                        pAdjustMethod = "fdr",
                        nShuffle      = 50,
                        chainFile     = NULL,
                        verbose       = FALSE)
    })
    ov_pmat <- sapply(enr_overlap, function(i) i$p.adjust)
    ov_n <- sapply(enr_overlap, function(i) i$N_OL/i$qLen)
    rownames(ov_n) <- rownames(ov_pmat) <- names(enr_overlap)
    melt_n <- melt(round(ov_n,4))
    melt_pmat <- melt(round(ov_pmat,4))
    melt_pmat$n <- melt_n$value
    
    gg_group[['peak_ov']] <- ggplot(data=melt_pmat, aes(x=Var1, y = Var2, color=value, size = n)) + 
      geom_point() + 
      scale_color_gradientn(colours = c('#fdd49e', '#ef6548', '#636363', '#636363'), # rev(viridis::inferno(6))[-c(1,5)], 
                            breaks = c(0, 0.001, 0.05, 0.5), 
                            trans='log') + 
      cowplot::theme_cowplot() + 
      theme(axis.line  = element_blank()) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      ylab('') +
      theme(axis.ticks = element_blank()) 
  } else {
    gg_group[['peak_ov']] <- ggplot(data.frame(1))
  }
  
  return(gg_group)
})


## Function and plotting commands
plotIt <- function(ggps, id, outfile, adj_hw='neither'){
  max_f <- max(sapply(ggps, function(i) length(levels(i[[1]]$data$.id))))
  if(adj_hw=='height'){
    adj_height <- (max_f*1.5) + (7-2)
    adj_width <- 7
  } else if(adj_hw=='width'){
    adj_width <- (max_f*1.5) + (7-2)
    adj_height <- 7
  } else {
    adj_width <- adj_height <- 7
  }
  
  pdf(outfile, height=adj_height, width=adj_width)
  lapply(ggps, function(gp){
    if(!any(class(gp[[id]]) == 'gg')){
      lapply(gp[[id]], function(gp_i){
        print(gp_i)
      })
    } else {
      print(gp[[id]])
    }
  })
  dev.off()
}
plotIt(ggps, 'tss_cov', snakemake@output[['tss_cov']]) # pdf("~/xfer/tss_cov.pdf")
plotIt(ggps, 'tss_cov_avg', snakemake@output[['tss_cov_avg']], adj_hw = 'height') # pdf("~/xfer/tss_cov_avg.pdf", height=adj_height) 
plotIt(ggps, 'tagmat_deep', snakemake@output[['tagmat_deep']], adj_hw = 'width') # pdf("~/xfer/tagmat_deep.pdf", height=adj_height) 
plotIt(ggps, 'peak_anno', snakemake@output[['peak_anno']], adj_hw = 'height') # pdf("~/xfer/peak_anno.pdf", height=adj_height) 
plotIt(ggps, 'functional_anno', snakemake@output[['functional_anno']], adj_hw = 'width') # pdf("~/xfer/functional_anno.pdf", height=adj_height) 
plotIt(ggps, 'peak_ov', snakemake@output[['peak_ov']]) # pdf("~/xfer/peak_ov.pdf", height=adj_height) 