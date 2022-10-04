## Build Transcription Factor BED file
library("ChIPseeker")
library("clusterProfiler")
library("enrichplot")
library("ggplot2")
library("dplyr")
library("reshape2")
library("GenomicFeatures")
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)


gbuild <- 'mm10'
human_txdb <- '/cluster/projects/mcgahalab/ref/genomes/human/hg38/GTF/genome.gtf'
mouse_txdb <- '/cluster/projects/mcgahalab/ref/genomes/mouse/GRCm38/GTF/genome.gtf'
motif_dir <- '/cluster/projects/mcgahalab/ref/tf/JASPAR2020'
motif_interest <- setNames('MA0051.1', 'IRF2')

########################
#### Reference Data ####
txdb <- makeTxDbFromGFF(file = human_txdb, format = "gtf")

# Create a reference map of ENSEMBL to SYMBOL
gbuild <- 'hg38'
if(gbuild %in% c('GRCh38', 'GRCh37', 'hg19', 'hg38')){
  if(any(gbuild %in% c('hg19', 'GRCh37'))) {
    library(BSgenome.Hsapiens.UCSC.hg37)
    bsgenome <- BSgenome.Hsapiens.UCSC.hg37
  } else if(any(gbuild %in% c('hg38', 'GRCh38'))){
    library(BSgenome.Hsapiens.UCSC.hg38)
    bsgenome <- BSgenome.Hsapiens.UCSC.hg38
  }
  
  library("org.Hs.eg.db")
  species <- 'org.Hs.eg.db'
  jaspar_species <- 'Homo sapiens'
  genome <- org.Hs.eg.db
} else if(gbuild %in% c('mm10', 'mm9', 'GRCm38')){
  if(gbuild %in% c('mm9')) {
    library(BSgenome.Mmusculus.UCSC.mm9)
    bsgenome <- BSgenome.Mmusculus.UCSC.mm9
  } else if(any(gbuild %in% c('mm10', 'GRCm38'))){
    library(BSgenome.Mmusculus.UCSC.mm10)
    bsgenome <- BSgenome.Mmusculus.UCSC.mm10
  }
  
  library("org.Mm.eg.db")
  species <- 'org.Mm.eg.db'
  jaspar_species <- 'Mus musculus'
  genome <- org.Mm.eg.db
}
txby <- keys(genome, 'ENSEMBL')
gene_ids <- mapIds(genome, keys=txby, column='SYMBOL',
                   keytype='ENSEMBL', multiVals="first")

######################################
#### Get all predicted IRF2 sites ####
## get unique genes
genes <- ChIPseeker:::getGene(txdb, 'gene')
tss <- ifelse(strand(genes) == "+", start(genes), end(genes))
pr <- GRanges(seqnames = seqnames(genes), ranges = IRanges(tss, tss), strand = strand(genes))
pr_idx <- which(duplicated(pr))
genes <- genes[-pr_idx,]


## Get Motif data from JASPAR2020
if(is.null(motif_interest)){
  pfm <- getMatrixSet(x = JASPAR2020,
                      opts = list(collection = "CORE", species=jaspar_species))
} else {
  pfm <- getMatrixSet(x = JASPAR2020,
                      opts = list(collection = "CORE", ID=motif_interest))
}


# get promoters
dtss <- 3000
chrs <- GRanges(seqinfo(bsgenome))
chrs <- keepStandardChromosomes(chrs, pruning.mode='coarse')
promoter <- getPromoters(TxDb=txdb, upstream=dtss, downstream=dtss)

# motif inference from JASPAR2020
motif_ix <- matchMotifs(pwms = pfm,
                        subject = chrs,
                        out = 'positions',
                        bg='genome',
                        p.cutoff = 0.00005, #0.0001,
                        genome = bsgenome)
motif <- sort(motif_ix[[1]])

motif_df <- as.data.frame(motif)[,-4]
motif_df[,4] <- "."
write.table(motif_df, file=file.path(motif_dir, paste0(names(motif_interest), "_", motif_interest, ".genome.bed")),
            col.names = F, row.names = F, sep="\t", quote = F)
saveRDS(motif, file=file.path(motif_dir, "rds", paste0(names(motif_interest), "_", motif_interest, ".genome.rds")))