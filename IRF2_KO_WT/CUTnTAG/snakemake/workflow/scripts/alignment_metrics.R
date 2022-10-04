library(dplyr)
args = commandArgs(trailingOnly=TRUE)

ids       <- strsplit(args[3], ",")[[1]]  # c('HPB63_MDSC_IgG,HPB63_MDSC_K27')
log_path  <- args[1]  # "."
out_path  <- args[2]  # "."
C <- 10000

###################
#### Functions ####
# Process each bowtie2 log file and extract the sequencing depth, 
# alignment rate, and mapped fragments
genAlignMetric <- function(file_path, build){
  if(!file.exists(file_path)){
    resdf <- data.frame(SequencingDepth=NA,
                        MappedFragNum=NA,
                        AlignmentRate=NA)
  } else {
    res <- read.table(file_path, header = FALSE, fill = TRUE, stringsAsFactors = FALSE)
    
    res$V1  <- as.numeric(gsub("%", "", as.character(res$V1)))
    resdf   <- data.frame(SequencingDepth = res$V1[1], 
                          MappedFragNum = res$V1[4] + res$V1[5], 
                          AlignmentRate = res$V1[6])
    
  }
  
  colnames(resdf)[2:3] <- paste0(colnames(resdf)[2:3], "_", build)
  resdf
}

##############
#### Main ####
summ <- t(sapply(ids, function(id){
  # Cycles through each input ID and looks for bam and spikein bam log from bowtie2
  # It assumes a sample name format of [sample]_[celltype]_[antibody].log or .spikein.log
  hist_info   <- strsplit(gsub("\\..*", "", basename(id)), "_")[[1]]
  align_summ = data.frame(Sample = hist_info[1], 
                          Cell_type = hist_info[2], 
                          Histone = hist_info[3])
  
  bam_path      <- file.path(log_path, paste0(id, ".log"))
  spikein_path  <- file.path(log_path, paste0(id, ".spikein.log"))
    
  # Combine sample ID with the alignment metrics from bowtie logs
  unlist(cbind(align_summ, 
               genAlignMetric(bam_path, 'hg38'),
               genAlignMetric(spikein_path, 'spikein')[,-1]))
}))
summ <- as.data.frame(summ)

# Calculate scaling factor as c/mapped_fragments_to_spikein_genome
summ$scaleFactor <- round(C/as.numeric(summ$MappedFragNum_spikein),2)

# Output table
write.table(summ, file=file.path(out_path, "alignment_metrics.tsv"),
            sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)