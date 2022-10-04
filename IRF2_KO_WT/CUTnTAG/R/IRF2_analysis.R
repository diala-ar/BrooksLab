library(lsa)
library(GenomicRanges)
library(dplyr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(ChIPseeker)
library(org.Mm.eg.db)
library(GenomicFeatures)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(msigdbr)

pdir <- '/cluster/projects/mcgahalab/data/mcgahalab/sabelo_irf2/weiguo_cutandtag_july29'
outdir <- file.path(pdir, "results", "manual", "filtered")
gtf_file <- '/cluster/projects/mcgahalab/ref/genomes/mouse/GRCm38/GTF/genome.gtf'
gbuild <- 'GRCm38'
species <- 'Mus musculus'
genome <- org.Mm.eg.db
goi <- c('Cd274', 'Pdcd1', 'Lag3', 'Ctla4', 'Tox', 'Tlr3')

ext_range <- 100  # BP to add to a peak to extend the peak search overlap
quantile_vals <- c(seq(0.1, 0.8, by=0.1), 
                   seq(0.81, 0.9, by=0.01),
                   seq(0.905, 1, by=0.005))


# A directory containing the output bed files from getPredictedTFSites.R
motif_dir <- file.path(pdir, "ref")

# Path to the seacr_single files
peaks_dir <- file.path(pdir, "results", "peaks", "seacr_single")

#######################
#### Functions Git ####
source("~/git/mini_projects/mini_functions/iterateMsigdb.R")

###################
#### Functions ####
# An expansion to the ChIPseeker::annotatePeak() function that will annotate
# all genes that overlap a TSS
annotateMultipleTSS <- function(peak, TxDb, level='gene'){
  # Copy annotatePeak() function for select features
  TxDb <- ChIPseeker:::loadTxDb(TxDb)
  if (level == "transcript") {
    features <- ChIPseeker:::getGene(TxDb, by = "transcript")
  }else {
    features <- ChIPseeker:::getGene(TxDb, by = "gene")
  }
  
  # subset for just the TSS start site
  features_tss <- resize(features, width=1) # faster
  
  # Find overlaps between peaks and TSS
  hit <- findOverlaps(peak, unstrand(features_tss)) %>%
    as.data.frame %>% 
    mutate(genes=features[subjectHits]$gene_id)
  
  # Append a comma-separated listed of all genes overlapping each TSS
  hit_qh <- unique(hit$queryHits)
  peak$multipleTSS <- NA
  peak[hit_qh]$multipleTSS <- sapply(split(hit$genes, hit$queryHits), paste, collapse=",")
  
  return(peak)
}

# Function to process batches of matrix-motifs in S3 and run a simple functino,
# then melt the data-structure into percentiles
meltPropRatio <- function(prop_l, fun, id){
  prop_ratio_trans <- sapply(prop_l, fun) %>% t()
  melt_prop_ratio <- prop_ratio_trans %>%
    as.data.frame() %>%
    rename_with(., ~ gsub("%", "", colnames(prop_ratio_trans))) %>%
    t() %>% melt() %>%
    setNames(., c('percentile', 'category', id))
  return(melt_prop_ratio)
}

# Reads in the peak bed files and formats it into GR objects
readPeak <- function(pathto, ext_range=0){
  grpeak <- setNames(read.table(pathto, sep="\t"),
                     c("chr", "start_raw", "end_raw", "AUC", "max", "range")) %>%
    mutate(start=start_raw - ext_range,
           end=end_raw + ext_range) %>%
    makeGRangesFromDataFrame(., keep.extra.columns = T)
  seqlevelsStyle(grpeak) <- 'UCSC'
  grpeak
}

# merge overlapping intervals if they overlap with a minimum length of one of the intervals
mergeOverlapping <- function(x, y, minfrac=0.6) {
  x <- granges(x)
  y <- granges(y)
  ov <- findOverlaps(x, y)
  xhits <- x[queryHits(ov)]
  yhits <- y[subjectHits(ov)]
  
  frac <- width(pintersect(xhits, yhits)) / pmin(width(xhits), width(yhits))
  merge <- frac >= minfrac
  c(reduce(c(xhits[merge], yhits[merge])),
    xhits[!merge], yhits[!merge],
    x[-queryHits(ov)], y[-subjectHits(ov)])
}

# Populate the target-control signal matrix with the normalized AUC-signal
# per catalogue peak ((AUC / width of peak) * scaling_factor)
populateCatalogue <- function(cnr_catalogue, ov_idx, gr_target, gr_ctrl){
  targ_ctrl_df <- matrix(nrow=length(cnr_catalogue), ncol=2) %>%
    as.data.frame() %>%
    setNames(., c("target", "control"))
  # Create a matrix of motif-significance per catalogue-peak
  motif_mat <- cnr_catalogue@elementMetadata[,names(ov_idx)]
  targ_ctrl_df <- as.data.frame(cbind(targ_ctrl_df, motif_mat))
  
  ## scale signal
  scale_signal <- function(signal, width, c=10) {
    round((signal / width) * c,3) %>%
      mean
  }
  
  ## target overlap
  ov <- findOverlaps(cnr_catalogue, gr_target)
  dup <- unique(queryHits(ov)[which(duplicated(queryHits(ov)))])
  ov_single <- ov[which(! queryHits(ov) %in% dup),]
  targ_ctrl_df[queryHits(ov_single),1] <- scale_signal(signal=sum(gr_target[subjectHits(ov_single)]$AUC), 
                                                       width=sum(width(gr_target[subjectHits(ov_single)])))
  ov_dup <- ov[which(queryHits(ov) %in% dup),]
  ov_dupl <- split(ov_dup, queryHits(ov_dup))
  targ_ctrl_df[unique(queryHits(ov_dup)),1] <- sapply(ov_dupl, function(ov_d){
    scale_signal(signal=gr_target[subjectHits(ov_d)]$AUC,
                 width=width(gr_target[subjectHits(ov_d)]))
  }) %>% as.numeric()
  
  ## control overlap
  ov <- findOverlaps(cnr_catalogue, gr_ctrl)
  dup <- unique(queryHits(ov)[which(duplicated(queryHits(ov)))])
  ov_single <- ov[which(! queryHits(ov) %in% dup),]
  targ_ctrl_df[queryHits(ov_single),2] <- scale_signal(signal=sum(gr_ctrl[subjectHits(ov_single)]$AUC), 
                                                       width=sum(width(gr_ctrl[subjectHits(ov_single)])))
  ov_dup <- ov[which(queryHits(ov) %in% dup),]
  ov_dupl <- split(ov_dup, queryHits(ov_dup))
  ov_targ_dupl_val <- sapply(ov_dupl, function(ov_d){
    scale_signal(signal=sum(gr_ctrl[subjectHits(ov_d)]$AUC),
                 width=sum(width(gr_ctrl[subjectHits(ov_d)])))
  }) %>% as.integer()
  targ_ctrl_df[unique(queryHits(ov_dup)),2]  <- ov_targ_dupl_val
  
  return(targ_ctrl_df)
}

##################
#### 0. Setup ####
# Import GTF file
txdb <- makeTxDbFromGFF(file = gtf_file, format = "gtf")
txby <- keys(genome, 'ENSEMBL')
ens2sym_ids <- mapIds(genome, keys=txby, column='SYMBOL',
                      keytype='ENSEMBL', multiVals="first")
sym2ens_ids <- setNames(names(ens2sym_ids), ens2sym_ids)
ens2entrez_ids <- mapIds(genome, keys=txby, column='ENTREZID',
                         keytype='ENSEMBL', multiVals="first")
entrez2ens_ids <- setNames(names(ens2entrez_ids), ens2entrez_ids)

# Read in IRF2 motif matrix
gr_motifs <- lapply(list.files(motif_dir, pattern=".bed"), function(f){
  motif <- setNames(read.table(file.path(motif_dir, f), sep="\t", check.names = F),
                    c("chr", "start", "end", "strand", "sig")) %>%
    makeGRangesFromDataFrame(., keep.extra.columns = T)
  seqlevelsStyle(motif) <- 'UCSC'
  motif
})
names(gr_motifs) <- gsub("_.*", "", list.files(motif_dir, pattern=".bed"))

# Generate Cut&Run peak catalogue - based off of Unfixed IRF2 
cnr_catalogue <- sapply(list.files(peaks_dir, pattern=".stringent.bed"), function(pfile){
  readPeak(file.path(peaks_dir, pfile),
           ext_range=ext_range)
})
cnr_catalogue_bkup <- cnr_catalogue
select_samples <- grep("IRF2", names(cnr_catalogue_bkup), value=T) %>%
  grep("Naive", ., invert = T, ignore.case = T, value=T) %>%
  grep("Act_fix", ., invert = T, ignore.case = T, value=T) %>%
  grep("-[123]", ., invert = T, ignore.case = T, value=T)
cnr_catalogue <- lapply(cnr_catalogue_bkup[select_samples], function(gr){
  # Set a min Signal_AUC threshold to consider a peak
  min_auc <- quantile(gr$AUC, 0.95)
  gr[gr$AUC > min_auc]
})
cnr_catalogue <- as(cnr_catalogue, "GRangesList") %>%
  unlist() %>%
  reduce()
seqlevelsStyle(cnr_catalogue) <- 'NCBI'
cnr_catalogue <- annotatePeak(peak=cnr_catalogue, TxDb=txdb, level='gene',
                              tssRegion=c(-3000, 3000), verbose=FALSE)@anno
cnr_catalogue <- annotateMultipleTSS(peak=cnr_catalogue, TxDb=txdb, level='gene')
cnr_catalogue$symbol <- ens2sym_ids[cnr_catalogue$geneId]
cnr_catalogue$entrez <- ens2entrez_ids[cnr_catalogue$geneId]
multi_tss_idx <- grep(",", cnr_catalogue$multipleTSS)
sym_ids <- sapply(cnr_catalogue[multi_tss_idx,]$multipleTSS, function(ens_ids){
  ens2sym_ids[strsplit(ens_ids, ",")[[1]]] %>% 
    paste(., collapse=",")
})
cnr_catalogue$symbol[multi_tss_idx] <- sym_ids
seqlevelsStyle(cnr_catalogue) <- 'UCSC'

# Annotate the peak-catalogue with motif overlap
catalogue_motifs <- sapply(gr_motifs, function(gr0){
  ov <- findOverlaps(gr0, cnr_catalogue)
  motif_bool <- rep(0, length(cnr_catalogue))
  motif_bool[subjectHits(ov)] <- round(gr0$sig[queryHits(ov)],2)
  return(motif_bool)
})
cnr_catalogue@elementMetadata <- cbind(cnr_catalogue@elementMetadata, 
                                       catalogue_motifs)

# Create a mapping between catalogue and IRF2 motifs
motifs_ov <- lapply(gr_motifs, findOverlaps, query=cnr_catalogue)
irf2_cat_ov <- motifs_ov$IRF2

# Create listing of TARGET - CONTROL pair of files
f_pairs <- list.files(peaks_dir, pattern=".stringent.bed")
f_grps <- gsub("_[a-zA-Z0-9]*\\..*", "", f_pairs) %>%
  gsub("_IRF2.*", "", .)
f_pairs <- split(f_pairs, f_grps)
f_pairs <- lapply(f_pairs, function(fp){
  data.frame("target"=grep("igg", fp, ignore.case = T, invert = T, value = T)) %>%
    mutate(control=grep("igg", fp, ignore.case = T, value=T))
}) %>%
  do.call(rbind, .)
rownames(f_pairs) <- gsub(".stringent.*", "", f_pairs$target)

f_pairs <- f_pairs %>% 
  filter(!grepl("-[123]", target))



###################################
#### 1. IRF2 Overlapping Peaks ####
dist <- sapply(list.files(peaks_dir, pattern=".stringent.bed"), function(pfile){
  grpeak <- readPeak(file.path(peaks_dir, pfile),
                     ext_range=ext_range)
  
  quantiles <- quantile(grpeak$AUC, quantile_vals)
  sapply(quantiles, function(q){
    grpeak_q <- grpeak[grpeak$AUC > q,]
    ov_idx <- findOverlaps(grpeak_q, gr_motifs$IRF2)
    length(unique(queryHits(ov_idx))) / length(grpeak_q)
  })
})

dist2 <- round(dist, 3)
colnames(dist2) <- gsub(".stringent.*", "", colnames(dist2))

mdist2 <- melt(dist2) %>%
  rename_with(., ~ c("signalAUC_percentile", "sample", "IRF2_proportion")) %>%
  mutate(signalAUC_percentile = gsub("%", "", signalAUC_percentile),
         group = gsub("^.*_", "", sample),
         marker = gsub("(.*_.*)_.*", "\\1", sample))
mdist2$signalAUC_percentile <- as.numeric(as.character(mdist2$signalAUC_percentile))
mdist2$group <- gsub("-44$", "", mdist2$group)
pdf("~/xfer/irf2_motif_peak.all.pdf", width = 10)
ggplot(mdist2, aes(x=signalAUC_percentile, IRF2_proportion, color=marker, group=sample)) +
  geom_line() + 
  facet_grid(rows=vars(group)) +
  theme_classic() +
  theme(axis.text.x=element_text(size=6, angle=45))
dev.off()

############################################
#### 2. Cosine similarity between peaks ####
# Fill in the Signal matrix for all instances of the Peaks catalogue
signal_mat <- sapply(list.files(peaks_dir, pattern=".stringent.bed"), function(pfile){
  gr0 <- readPeak(file.path(peaks_dir, pfile), ext_range = 0)
  
  setNames(read.table(file.path(peaks_dir, pfile), sep="\t"),
                  c("chr", "start", "end", "AUC", "max", "range")) %>%
    makeGRangesFromDataFrame(., keep.extra.columns = T)
  signal_vec <- matrix(nrow=length(cnr_catalogue), ncol=1)
  
  ov <- findOverlaps(cnr_catalogue, gr0)
  signal_vec[queryHits(ov),] <- round(gr0[subjectHits(ov)]$AUC,2)
  return(signal_vec)
})

# Generate a cosine similarity matrix between the catalogue
jaccard <- function(a, b) {
  sum(!is.na(a) & !is.na(b)) / sum(!is.na(a) | !is.na(b))
}

metric <- 'cosine'
signal_mat[is.na(signal_mat)] <- 0
colnames(signal_mat) <- gsub(".stringent.*", "", colnames(signal_mat))
sim_mats <- lapply(quantile_vals, function(min_percentile){
  sim_mat <- apply(signal_mat, 2, function(i){
    i[which(i < quantile(i[i!=0], min_percentile))] <- switch(metric, 
                                                              jaccard=NA,
                                                              cosine=0)
    apply(signal_mat, 2, function(j){
      j[which(j < quantile(j[j!=0], min_percentile))] <- switch(metric, 
                                                                jaccard=NA,
                                                                cosine=0)
      if(metric == 'cosine'){
        lsa::cosine(i,j)
      } else {
        jaccard(i,j)
      }
    })
  })
  return(sim_mat)
})
names(sim_mats) <- quantile_vals
  
melt_sims <- melt(sim_mats) %>%
  mutate(L1 = gsub("%", "", L1),
         target_j = gsub("^.*_", "", Var2),
         sample = gsub("(.*_.*)_.*", "\\1", Var2)) %>%
  rename_with(., ~ c("sample_i", "sample_j", "cosine_similarity", 
                     "percentile", "target", "sample")) %>% 
  filter(., as.numeric(cosine_similarity) < 0.99)
melt_sims$percentile <- as.numeric(as.character(melt_sims$percentile)) 
pdf("~/xfer/irf2_sim_matrix.pdf", height = 12)
ggplot(melt_sims, aes(x=percentile, y=cosine_similarity, 
                      group=sample_j, col=target, linetype=sample)) +
  facet_grid(rows=vars(sample_i), scales="free") +
  ylim(0,0.5) +
  geom_line()
dev.off()

#######################################################
#### 3. Target/Control overlap and Motif breakdown ####
# cnr_catalogue   # peak-catalogue across all samples
# grirf2          # interval catalogue for IRF2 motifs

dir.create(file.path(outdir, "rds"), showWarnings = F)
irf2_props <- apply(f_pairs, 1, function(pfiles){
  id <- paste(gsub(".stringent.*", "", pfiles), collapse="_")
  targctrl_rds <- paste0("targctrl.", id, ".rds")
  rds <- paste0("paired.", id, ".rds")
  print(paste0(id, "..."))
  if(file.exists(file.path(outdir, "rds", rds))){
    list_all <- readRDS(file.path(outdir, "rds", rds))
    return(list_all)
  }
  
  gr_target <- readPeak(file.path(peaks_dir, pfiles['target']), ext_range = 0)
  gr_ctrl <- readPeak(file.path(peaks_dir, pfiles['control']), ext_range = 0)
  
  # Populate the target-control signal matrix
  if(file.exists(file.path(outdir, "rds", targctrl_rds))){
    print("Reading in existing target-cotnrol signal matrix...")
    targ_ctrl_df <- readRDS(file.path(outdir, "rds", targctrl_rds))
  } else {
    print("Populating target-control signal matrix...")
    targ_ctrl_df <- populateCatalogue(cnr_catalogue, motifs_ov, gr_target, gr_ctrl)
    saveRDS(targ_ctrl_df, file.path(outdir, "rds", targctrl_rds))
  }
  
  ## Create percentile-cutoffs for each sample to look at only peaks exceeding
  # a set signal value
  quantiles <- apply(targ_ctrl_df[,c('target', 'control')], 2, quantile, 
                     probs=quantile_vals[-length(quantile_vals)], na.rm=T)
  irf2_prop <- apply(quantiles, 1, function(qs){
    # Remove peaks that don't meet the threshold
    targ_ctrl_mat <- as.matrix(targ_ctrl_df[,c(1,2)])
    thresh_mat <- matrix(qs, byrow = T, 
                         nrow=nrow(targ_ctrl_mat), ncol=ncol(targ_ctrl_mat))
    thresh_mat_lowCtrl <- thresh_mat
    targ_ctrl_mat[targ_ctrl_mat < thresh_mat] <- NA
    
    # Keep low signal control peaks to compare against high-signal target peaks
    targ_ctrl_mat_lowCtrl <-  as.matrix(targ_ctrl_df[,c(1,2)])
    thresh_mat_lowCtrl[,2] <- 0.1
    targ_ctrl_mat_lowCtrl[targ_ctrl_mat_lowCtrl < thresh_mat_lowCtrl] <- NA
    
    # Subset for peaks that are found only in target, or target+control
    keep_idx <- rowSums(is.na(targ_ctrl_mat))
    keep_idx <- setNames(split(c(1:nrow(targ_ctrl_mat)), f=keep_idx),
                         c('bothHigh', 'one', 'none'))
    
    # Subset for peaks that are found both in thresholded-target and no-threshold control
    keep_idx_lowCtrl <- rowSums(is.na(targ_ctrl_mat_lowCtrl))
    keep_idx_lowCtrl <- setNames(split(c(1:nrow(targ_ctrl_mat)), f=keep_idx_lowCtrl),
                                 c('bothLow', 'one', 'none'))['bothLow']
    
    # Subset for peaks that are found in target-only or control-only
    ss_keep_idx <- is.na(targ_ctrl_mat[keep_idx[['one']],1])
    ss_keep_idx <- setNames(split(keep_idx[['one']], ss_keep_idx),
                            c('trgPk', 'ctrlPk'))
    keep_idx2 <- c(keep_idx['bothHigh'], keep_idx_lowCtrl['bothLow'],
                   ss_keep_idx[c('trgPk', 'ctrlPk')])
    
    # Calculate proportion of motif+ peaks within each subgroup
    sapply(keep_idx2, function(kidx){
      colSums(targ_ctrl_df[kidx,-c(1,2)] != 0) / length(kidx)
      # sum(!is.na(targ_ctrl_df$irf2[kidx])) / length(kidx)
    }) %>% 
      as.data.frame() %>%
      unlist()
  })
  irf2_prop_ratio <- apply(irf2_prop, 2, function(i) i / irf2_prop[,1]) %>%
    as.data.frame() %>%
    split(., f=gsub("[0-9]*", "", rownames(.)))
  colSd <- function(i) apply(i, 2, sd, na.rm=T)
  prop_ratio_melt <- irf2_prop_ratio %>%
    meltPropRatio(., colMeans, 'motif_mean') %>%
    full_join(irf2_prop_ratio %>%
                meltPropRatio(., colSd, 'motif_sd')) %>%
    mutate(log2_motif_mean = log2(motif_mean),
           log2_motif_low = log2(motif_mean - (2*motif_sd)),
           log2_motif_hi = log2(motif_mean + (2*motif_sd)),
           percentile = percentile / 100)
  
  targ_ctrl_df <- targ_ctrl_df %>%
    mutate(quantile_target = round(ecdf(target)(target),2),
           quantile_control = round(ecdf(control)(control),2),
           gene = cnr_catalogue$symbol,
           annotation = cnr_catalogue$annotation,
           loc=paste0(seqnames(cnr_catalogue), ":",
                      start(cnr_catalogue), "-",
                      end(cnr_catalogue)))
  
  gg <- ggplot(prop_ratio_melt, aes(x=percentile, y=log2_motif_mean, 
                                    group=category, col=category, fill=category)) + 
    geom_line() +
    geom_ribbon(aes(ymin=log2_motif_low, ymax=log2_motif_hi), alpha=0.1, colour = NA) +
    ggtitle(paste(gsub(".stringent.*", "", pfiles), collapse=" - ")) + 
    theme_classic() + ylim(-3, 3) + xlim(0,1)

  if(!is.null(goi)){
    idx <- which((targ_ctrl_df$gene %in% goi) &
                   grepl("Promoter.*1kb", targ_ctrl_df$annotation))
    goi_df <- targ_ctrl_df[idx,,drop=F] %>%
      filter(., !is.na(target)) %>%
      melt() %>%
      filter(., grepl("quantile", variable)) %>%
      mutate(variable = gsub("quantile_", "", variable))
    goi_df$variable <- factor(goi_df$variable, levels=c("target", "control"))
    gg_point <- ggplot(goi_df, aes(x=value, y=variable, col=variable, fill=variable)) +
      facet_grid(rows=vars(gene), scales='free', switch='y') +
      geom_point(alpha=0.5) +
      theme_classic() + xlim(0,1) +
      xlab("percentile") +
      theme(strip.text.y.left = element_text(angle = 0),
            strip.placement = 'outside',
            axis.text.y = element_blank(),
            axis.title.y = element_blank(), 
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            panel.background = element_rect(fill = 'grey90'),
            strip.background = element_rect(colour="white", fill="white"),
            panel.spacing.y = unit(0.1, "line"))
    gg_all <- plot_grid(gg + xlab("") + theme(axis.text.x = element_blank()), 
                        gg_point, ncol=1, rel_heights=c(4,1), 
                        align='v', axis='lr')
  }
  
  list_all <- list("gg"=gg, "gg_point"=gg_point, "gg_all"=gg_all, 
                   "mat"=irf2_prop, "df"=targ_ctrl_df, goi=goi)
  saveRDS(list_all, file=file.path(outdir, "rds", rds))

  return(list_all)
})

pdf("~/xfer/target_control_irf2.pdf"); 
lapply(irf2_props, function(id){ id$gg_all})
dev.off()

#####################################################################
#### 4. Target/Control identifying thresholds for IRF2 detection ####
# cnr_catalogue   # Sec2: peak-catalogue across all samples
# grirf2          # Sec1: interval catalogue for IRF2 motifs
# irf2_cat_ov     # Sec3: Overlap between Peak-catalogue and IRF2 motifs
# f_pairs         # Sec3: list of all target-control sample pairs
irf2_props_mat <- apply(f_pairs, 1, function(pfiles){
  id <- paste(gsub(".stringent.*", "", pfiles), collapse="_")
  targctrl_rds <- paste0("targctrl.", id, ".rds")
  rds <- paste0("propmat.", id, ".rds")
  print(paste0(id, "..."))
  if(file.exists(file.path(outdir, "rds", rds))){
    list_all <- readRDS(file.path(outdir, "rds", rds))
    return(list_all)
  }
  
  gr_target <- readPeak(file.path(peaks_dir, pfiles['target']), ext_range = 0)
  gr_ctrl <- readPeak(file.path(peaks_dir, pfiles['control']), ext_range = 0)
  
  # Populate the target-control signal matrix
  if(file.exists(file.path(outdir, "rds", targctrl_rds))){
    targ_ctrl_df <- readRDS(file=file.path(outdir, "rds", targctrl_rds))
  } else {
    targ_ctrl_df <- populateCatalogue(cnr_catalogue, motifs_ov, gr_target, gr_ctrl)
    saveRDS(targ_ctrl_df, file=file.path(outdir, "rds", targctrl_rds))
  }
  
  # Create index of peaks where both Target and Control overlap
  tc_ov <- as.matrix(targ_ctrl_df[,1:2]) %>%
    is.na() %>%
    rowSums()
  tc_ov_idx <- split(c(1:nrow(targ_ctrl_df)), f=tc_ov)[['0']]
  
  ## Create percentile-cutoffs for each sample to look at only peaks exceeding
  # a set signal value
  quantile_vals <- seq(0.01, 1, by=0.01)
  quantiles <- apply(targ_ctrl_df[,c('target', 'control')], 2, quantile, 
                     probs=quantile_vals[-length(quantile_vals)], na.rm=T) %>%
    as.data.frame()
  # Establish base distribution
  targ <- targ_ctrl_df$target[tc_ov_idx]
  ctrl <- targ_ctrl_df$control[tc_ov_idx]
  base_prop <- apply(quantiles[seq(1,31, by=5),], 1, function(q){
    targ[targ < q[1]] <- NA
    ctrl[ctrl < q[2]] <- NA
    both_idx <- which(!is.na(targ) & !is.na(ctrl))
    colSums(targ_ctrl_df[tc_ov_idx,][both_idx,-c(1,2)] != 0) / length(both_idx)
  })
  
  
  tstat <- function(x, x0=0){
    s <- (mean(x) - x0) / (sd(x) / sqrt(length(x)))
    if(is.nan(s)) s<- 0
    s
  }
  prop_mat <- sapply(setNames(quantiles$target, rownames(quantiles)), function(q_t){
    targ[targ < q_t] <- NA
    sapply(setNames(quantiles$control, rownames(quantiles)), function(q_c){
      ctrl[ctrl < q_c] <- NA
      both_idx <- which(!is.na(targ) & !is.na(ctrl))
      tc_prop <- colSums(targ_ctrl_df[tc_ov_idx,][both_idx,-c(1,2)] != 0) / length(both_idx)
      base_dist <- apply(base_prop, 2, function(i) i - base_prop) %>%
        as.numeric()
      base_dist <- base_dist[base_dist!=0]
      
      tval <- apply(base_prop, 2, function(i) tc_prop - i) %>%
        as.numeric() %>% 
        t.test(., base_dist)
      return(tval$statistic)
    })
  })
  
  list_all <- list("matrix"=prop_mat, "quantiles"=quantiles)
  saveRDS(list_all, file=file.path(outdir, "rds", rds))
  return(list_all)
})
names(irf2_props_mat) <- apply(f_pairs, 1, paste, collapse=" - ") %>%
  gsub(".stringent[a-zA-Z0-9.]*", "", .)

pdf("~/xfer/irf2_prop_matrix.pdf")
rng <- c(-1, 25)
lapply(names(irf2_props_mat), function(id){
  melt(irf2_props_mat[[id]]$matrix) %>%
    mutate(Var1 = gsub("%", "", Var1) %>% gsub(".t", "", .) %>% as.numeric(),
           Var2 = gsub("%", "", Var2) %>% as.numeric(),
           value = value %>% 
             replace(., . > max(rng), max(rng)) %>%
             replace(., . < min(rng), min(rng))) %>%
  ggplot(., aes(Var1, Var2)) +
    # geom_raster(aes(fill = value), interpolate = FALSE) +
    geom_tile(aes(fill=value)) +
    ggtitle(id) + ylab("Signal percentile threshold (target)") + 
    xlab("Signal percentile threshold (control)") +
    theme(axis.text.x = element_text(angle=45, size=6),
          axis.text.y = element_text(size=6)) +
    scale_fill_gradientn(colours=c("black", "grey", "yellow", "red")) 
})
dev.off()

#####################################
#### 5. Refining confident peaks ####
write_tables <- FALSE
thresholds <- list("Act_unfixed_IRF2"=c("target"=0.69, "control"=0.86),
                   "Act_unfixed_H3K4me3"=c("target"=0.54, "control"=0.87),
                   "CnT_Act_IRF2-44"=c('target'=0.9, 'control'=0.85),
                   # "AC_Fx_IRF2-44"=c('target'=, 'control'=),
                   "AC_uF_IRF2-44"=c('target'=0.81, 'control'=0.82))#,
                   # "AC_Fx_IRF2-58"=c('target'=, 'control'=),
                   # "AC_uF_IRF2-58"=c('target'=, 'control'=))
# thresholds <- setNames(c(0.985, 0.985, 0.915), 
#                        paste0("Act_unfixed_", c("IRF2", "H3K4me3", "IgG"), ".stringent.bed"))
grs_anno <- apply(f_pairs[names(thresholds)[4],,drop=F], 1, function(pfiles){
  id <- paste(gsub(".stringent.*", "", pfiles), collapse="_")
  targctrl_rds <- paste0("targctrl2.", id, ".rds")  # targctrl2 is equivalent to catalogue of IRF2 peaks, targctrl is all IRF2+H3K4me3+IgG
  print(paste0(id, "..."))
  
  gr_target <- readPeak(file.path(peaks_dir, pfiles['target']), ext_range = 0)
  gr_ctrl <- readPeak(file.path(peaks_dir, pfiles['control']), ext_range = 0)
  target_id <- gsub("\\..*", "", pfiles['target'])
  targ_thresh <- thresholds[[target_id]]['target']
  ctrl_thresh <- thresholds[[target_id]]['control']
  
  # Populate the target-control signal matrix
  if(file.exists(file.path(outdir, "rds", targctrl_rds))){
    targ_ctrl_df <- readRDS(file.path(outdir, "rds", targctrl_rds))
  } else {
    targ_ctrl_df <- populateCatalogue(cnr_catalogue, motifs_ov, gr_target, gr_ctrl)
    saveRDS(targ_ctrl_df, file.path(outdir, "rds", targctrl_rds))
  }
  targ_ctrl_df$t_percentile <- ecdf(targ_ctrl_df$target)(targ_ctrl_df$target)
  targ_ctrl_df$c_percentile <- ecdf(targ_ctrl_df$control)(targ_ctrl_df$control)
  
  # Subset for peaks that are found only in target, or target+control
  targ_ctrl_mat <- as.matrix(targ_ctrl_df[,c(1,2)])
  keep_idx <- rowSums(is.na(targ_ctrl_mat))
  keep_idx <- setNames(split(c(1:nrow(targ_ctrl_mat)), f=keep_idx),
                       c('both', 'one', 'none'))
  
  # Remove peaks that don't meet the threshold
  thresh_mat <- matrix(c(quantile(targ_ctrl_df$target, targ_thresh, na.rm=T),
                         quantile(targ_ctrl_df$control, ctrl_thresh, na.rm=T)), 
                       byrow = T, nrow=nrow(targ_ctrl_mat), ncol=ncol(targ_ctrl_mat))
  targ_ctrl_mat[targ_ctrl_mat[,1] < thresh_mat[,1],1] <- NA
  targ_ctrl_mat[targ_ctrl_mat[,2] >= thresh_mat[,2],2] <- NA

  # Identify peaks that are found  in targetOnly
  idx <- list('both'=intersect(which(rowSums(!is.na(targ_ctrl_mat)) == 2), keep_idx[['both']]), 
              'trgPk'=intersect(which(!is.na(targ_ctrl_mat[,1])), keep_idx[['one']]) )
  
  ## Annotate the peaks with their gene and genomic annotation
  grs <- lapply(names(idx), function(idx_id){
    idx_i <- idx[[idx_id]]
    gr0 <- cnr_catalogue[idx_i]
    gr0@elementMetadata <- cbind(gr0@elementMetadata,
                                 as(setNames(targ_ctrl_df[idx_i, c('target', 't_percentile',
                                                                   'control', 'c_percentile')],
                                       c('targetSignal', 'targetPercentile', 
                                         'ctrlSignal', 'ctrlPercentile')), 
                                    "DataFrame"))
    gr0$classification <- idx_id
    return(gr0)
  })
  names(grs) <- names(idx)
  
  return(grs)
})
grs_anno[[1]]$trgPk[grep("Ctla4", grs_anno[[1]]$trgPk$symbol, ignore.case=T),]
sapply(grs_anno[[1]], function(gr) gr[grep("Promoter", gr$annotation),] %>% length) %>% sum
saveRDS(grs_anno, file=file.path(outdir, "grs_anno.rds"))

sapply(grs_anno[[1]], function(i) i[grep("Promoter", i$annotation),] %>% length) %>% sum
# Writing the annotated peaks to a file
if(write_tables){
  lapply(names(grs_anno), function(id){
    grl <- unlist(grs_anno[[id]], recursive=F)
    lapply(names(grl), function(id_i){
      as.data.frame(grl[[id_i]]) %>%
        mutate(annotation=gsub(",", "-", annotation)) %>%
        write.table(., file=paste0("~/xfer/", id, ".", gsub("\\.", "_", id_i), ".tsv"),
                    sep="\t", row.names = F, col.names = T, quote = F)
      paste0(id, ".", gsub("\\.", "_", id_i), ".tsv")
    })
  })
}



###########################################################
#### 5.b) Get List of genes within the confident peaks ####
write_tables <- FALSE
grs_anno <- readRDS(file.path(outdir, "grs_anno.rds"))
id='AC_uF_IRF2-44'

.getAllGenes <- function(grl, motif=NULL, direction='greater'){
  x <- unlist(sapply(grl, function(gr0) {
    if(!is.null(motif)){
      gr0 <- switch(direction,
                    equal=gr0[which(gr0@elementMetadata[,motif]==0),],
                    greater=gr0[which(gr0@elementMetadata[,motif]>0),])
    } 
    gr0[grepl("Promoter", gr0$annotation),]$symbol
  }))
  x <- x[!is.na(x)]
  as.character(unlist(strsplit(x, split=",")))
}
genes <- list()
genes[['all']] <- .getAllGenes(grs_anno[[id]], motif=NULL)
genes[['irf2']] <- .getAllGenes(grs_anno[[id]], motif='IRF2', direction='greater')
genes[['non_irf2']] <- .getAllGenes(grs_anno[[id]], motif='IRF2', direction='equal')


# Write out the confident genes
if(write_tables){
  write.table(genes$all, file=file.path("~/xfer", "irf2_genes.csv"), sep=",",
              col.names = F, row.names = F, quote = F)
}


# Write out a bed file for the peaks
bedfile <- lapply(grs_anno[[id]], function(gr0){
  gr0[grep("Promoter", gr0$annotation),] %>%
    as.data.frame %>% 
    dplyr::select(seqnames, start, end)
}) %>%
  do.call(rbind, .) %>% 
  arrange(seqnames, start) %>%
  mutate(seqnames=gsub("chr", "", seqnames))
if(write_tables){
  write.table(bedfile, file=file.path(outdir, paste0(id, ".bed")),
              sep="\t", col.names = F, row.names = F, quote = F)
}


###########################################
#### 5.c) Over representation analysis ####
write_tables <- FALSE
grs_anno <- readRDS(file.path(outdir, "grs_anno.rds"))
genome_gse <- org.Mm.eg.db
txby <- keys(genome_gse, 'SYMBOL')
sym2entrez_ids <- mapIds(genome_gse, keys=txby, column='ENTREZID',
                         keytype='SYMBOL', multiVals="first")

species="Mus musculus"
all_gene_sets = msigdbr(species = species)
genesets <- grep("exhaust", all_gene_sets$gs_name, ignore.case = T, value = T) %>%
  sort %>% unique

oras <- lapply(genes, function(x0){
  iterateMsigdb(species='Mus musculus', fun=oraFun, 
                msig_lvls=list("H"=list(NULL),
                               "C7"=list("IMMUNESIGDB")), #list("custom"=msig_ds),
                entrez_genes=sym2entrez_ids[x0])
})

.getGenesetsAndMelt <- function(oras2){
  .reshape_ora <- function(x){
    x %>% 
      mutate(logq=(-1*log(qvalue))) %>%
      select(ID, Count, logq) %>%
      mutate(Group=gsub("_.*", "", ID),
             ID=gsub("^.*?_", "", ID))
  }
  hallmark_ora <- oras2$H %>% .reshape_ora
  immune_ora <- oras2$C7[which(oras2$C7$ID %in% genesets),] %>% 
    .reshape_ora %>% mutate(Group='IMMUNESIGDB')
  
  melt_ora <- as.data.frame(rbind(hallmark_ora, immune_ora))
  melt_ora$ID <- factor(melt_ora$ID, levels=melt_ora$ID)
  return(melt_ora)
}

oras_all <- lapply(oras$all, summarizeOra, keep_genes=FALSE,  qcutoff=0.05)
oras_irf2 <- lapply(oras$irf2, summarizeOra, keep_genes=FALSE,  qcutoff=1)
oras_nonirf2 <- lapply(oras$non_irf2, summarizeOra, keep_genes=FALSE,  qcutoff=1)
melt_ora <- .getGenesetsAndMelt(oras_all)  %>% mutate(subset="All")
melt_ora_irf2 <- .getGenesetsAndMelt(oras_irf2) %>% 
  mutate(subset="IRF2") %>%
  filter(ID %in% melt_ora$ID)
melt_ora_nonirf2 <- .getGenesetsAndMelt(oras_nonirf2)  %>% 
  mutate(subset="nonIRF2") %>%
  filter(ID %in% melt_ora$ID)

melt_ora <- rbind(melt_ora, melt_ora_irf2, melt_ora_nonirf2) %>%
  as.data.frame %>%
  mutate(sig=(logq > (-1*log(0.05))))

pdf("~/xfer/ora_pathways.pdf", height = 12, width = 6)
ggplot(melt_ora, aes(x=ID, y=logq, size=Count)) +
  geom_segment(aes(x = ID, xend = ID, y = 0, yend = logq, size=1)) +
  geom_point(aes(color=sig)) + 
  scale_color_manual(values=c('TRUE'='#1a9850', 'FALSE'='#4d4d4d')) +
  facet_grid(subset~Group, scales='free', space = 'free') +
  ylab("-log(q-value)") + xlab("") + ylim(0,50) +
  scale_size(range = c(0.5, 4)) +
  theme_classic() +
  theme(axis.text.x=element_text(angle=90, hjust = 1),
        strip.background = element_blank())
dev.off()

pdf("~/xfer/ora_pathways.paper.pdf", height = 2.5, width = 8)
melt_ora_paper <- melt_ora %>%
  filter(subset =='All',
         Group == 'HALLMARK' | grepl("DAY6_EFFECTOR_VS_DAY30", ID)) %>%
  arrange(logq) %>%
  mutate(ID = factor(as.character(ID), levels=as.character(ID)))
ggplot(melt_ora_paper, aes(x=ID, y=logq, size=Count)) +
  geom_segment(aes(x = ID, xend = ID, y = 0, yend = logq, size=1)) +
  geom_point() + 
  ylab("-log(q-value)") + xlab("") + ylim(0,40) +
  scale_size(range = c(0.5, 4)) +
  theme_classic() +
  coord_flip() +
  theme(strip.background = element_blank())
dev.off()
