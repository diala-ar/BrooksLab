
#Disable error messaging in v1.0, prevent early exit.
options(warn = - 1)

library(Seurat)
library(data.table)


loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

######################
######################


#%%%%%%%%%%%%%%%%%%%%%%%#
## Read Data from Java ##
#%%%%%%%%%%%%%%%%%%%%%%%#

runID <- ""

model <- "None"


inputString = snakemake@input[['seur']]
extenShaun  = strsplit(basename(inputString), split="\\.")[[1]][2]
out_dir     = snakemake@params[['out_dir']]

# Read CSV from SeqGeq
if (extenShaun == "rda") {
  import = loadRData(inputString)
} else if (extenShaun == "rds") {
  import <- readRDS(inputString)
} else if (extenShaun == "Robj"){
  import <- get(load(inputString, verbose=TRUE))
} else {
  writeLines("Could not load the data file specified. Make sure output from Seurat is in RDA or RDS format.")
}

## Deal with data like the Sade Feldman example:
writeLines("You've got a Seurat dataset")
## TODO - What's the cellID going to be?

###########################################
######## What is the raw Matrix?: #########

isNorm <- F
hasAb <- F
if (!is.null(import@assays$ADT)) { # if there is an ADT assay
    AB_counts = import@assays$ADT@counts
    AB_norm   = import@assays$ADT@data
    AB_counts = as.matrix(AB_counts)
    AB_norm   = as.matrix(AB_norm)
    hasAb     = T
}



if ( !is.null(import@assays[["RNA"]]) ) {
  dataMatrix <- import@assays[["RNA"]]@counts
  DM <- as.matrix(dataMatrix)
  if(hasAb){
    DM_counts <- cbind(t(AB_counts), t(DM))
  } else {
    DM_counts <- t(DM)
  }
  DMRows <- as.data.frame(rownames(DM_counts))
  colnames(DMRows) <- "CellId"
  DM_counts <- cbind(DMRows, DM_counts)

  dataMatrix <- import@assays[["RNA"]]@data
  DM <- as.matrix(dataMatrix)
  if(hasAb){
    DM_norm <- cbind(t(AB_norm), t(DM))
  } else {
    DM_norm <- t(DM)
  }
  DMRows <- as.data.frame(rownames(DM_norm))
  colnames(DMRows) <- "CellId"
  DM_norm <- cbind(DMRows, DM_norm)
}


###########################################
######## What are the DPs?: #########

DP <- c()
if ( !is.null(import@reductions) ) {
    dimRedux <- import@reductions

    if (!is.null(import@meta.data$orig.ident)){
        SampleId <- as.factor( import@meta.data$orig.ident )
        # SampleId <- import@meta.data$orig.ident
        DP <- cbind(DP, SampleId)
    }

    ###PCA###
    if ( exists("pca", where=dimRedux) ) {
        PC.DR <- dimRedux[["pca"]]@cell.embeddings
        pc.dr <- as.matrix(PC.DR)
        DP <- cbind(DP, pc.dr)
    }

    ###UMAP###
    if ( exists("wnn.umap", where=dimRedux) ) {
        U.DR <- dimRedux[["wnn.umap"]]@cell.embeddings
        u.dr <- as.matrix(U.DR)
        DP <- cbind(DP, u.dr)
    }
    if ( exists("umap", where=dimRedux) ) {
      U.DR <- dimRedux[["umap"]]@cell.embeddings
      u.dr <- as.matrix(U.DR)
      DP <- cbind(DP, u.dr)
    }
}

###CLUSTERS###
resol = snakemake@params[['resolution']]
cluster_by = ifelse('ADT' %in% Assays(import), paste0('wsnn_res.', resol), paste0('RNA_snn_res.', resol))
if ( !is.null(import@meta.data[, cluster_by]) ) {
    CLUST <- import@meta.data[, cluster_by]
    clust <- as.matrix(CLUST)
    colnames(clust) <- cluster_by
    DP <- cbind(DP, clust)
}

# write files
if (!is.null(DP)) {
  DM_counts = cbind(DM_counts, DP)
  fwrite(DM_counts, file=snakemake@output[['seur_counts']], na="0", quote=F, row.names=F)
  
  DM_norm = cbind(DM_norm, DP)
  fwrite(DM_norm, file=snakemake@output[['seur_norm_counts']], na="0", quote=F, row.names=F)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%#
##     Total Meta Data     ##
#%%%%%%%%%%%%%%%%%%%%%%%%%%%#

if( !is.null(import@meta.data) ){
  
  df <- import@meta.data
  fwrite(df, file=snakemake@output[['seur_metadata']], na="0", quote=F, row.names=F)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%#
## Export differentially expressed genes
#%%%%%%%%%%%%%%%%%%%%%%%%%%%#
degs = read.csv(snakemake@input[['degs']])
clusters_1 = unique(degs$cluster_1)
clusters_2 = unique(degs$cluster_2)
clusters   = paste0(clusters_1, '_vs_', clusters_2)
for (i in 1:length(clusters)) {
  cluster_1 = clusters_1[i]
  this_clust_degs = degs[degs$cluster_1==cluster_1, c('symbol', 'avg_log2FC', 'p_val_adj', "p_val")]
  df_header = c('[Genes]', 'Log2 Fold-Change', 'q-Value', 'P-value')
  this_clust_degs = rbind(df_header, this_clust_degs)
  this_clust_degs = rbind(c('[GeneSet]', '', '', ''), this_clust_degs)
  fwrite(this_clust_degs, file=paste0(out_dir,"degs_", clusters[i],".csv"), na="0", quote=F, row.names=F, col.names=F)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%#
## Export cluster markers
#%%%%%%%%%%%%%%%%%%%%%%%%%%%#
rna_markers = fread(file=snakemake@input[['rna_markers']])
if ('ADT' %in% Assays(import)) {
  adt_markers = fread(file=snakemake@input[['adt_markers']])
  all_markers = as.data.frame(rbind(adt_markers, rna_markers))
} else {
  all_markers = as.data.frame(rna_markers)
}  
if ('resolution' %in% names(all_markers)) {
  all_markers = all_markers[all_markers$resolution==resol, ]
}

clusters = sort(unique(import@meta.data[, cluster_by]))
for (i in 1:length(clusters)) {
  clust = clusters[i]
  this_clust_markers = all_markers[all_markers$cluster==clust, c('gene', 'avg_log2FC', 'p_val_adj')]
  this_clust_markers = this_clust_markers[this_clust_markers$p_val_adj<0.05, ]
  df_header = c('[Genes]', 'Log2 Fold-Change', 'q-Value')
  this_clust_markers = rbind(df_header, this_clust_markers)
  this_clust_markers = rbind(c('[GeneSet]', '', ''), this_clust_markers)
  fwrite(this_clust_markers, file=paste0(out_dir,"C", clusters[i],"_markers.csv"), na="0", quote=F, row.names=F, col.names=F)
}
