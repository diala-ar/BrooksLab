# identify doublets
library(Seurat)
library(ggplot2)
library(dplyr)
library(Chord)

proj        = unique(read.csv(snakemake@params[['conds']])$project)
seur        = readRDS(file=snakemake@input[['seur']])
doubletrate = ncol(seur) * 0.009 / 1000  # there are 0.9% doublets per 1000 cells
doubletrate
proj
seur


system('mkdir -vp results/Chord/')
setwd('results/Chord/')
getwd()
chord(seu=seur, doubletrate=doubletrate, overkill=T, outname=paste0(proj, '_'))
load('seu.robj')
doublets = read.csv(paste0(proj, '__doublet.csv'))$x
setwd('../../')
getwd()


# # re-normalize, find variable features, and scaledata after removing doublets
meta_data = seur@meta.data
all(row.names(seur@meta.data)==row.names(seu@meta.data))
meta_data = cbind(meta_data, seu@meta.data$chord)
meta_data$doublet = sapply(row.names(meta_data), function(x) ifelse(x %in% doublets, T, F))
seur@meta.data = meta_data
head(seur, 3)
rm(seu); gc()

seur = subset(seur, doublet==F)
saveRDS(seur, file=snakemake@output[['seur']])
