library(dplyr)
library(Seurat)
library(ggplot2)


conds = read.csv(snakemake@params[['conds']], row.names=NULL, header=T)$conditions
seur  = readRDS(snakemake@input[['seur']])

# load cell cycle genes
load(snakemake@params[['cellcycle_genes']])
length(s.genes)
length(g2m.genes)

seur = CellCycleScoring(seur, g2m.features=g2m.genes, s.features=s.genes, set.ident=TRUE)
head(seur[[]])

seur = NormalizeData(seur) %>%
  FindVariableFeatures() %>%
  ScaleData(features=row.names(seur)) %>%
  RunPCA()

# Normalize AB-seq data
seur <-  NormalizeData(seur,
                       assay="ADT",
                       normalization.method = "CLR") %>%  
  ScaleData(assay="ADT") %>% 
  RunPCA(reduction.name='apca')

seur = RunPCA(seur, features=c(s.genes, g2m.genes))

Idents(seur) = 'orig.ident'
p1 = DimPlot(seur, reduction='pca', split.by='Phase') + ggtitle('Before cell cycle regression')


seur$CC.Difference <- seur$S.Score - seur$G2M.Score
seur = ScaleData(seur, vars.to.regress=c("pc_mito", "CC.Difference"), features=rownames(seur))

seur = RunPCA(seur, features=c(s.genes, g2m.genes))
Idents(seur) = 'orig.ident'
p2 = DimPlot(seur, reduction='pca', split.by='Phase') + ggtitle('After cell cycle regression')

pdf(snakemake@output[['plot']], width=10, height=4)
p1 + p2
dev.off()

saveRDS(seur, file=snakemake@output[['seur']])


