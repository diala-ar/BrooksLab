library(Seurat)

seur = readRDS(file=snakemake@input[['seur']])

p1 = ElbowPlot(seur, ndims=50, reduction='pca')
if ('ADT' %in% Assays(seur)) {
  p2 = ElbowPlot(seur, ndims=nrow(GetAssayData(seur, assay='ADT')), reduction='apca')
  pdf(snakemake@output[['elbow_plots']], width=10)
  p1 + p2
} else {
  pdf(snakemake@output[['elbow_plots']])
  print(p1)
}
dev.off()
