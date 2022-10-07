scranDB<-function(seu=seu){
  require(scran)
  dbl<-scDblFinder::computeDoubletDensity(seu@assays$RNA@counts)
  seu$scran<-dbl
  return(seu)
}
