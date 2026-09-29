## Run SCENIC in R from a Seurat object

library(Seurat)
library(SCENIC)
library(data.table)
library(Matrix)
library(AUCell)
library(RcisTarget)
library(reticulate)
library(ggplot2)

renv::use_python(type='virtualenv', name='renv/py_env')
arboreto = import('arboreto.algo', convert=F)
distributed = import('distributed', convert=F)

dir_seur   = file.path('results', 'seurat_object')

## directory containing cisTarget feather databases
dir_resources = file.path('/', 'Users', 'dabdrabb', 'projects', 'data', 'cisTarget_databases')

org = 'mgi' # mouse
nCores = 4
set.seed(123)


## 1. Load and filter Seurat object
seurs_ls = readRDS(file.path(dir_seur, '04_integrated_seurs_ls_n_clustered_CD8_T_cells.rds'))
samples = names(seurs_ls)

### Initialize settings
data('motifAnnotations_mgi_v9', package='RcisTarget')
motifAnnotations_mgi = motifAnnotations_mgi_v9

old_wd = getwd()

for (sample_i in samples) {
  dir_scenic = file.path('results', paste0('SCENIC_', sample_i))
  dir.create(dir_scenic)
  setwd(dir_scenic)
  dir_init   = 'init'
  dir.create(dir_init)

  scenicOptions = initializeScenic(org=org, dbDir=dir_resources, nCores=4, 
                                   dbs=c('500bp'='mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather',
                                         '10kb'='mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather'))
  saveRDS(scenicOptions, file=file.path(dir_init, "scenicOptions.Rds")) 
  
  
  seur = seurs_ls[[sample_i]]
  expr = GetAssayData(seur, assay='RNA', slot='counts')
  expr = as.matrix(expr)
  message('Initial matrix dimensions:')
  print(dim(expr))
  
  genesKept <- geneFiltering(as.matrix(expr), scenicOptions=scenicOptions,
                             minCountsPerGene = 3 * .01 * ncol(expr),
                             minSamples = ncol(expr) * .01)
  # Maximum value in the expression matrix: 5289
  # Ratio of detected vs non-detected: 0.32
  # Number of counts (in the dataset units) per gene:
  #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 0     203    1804    9360    5664 2644615 
  # Number of cells in which each gene is detected:
  #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 0     184    1522    2464    3844   10101 
  # 
  # Number of genes left after applying the following filters (sequential):
  #   9815	genes with counts per gene > 303.03
  # 9804	genes detected in more than 101.01 cells
  # Using the column 'features' as feature index for the ranking database.
  # 9288	genes available in RcisTarget database
  # Gene list saved in int/1.1_genesKept.Rds


  expr_filtered = expr[genesKept, , drop=F]
  dim(expr_filtered)
  # [1]  9288 10101  D8
  # [1] 6636 8423    D50
  # Get the TF list corresponding to mouse SCENIC setup
  tf_names = intersect(getDbTfs(scenicOptions), rownames(expr_filtered))
  length(tf_names)
  # [1] 848  D8
  # [2] 598  D50
  
  runCorrelation(expr_filtered, scenicOptions)
  expr_grn = as(expr_filtered, 'dgCMatrix')
  expr_grn@x = log2(expr_grn@x + 1)
  rm(expr); gc()
  
  # Run GRNboost
  expr_py = r_to_py(t(expr_grn), convert=F)
  gene_names_py = r_to_py(as.list(rownames(expr_grn)), convert=F)
  tf_names_py = r_to_py(as.list(tf_names), convert=F)
 
  # create Dask client
  cluster = distributed$LocalCluster(n_workers = 4L,
                                     threads_per_worker=1L,
                                     processes = F)
  client = distributed$Client(cluster)
  adj_py = arboreto$grnboost2(expression_data = expr_py,
                              gene_names = gene_names_py,
                              tf_names = tf_names_py,
                              client_or_address = client,
                              seed = 123L,
                              verbose = T)
  adj_file = file.path(dir_init, 'grnboost_adjacencies.tsv')
  adj_py$to_csv(adj_file, sep='\t', index=F)
  
  client$close()
  cluster$close()
  
  adj = data.table::fread(adj_file, data.table=F)
  colnames(adj) = c('TF', 'Target', 'weight')
  adj = adj[order(adj$weight, decreasing=T), ]
  file.remove(adj_file)
  
  
  scenicOptions = runSCENIC_1_coexNetwork2modules(scenicOptions, linkList=adj)
  scenicOptions = runSCENIC_2_createRegulons(scenicOptions)
  regulons = loadInt(scenicOptions, 'regulons')

  gc()
  expr = GetAssayData(seur, assay='RNA', slot='counts')
  expr_log = as(expr, 'dgCMatrix')
  expr_log@x = log2(expr_log@x + 1)
  scenicOptions = runSCENIC_3_scoreCells(scenicOptions, exprMat=expr_log)
  regulonAUC = loadInt(scenicOptions, 'aucell_regulonAUC')
  aucMat = AUCell::getAUC(regulonAUC)
  dim(aucMat)
  seur[['SCENIC']] = CreateAssayObject(data=aucMat)
  seurs_ls[[sample_i]] = seur
  
  saveRDS(scenicOptions, file.path(dir_init, 'scenicOptions.rds'))
  saveRDS(adj,  file.path(dir_init, 'grnboost_adjacencies.rds'))
  saveRDS(regulons,  file.path(dir_init, 'regulons.rds'))
  
  setwd(old_wd)
  saveRDS(seurs_ls, file.path(dir_seur, '05_seurs_ls_with_SCENIC_AUC_scores.rds'))
  
  
  #### generate the auc_stats excel file of differentially activated/inhibited regulons
  # to reload data after relaunching R
  # scenicOptions = readRDS('scenicOptions.rds')
  # adj = readRDS('int/grnboost_adjacencies.rds')
  # regulons = loadInt(scenicOptions, 'regulons')
  # regulonAUC = loadInt(scenicOptions, 'aucell_regulonAUC')
  # aucMat = AUCell::getAUC(regulonAUC)
  if (sample_i == 'D8') {
    tbet_reg = grep('Tbx21', names(regulons), value=T)
    # [1] "Tbx21_extended"
    scenic_tbet_genes = regulons[[tbet_reg]]
    # [1] "Fam117a" "Pycard"  "Rara"    "Tbx21"   "Zeb2"  in D8
    VlnPlot(seur, features=scenic_genes_tbet, group.by='sample', pt.size=0, combine=F) |>
      lapply(function(p) {
        p + geom_boxplot(width=0.12, outlier.shape=NA)}) |>
      patchwork::wrap_plots()
  }
  
  
  Idents(seur) = 'seurat_clusters'
  clusters = c(levels(seur$seurat_clusters), 'Bulk')
  auc_stats = lapply(clusters, function(clust) {
    auc_stats_clust = do.call(rbind, lapply(rownames(aucMat), function(reg) {
      if (clust == 'Bulk') {
        x1 = aucMat[reg, seur$sample == paste0('ON_', sample_i)]
        x2 = aucMat[reg, seur$sample == paste0('OFF_', sample_i)]
      } else {
        x1 = aucMat[reg, seur$seurat_clusters==clust & seur$sample == paste0('ON_', sample_i)]
        x2 = aucMat[reg, seur$seurat_clusters==clust & seur$sample == paste0('OFF_', sample_i)]
      }
      x = c(x1, x2)
      
      df = data.frame(regulon = reg,
                      mean_1 = round(mean(x1), 3),
                      mean_2 = round(mean(x2), 3),
                      delta_AUC = round(mean(x1) - mean(x2), 4),
                      SD = round(sd(x), 3),
                      delta_SD = round((mean(x1) - mean(x2)) / sd(x), 2),
                      p_val = wilcox.test(x1, x2)$p.value,
                      cluster = clust)
    })) %>%
      dplyr::arrange(desc(delta_SD)) %>% 
      dplyr::mutate(adj_p_val = p.adjust(p_val, method='fdr')) %>% 
      dplyr::filter(adj_p_val < 0.05)
    if (nrow(auc_stats_clust) > 0) {
      names(auc_stats_clust)[2:3] = paste0(c('mean_ON_', 'mean_OFF_'), sample_i)
    }
    auc_stats_clust
  })
  names(auc_stats) = c(paste0('C', clusters[-6]), 'Bulk')
  openxlsx::write.xlsx(auc_stats, file.path(dir_scenic, paste0('SCENIC_regulon_AUC_ON_vs_OFF_', sample_i, '.xlsx')))
}


