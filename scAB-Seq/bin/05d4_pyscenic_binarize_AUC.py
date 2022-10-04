import sys
#import binarization as fl
#fl.binarize
import pandas as pd

cond = sys.argv[1]

exec(open('bin/binarization.py').read())

AUC_MTX_FNAME = "results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/pyscenic_100_runs/aucs_" + cond + ".csv" # Gene expression as (cell, gene) - matrix.

auc_matrix = pd.read_csv(AUC_MTX_FNAME, sep='\,', header=0, index_col=0, engine='python')
auc_matrix.shape

bin = binarize(auc_matrix, num_workers=10)

auc_bin_pd = bin[0]
auc_thresholds = bin[1]

auc_bin_pd.to_csv('results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/pyscenic_100_runs/binary_aucs_' + cond + '.csv')
auc_thresholds.to_csv('results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/pyscenic_100_runs/aucs_thresholds_' + cond + '.csv')
