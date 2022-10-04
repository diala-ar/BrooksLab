import sys
import pandas as pd
from pyscenic.genesig import GeneSignature
from pyscenic.aucell import aucell, derive_auc_threshold, create_rankings

cond = sys.argv[1]

# A module from GeneSigDB (C6)
GMT_FNAME = 'results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/pyscenic_100_runs/high_confidence_reg_' + cond + '.gmt'
# An expression matrix from GEO
EXPRESSION_MTX_FNAME = "results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/pyscenic_100_runs/Seur_" + cond + "_counts.tsv" # Gene expression as (cell, gene) - matrix.

signatures = GeneSignature.from_gmt(GMT_FNAME, field_separator='\t', gene_separator='\t')
len(signatures)


ex_matrix = pd.read_csv(EXPRESSION_MTX_FNAME, sep='\t', header=0, index_col=0).T
ex_matrix.shape

percentiles = derive_auc_threshold(ex_matrix)
percentiles

aucs_mtx = aucell(ex_matrix, signatures, auc_threshold=percentiles[0.01], num_workers=8)
aucs_mtx.head()

aucs_mtx.to_csv('results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/pyscenic_100_runs/aucs_' + cond + '.csv')
