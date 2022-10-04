#!/bin/bash

cd /cluster/projects/mcgahalab/data/brookslab/sabelo/AB_seq_v2/

conds=('Tum_CD8' 'Tum_CD8_KO.C3_WT.C0')
occurrence_cutoff=0.8

for cond in ${conds[@]}; do
  module load R/4.1.0
  Rscript bin/05d2_pyscenic_aggregate_100_runs.R ${cond} ${occurrence_cutoff}

  python3 bin/05d3_pyscenic_get_aucell.py ${cond}
  python3 bin/05d4_pyscenic_binarize_AUC.py ${cond}
done

# salloc --partition=himem -c 10 -t 10:0:0 --mem 60G
# cd /cluster/projects/mcgahalab/data/brookslab/sabelo/AB_seq_v2/bin/
# conda activate pyscenic
# sh 05d1_pyscenic_aggregate_100_runs_n_get_aucell.sh
