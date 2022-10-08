#!/bin/bash
#SBATCH -J pyscenic
#SBATCH -t 2-00:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -p himem #veryhimem
#SBATCH --mem=60gb #180gb
#SBATCH --array=1-100
#SBATCH --output=%x-%j.out
#SBATCH --mail-user=diala.abdrabbo@uhnresearch.ca


cd /cluster/projects/mcgahalab/data/brookslab/sabelo/AB_seq_v2/
conds=('Tum_CD8' 'Tum_CD8_KO.C3_WT.C0')

for cond in ${conds[@]}; do

  adj_name=results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/pyscenic_100_runs/adj_${cond}_run_${SLURM_ARRAY_TASK_ID}.tsv
  reg_name=results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/pyscenic_100_runs/reg_${cond}_run_${SLURM_ARRAY_TASK_ID}.csv

  echo ***${adj_name}***
  arboreto_with_multiprocessing.py \
      results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/pyscenic_100_runs/Seur_${cond}.loom \
      /cluster/home/dabdrabb/projects/data/mm_mgi_tfs.txt \
      --method grnboost2 \
      --output ${adj_name} \
      --num_workers 24

  echo ***${reg_name}***
  pyscenic ctx \
      ${adj_name} \
      /cluster/home/dabdrabb/projects/data/cisTarget_databases/mm10__*.feather \
      --annotations_fname /cluster/home/dabdrabb/projects/data/motifs-v9-nr.mgi-m0.001-o0.0.tbl \
      --expression_mtx_fname results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/pyscenic_100_runs/Seur_${cond}.loom \
      --mode "dask_multiprocessing" \
      --output ${reg_name} \
      --num_workers 24
done


# to run pyscenic scheduler on the server, log to a computing node on the cluster
# cd /cluster/projects/mcgahalab/data/brookslab/sabelo/AB_seq_v2/bin/
# conda activate pyscenic
# sbatch 05c_pyscenic_scheduler_100_runs.sh
