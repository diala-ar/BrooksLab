# scAB-Seq Analyses


## About
Analyses of the scAB-Seq data of IRF2 knockout vs wild type tumour cells (IRF2_KO vs WT).

## Prerequisites:
The below files are needed to create the seurat object:
* peaks.bed
* singlecell.csv
* fragments.tsv.gz

Those files are compressed as a tar.gz file and submitted to GEO (Accession number GSE199177).

## Code description
Code should be run in the order in which it was numbered
* `01_preprocess_scATAC-seq.R`: Read in peaks data, perform quality control to filter out poor-quality cells and doublets, integrate Cre samples apart and KO and WT apart, annotate cell types using geneactivity,
and finally remove non-CD8 T cells
* `02_cluster_n_da_peaks.R`: cluster cells by peaks counts called using MACS2 and identify peaks differentially accessible between CD8_cre_pos vs CD8_cre_neg and IRF2_KO vs WT.
* `03_genes_with_IRF2_peaks.R`: Identify IRF2 motifs accessible in genes promoter of the mouse genome.
