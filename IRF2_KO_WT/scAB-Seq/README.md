# scAB-Seq Analyses


## About
Analyses of the scAB-Seq data of IRF2 knockout vs wild type tumour cells (IRF2_KO vs WT).

## Prerequisites:
To reproduce the scAB-Seq results, two types of code need to be run: A conda workflow and some individual R scripts. To run the conda workflow:
1. Install mambaforge from [conda forge website](https://github.com/conda-forge/miniforge) (look for the mambaforge section in this link).
2. Install snakemake.
2. Unzip the raw counts matrix file called `tagged_cells_raw_counts.csv.zip`.
2. Build the environment of this workflow:
  * Activate the snakemake environment: `conda activate snakemake`
  * Run snakemake with the first line of command just to build the environment: `sh run_snakemake.sh`
  * Check the name of the create conda environment; it can be found in
  * Activate the create conda environment: `conda activate name_of_created_environment`

Rscript -e "install.packages('adabag', version=4.2 , repos='https://cloud.r-project.org', dependencies=F)"
Rscript -e "devtools::install_github('chris-mcginnis-ucsf/DoubletFinder', dependencies=F)"
Rscript -e "devtools::install_github('aertslab/SCopeLoomR', dependencies=F)"
Rscript -e "devtools::install_github('aertslab/SCENIC', dependencies=F)"
Rscript -e "install.packages('/Users/dabdrabb/projects/updated_packages/Chord-main', dependencies=F, type='source', repos=NULL)"
Rscript -e "devtools::install_github('siyao-liu/MultiK', dependencies=F)"

conda deactivate

continue running the rest of your snakemake


cd /Users/dabdrabb/projects/sabelo/cytof
sh run_snakemake.sh

in order to make the interactive plots, we need to add the path of pando to the path variable. Add the below to the top of the .zshrc file.
export PATH=$PATH:/Applications/RStudio.app/Contents/MacOS/pandoc



Fastq files were processed by the BD Rhapsody™ WTA Analysis Pipeline, available on SevenBridges server, run with default parameters, to generate the read counts matrices. The raw counts matrix file called `tagged_cells_raw_counts.csv.zip` is needed for the AB-Seq workflow and can be found in the [Seurat_AB_seq]() folder. Unzip this file and run the workflow by running `snakemake -F --cores 2 --use-conda`
* peaks.bed
* singlecell.csv
* fragments.tsv.gz

Those files are compressed as a tar.gz file and submitted to the Gene Expression Omnibus (GEO accession number: GSE199177).

## Code description
A Snakemake workflow is provided to analyze all scAB-Seq related data, except for the pySCENIC part. pySCENIC code is provided in the bin folder. Another script that

Code should be run in the order in which it was numbered
* `01_preprocess_scATAC-seq.R`: Read in peaks data, perform quality control to filter out poor-quality cells and doublets, integrate Cre samples apart and KO and WT apart, annotate cell types using geneactivity,
and finally remove non-CD8 T cells
* `02_cluster_n_da_peaks.R`: cluster cells by peaks counts called using MACS2 and identify peaks differentially accessible between CD8_cre_pos vs CD8_cre_neg and IRF2_KO vs WT.
* `03_genes_with_IRF2_peaks.R`: Identify IRF2 motifs accessible in genes promoter of the mouse genome.
