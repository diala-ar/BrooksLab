# scAB-Seq Analyses


## About
Analyses of the scAB-Seq data of IRF2 knockout vs wild type tumour cells (IRF2_KO vs WT).

## Prerequisites:
To reproduce the scAB-Seq results, two types of code need to be run: A snakemake workflow to analyze the scAB-Seq data and some individual R scripts to perform the pySCENIC analyses.

### Snakemake workflow
To run the snakemake workflow:
1. Clone the github repository and unzip the downloaded file.
2. Install mambaforge from [conda forge website](https://github.com/conda-forge/miniforge) (look for the mambaforge section in this link).
3. Install snakemake by running the following in the terminal: `conda install -n base -c conda-forge mamba`.
4. Unzip the raw counts matrix file `BrooksLab-main/IRF2_KO_WT/scAB-Seq/resources/Seurat_AB_seq/tagged_cells_raw_counts.csv.zip`.
5. Build the environment of snakemake workflow:
    * Open the file `BrooksLab/IRF2_KO_WT/scAB-Seq/workflow/Snakefile` and change the `workdir` variable to the full path of the `scAB-Seq` folder on your computer.
    * To activate the snakemake environment, type the following in the terminal: `conda activate snakemake`
    * To build the environment, run snakemake : `sh run_snakemake.sh`.
    * Check the name of the created conda environment; it is the name of the folder in the `BrooksLab-main/IRF2_KO_WT/scAB-Seq/.snakemake/conda/` folder. Let's suppose it is called name_of_created_environment.
6. change the `condaprefix` variable in the `BrooksLab-main/IRF2_KO_WT/scAB-Seq/run_snakemake.sh` file to the full path of the `BrooksLab-main/IRF2_KO_WT/scAB-Seq/.snakemake/conda/` folder on your computer. This variable contains the path of the created conda environment.
7. Activate the created conda environment by typing the following in the terminal: `conda activate name_of_created_environment`.
8. Install some R packages to the created conda environment by running the following in the specified order in the terminal:

  `Rscript -e "install.packages('adabag', version=4.2 , repos='https://cloud.r-project.org', dependencies=F)"`

  `Rscript -e "devtools::install_github('chris-mcginnis-ucsf/DoubletFinder', dependencies=F)"`

  `Rscript -e "devtools::install_github('aertslab/SCopeLoomR', dependencies=F)"`

  `Rscript -e "devtools::install_github('aertslab/SCENIC', dependencies=F)"`

  `# We updated the Chord package to prevent an error issued while running the scranDB function`
  `# We provide the updated version of the Chord package, and you can install it by replacing `chord_path` by the full path of the "BrooksLab-main/IRF2_KO_WT/scAB-Seq/updated_chord_package" folder on your computer.`

  `Rscript -e "install.packages('chord_path', dependencies=F, type='source', repos=NULL)"`
  
  `conda deactivate`
9. To continue running the rest of the snakemake workfow, open the file `BrooksLab-main/IRF2_KO_WT/scAB-Seq/workfow/rules/Seurat_AB_seq.smk` and uncomment code from lines 14 to lines 52.
10. To run the snakemake workflow and generate results, run: `sh run_snakemake.sh`.
N.B.: In order to make the interactive plots, add the path of pando to the path variable. Add the following to the top of the .zshrc file: `export PATH=$PATH:/Applications/RStudio.app/Contents/MacOS/pandoc`/

### pySCENIC code

## Code description
A Snakemake workflow is provided to analyze all scAB-Seq related data, except for the pySCENIC part. pySCENIC code is provided in the bin folder. Another script that

Code should be run in the order in which it was numbered
