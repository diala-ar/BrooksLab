# scAB-Seq Analyses


## About
Analyses of the scAB-Seq data of IRF2 knockout vs wild type tumour cells (IRF2_KO vs WT). To reproduce the scAB-Seq results, two types of code need to be run: A snakemake workflow to analyze the scAB-Seq data and a bunch of individual R and sh scripts that needs to be run in a specific order to perform the pySCENIC analyses.

## How to run the snakemake workflow
To run the snakemake workflow, follow the below steps:
1. Install mambaforge from [conda forge website](https://github.com/conda-forge/miniforge) (look for the mambaforge section in this link).
2. Install snakemake by running the following in the terminal: `conda install -n base -c conda-forge mamba`.
3. Unzip the raw counts matrix file `BrooksLab-main/IRF2_KO_WT/scAB-Seq/resources/Seurat_AB_seq/tagged_cells_raw_counts.csv.zip`.
4. Build the environment of snakemake workflow:
    * Open the file `BrooksLab/IRF2_KO_WT/scAB-Seq/workflow/Snakefile` and change the `workdir` variable to the full path of the `scAB-Seq` folder on your computer.
    * To activate the snakemake environment, type the following in the terminal: `conda activate snakemake`
    * To build the conda environment of the snakemake workflow, we will run the snakemake workflow with just the first line of command, to do so run : `sh run_snakemake.sh`.
    * Check the name of the created conda environment; it is the name of the folder in the `BrooksLab-main/IRF2_KO_WT/scAB-Seq/.snakemake/conda/` folder. Let's suppose it is called name_of_created_environment.
5. change the `condaprefix` variable in the `BrooksLab-main/IRF2_KO_WT/scAB-Seq/run_snakemake.sh` file to the full path of the `BrooksLab-main/IRF2_KO_WT/scAB-Seq/.snakemake/conda/` folder on your computer. This variable contains the path of the created conda environment.
6. Activate the created conda environment by running: `conda activate name_of_created_environment`.
7. Install some R packages to the created conda environment by running the following in the specified order in the terminal:
```
    Rscript -e "install.packages('adabag', version=4.2 , repos='https://cloud.r-project.org', dependencies=F)"
    Rscript -e "devtools::install_github('chris-mcginnis-ucsf/DoubletFinder', dependencies=F)"
    Rscript -e "devtools::install_github('aertslab/SCopeLoomR', dependencies=F)"
    Rscript -e "devtools::install_github('aertslab/SCENIC', dependencies=F)"
    # We updated the Chord package to prevent an error issued while running the scranDB function.
    # We provide the updated version of the Chord package, and you can install it by replacing `chord_path` by the full path of the "BrooksLab-main/IRF2_KO_WT/scAB-Seq/updated_chord_package" folder on your computer.
    Rscript -e "install.packages('chord_path', dependencies=F, type='source', repos=NULL)"`
    conda deactivate
```
8. Now that the conda environment is created, we can run the rest of the workflow. To do so, open the file `BrooksLab-main/IRF2_KO_WT/scAB-Seq/workfow/rules/Seurat_AB_seq.smk` and uncomment code from lines 14 to lines 52.
9. To run the snakemake workflow and generate results, run: `sh run_snakemake.sh`.
N.B.: In order to make the interactive volcano plots, add the path of pando to the path variable. For Mac users, then add the following to the top of the .zshrc file: `export PATH=$PATH:/Applications/RStudio.app/Contents/MacOS/pandoc`/

### How to run the pySCENIC code
Scripts that will reproduce pySCENIC results needs to be run a cluster.

1. Create a pyscenic environment on the cluster:
```
mamba create -y -n pyscenic python=3.7
conda activate pyscenic
mamba install -y numpy
mamba install -y -c anaconda cytoolz
pip3 install pyscenic
mamba install sklearn
```
2. Run pySCENIC Scripts on a computing node on the cluster
```
Rscript 05a_pyscenic_prepare_seurat_objects.R
Rscript 05b_pyscenic_convert_seur_obj_to_loom.R
# to run the pySCENIC scheduler, make the bin folder your current directory on the cluster using the cd command.
sbatch 05c_pyscenic_scheduler_100_runs.sh
sh 05d1_pyscenic_aggregate_100_runs_n_get_aucell.sh
Rscript 05e_pyscenic_visualize_pyscenic_results.R
```

## Code description
### Snakemake workflow rules
1. `Seurat_AB_seq.smk` runs the below R scripts in the listed order:
    * `seur_create_multimodal_seurat_object.R` Creates a multimodal seurat object containing proteins expression in the ADT assay and genes expression in the RNA assay
    * `seur_quality_control_plots.R` Generates quality control plots.
    * `seur_filter_cells_n_genes.R` Filters out poor-quality cells and doublets, and filters out genes that are expressed in less than 10 cells.  
    * `seur_remove_doublets_Chord` Removes remaining doublets identified with the Chord R package
    * `seur_regress_out_cell_cycle_n_mito_prct.R` Regresses out cell cycle genes and mito chondrial genes to prevent any biased results in subsequent analyses.
    * `seur_elbow_plots.R` Generates an elbow plot to help choose how many principle components need to be used for clustering later.
    * `seur_cluster_cells.R` Clusters or reclusters cells after gating on cells (depending on the recluster parameter in each rule). Clustering uses the FindMultiModalNeighbors method to benefit from both protein and gene expression data, in case one of these modalities does not provide enough data and confidence to cluster a given cell. Clustering is performed with various resolutions to help decide which resolution is best to choose for subsequent analyses. Also generates different UMAPs and heatmaps helping specifying the cell type of each cluster in all resolutions.
    * `seur_gate_on_clusters_at_resol.R` gates on cells in specific cluter(s) given a predefined resolution.
    * `seur_visualize_umap_n_heatmap_for_one_resol.R` Generates different different UMAPs and heatmaps helping specifying the cell type of each cluster in all resolutions to help identify cell type of each cluster for a given resolution.
    * `seur_identify_DEGs.R` Identifies differentially expressed genes (DEGs) between the 2 conditions within each cluster and between total cells of the 2 conditions.
    * `seur_KO.C3_vs_WT.C0_analyses.R` Identifies DEGs between cluter 3 of KO cells and cluster 0 of WT cells and generate heatmaps of those DEGs or a subset of those DEGs
    * `export_seur_data_to_SeqGeq.R` Exports data of the Seurat object to files that can be imported into the SeqGeq app.
    * `seur_visualize_DEGs.R` Generates various types of plots to visualize DEGs: static and interactive volcano plots, heatmaps of their gene expression and some plots for a subset of genes of interest.
2. `GSEA.smk` permits to run the `GSEA_of_custom_gs.R` script to perform a GSEA analyses on, only, cells of cluster 3 in KO vs cluster 0 in WT for 2 genesets ('EFFECTOR_VS_EXHAUSTED_CD8_TCELL' taken from C7 misgdb and 'tox_degs' representing DEGs between Tox KO cells vs WT cells taken from Table S2 from Khan et al., 2019, Nature, PMID: 31207603).

### pySCENIC scripts and bash files:
* `05a_pyscenic_prepare_seurat_objects.R` Preprocesses data: 1)Rename some meta_data columns to prepare data for pyscenic analyses. 2)Subset data to include only clusters of interest (cluster 3 from KO cells and cluster 0 from WT cells).
* `05b_pyscenic_convert_seur_obj_to_loom.R` Converts RDS files containing seurat objects to loom file which is needed to run the arboreto_with_multiprocessing function and generates tsv files that are needed to run the aucell function.
* `05c_pyscenic_scheduler_100_runs.sh` Runs the pySCENIC implementation of the SCENIC algorithm 100 times.
* `05d1_pyscenic_aggregate_100_runs_n_get_aucell.sh` aggregate results of the 100 runs and generates AUC scores called also regulon activity scores (RASs).
* `05d2_pyscenic_aggregate_100_runs.R` Aggregates results of the 100 pySCENIC runs.
* `05d3_pyscenic_get_aucell.py` Generates AUC scores or RAS.
* `05d4_pyscenic_binarize_AUC.py` Binarizes AUC scores using predefined cutoffs (for each regulon, a cutoff is defined based on a histogram of AUC scores of this regulon).
* `05e_pyscenic_visualize_pyscenic_results.R` Generates dotplots of AUC scores.
* `binarization.py` Called by `05d4_pyscenic_binarize_AUC.py`
* `prepare_data_for_GSEA_app.R` Exports files used as input to generate the enrichment map visualization in the GSEA application.
* `05_pyscenic_functions_library.R` Contains a library of functions called by some scripts called above.
N.B.: for more information about running and aggregating pySCENIC results, see Methods section of the article.
