# Rbpj_ON_OFF project
This project analyses scRNA-Seq data and scTCR-seq data

## Abstract
Acute viral infection generates heterogeneous effector CD8 T cells, most of which are eliminated after viral clearance, while a minority persist as long-lived memory. Although Notch signaling is known to promote terminal effector cell differentiation, its role in memory T cell development remains unclear. Using complementary Rbpj-inducible mouse models with temporal and cell-restricted control of canonical Notch signaling, we examined how Notch activity programs antiviral CD8 T cell immunity. Across systemic and respiratory viral infections, temporal RBPJ-deficiency redirected virus-specific CD8 T cells away from CX3CR1high effector differentiation toward TCF1high/CD62L+ stem-like and tissue-resident–biased states. Paired scRNA-seq/scTCR-seq revealed that RBPJ-dependent Notch signaling maintained effector-associated transcriptional programs and favored expansion of effector-biased TCR clones, whereas RBPJ deficiency shifted both transcriptional and clonal landscapes toward memory-precursor and tissue-resident programs. This imprint persisted into memory, where RBPJ deficiency reduced the size of the virus-specific pool, depleted CX3CR1high TEM, enriched TCM/TRM subsets, and impaired polyfunctional recall responses. Competitive co-transfer and mixed bone marrow chimera experiments established that this reprogramming occurred cell-intrinsically during antigen priming. Despite augmented mitochondrial mass and respiratory capacity, RBPJ-deficient memory cells mediated weaker secondary protection, indicating that metabolic fitness alone is insufficient for CD8 T cell memory recall responses without RBPJ-dependent effector imprinting.

To reproduce the results of this project, first restore the R environment using the `renv.lock` file. Then, execute the R scripts located in the `bin/` folder in the numerical order they are listed in.

### 00_functions.R
This file contains all the functions that were used in the other R scripts.

### 01_scRNA_analyses_with_regressing_out_TCR_from_D8.R
The script is used to:
* create seurat object, identify doublets, poor-quality and dead cells and filter them out.
* preprocess samples, integrate them, cluster cells in to main cell types and re-cluster specific subpopulation.
* generate different illustrations and excel files that help in annotating cell types.
* perform DEGs analyses.
* perform GSEA analyses.

### 02_visualizations.R
This script contains different types of visualizations for scRNA-Seq data results (DEGs and GSEA) and for scTCR-Seq data (TCR genes, GSEA of top 3 clones in specific groups).

### 03_assess_clonality_using_scRepertoire.R
This script is used to:
* read in contig annotations of each sample and combine the list of T cell receptor contigs into clones.
* identify the top 3 TCRs in each group based on their abundance
* perform DEGs and GSEA analyses on the top 3 clones of each group.
* generate different visualizations.

### 04_recluster_clust_0_of_D8_samples.R
This script recluster cluster 0 of Day 8 samples in order to identify clusters that are associated with Rbpj knock out. 
It repeats the scRNA-Seq and the scTCR-Seq analyses on this specific subpopulation
