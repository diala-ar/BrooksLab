<h3 align="center">IRF2 CUT&TAG Analysis</h3>

---

## Table of Contents
- [About](#about)
- [Prerequisite Files](#prerequisites)
- [Code Description](#code)
- [Analyses](#analyses)

## About <a name = "about"></a>
Analysis of the IRF2 CUT&TAG data from unfixed activated T-cells.


## Prerequisites <a name = "prerequisites"></a>
The following reference files were obtained and are needed in the R-script analysis

**GTF File**
```
branch=''
release='102'
fmt='gtf' #gff3
species='mus_musculus'
species_cap='Mus_musculus'
build='GRCm38'
flavor=""
suffix="gtf.gz" # gff3.gz
datatype="dna"

url="ftp://ftp.ensembl.org/pub/${branch}release-${release}/${fmt}/${species}/${species_cap}.${build}.${release}.${flavor}${suffix}"
````

## Code Description <a name = "code"></a>

 * `analysis/getPredictedTFSites.R`: Gets the predicted TF binding sites for a given JASPAR2020 TF
 * `analysis/iterateMsigdb.R`: Wrapper function to iterate through msigdb levels to execute GSEA
 * `analysis/IRF2_analysis.R`: Main code used to identify significant IRF2 peaks and significance
 
## Analyzing the IRF2 CUT&TAG data <a name = "analyses"></a>
### a) Building a reference of predicted motif sites
The code in getPredictedTFSites.R can be executed to generate human (GRCh38) or
mouse (GRCm38) based motifs from the JASPAR2020 database using the corresponding
R package

This piece of code requires a hardcoded path to the corresponding GTF file

```
Rscript analysis/getPredictedTFSites.R

# Gets the predicted TF binding sites for a given JASPAR2020 TF
# *  Uses the entire genome with the genome as a background model
# *  set to use a p-threshold of 0.00005
# * Returns the bed file for IGV visualization
# * Returns the RDS file to load in for other analysis
# * Currently hardcoded the IRF2 TF
```

### b) Identifying SEACR Peaks
The main code can be found in `analysis/IRF2_analysis.R` and was run interactively to
perform the following tasks:

1) Create a catalogue of all potential IRF2 Peaks. If multiple IRF2 samples
are provided, it will create the union of all IRF2 peaks that exceed the 95th
percentile from SEACR_single sample mode. It will then subset all the IRF2
samples for the catalogue peaks to standardize across all samples.

2) Using the IRF2 (And other) motifs created from `getPredictedTFSites.R`, it will
annotate the peaks for each sample using any overlap with these motifs.

3) Execute some QC metrics to look for the cosine similarity between samples
and calculate the overlap between the motifs and the peaks.

4) Identify a threshold for what is considered a peak found only in IRF2 sample,
what is found in both IRF2 and IgG, and what peaks should be removed.

5) Extracts all the genes found within the confident peaks

6) Perform an Over-Representation Analysis using the msigdbr on the
HALLMARK and IMMUNE datasets for the genes found within the IRF2 confident Peaks
