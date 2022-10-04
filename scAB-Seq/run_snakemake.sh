#! /bin/bash

condaprefix='/Users/dabdrabb/projects/sabelo/AB_seq_v2/.snakemake/conda'

Snakemake \
--cores 4 \
--use-conda \
--conda-prefix ${condaprefix}


# to run the workflow run: 'snakemake -F --cores 2 --use-conda'
