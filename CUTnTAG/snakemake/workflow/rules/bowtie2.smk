rule bowtie2:
  input:
    sample=["results/trim/trimgalore/{sample}_{celltype}_{antibody}.1_val_1.fq.gz",
            "results/trim/trimgalore/{sample}_{celltype}_{antibody}.2_val_2.fq.gz"]
  output:
    "results/alignment/bam/{sample}_{celltype}_{antibody}.bam",
  log:
    "logs/bowtie2/{sample}_{celltype}_{antibody}.log"
  params:
    # https://github.com/nf-core/cutandrun/blob/049f57c17399579cb62dbd62b66eb21eac180ce7/nextflow.config
    # https://yezhengstat.github.io/CUTTag_tutorial/#312_Alignment_to_spike-in_genome_for_spike-in_calibration_[optionalrecommended]
    index=config['params']['bowtie2']['genome'],  # prefix of reference genome index (built with bowtie2-build)
    extra=combine_args(config["params"]["bowtie2"]["extra"]),  # optional parameters
  threads: 8  # Use at least two threads
  wrapper:
    # https://github.com/nf-core/cutandrun/blob/36c9c4fea4db8f1b50729f292c32bbd939eb9818/modules/nf-core/software/bowtie2/align/main.nf
    "v0.75.0/bio/bowtie2/align"

rule bowtie2_spikein:
  input:
    sample=["results/trim/trimgalore/{sample}_{celltype}_{antibody}.1_val_1.fq.gz",
            "results/trim/trimgalore/{sample}_{celltype}_{antibody}.2_val_2.fq.gz"]
  output:
    "results/alignment/spikein/{sample}_{celltype}_{antibody}.spikein.bam"
  log:
    "logs/bowtie2/{sample}_{celltype}_{antibody}.spikein.log"
  params:
    # https://github.com/nf-core/cutandrun/blob/049f57c17399579cb62dbd62b66eb21eac180ce7/nextflow.config
    # https://yezhengstat.github.io/CUTTag_tutorial/#312_Alignment_to_spike-in_genome_for_spike-in_calibration_[optionalrecommended]
    index=config['params']['bowtie2']['spikein'],  # prefix of reference genome index (built with bowtie2-build)
    extra=" ".join([combine_args(config["params"]["bowtie2"]["extra"]),  # optional parameters
                    combine_args(config["params"]["bowtie2"]["extraspikein"])]),
  threads: 8  # Use at least two threads
  wrapper:
    # https://github.com/nf-core/cutandrun/blob/36c9c4fea4db8f1b50729f292c32bbd939eb9818/modules/nf-core/software/bowtie2/align/main.nf
    "v0.75.0/bio/bowtie2/align"

rule summarize_alignment:
  input:
    spikein=expand("results/alignment/spikein/{sid}.spikein.bam", sid=["_".join(i) for i in units.index]),
    sample=expand("results/alignment/bam/{sid}.bam", sid=["_".join(i) for i in units.index]),
  output:
    'results/alignment/qc/alignment_metrics.tsv'
  conda:
    "../envs/r.yaml"
  params:
    logdir='logs/bowtie2/',
    samples=",".join(expand("{sid}", sid=["_".join(i) for i in units.index])),
    outdir='results/alignment/qc',
  log:
    "logs/bowtie/summary_metrics.txt"
  shell:
    "Rscript scripts/alignment_metrics.R "
    "{params.logdir} {params.outdir} {params.samples}"
