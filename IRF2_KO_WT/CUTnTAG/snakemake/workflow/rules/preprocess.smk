rule samtools_mapq_filter:
  input:
    "results/alignment/bam/{sample}_{celltype}_{antibody}.bam"
  output:
    temp("results/alignment/bam/{sample}_{celltype}_{antibody}.filtered.bam"),
  log:
    "logs/samtools/qfilter_intermed/{sample}_{celltype}_{antibody}.log"
  params:
    extra=combine_args(config["params"]["samtools"]["qfilter_intermed"])
  wrapper:
    "v0.75.0/bio/samtools/view"

rule samtools_mapq_index:
  input:
    "results/alignment/bam/{sample}_{celltype}_{antibody}.filtered.bam",
  output:
    "results/alignment/bam/{sample}_{celltype}_{antibody}.filtered.bam.bai",
  log:
    "logs/samtools/qfilter_intermed/{sample}_{celltype}_{antibody}.index.log"
  params:
    "" # optional params string
  threads:  # Samtools takes additional threads through its option -@
    4     # This value - 1 will be sent to -@
  wrapper:
    "v0.75.0/bio/samtools/index"

rule samtools_sort:
  input:
    "results/alignment/bam/{sample}_{celltype}_{antibody}.filtered.bam",
  output:
    "results/alignment/bam/{sample}_{celltype}_{antibody}.sorted.bam",
  log:
    "logs/samtools/sort/{sample}_{celltype}_{antibody}.log"
  params:
    extra = "-m 4G",
  threads:  # Samtools takes additional threads through its option -@
    4     # This value - 1 will be sent to -@
  wrapper:
    "v0.75.0/bio/samtools/sort"

rule samtools_sort_index:
  input:
    "results/alignment/bam/{sample}_{celltype}_{antibody}.sorted.bam",
  output:
    "results/alignment/bam/{sample}_{celltype}_{antibody}.sorted.bam.bai",
  log:
    "logs/samtools/sort/{sample}_{celltype}_{antibody}.index.log"
  params:
    "" # optional params string
  threads:  # Samtools takes additional threads through its option -@
    4     # This value - 1 will be sent to -@
  wrapper:
    "v0.75.0/bio/samtools/index"

rule samtools_spikein_sort:
  input:
    "results/alignment/spikein/{sample}_{celltype}_{antibody}.spikein.bam",
  output:
    "results/alignment/spikein/{sample}_{celltype}_{antibody}.spikein.sorted.bam",
  log:
    "logs/samtools/sort/{sample}_{celltype}_{antibody}.spikein.log"
  params:
    extra = "-m 4G",
  threads:  # Samtools takes additional threads through its option -@
    4     # This value - 1 will be sent to -@
  wrapper:
    "v0.75.0/bio/samtools/sort"

rule samtools_spikein_sort_index:
  input:
    "results/alignment/spikein/{sample}_{celltype}_{antibody}.spikein.sorted.bam",
  output:
    "results/alignment/spikein/{sample}_{celltype}_{antibody}.spikein.sorted.bam.bai",
  log:
    "logs/samtools/sort/{sample}_{celltype}_{antibody}.spikein.index.log"
  params:
    "" # optional params string
  threads:  # Samtools takes additional threads through its option -@
    4     # This value - 1 will be sent to -@
  wrapper:
    "v0.75.0/bio/samtools/index"

rule rm_unused_chr:
  input:
    bam="results/alignment/bam/{sample}_{celltype}_{antibody}.sorted.bam"
  output:
    bam=temp("results/alignment/dedup/{sample}_{celltype}_{antibody}.tmp.bam"),
    idx=temp("results/alignment/dedup/{sample}_{celltype}_{antibody}.tmp.bam.bai"),
  log:
    "logs/samtools/view/{sample}_{celltype}_{antibody}.rm_unused_chr.log",
  conda:
    "../envs/bamtools.yaml",
  params:
    chrs=config["params"]["samtools"]["chrs"],  # optional params stringparams:
  shell:
    """
    samtools index {input.bam}; \
    samtools view -b {input.bam} {params.chrs} | samtools sort > {output.bam}; \
    samtools index {output.bam};
    """

rule mark_duplicates:
  input:
    "results/alignment/dedup/{sample}_{celltype}_{antibody}.tmp.bam"
  output:
    bam="results/alignment/dedup/{sample}_{celltype}_{antibody}.bam",
    metrics="results/alignment/dedup/{sample}_{celltype}_{antibody}.metrics.txt"
  log:
    "logs/picard/dedup/{sample}_{celltype}_{antibody}.log"
  params:
    lambda w: "{common_params} {rm_duplicates}".format(
      common_params="{}".format(combine_args(config["params"]["picard"]["mark_duplicates"]["args"])),
      rm_duplicates="{}".format(combine_args(config["params"]["picard"]["mark_duplicates"]["control"])) if is_control else "{}".format(combine_args(config["params"]["picard"]["mark_duplicates"]["target"]))
      )
  wrapper:
    "v0.75.0/bio/picard/markduplicates"

rule samtools_dedup_index:
  input:
    "results/alignment/dedup/{sample}_{celltype}_{antibody}.bam",
  output:
    "results/alignment/dedup/{sample}_{celltype}_{antibody}.bam.bai",
  log:
    "logs/samtools/dedup/{sample}_{celltype}_{antibody}.index.log"
  params:
    "" # optional params string
  threads:  # Samtools takes additional threads through its option -@
    4     # This value - 1 will be sent to -@
  wrapper:
    "v0.75.0/bio/samtools/index"
