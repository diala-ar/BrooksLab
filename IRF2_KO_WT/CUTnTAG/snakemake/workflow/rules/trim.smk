rule trimmomatic_pe:
  input:
    get_fastqs
  output:
    r1="results/trim/trimmomatic/{sample}_{celltype}_{antibody}.1.fastq.gz",
    r2="results/trim/trimmomatic/{sample}_{celltype}_{antibody}.2.fastq.gz",
    # reads where trimming entirely removed the mate
    r1_unpaired="results/trim/trimmomatic/{sample}_{celltype}_{antibody}.1.unpaired.fastq.gz",
    r2_unpaired="results/trim/trimmomatic/{sample}_{celltype}_{antibody}.2.unpaired.fastq.gz"
  log:
    "logs/trimmomatic/{sample}_{celltype}_{antibody}.log"
  params:
    # list of trimmers (see manual)
    trimmer=["TRAILING:20", "LEADING:20", "SLIDINGWINDOW:4:15", "MINLEN:25",
        "ILLUMINACLIP:data/ref/Truseq3.PE.fa:2:15:4:4:true"],
    # optional parameters
    extra="-phred33",
    compression_level="-9"
  threads: 8
  resources:
    mem_mb=1024
  wrapper:
    "v0.75.0/bio/trimmomatic/pe"

rule trim_galore_pe:
  input:
    get_fastqs
  output:
    "results/trim/trimgalore/{sample}_{celltype}_{antibody}.1_val_1.fq.gz",
    "results/trim/trimgalore/{sample}_{celltype}_{antibody}.1.fastq.gz_trimming_report.txt",
    "results/trim/trimgalore/{sample}_{celltype}_{antibody}.2_val_2.fq.gz",
    "results/trim/trimgalore/{sample}_{celltype}_{antibody}.2.fastq.gz_trimming_report.txt"
  params:
    # https://github.com/nf-core/cutandrun/blob/049f57c17399579cb62dbd62b66eb21eac180ce7/nextflow.config
    extra=combine_args(config["params"]["trimgalore"]),
  log:
    "logs/trim_galore/{sample}_{celltype}_{antibody}.log"
  wrapper:
    # https://github.com/nf-core/cutandrun/blob/36c9c4fea4db8f1b50729f292c32bbd939eb9818/modules/nf-core/software/trimgalore/main.nf
    # https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md
    "v0.75.0/bio/trim_galore/pe"
