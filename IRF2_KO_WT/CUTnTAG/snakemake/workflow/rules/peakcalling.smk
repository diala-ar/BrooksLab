rule seacr_stringent:
  input:
    sample="results/alignment/bed/{sample}_{celltype}_{antibody}.bedGraph",
    control=lambda w: expand("results/alignment/bed/{ctrl}.bedGraph", ctrl=get_sample_control(w)),
  output:
    #"results/peaks/seacr/{sample}_{celltype}_{antibody}.seacr.txt",
    "results/peaks/seacr/{sample}_{celltype}_{antibody}.stringent.bed",
  log:
    "logs/seacr/{sample}_{celltype}_{antibody}.seacr.log"
  conda:
    "../envs/seacr.yaml"
  params:
    outid="results/peaks/seacr/{sample}_{celltype}_{antibody}",
    control=has_a_control,
    threshold=config['params']['seacr']['signal_threshold'],
    norm=config['params']['seacr']['normalization'],
    stringency=config['params']['seacr']['stringency']
  shell:
    # SEACR_1.3.sh CTCF_DE_chr1_100Mb.bedgraph.txt IgG_DE_chr1_100Mb.bedgraph.txt norm stringent CTCF_DE_chr1_100Mb
    # SEACR_1.3.sh CTCF_DE_chr1_100Mb.bedgraph.txt 0.05 norm stringent CTCF_DE_chr1_thresh
    "echo '{input.sample} - {input.control}'; "
    "if [ '{params.control}' == 'True' ]; then "
        "echo 'SEACR: ctrl bedGraph ' > {log}; "
        "SEACR_1.3.sh {input.sample} {input.control} "
        "{params.norm} stringent {params.outid}; "
    "else "
        "echo 'SEACR: threshold ({params.threshold}) ' > {log}; "
        "SEACR_1.3.sh {input.sample} {params.threshold} "
        "{params.norm} stringent {params.outid}; "
    "fi "

rule seacr_relaxed:
  input:
    sample="results/alignment/bed/{sample}_{celltype}_{antibody}.bedGraph",
    control=lambda w: expand("results/alignment/bed/{ctrl}.bedGraph", ctrl=get_sample_control(w)),
  output:
    #"results/peaks/seacr/{sample}_{celltype}_{antibody}.seacr.txt",
    "results/peaks/seacr/{sample}_{celltype}_{antibody}.relaxed.bed",
  log:
    "logs/seacr/{sample}_{celltype}_{antibody}.seacr.log"
  conda:
    "../envs/seacr.yaml"
  params:
    outid="results/peaks/seacr/{sample}_{celltype}_{antibody}",
    control=has_a_control,
    threshold=config['params']['seacr']['signal_threshold'],
    norm=config['params']['seacr']['normalization'],
    stringency=config['params']['seacr']['stringency']
  shell:
    # SEACR_1.3.sh CTCF_DE_chr1_100Mb.bedgraph.txt IgG_DE_chr1_100Mb.bedgraph.txt norm stringent CTCF_DE_chr1_100Mb
    # SEACR_1.3.sh CTCF_DE_chr1_100Mb.bedgraph.txt 0.05 norm stringent CTCF_DE_chr1_thresh
    "echo '{input.sample} - {input.control}'; "
    "if [ '{params.control}' == 'True' ]; then "
        "echo 'SEACR: ctrl bedGraph ' > {log}; "
        "SEACR_1.3.sh {input.sample} {input.control} "
        "{params.norm} relaxed {params.outid}; "
    "else "
        "echo 'SEACR: threshold ({params.threshold}) ' > {log}; "
        "SEACR_1.3.sh {input.sample} {params.threshold} "
        "{params.norm} relaxed {params.outid}; "
    "fi "

rule seacr_single:
  input:
    sample="results/alignment/bed/{sample}_{celltype}_{antibody}.bedGraph",
  output:
    "results/peaks/seacr_single/{sample}_{celltype}_{antibody}." + config['params']['seacr']['stringency'] + ".bed",
  log:
    "logs/seacr/{sample}_{celltype}_{antibody}.seacr_single.log"
  conda:
    "../envs/seacr.yaml"
  params:
    outid="results/peaks/seacr_single/{sample}_{celltype}_{antibody}",
    threshold=config['params']['seacr']['signal_threshold'],
    norm=config['params']['seacr']['normalization'],
    stringency=config['params']['seacr']['stringency']
  shell:
    # SEACR_1.3.sh CTCF_DE_chr1_100Mb.bedgraph.txt IgG_DE_chr1_100Mb.bedgraph.txt norm stringent CTCF_DE_chr1_100Mb
    # SEACR_1.3.sh CTCF_DE_chr1_100Mb.bedgraph.txt 0.05 norm stringent CTCF_DE_chr1_thresh
    "echo '{input.sample} - no control mode'; "
    "echo 'SEACR: threshold ({params.threshold}) ' > {log}; "
    "SEACR_1.3.sh {input.sample} {params.threshold} "
    "{params.norm} {params.stringency} {params.outid}; "

rule macs2_narrow:
  input:
    treatment="results/alignment/dedup/{sample}_{celltype}_{antibody}.bam",   # required: treatment sample(s)
    control=lambda w: expand("results/alignment/dedup/{ctrl}.bam",
                             ctrl=get_sample_control(w, retself=False))  if get_sample_control(w, retself=False) else ""     # optional: control sample(s)
  output:
    #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/macs2/callpeak.html.
    multiext("macs2/narrow/{sample}_{celltype}_{antibody}",
      "_peaks.xls",   ### required
      "_control_lambda.bdg",
      "_treat_pileup.bdg",
      "_peaks.narrowPeak",
      "_summits.bed"
      )
  log:
    "logs/macs2/{sample}_{celltype}_{antibody}.narrow_callpeak.log"
  params:
    #  --SPMR ask MACS2 to generate pileup signal file of 'fragment pileup per million reads' in bedGraph format.
    #  -g lets MACS2 consider 'hm' (human) or 'mm' (mouse) genome as background
    #  --nomodel will bypass building the shifting model (fragment size estimation)
    #  --extsize sets specifies fragment size related to binding region for your transcription factor for pileup of sequencing reads
    "-f BAM -g {} -i {} --nomodel --extsize {} --SPMR -q {}".format(config["common"]["species"],
      "{wildcards.sample}_{wildcards.celltype}_{wildcards.antibody}",
      config["params"]["macs2"]["extsize"],
      config['params']['macs2']['q'])
  wrapper:
    "0.77.0/bio/macs2/callpeak"

rule macs2_broad:
  input:
    treatment="results/alignment/dedup/{sample}_{celltype}_{antibody}.bam",   # required: treatment sample(s)
    control=lambda w: expand("results/alignment/dedup/{ctrl}.bam",
                             ctrl=get_sample_control(w, retself=False)) if get_sample_control(w, retself=False) else ''     # optional: control sample(s)
  output:
    #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/macs2/callpeak.html.
    multiext("macs2/broad/{sample}_{celltype}_{antibody}",
      "_peaks.xls",   ### required
      "_peaks.broadPeak",
      "_summits.gappedPeak"
      )
  log:
    "logs/macs2/{sample}_{celltype}_{antibody}.broad_callpeak.log"
  params:
    #  --SPMR ask MACS2 to generate pileup signal file of 'fragment pileup per million reads' in bedGraph format.
    #  -g lets MACS2 consider 'hm' (human) or 'mm' (mouse) genome as background
    #  --nomodel will bypass building the shifting model (fragment size estimation)
    #  --extsize sets specifies fragment size related to binding region for your transcription factor for pileup of sequencing reads
    "-f BAM -g {} -i {} --nomodel --extsize {} --SPMR -q {}".format(config["common"]["species"],
      "{wildcards.sample}_{wildcards.celltype}_{wildcards.antibody}",
      config["params"]["macs2"]["extsize"],
      config['params']['macs2']['q'])
  wrapper:
    "0.77.0/bio/macs2/callpeak"

rule chipseeker:
  input:
    get_all_peak_files,
  output:
    tss_cov='results/peaks/chipseeker/tss_cov.pdf',
    tss_cov_avg='results/peaks/chipseeker/tss_cov_avg.pdf',
    tagmat_deep='results/peaks/chipseeker/tagmat_deep.pdf',
    peak_anno='results/peaks/chipseeker/peak_anno.pdf',
    functional_anno="results/peaks/chipseeker/functional_anno.pdf",
    peak_ov='results/peaks/chipseeker/peak_ov.pdf',
  log:
    "logs/peak_analysis/chipseeker.log"
  conda:
    "../envs/chipseeker.yaml"
  params:
    files=lambda w, input: ','.join(input),
    build=config['common']['build'],
    gtf=config['common']['gtf'],
    promoter_range=config['peakanalyze']['promoter_range'],
  script:
    "../scripts/chipseeker.R"

rule homer_annotatepeaks:
  input:
    peaks="results/peaks/seacr/{sample}_{celltype}_{antibody}." + config['params']['seacr']['stringency'] + ".bed",
    genome=config['common']['genome'],
    gtf=config['common']['gtf']
  output:
    annotations="results/homer/annotate_peaks/{sample}_{celltype}_{antibody}.seacr_peaks.annotatePeaks.txt"
  threads:
    2
  params:
    mode="",
    extra="-gid"
  log:
    "logs/homer/annotate_peaks/{sample}_{celltype}_{antibody}.log"
  wrapper:
    "v0.75.0//bio/homer/annotatePeaks"

'''
rule plot_homer_annotatepeaks:
  input:
    get_plot_homer_annotatepeaks_input()
  output:  #ToDo: add description to report caption
    summmary="results/homer/plots/plot_{peak}_annotatepeaks_summary.txt",
    plot=report("results/homer/plots/plot_{peak}_annotatepeaks.pdf", caption="../report/plot_annotatepeaks_homer.rst", category="CallPeaks")
  params:
    input = lambda wc, input: ','.join(input),
    sample_control_combinations = ','.join(get_sample_control_peak_combinations_list())
  log:
    "logs/homer/plot_{peak}_annotatepeaks.log"
  conda:
    "../envs/plot_macs_annot.yaml"
  shell:
    "Rscript ../workflow/scripts/plot_homer_annotatepeaks.R "
    "-i {params.input} "
    "-s {params.sample_control_combinations}  "
    "-o {output.plot} "
    "-p {output.summmary} 2> {log}"

rule memesuite_jaspar:
  input:

  output:
  params:
  log:
    "logs/memesuite/{sample}_{celltype}_{antibody}.jaspar2.log"
  conda:
    "../envs/peakanalyze.yaml"
  shell:
    "centrimo "
'''

#module load meme/5.1.0
#centrimo Catalogue-CD8pos-200.fa /mnt/work1/users/pughlab/references/MEMESuite/motif_databases/JASPAR/JASPAR_CORE_2016.meme \
#--oc MEME_KOMBAT-CD8pos_JASPAR2 \
#--neg  /mnt/work1/users/pughlab/projects/KOMBAT/atacseq_analysis/Sorted_ATACSeq_hg19/Peaks_KOM-Immune_BC/ATAC_ERposPRpos/RoadmapEncode_Catalogue_200.fa \
#--use-pvalues --seqlen 200
