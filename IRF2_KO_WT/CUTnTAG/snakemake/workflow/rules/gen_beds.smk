rule bam_to_bed:
  input:
    bam="results/alignment/dedup/{sample}_{celltype}_{antibody}.bam",
    idx="results/alignment/dedup/{sample}_{celltype}_{antibody}.bam.bai",
  output:
    bed=temp("results/alignment/bed/{sample}_{celltype}_{antibody}.bed"),
    clean=temp("results/alignment/bed/{sample}_{celltype}_{antibody}.clean.bed"),
    frags="results/alignment/bed/{sample}_{celltype}_{antibody}.frags",
  log:
    "logs/bed/{sample}_{celltype}_{antibody}.bed.log"
  conda:
    "../envs/bamtools.yaml"
  params:
    bedtools=config["params"]['bam_to_bed']['bedtools'],
    awk=config["params"]['bam_to_bed']['awk'],
    cut=config["params"]['bam_to_bed']['cut']['cut_arg'],
    sort=config["params"]['bam_to_bed']['cut']['sort_arg'],
  shell:
    # [.bed]: Returns the intervals for R1 and R2
    "bedtools bamtobed -i {input.bam} {params.bedtools} > {output.bed}; "
    # [.clean.bed]: Returns rows where both R1/2 are on same chromosome and span < 1000bp
    "awk {params.awk} {output.bed} > {output.clean}; "
    # [.frags]: Returns the (Chr, R1_start, R2_end)
    "cut {params.cut} {output.clean} | sort {params.sort} > {output.frags}; "

rule samtools_stats:
  input:
    "results/alignment/dedup/{sample}_{celltype}_{antibody}.bam",
  output:
    "results/alignment/dedup/stats/{sample}_{celltype}_{antibody}.stats",
  params:
    extra=""
  log:
    "logs/samtools/stats/{sample}_{celltype}_{antibody}.stats.log"
  wrapper:
    "v0.75.0/bio/samtools/stats"

rule bam_to_bedgraph:
  input:
    spikein="results/alignment/spikein/{sample}_{celltype}_{antibody}.spikein.sorted.bam",
    bam="results/alignment/dedup/{sample}_{celltype}_{antibody}.bam",
    sampleidx="results/alignment/dedup/{sample}_{celltype}_{antibody}.bam.bai",
    stats="results/alignment/dedup/stats/{sample}_{celltype}_{antibody}.stats",
  output:
    bedgraph="results/alignment/bed/{sample}_{celltype}_{antibody}.unfiltered.bedGraph",
  log:
    "logs/bed/{sample}_{celltype}_{antibody}.bedGraph.log",
  conda:
    "../envs/bamtools.yaml",
  params:
    constant=config['params']['bedgraph']['constant'],
  shell:
    #"spikein_depth=`samtools view -F 0x04 {input.spikein} | wc -l`; "
    #"seqDepth=$((spikein_depth/2)); "
    #'scale_factor=`echo "{params.constant} / $seqDepth" | bc -l`; '
    'reads=$(grep "mapped and paired" {input.stats}  | sed "s/^.*:\\\s*//" | sed "s/\\\s.*//"); '
    "echo 'Reads = '$reads; "
    "scale_factor=$(echo $reads | perl -ne 'print {params.constant}/($_/2);'); "
    "echo 'Scale factor = '$scale_factor; "
    "bedtools genomecov -ibam {input.bam} -bg -scale $scale_factor | "
    "bedtools sort > {output.bedgraph} "

rule emit_sizes:
  input:
    config["common"]["index"],
  output:
    "resources/" + config['common']['build'] + ".sizes",
  log:
    "logs/resources/chrsizes.txt"
  shell:
    "cut -f 1,2 {input} > {output}"

rule bedclip:
  input:
    bed="results/alignment/bed/{sample}_{celltype}_{antibody}.unfiltered.bedGraph",
    sizes=rules.emit_sizes.output,
  output:
    "results/alignment/bed/{sample}_{celltype}_{antibody}.clip.bedGraph",
  log:
    "logs/bed/{sample}_{celltype}_{antibody}.bedclip.log"
  conda:
    "../envs/bamtools.yaml"
  params:
    ""
  shell:
    "bedClip {input.bed} {input.sizes} {output}"


rule blacklist_filter:
  input:
    left="results/alignment/bed/{sample}_{celltype}_{antibody}.clip.bedGraph",
    right=config['common']['blacklist'],
  output:
    "results/alignment/bed/{sample}_{celltype}_{antibody}.bedGraph",
  log:
    "logs/filter/{sample}_{celltype}_{antibody}.blacklist.log"
  params:
    ## Add optional parameters
    extra="-v" ## exclude all blacklist regions from the bed file
  wrapper:
    "v0.75.0/bio/bedtools/intersect"

rule bedgraph_to_bigwig:
  input:
    bed="results/alignment/bed/{sample}_{celltype}_{antibody}.bedGraph",
    sizes=rules.emit_sizes.output,
  output:
    "results/alignment/bed/{sample}_{celltype}_{antibody}.bigWig",
  log:
    "logs/bed/{sample}_{celltype}_{antibody}.bedwig.log"
  conda:
    "../envs/bamtools.yaml"
  params:
    ""
  shell:
    "bedGraphToBigWig {input.bed} {input.sizes} {output}"

rule bam_to_bigwig:
  input:
    spikein="results/alignment/spikein/{sample}_{celltype}_{antibody}.spikein.sorted.bam",
    bam="results/alignment/dedup/{sample}_{celltype}_{antibody}.bam",
    sampleidx="results/alignment/dedup/{sample}_{celltype}_{antibody}.bam.bai",
  output:
    "results/alignment/bed/{sample}_{celltype}_{antibody}.deeptools.bigWig",
  log:
    "logs/bed/{sample}_{celltype}_{antibody}.bamcoverage.log"
  conda:
    "../envs/bamtools.yaml"
  params:
    constant=config['params']['bedgraph']['constant'],
  shell:
    "spikein_depth=`samtools view -F 0x04 {input.spikein} | wc -l`; "
    "seqDepth=$((spikein_depth/2)); "
    'scale_factor=`echo "{params.constant} / $seqDepth" | bc -l`; '
    "echo 'Scale factor = '$scale_factor; "
    "bamCoverage -b {input.bam} --scale $scale_factor --binSize 10 -o {output}"

#rule calculate_scale_factor:
#  input:
#  output:
#  log:
#  params:
#  shell:
